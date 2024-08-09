#!/usr/bin/env python
import fileinput
import argparse
import subprocess
import tempfile
from pathlib import Path
from typing import Iterable
from collections import defaultdict

import pysam
from dark.bowtie2 import Bowtie2

from filter_spurious_reads import filter_spurious_reads


def makeParser() -> argparse.ArgumentParser:
    _parser = argparse.ArgumentParser(
        description="Reduce duplicate reads by calling a consensus."
    )
    _parser.add_argument("samfile", metavar="FILE.sam", help="The SAM file to read.")
    _parser.add_argument(
        "reference",
        metavar="REFERENCE.fa",
        help="A FASTA file containing the reference sequence.",
    )
    _parser.add_argument(
        "-t",
        "--target",
        metavar="TARGETDIR",
        help="The directory to write the consensus sequences to.",
    )
    _parser.add_argument(
        "--cluster-size",
        metavar="CLUSTERSIZE",
        help="The maximum distance from a read to the center of a cluster to be included in the cluster.",
        type=int,
        default=10,
        dest="cluster_size",
    )
    _parser.add_argument(
        "--min-n-reads",
        metavar="MIN_N_READS",
        help="The minimum number of reads in a cluster to be included in filtered output.",
        type=int,
        default=3,
        dest="min_n_reads",
    )
    _parser.add_argument(
        "--max-magnitude-difference",
        metavar="MAXMAGNITUDE",
        help="The maximum difference in magnitude between the most frequent and the second most frequent fragment to be included in filtered output.",
        type=float,
        default=1,
        dest="max_magnitude_difference",
    )
    _parser.add_argument(
        "-F",
        "--no-filtered-output",
        help="Use this option to not generate output with reads filtered.",
        default=True,
        dest="filtered_output",
        action="store_false",
    )
    _parser.add_argument(
        "-U",
        "--no-unfiltered-output",
        default=True,
        dest="unfiltered_output",
        action="store_false",
        help="Use this option to not generate output for all reads.",
    )
    return _parser


def sort_reads_by_position(
    alignment_iterable: Iterable[pysam.AlignedSegment],
) -> dict[tuple[int, int], list[pysam.AlignedSegment]]:
    """Takes an iterable returning AlignedSegment (such as an AlignmentFile object)
    and returns a dict grouping the reads by start and end position."""
    reads_dict = defaultdict(list)
    for read in alignment_iterable:
        # pysam returns 0-based coordinates, but the standard uses 1-based coordinates
        reads_dict[(read.pos + 1, read.aend + 1)].append(read)  # pyright: ignore
    return reads_dict


def summarize_reads_dict(reads_dict: dict) -> None:
    """Takes a dict of reads and summarizes it."""
    for (start, end), reads in reads_dict.items():
        print(f"Number of reads at {start}-{end}: {len(reads)}")


def deduplicate_samfile(
    samfile_name: str,
    target_dir: Path,
    filtered_output: bool,
    unfiltered_output: bool,
    cluster_size: int,
    min_n_reads: int,
    max_magnitude_difference: float,
) -> tuple[list[Path], list[Path]]:
    """
    Deduplicate the reads in a SAM file by calling consensuses on reads that start
    and end at the same position. Such collections of reads are called fragments.

    - filtered_output: if True, fragments will be filtered by a minimum number of
    reads and clustered by combining fragments that are within cluster_size of each
    other. The fragments will then be filtered by the maximum difference in magnitude
    of the number of reads in the fragment.

    - unfiltered_output: if True, all fragments will be deduplicated.
    """
    with (pysam.AlignmentFile(samfile_name, "r") as samfile,):
        fastq_dir = Path(target_dir, "fastq")
        fastq_dir.mkdir(parents=True, exist_ok=True)

        # get all the reads in the samfile
        unfiltered_reads_dict = sort_reads_by_position(samfile)
        filtered_file_list = []
        unfiltered_file_list = []
        if filtered_output:
            filtered_reads_list = filter_spurious_reads(
                unfiltered_reads_dict,
                cluster_size,
                min_n_reads,
                max_magnitude_difference,
            )
            filtered_reads_dict = sort_reads_by_position(filtered_reads_list)
            if not unfiltered_output:
                # If we only want the filtered output, we only need to loop
                # over the filtered reads.
                reads_dict = filtered_reads_dict
        if unfiltered_output:
            reads_dict = unfiltered_reads_dict

        # For each fragment in reads_dict: write the fragment to a file.
        # If there is more than one read for that fragment use `samtools`
        # to create a consensus for that fragment.
        with tempfile.TemporaryDirectory() as tmp_bam_dir:
            for (start, end), reads in reads_dict.items():  # pyright: ignore
                # Fragments are identified by their start and end position and
                # the number of reads in the fragment.
                fragment_name = f"{len(reads)}-reads-at-{start}-{end}"

                # BAM file to write the consensus for the fragment to.
                tmp_bamfilename = Path(
                    tmp_bam_dir, fragment_name+".bam"
                )
                # Output FASTQ file for the consensus.
                fastq_filename = fastq_dir / (fragment_name+".fastq")

                # Write the reads to the temporary BAM file and count them.
                # Since the next step expects a file this file is
                # created regardless of whether there are no reads.
                with pysam.AlignmentFile(
                    str(tmp_bamfilename), "w", header=samfile.header
                ) as tmp_bamfile:
                    n_reads = 0
                    for read in reads:
                        tmp_bamfile.write(read)
                        n_reads += 1

                if n_reads > 1:
                    # If there is more than one read in the fragment, create a consensus.
                    # Start by sorting the reads by position.
                    pysam.sort(str(tmp_bamfilename), "-o", str(tmp_bamfilename))
                    # Create a consensus sequence from the reads in the fragment.
                    # the `-l 0` option tells pysam that there is no maximum
                    # line length, i.e. the sequence should be written on one line.
                    # The `-f fastq` option tells pysam to write the consensus
                    # in FASTQ format.
                    pysam.consensus(
                        str(tmp_bamfilename),
                        "-l",
                        "0",
                        "-f",
                        "fastq",
                        "-o",
                        str(fastq_filename),
                    )

                    # Because pysam uses the same name for all sequences we need to
                    # replace the read name with the fragment name.
                    # Note that we can do this in place and using this very
                    # basic parsing-method because the BAM file contains only one
                    # read, so the name is guaranteed to be in the first line.
                    with fileinput.input(files=fastq_filename, inplace=True) as fh:
                        for i, line in enumerate(fh):
                            if i == 0:
                                line = "@" + fragment_name + "\n"
                            print(line, end="")
                else:
                    # If there is only one read in the fragment, we only need to
                    # convert from BAM to FASTQ.
                    pysam.fastq(
                        str(tmp_bamfilename),
                        "-0",
                        str(fastq_filename),
                        "-o",
                        str(fastq_filename),
                        )
                if (
                    filtered_output and (start, end) in filtered_reads_dict  #pyright: ignore
                ):
                    filtered_file_list.append(fastq_filename)
                if unfiltered_output:
                    unfiltered_file_list.append(fastq_filename)

    return filtered_file_list, unfiltered_file_list


def align_with_bowtie(
    reference_file: Path,
    sequences_file: Path,
    target_dir: Path,
    output_file_name: str | Path = "dedup.sam",
    bt2args="--xeq --local --very-sensitive-local -N 1 --no-unal",
) -> Path:
    """
    Align the sequences in the FASTQ file to the reference using Bowtie2 and
    write the output to the `output_file_name` in the target directory.
    """
    bt = Bowtie2()
    bt.buildIndex(str(reference_file.absolute()))
    bt.align(fastq1=str(sequences_file.absolute()), bowtie2Args=bt2args)

    output_file = target_dir / output_file_name
    subprocess.run(["mv", bt.outputFile(), str(output_file.absolute())], check=True)
    return output_file


def concatenate_files(
    target_dir: Path,
    list_of_files: Iterable[str | Path],
    output_file_name: str | Path,
) -> Path:
    """
    Concatenate the content of the files in `list_of_files` and write the output to
    `output_file_name` in the target directory.
    Return the path to the output file.
    """
    output_file = target_dir / output_file_name
    with open(output_file, "w") as output_fh:
        for filename in list_of_files:
            with open(filename, "r") as input_fh:
                output_fh.write(input_fh.read())
    return output_file


def main():
    parser = makeParser()
    args = parser.parse_args()
    reference_file = Path(args.reference)
    # Set target directory either to the given argument or to cwd. Create it if necessary.
    if args.target is None:
        target_dir = Path.cwd()
    else:
        target_dir = Path(args.target)
        if not target_dir.exists():
            print(f"Creating target directory {target_dir}")
            target_dir.mkdir(parents=True)

    # Deduplicate the reads in the SAM file using consensus deduplication.
    # The deduplicated reads are written to the fastq dir in the target directory.
    filtered_list, unfiltered_list = deduplicate_samfile(
        args.samfile,
        target_dir,
        args.filtered_output,
        args.unfiltered_output,
        args.cluster_size,
        args.min_n_reads,
        args.max_magnitude_difference,
    )

    # Next we concatenate the deduplicated reads and align them to the reference.
    if filtered_list:
        sequences_file = concatenate_files(
            target_dir, filtered_list, "filtered_deduplicated_sequences.fq"
        )
        align_with_bowtie(
            reference_file,
            sequences_file,
            target_dir,
            "filtered_deduplicated.sam",
            "--xeq --local --very-sensitive-local -N 1 --no-unal",
        )
    if unfiltered_list:
        sequences_file = concatenate_files(
            target_dir, unfiltered_list, "all_deduplicated_sequences.fq"
        )
        align_with_bowtie(
            reference_file,
            sequences_file,
            target_dir,
            "deduplicated.sam",
            "--xeq --local --very-sensitive-local -N 1 --no-unal",
        )


if __name__ == "__main__":
    main()
