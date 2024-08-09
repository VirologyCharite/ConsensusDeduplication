# Co(nsensus)De(duplication)
This repo can be used to deduplicate reads by consensus calling.
The resulting "reads" are then aligned against a reference sequence.

The consensus is called using the pysam (samtools) consensus caller and the alignment is done using bowtie2.

The reads must be in a sam file and the reference should be in a fasta file.
Note that bowtie2 does not like Uracil, so if you have a RNA reference you must substitute Ts for all the Us.

The script will also produce a "filtered" output which aims to remove spurious reads (e.g. single digit number of reads from a fragment that is only a few nucleotides short or longer than a fragment with thousands of reads).
This filtering is done by first throwing out fragment that have too few reads, then clustering the fragment by their start and end position and then for each cluster throwing out all the fragments that have an order of 10 fewer reads than the most numerous fragment.
(A fragment in this context means all the reads having the same start an end position.)

# Requirements
Requires python version 3.10, samtools (tested with 1.19.2, using htslib 1.19.1) and bowtie2 (tested with version 2.5.1) and requirements.txt.

# Usage
Evertyhing is done using the `sort_reads.py` script.
```shell
$ ./sort_reads.py -h                                                                                                                                             [9:33:02]
usage: sort_reads.py [-h] [-t TARGETDIR] [--cluster-size CLUSTERSIZE] [--min-n-reads MIN_N_READS] [--max-magnitude-difference MAXMAGNITUDE] [-F] [-U] FILE.sam REFERENCE.fa

Reduce duplicate reads by calling a consensus.

positional arguments:
  FILE.sam              The SAM file to read.
  REFERENCE.fa          A FASTA file containing the reference sequence.

options:
  -h, --help            show this help message and exit
  -t TARGETDIR, --target TARGETDIR
                        The directory to write the consensus sequences to.
  --cluster-size CLUSTERSIZE
                        The maximum distance from a read to the center of a cluster to be included in the cluster.
  --min-n-reads MIN_N_READS
                        The minimum number of reads in a cluster to be included in filtered output.
  --max-magnitude-difference MAXMAGNITUDE
                        The maximum difference in magnitude between the most frequent and the second most frequent fragment to be included in filtered output.
  -F, --no-filtered-output
                        Use this option to not generate output with reads filtered.
  -U, --no-unfiltered-output
                        Use this option to not generate output for all reads.

```
