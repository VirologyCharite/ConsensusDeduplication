from collections import defaultdict
import numpy as np
import pysam
import sys

sys.setrecursionlimit(10**8)

# Some clarification on the types on this module:
# pysam.AlignedSegment is the type of the reads,
# a list of pysam.AlignedSegment is a fragment,
# a list of fragments (i.e. a list of list of pysam.AlignedSegment) is a cluster.
# When we cluster the reads we generate a list of clusters, which is a list of
# list of list of pysam.AlignedSegment.


def cluster_reads(
    reads_dict: dict, cluster_size: int
) -> list[list[list[pysam.AlignedSegment]]]:
    cluster_list: list[list[list[pysam.AlignedSegment]]] = []
    if len(reads_dict) == 0:
        return cluster_list
    positions = list(reads_dict.keys())
    positions.sort()
    positions_array = np.array(positions)
    X = positions_array - positions_array[:, np.newaxis]
    # X is a #reads x #reads x 2 array of differences between read start and end positions
    # i.e. the following holds true
    # for i in range(len(positions_array)):
    #     assert np.all(np.isclose(X[i], positions_array - positions_array[i]))
    # Next we convert these to distances between the reads in the taxi-cab metric
    distances = np.abs(np.sum(X, axis=2))
    # And then find the total distance from each read to all other reads
    total_distances = np.sum(distances, axis=1)
    # And finally find the read with the smallest total distance
    center = np.argmin(total_distances)
    # Next we find all the reads that are within cluster_size of center
    # and pop them
    reads_to_be_merged = []
    for start, end in positions_array[np.where(distances[center] < cluster_size)]:
        # print(start,end, len(reads_dict[start,end]))
        reads_to_be_merged.append(reads_dict.pop((start, end)))
    # print()
    cluster_list.append(reads_to_be_merged)
    cluster_list.extend(cluster_reads(reads_dict, cluster_size))
    return cluster_list


def sort_by_magnitude(
    cluster: list[list[pysam.AlignedSegment]], max_difference: float = 1
) -> list[pysam.AlignedSegment]:
    """Returns those fragment from a cluster of fragment that are within
    max_difference of the largest fragment in terms of order of magnitude
    (i.e. log base 10)."""
    # Determine the magnitude (log_10) of the number of reads for each fragment
    magnitudes = [(np.log10(len(reads)), reads) for reads in cluster]
    magnitudes = sorted(magnitudes, key=lambda x: x[0], reverse=True)
    output = magnitudes[0][1]
    for magnitude, reads in magnitudes[1:]:
        if magnitudes[0][0] - magnitude < max_difference:
            output.extend(reads)
    return output


def filter_spurious_reads(
    reads_dict: dict[tuple[int, int], list[pysam.AlignedSegment]],
    cluster_size: int,
    min_n_reads: int,
    max_magnitude_difference: float,
) -> list[pysam.AlignedSegment]:
    filtered_reads_dict = defaultdict(list)
    # Throw away fragment with less than min_n_reads reads
    for (start, stop), reads in reads_dict.items():
        if len(reads) > min_n_reads:
            filtered_reads_dict[(start, stop)] = reads
    cluster_list = cluster_reads(filtered_reads_dict, cluster_size)
    # Filter by order of magnitude and split up clusters
    filtered_cluster_list = []
    for cluster in cluster_list:
        filtered_cluster_list.extend(
            sort_by_magnitude(cluster, max_magnitude_difference)
        )
    return filtered_cluster_list


def summarize_cluster_list(cluster_list: list) -> None:
    for cluster in cluster_list:
        print(f"There are {len(cluster)} fragments in this cluster.")
        for reads in cluster:
            print(f"\tThere are {len(reads)} reads in this fragment.")
            print(f"\tAt position {reads[0].reference_start}-{reads[0].reference_end}.")
        print()
