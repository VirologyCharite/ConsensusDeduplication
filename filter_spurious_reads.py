from collections import defaultdict
import numpy as np
import pysam

# Some clarification on the types on this module:
# pysam.AlignedSegment is the type of the reads,
# a list of pysam.AlignedSegment is a fragment,
# a list of fragment (ie a list of list of pysam.AlignedSegment) is a cluster.
# When we cluster the reads we generate a list of clusters, which is a list of list of list of pysam.AlignedSegment.


class Fragment:
    def __init__(self, start: int, end: int, reads: list[pysam.AlignedSegment]):
        self.start = start
        self.end = end
        self.reads = reads

    def __len__(self):
        return len(self.reads)


def cluster_reads(reads_dict: dict, cluster_size: int) -> list[list[Fragment]]:
    cluster_list: list[list[Fragment]] = []
    """
    Clusters reads based on their start and end positions. The reads are
    popped from the reads_dict as they are clustered.
    The centers are found iteratively by finding the read with the smallest
    total distance to all other reads in the cluster.
    The cluster_size is the maximum L1 distance between reads in a cluster and its center.
    """
    while reads_dict:
        positions = np.array(sorted(reads_dict))
        X = positions - positions[:, np.newaxis]
        # X is a #reads x #reads x 2 array of differences between read start and
        # end positions. I.e. the following holds true
        # for i in range(len(positions)):
        #     assert np.all(np.isclose(X[i], positions - positions[i]))
        # Next we convert these to distances between the reads in the taxi-cab metric (L1).
        distances = np.sum(np.abs(X), axis=2)
        # And then find the total distance from each read to all other reads.
        total_distances = np.sum(distances, axis=1)
        # And finally find the read with the smallest total distance.
        center = np.argmin(total_distances)
        # Next we find all the reads that are within cluster_size of center
        # and pop them.
        reads_to_be_merged = []
        for start, end in positions[np.where(distances[center] < cluster_size)]:
            reads_to_be_merged.append(Fragment(start, end, reads_dict.pop((start, end))))
        cluster_list.append(reads_to_be_merged)
    return cluster_list


def sort_by_magnitude(
    cluster: list[Fragment], max_difference: float = 1
) -> list[pysam.AlignedSegment]:
    """
    Combines and returns in a list all the reads from those fragments in a
    cluster that, in terms of read count, are within max_difference of the
    largest fragment in terms of order of magnitude (i.e. log base 10).
    """
    # Determine the magnitude (log_10) of the number of reads for each fragment
    # and sort the fragments by magnitude.
    magnitudes = [(np.log10(len(f)), f.reads) for f in cluster]
    magnitudes = sorted(magnitudes, key=lambda x: x[0], reverse=True)

    # Extend the largest fragment with the reads of all smaller fragments
    # that are within max_difference of it.
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
    """
    Takes a dict of fragments and filters out the spurious reads.

    First, fragments with less than min_n_reads reads are discarded.
    Then, the fragments are clustered by combining fragments that are within
    cluster_size of each other as measured by the L1 distance between their
    start and end positions.
    Within each cluster, the fragments are then filtered by order of magnitude
    of their reads, with those fragments that are negligible being discarded.

    The reads from the remaining fragments are then returned in a list.
    """
    filtered_reads_dict = defaultdict(list)
    # Throw away fragments with less than min_n_reads reads.
    for (start, stop), fragment in reads_dict.items():
        if len(fragment) > min_n_reads:
            filtered_reads_dict[(start, stop)] = fragment

    all_clusters = cluster_reads(filtered_reads_dict, cluster_size)

    # For each cluster, filter by read count order of magnitude.
    # Join the reads of the remaining fragments from all the clusters in a list.
    filtered_reads_list = []
    for cluster in all_clusters:
        filtered_reads_list.extend(
            sort_by_magnitude(cluster, max_magnitude_difference)
        )
    return filtered_reads_list


def summarize_cluster_list(cluster_list: list) -> None:
    for cluster in cluster_list:
        print(f"There are {len(cluster)} fragments in this cluster.")
        for reads in cluster:
            print(f"\tThere are {len(reads)} reads in this fragment.")
            print(f"\tAt position {reads[0].reference_start}-{reads[0].reference_end}.")
        print()
