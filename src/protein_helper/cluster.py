from collections import namedtuple, Generator
import os
from typing import List

import matplotlib.pyplot as plt

from protein_helper import cluster_tools

CdhitCluster = namedtuple('CdhitCluster', ['cluster_id', 'proteins', 'representative'])


def iter_cdhit_clusters(clstr_handle):
    """
    Returns a Generator that yields a CdhitCluster namedtuple.
    """
    def _sequence_name(cstr_line):
        return cstr_line.split(">")[1].rsplit("...", 1)[0]

    def _is_representative(cstr_line):
        return cstr_line[-1] == '*'

    def _cluster_tuple(cluster_id, cluster_lines):
        return CdhitCluster(
            cluster_id=cluster_id,
            proteins=[_sequence_name(line) for line in cluster_lines],
            representative=next(
                _sequence_name(line) for line in cluster_lines if not _is_representative(line))
        )

    cluster_lines = []
    cluster_id = None
    for line in clstr_handle:
        if line.startswith('>Cluster'):
            if cluster_lines:
                yield _cluster_tuple(cluster_id, cluster_lines)
            cluster_lines = []
            cluster_id = line.strip('>\n')
        else:
            cluster_lines.append(line)
    if cluster_lines:
        yield _cluster_tuple(cluster_id, cluster_lines)


def get_cdhit_clusters(
        input_fasta: str,
        percent_identity: float,
        length_difference_cutoff: float = None,
        min_alignment_coverage: float = None,
        percent_identity_suffix: bool = False,
        output_dir: str = None,
) -> Generator:
    """Runs cd-hit when clstr file isn't present and parses cd-hit output from the cstr file

    Args:
        fasta: Input fasta file
        percent_identity: Minimum percent identity for edge inclusion.
        length_difference_cutoff: Sequences need to be at least this percent length of the
            representative sequence.
        min_alignment_coverage: Alignment must cover at least this percent of both sequences.
        percent_identity_suffix: Include percent identity in the filenames of the cdhit output.
        output_dir: If provided, cd-hit files will be read and written from this directory

    Returns:
        A Generator that yields CdhitClusters

    Raises:
        CalledProcessError: If the subprocess running the cd-hit program can not complete
            successfully.
    """
    clstr_filepath = get_cluster_filepath(
        input_fasta=input_fasta, percent_identity=percent_identity, output_dir=output_dir
    )

    if not os.path.exists(clstr_filepath):
        cluster_tools.cdhit(
            input_fasta=input_fasta,
            percent_identity=percent_identity,
            length_difference_cutoff=length_difference_cutoff,
            min_alignment_coverage=min_alignment_coverage,
            percent_identity_suffix=percent_identity_suffix,
            output_dir=output_dir,
        )

    yield from iter_cdhit_clusters(open(clstr_filepath))


def get_cluster_filepath(
        input_fasta: str,
        percent_identity: float,
        output_dir: str = None
) -> str:
    # TODO: generalize this function by having it accept a parameter percent_identity_suffix, so
    #  that this can be used to interpolate the filename when the percend identity suffix is desired
    if output_dir is not None:
        root = os.path.join(output_dir, os.path.splitext(os.path.basename(input_fasta))[0])
    else:
        root = os.path.splitext(input_fasta)[0]
    return f'{root}{percent_identity}.clstr'


def get_cdhit_cluster_sizes(
        input_fasta: str,
        start_percent_identity: int = 40,
        end_percent_identity: int = 100,
        step: int = 5,
        length_difference_cutoff: float = 0.1,
        min_alignment_coverage: float = 0.6,
        output_dir: str = None
) -> List[tuple]:
    """
    TODO: Finish docstring
    """
    return [
        (
            pid,
            len(list(get_cdhit_clusters(
                input_fasta=input_fasta,
                percent_identity=pid/100,
                length_difference_cutoff=length_difference_cutoff,
                min_alignment_coverage=min_alignment_coverage,
                percent_identity_suffix=True,
                output_dir=output_dir,)))
        )
        for pid in range(start_percent_identity, end_percent_identity + 1, step)
    ]


def generate_cdhit_cluster_number_plot(
        cluster_counts: List[tuple],
        output_png: str,
) -> None:
    """Generate a plot of percent identity on the x-axis versus cluster count on the y-axus

    Args:
        cluster_counts: A list of tuples. The first element is the perecent identity. Teh second
            element is the count clusters at that percent identity.
        output_png: Path to output plot.
    """
    fig, ax = plt.subplots()
    percent_identity, cluster_sizes = zip(*cluster_counts)
    ax.plot(percent_identity, cluster_sizes, marker='o', markersize=4)
    ax.set(xlabel='percent identity', ylabel='number of cd-hit cluster',
           title='cd-hit cluster size by percent identity')
    ax.grid()
    fig.savefig(output_png)
