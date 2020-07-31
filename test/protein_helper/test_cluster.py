from importlib.resources import path

from protein_helper.cluster import get_cdhit_cluster_sizes
from test.fixtures import cdhit_cluster_sizes


def test_get_cdhit_cluster_sizes(tmp_path):
    expected_size_tups = [
        (40, 123), (45, 197), (50, 280), (55, 368), (60, 470), (65, 552), (70, 666), (75, 804),
        (80, 959), (85, 1119), (90, 1311), (95, 1595), (100, 2222),
    ]
    with path(cdhit_cluster_sizes, "input.fa") as fasta:
        size_tups = get_cdhit_cluster_sizes(input_fasta=fasta, output_dir=tmp_path)
    assert expected_size_tups == size_tups
