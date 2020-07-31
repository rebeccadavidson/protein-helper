from importlib.resources import path
import os

from protein_helper.visualization import generate_network
from test.fixtures import blastp_all_by_all


def test_generate_cytoscape_network(tmp_path):
    with path(blastp_all_by_all, "input.fasta") as input_fasta:
        fasta_prefix = os.path.basename(os.path.splitext(input_fasta)[0])
        cytoscape_network_path = os.path.join(tmp_path, f'{fasta_prefix}.cyjs')
        output_plot_path = os.path.join(tmp_path, f'{fasta_prefix}.png')
        generate_network(
            fasta=input_fasta,
            cytoscape_network_path=cytoscape_network_path,
            output_plot_path=output_plot_path,
            temp_dir=tmp_path
        )


def test_generate_cytoscape_network_minimum_percent_identity_real_data(tmp_path):
    with path(blastp_all_by_all, "PF00135_seed.fasta") as input_fasta:
        fasta_prefix = os.path.basename(os.path.splitext(input_fasta)[0])
        cytoscape_network_path = os.path.join(tmp_path, f'{fasta_prefix}.cyjs')
        output_plot_path = os.path.join(tmp_path, f'{fasta_prefix}.png')
        generate_network(
            fasta=input_fasta,
            cytoscape_network_path=cytoscape_network_path,
            output_plot_path=output_plot_path,
            minimum_percent_identity=42,
            temp_dir=tmp_path
        )


def test_generate_cytoscape_network_real_data(tmp_path):
    with path(blastp_all_by_all, "PF00135_seed.fasta") as input_fasta:
        fasta_prefix = os.path.basename(os.path.splitext(input_fasta)[0])
        cytoscape_network_path = os.path.join(tmp_path, f'{fasta_prefix}.cyjs')
        output_plot_path = os.path.join(tmp_path, f'{fasta_prefix}.png')
        generate_network(
            fasta=input_fasta,
            cytoscape_network_path=cytoscape_network_path,
            output_plot_path=output_plot_path,
            temp_dir=tmp_path
        )
