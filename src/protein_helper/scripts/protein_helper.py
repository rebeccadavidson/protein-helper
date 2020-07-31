import click
from Bio import SearchIO
from click import Path, File

from protein_helper import (
    cluster,
    hmm_utils,
    utils,
    visualization,
)


@click.group()
def cli():
    pass


@cli.command('network')
@click.option(
    '--input-fasta',
    type=Path(exists=True, file_okay=True, dir_okay=True, readable=True, resolve_path=True,),
    required=True,
    help="Input fasta file")
@click.option(
    '--cytoscape-cyjs',
    type=Path(exists=False, file_okay=True, dir_okay=True, readable=True, resolve_path=True,),
    help="Output Cytoscape cys file")
@click.option(
    '--mcl-format',
    type=Path(exists=False, file_okay=True, dir_okay=True, readable=True, resolve_path=True,),
    help="Output file to use for mcl clustering input")
@click.option(
    '--mcl-edge-type',
    type=click.Choice(['percent_identity', 'bitscore', 'evalue'], case_sensitive=False),
    help="edge type to use for MCL clustering file")
@click.option(
    '--output-plot',
    type=Path(exists=False, file_okay=True, dir_okay=True, readable=True, resolve_path=True,),
    help="If provided, an output plot with a networx rendering of the graph will be provided as a"
         "png.")
@click.option(
    '--min-percent-identity',
    type=int, default=0,
    help="Edge type to use for MCL clustering file.")
@click.option(
    '--temp-dir',
    type=Path(exists=True, file_okay=True, dir_okay=True, readable=True, resolve_path=True,),
    help="Temp directory to place cd-hit files in")
@click.option(
    '--threads',
    type=int,
    help="Number of threads to use")
def generate_network(
    input_fasta: str,
    cytoscape_cyjs: str,
    mcl_format: str,
    mcl_edge_type: str,
    min_percent_identity: int,
    temp_dir: str,
    threads: int,
    output_plot: str = None,
) -> None:
    visualization.generate_network(
            fasta=input_fasta,
            cytoscape_network_path=cytoscape_cyjs,
            mcl_format_filepath=mcl_format,
            mcl_edge_type=mcl_edge_type,
            output_plot_path=output_plot,
            minimum_percent_identity=min_percent_identity,
            temp_dir=temp_dir,
            threads=threads,
    )


@cli.command('sequences')
@click.option(
    '--hmm-search-tab',
    type=File(mode='r', encoding=None, errors='strict', lazy=None, atomic=False),
    required=True,
    help="hmm search tab file from which a fasta will be generate from the hits")
@click.argument(
    "sequence_db_fasta",
    type=Path(exists=True, file_okay=True, dir_okay=True, readable=True, resolve_path=True,))
@click.argument(
    "output_fasta",
    type=Path(exists=False, file_okay=False, dir_okay=True, writable=True, resolve_path=True,))
def get_sequences(hmm_search_tab, sequence_db_fasta, output_fasta):
    hits = hmm_utils.iter_hits(SearchIO.parse(hmm_search_tab, 'hmmer3-tab'))
    seq_ids = [hit.target_protein for hit in hits]
    seq_records = utils.get_records_from_sequence_database(
        sequence_db_fasta=sequence_db_fasta,
        identifiers=seq_ids)
    utils.write_fasta(
        sequence_records=seq_records,
        output_fasta_path=output_fasta)


@cli.command('sample')
@click.option(
    '--input-fasta',
    type=Path(exists=True, file_okay=True, dir_okay=True, readable=True, resolve_path=True,),
    required=True,
    help="Input fasta file")
@click.option(
    '--percent-identity',
    type=float, required=True,
    help="Percent identity to run clustering at.")
@click.option(
    '--length-difference-cutoff',
    type=float, default=0.1,
    help="Sequences need to be at least this percent length of the representative sequence.")
@click.option(
    '--min-align-coverage',
    type=float, default=0.6,
    help="Alignment must cover at least this percent of both sequences.")
@click.option(
    '--percent-identity-suffix',
    type=bool, default=True,
    help="Include percent identity in the filenames of the cdhit output.")
@click.option(
    '--output', type=File(mode='w', encoding=None, errors='strict', lazy=None, atomic=False),
    required=True,
    help="File to write identifiers of sampled sequences"
    )
@click.option(
    '--output-dir',
    type=Path(exists=False, file_okay=True, dir_okay=True, readable=True, resolve_path=True,),
    required=False,
    help="Output dir to write cdhit files. When not provided will be written to the same dir as the"
         "input fasta")
def get_clusters(
        input_fasta: str,
        percent_identity: float,
        length_difference_cutoff: float,
        min_align_coverage: float,
        percent_identity_suffix: bool,
        output: str,
        output_dir: str,
) -> None:
    clusters = cluster.get_cdhit_clusters(
        input_fasta=input_fasta,
        percent_identity=percent_identity,
        length_difference_cutoff=length_difference_cutoff,
        min_alignment_coverage=min_align_coverage,
        percent_identity_suffix=percent_identity_suffix,
        output_dir=output_dir,
    )
    output.writelines([
        f'{c.representative}\n'
        for c in clusters
    ])


@cli.command('plot')
@click.option(
    '--input-fasta',
    type=Path(exists=True, file_okay=True, dir_okay=True, readable=True, resolve_path=True,),
    required=True,
    help="Input fasta file")
@click.option(
    '--output-png',
    type=Path(exists=False, file_okay=True, dir_okay=True, readable=True, resolve_path=True,),
    required=True,
    help="Output png of plot.")
@click.option(
    '--start-percent-identity',
    type=int, default=40,
    help="Beginning of range of percent identity to run clusterings at.")
@click.option(
    '--end-percent-identity',
    type=int, default=100,
    help="End of range of percent identity to run clusterings at.")
@click.option(
    '--step',
    type=int, default=5,
    help="Step to calculate percent identity at and mark on plot.")
@click.option(
    '--length-difference-cutoff',
    type=float, default=0.1,
    help="Sequences need to be at least this percent length of the representative sequence.")
@click.option(
    '--min-align-coverage',
    type=float, default=0.6,
    help="Alignment must cover at least this percent of both sequences.")
@click.option(
    '--output-dir',
    type=Path(exists=False, file_okay=True, dir_okay=True, readable=True, resolve_path=True,),
    required=False,
    help="Output dir to write cdhit files. When not provided will be written to the same dir as the"
         "input fasta")
def generate_cluster_number_plot(
        input_fasta: str,
        output_png: str,
        start_percent_identity: int,
        end_percent_identity: int,
        step: int,
        length_difference_cutoff: float,
        min_align_coverage: float,
        output_dir: str,
) -> None:
    cluster_count_tups = cluster.get_cdhit_cluster_sizes(
        input_fasta=input_fasta,
        start_percent_identity=start_percent_identity,
        end_percent_identity=end_percent_identity,
        step=step,
        length_difference_cutoff=length_difference_cutoff,
        min_alignment_coverage=min_align_coverage,
        output_dir=output_dir,
    )
    cluster.generate_cdhit_cluster_number_plot(
        cluster_counts=cluster_count_tups,
        output_png=output_png,
    )
