import json

import matplotlib.pyplot as plt
import networkx

from protein_helper.align import run_blastp_all_by_all


def generate_network(
    fasta: str,
    cytoscape_network_path: str,
    mcl_format_filepath: str = None,
    mcl_edge_type: str = None,
    output_plot_path: str = None,
    minimum_percent_identity: int = 0,
    temp_dir: str = None,
    threads: int = None,
) -> None:
    """
    TODO: Finish docstring
    """
    edge_types = ['percent_identity', 'evalue', 'bitscore']
    if mcl_edge_type is not None and mcl_edge_type not in edge_types:
        raise ValueError(f'mcl_edge_type must be one of {", ".join(edge_types)}')  # TODO: add test.
    hits = run_blastp_all_by_all(
        fasta=fasta,
        work_dir=temp_dir,
        percent_identity=minimum_percent_identity,
        threads=threads,
    )
    g = networkx.Graph()
    g.add_edges_from([
        (
            hit.query,
            hit.target,
            {
                'percent_identity': hit.percent_identity,
                'evalue': hit.evalue,
                'bitscore': hit.bitscore
            }
        )
        for hit in hits
    ])

    if output_plot_path is not None:
        networkx.draw(g)
        plt.savefig(output_plot_path)

    if mcl_format_filepath is not None:
        edge_weight = networkx.get_edge_attributes(g, mcl_edge_type)
        print(edge_weight)
        with open(mcl_format_filepath, 'w') as out_mcl:
            out_mcl.writelines([
                f"{e[0]}\t{e[1]}\t{edge_weight[e]}\n"
                for e in g.edges
            ])

    cyjs_json = networkx.readwrite.json_graph.cytoscape_data(g)

    with open(cytoscape_network_path, 'w') as output_network:
        json.dump(cyjs_json, output_network)
