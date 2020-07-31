from collections import namedtuple
import os
from typing import (
    Generator,
    List,
)

from protein_helper.alignment_tools import (
    blastp,
    make_database,
)

Hit = namedtuple('Hit', ['query', 'target', 'percent_identity', 'evalue', 'bitscore'])


def run_blastp_all_by_all(
    fasta: str,
    percent_identity: int = 0,
    work_dir: str = None,
    threads: int = None,
) -> List[Hit]:
    """Runs a blastp all by all using Diamond.

    Args:
        fasta: Input fasta file
        percent_identity: Minimum percent identity for edge inclusion.
        work_dir: If provided diamond database and all outputfiles will be written to this dir
        threads: Number of threads to be used for blastp program

    Returns:
        A list of Hits

    Raises:
        CalledProcessError: If the subprocess running the Diamond program can not complete
            successfully.
    """
    fasta_prefix = os.path.basename(os.path.splitext(fasta)[0])

    if work_dir:
        database_path = os.path.join(work_dir, f'{fasta_prefix}.dmnd')
        diamond_out = os.path.join(work_dir, f'{fasta_prefix}.diamond_out.tab')
    else:
        fasta_root = os.path.splitext(fasta)[0]
        database_path = f'{fasta_root}.dmnd'
        diamond_out = f'{fasta_root}.diamond_out.tab'

    make_database(database_name=database_path, fasta=fasta)
    return(
        list(run_blastp(
            database=database_path,
            output_tabfile=diamond_out,
            query_fasta=fasta,
            percent_identity=percent_identity,
            threads=threads,
        )))


def hits(blastp_tabfile):
    """
    Yields Hits from a blastp tabfile

    Args:
        blastp_tabfile: Fielhandle for blastp tab results

    Returns:
        A Generator that yields a Hit
    """
    for row in blastp_tabfile:
        tokens = row.split()
        yield(Hit(query=tokens[0], target=tokens[1], percent_identity=float(tokens[2]),
                  evalue=tokens[10], bitscore=tokens[11]))


def run_blastp(
    database: str,
    output_tabfile: str,
    query_fasta: str,
    percent_identity: int,
    threads: int = None,
) -> Generator[Hit, any, None]:
    """Run Diamond blastp and parse results

    Args:
        database: Full path to diamond formatted sequence database
        output_tabfile: Full path to file to write output tabfile
        query_fasta: Full path to fasta file of query sequence(s)
        work_dir: If provided diamond database and all outputfiles will be written to this dir
        threads: Number of threads to be used for blastp program

    Returns:
        A Generator that yields a Hit

    Raises:
        CalledProcessError: If the subprocess running the Diamond program can not complete
            successfully.
    """
    blastp(
        database=database,
        output_tabfile=output_tabfile,
        query_fasta=query_fasta,
        percent_identity=percent_identity,
        threads=threads,
    )
    with open(output_tabfile) as tabfile:
        yield from hits(tabfile)
