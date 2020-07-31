import subprocess


def make_database(
        database_name: str,
        fasta: str,
):
    """Creates a diamond formatted sequence database.

    Args:
        database_name: Full path to output diamond formatted sequence database
        fasta: Full path to fasta file of sequence to create database from

    Raises:
        CalledProcessError: If the subprocess running the Diamond program can not complete
            successfully.
    """
    subprocess.check_call(
        [
            'diamond', 'makedb',
            '--db', database_name,
            '--in', fasta,
        ]
    )


def blastp(
    database: str,
    output_tabfile: str,
    query_fasta: str,
    percent_identity: int,
    threads: int = None,
) -> None:
    """Runs protein alignments against a reference database using Diamond.

    Args:
        database: Full path to diamond formatted sequence database
        output_tabfile: Full path to file to write output tabfile
        query_fasta: Full path to fasta file of query sequence(s)
        threads: Number of threads to be used for blastp program

    Raises:
        CalledProcessError: If the subprocess running the Diamond program can not complete
            successfully.
    """
    params = [
        'diamond', 'blastp',
        '--db', database,
        '--out', output_tabfile,
        '--outfmt', '6',
        '--query', query_fasta,
        '--max-target-seqs', '0',
        '--max-hsps', '1',
        '--id', str(percent_identity),
        '--more-sensitive',
        '--no-self-hits',
    ]
    if threads is not None:
        params.extend(['--threads', str(threads)])
    subprocess.check_call(params)
