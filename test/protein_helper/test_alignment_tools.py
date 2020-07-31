from unittest import mock

from protein_helper.alignment_tools import blastp


class TestBlastp:

    def test_blastp(self):
        database = 'database.dmnd'
        output_tabfile = 'out.tab'
        query_fasta = 'query.fa'
        percent_identity = '30'
        with mock.patch('subprocess.check_call') as check_call_patch:
            blastp(
                database=database,
                output_tabfile=output_tabfile,
                query_fasta=query_fasta,
                percent_identity=percent_identity,
            )
        check_call_patch.assert_called_once_with([
            'diamond', 'blastp',
            '--db', database,
            '--out', output_tabfile,
            '--outfmt', '6',
            '--query', query_fasta,
            '--max-target-seqs', '0',  # Get all of the alignments.
            '--max-hsps', '1',
            '--id', str(percent_identity),
            '--more-sensitive',
            '--no-self-hits',
        ])

    def test_blastp_more_threads(self):
        database = 'database.dmnd'
        output_tabfile = 'out.tab'
        query_fasta = 'query.fa'
        percent_identity = '30'
        threads = 2
        with mock.patch('subprocess.check_call') as check_call_patch:
            blastp(
                database=database,
                output_tabfile=output_tabfile,
                query_fasta=query_fasta,
                percent_identity=percent_identity,
                threads=2,
            )
        check_call_patch.assert_called_once_with([
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
            '--threads', str(threads),
        ])
