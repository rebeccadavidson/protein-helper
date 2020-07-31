from importlib.resources import path
import os

import pytest

from protein_helper.align import (
    Hit,
    hits,
    run_blastp,
    run_blastp_all_by_all,
)
from test.fixtures import (
    blastp,
    blastp_all_by_all,
    parse_blastp,
)


@pytest.fixture
def expected_hits():
    return [
        Hit(query='seq1', target='seq2', percent_identity=95.0, evalue='1.4e-28', bitscore='105.9'),
        Hit(query='seq2', target='seq1', percent_identity=95.0, evalue='2.4e-28', bitscore='105.1'),
    ]


@pytest.fixture
def expected_hits_real_data():
    return [
        Hit(query='EST3A_MOUSE', target='EST1_PIG', percent_identity=42.6, evalue='3.6e-120',
            bitscore='416.4'),
        Hit(query='EST3A_MOUSE', target='H0VHN0_CAVPO', percent_identity=41.6, evalue='4.6e-115',
            bitscore='399.4'),
        Hit(query='H0VHN0_CAVPO', target='EST1_PIG', percent_identity=72.7, evalue='3.1e-236',
            bitscore='802.0'),
        Hit(query='H0VHN0_CAVPO', target='EST3A_MOUSE', percent_identity=41.4, evalue='1.7e-117',
            bitscore='407.5'),
        Hit(query='EST1_PIG', target='H0VHN0_CAVPO', percent_identity=72.7, evalue='3.8e-234',
            bitscore='795.0'),
        Hit(query='EST1_PIG', target='EST3A_MOUSE', percent_identity=42.6, evalue='1.3e-120',
            bitscore='417.9'),
    ]


def test_run_all_by_all(tmp_path, expected_hits):
    with path(blastp_all_by_all, "input.fasta") as fasta:
        hits = run_blastp_all_by_all(fasta=fasta, work_dir=tmp_path)
    assert expected_hits == hits


def test_run_all_by_all_real_data(tmp_path, expected_hits_real_data):
    with path(blastp_all_by_all, "PF00135_seed.fasta") as fasta:
        hits = list(run_blastp_all_by_all(fasta=fasta, work_dir=tmp_path))
    assert expected_hits_real_data == hits


def test_blastp(tmp_path, expected_hits):
    with path(blastp, "input.dmnd") as db, path(blastp, "input.fasta") as query_fasta:
        fasta_prefix = os.path.basename(os.path.splitext(query_fasta)[0])
        diamond_out = os.path.join(tmp_path, f'{fasta_prefix}.diamond_out.tab')
        hits = list(run_blastp(
            database=db,
            output_tabfile=diamond_out,
            query_fasta=query_fasta,
            percent_identity=0,
        ))
    assert expected_hits == hits


def test_blastp_real_data(tmp_path, expected_hits_real_data):
    with path(blastp, "PF00135_seed.dmnd") as db, path(blastp, "PF00135_seed.fasta") as query_fasta:
        fasta_prefix = os.path.basename(os.path.splitext(query_fasta)[0])
        diamond_out = os.path.join(tmp_path, f'{fasta_prefix}.diamond_out.tab')
        hits = list(run_blastp(
            database=db,
            output_tabfile=diamond_out,
            query_fasta=query_fasta,
            percent_identity=0,
        ))
    assert expected_hits_real_data == hits


def test_parse_blast(expected_hits):
    with path(parse_blastp, "input.diamond_out.tab") as tabfile, open(tabfile) as diamond_tab:
        hits_ = list(hits(diamond_tab))
    assert expected_hits == hits_


def test_parse_blast_real(expected_hits_real_data):
    with path(parse_blastp, "PF00135_seed.diamond_out.tab") as tabfile, \
            open(tabfile) as diamond_tab:
        hits_ = list(hits(diamond_tab))
    assert expected_hits_real_data == hits_
