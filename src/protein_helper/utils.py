from typing import (
    Generator,
    List,
)

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def get_records_from_sequence_database(
        sequence_db_fasta: str,
        identifiers: list,
) -> Generator[SeqRecord, any, None]:
    record_dict = SeqIO.index(sequence_db_fasta, "fasta")
    for id_ in identifiers:
        yield record_dict[id_]


def write_fasta(
        sequence_records: List[SeqRecord],
        output_fasta_path: str,
) -> None:
    with open(output_fasta_path, "w") as handle:
        for r in sequence_records:
            handle.write(r.format("fasta"))
        # SeqIO.write(sequence_records, handle, "fasta")
