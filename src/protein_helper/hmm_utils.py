from collections import namedtuple
from typing import List

from Bio.SearchIO import QueryResult


HmmHit = namedtuple('HmmHit', ['query_hmm', 'target_protein', 'evalue', 'bitscore'])


def iter_hits(
        results: List[QueryResult]
) -> None:
    """
    Returns:
        A Generator that yields the highest bitscore hit for each QueryResult
    """
    for r in results:
        for hit in r.hits:
            yield(HmmHit(query_hmm=r.id, target_protein=hit.id, evalue=hit.evalue,
                         bitscore=hit.bitscore))
