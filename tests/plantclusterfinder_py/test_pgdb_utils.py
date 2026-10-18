import pytest

from tests.plantclusterfinder_py.assets import PGDB_DIR


def test_read_in_pgdb_flat():
    from plantclusterfinder_py.pgdb_utils import read_in_pgdb_flat

    genes_all = read_in_pgdb_flat(
        PGDB_DIR / "genes.dat",
        ("UNIQUE-ID", "ACCESSION-1", "PRODUCT", "COMMON-NAME"),
        unique_keys=("UNIQUE-ID",),
    )
    assert len(genes_all) == 2053

    pubs_df = read_in_pgdb_flat(
        PGDB_DIR / "pubs.dat",
        ("UNIQUE-ID", "AUTHORS", "PUBMED-ID", "TITLE", "SOURCE"),
        unique_keys=("UNIQUE-ID",),
    )
    assert len(pubs_df) == 2872


@pytest.mark.benchmark
def test_match_pgdb_flat():
    from plantclusterfinder_py.pgdb_utils import match_pgdb_flat, read_in_pgdb_flat

    genes_all = read_in_pgdb_flat(
        PGDB_DIR / "genes.dat",
        ("UNIQUE-ID", "ACCESSION-1", "PRODUCT", "COMMON-NAME"),
        unique_keys=("UNIQUE-ID",),
    )
    # read in the PGDB file expecting the "BOGUS" column (not found)
    proteins = read_in_pgdb_flat(
        PGDB_DIR / "proteins.dat",
        ("UNIQUE-ID", "ACCESSION-1", "GENE", "COMMON-NAME", "BOGUS"),
        unique_keys=("UNIQUE-ID",),
    )
    # gene UNIQUE-ID should join on protein GENE column
    genes_and_prots = match_pgdb_flat(genes_all, proteins, join_on=("UNIQUE-ID", "GENE"))
    assert len(genes_and_prots) == 2053
    # gene UNIQUE-ID should **not** join on protein "BOGUS" column
    genes_and_nulls = match_pgdb_flat(genes_all, proteins, join_on=("UNIQUE-ID", "BOGUS"))
    assert len(genes_and_nulls) == 0
    proteins_and_bogus = match_pgdb_flat(proteins, proteins, join_on=("BOGUS", "UNIQUE-ID"))
    assert len(proteins_and_bogus) == 0
    with pytest.raises(Exception, match="Column-matching Mismatch") as _:
        match_pgdb_flat(proteins, proteins, join_on=("NOT-COLUMN", "UNIQUE-ID"))
    with pytest.raises(Exception, match="Column-matching Mismatch") as _:
        match_pgdb_flat(proteins, proteins, join_on=("UNIQUE-ID", "NOT-COLUMN"))
