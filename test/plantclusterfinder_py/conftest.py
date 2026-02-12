from pathlib import Path

import pytest


@pytest.fixture
def genelocation_obj():
    """
    Initialized GeneLocation object for testing
    """
    from plantclusterfinder_py.genelocation import GeneLoc

    return GeneLoc(
        file_path=Path(
            "csubellipsoidea/glof_CsubellipsoideaC_169_227_v2.0.gene.gff3.txt"
        )
    )


@pytest.fixture
def t_gene_loc():
    """
    Initialized GeneLocation object with simple test data for testing
    """
    from plantclusterfinder_py.genelocation import GeneLoc

    return GeneLoc(file_path=Path("test/data/gene_location_test.tsv"))


@pytest.fixture(scope="session")
def genelocation_transformed_obj():
    """
    Initialized and Transformed GeneLocation object for testing
    This shouldn't need further transformation,
    so it is only created once per session
    """
    from plantclusterfinder_py.cluster import (
        map_pgdb_gene_ids_to_conversion_file,
        map_pgdb_reactions_to_genes,
        read_in_conversion_file,
        read_in_pgdb_flat,
        read_MCL_clustering_file,
    )
    from plantclusterfinder_py.genelocation import GeneLoc

    gene_location = GeneLoc(
        file_path=Path(
            "csubellipsoidea/glof_CsubellipsoideaC_169_227_v2.0.gene.gff3.txt"
        )
    )
    conversion = read_in_conversion_file(
        Path(
            "./csubellipsoidea/gtpf_CsubellipsoideaC_169_227_v2.0.annotation_info.txt.txt"
        )
    )
    gene_location.remove_protein_transcript_info(conversion)
    gene_location.consolidate_gene_info()
    gene_location.sort_genes()
    pgdb_dir = Path("./csubellipsoidea/pgdb/csubellipsoideacyc/1.0/data")
    genes_all = read_in_pgdb_flat(
        pgdb_dir / "genes.dat",
        ("UNIQUE-ID", "ACCESSION-1"),
        unique_keys=("UNIQUE-ID", "ACCESSION-1"),
    )
    gene_info_from_genes_dat = ["ACCESSION-1"]
    ids_to_map = "TP"  # would be from args
    pgdb_conversions_df = map_pgdb_gene_ids_to_conversion_file(
        conversion, genes_all, gene_info_from_genes_dat, ids_to_map
    )
    gene2rxn = map_pgdb_reactions_to_genes(pgdb_dir)
    gene_location.annotate_enzyme_info(pgdb_conversions_df, gene2rxn)
    mcl_cluster_data = read_MCL_clustering_file(
        Path(
            "./csubellipsoidea/MCL_clustering_expected_results/dump.CsubellipsoideaC_169_227_v2.0.protein.pcf13.mci.I20"
        ),
        Path("./csubellipsoidea/CsubellipsoideaC_169_227_v2.0.protein.pcf13.fa"),
        conversion,
        False,  # do not remove proteins
    )
    gene_location.compute_local_duplication(
        mcl_cluster_data,
        100_000,  # kb for local duplication
        10,  # genes for local duplication
    )
    return gene_location
