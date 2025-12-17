import random
from pathlib import Path

import polars as pl
import pytest


def test_init_genelocation():
    from plantclusterfinder_py.genelocation import GeneLoc

    gene_location = GeneLoc(
        file_path=Path(
            "csubellipsoidea/glof_CsubellipsoideaC_169_227_v2.0.gene.gff3.txt"
        )
    )

    assert 5 == len(gene_location.data.columns)
    assert 9_629 == len(gene_location)


@pytest.fixture
def genelocation_obj():
    from plantclusterfinder_py.genelocation import GeneLoc

    return GeneLoc(
        file_path=Path(
            "csubellipsoidea/glof_CsubellipsoideaC_169_227_v2.0.gene.gff3.txt"
        )
    )


@pytest.fixture
def t_gene_loc():
    from plantclusterfinder_py.genelocation import GeneLoc

    return GeneLoc(file_path=Path("test/data/gene_location_test.tsv"))


def test_remove_protein_transcript_info_from_genepositionfile(
    genelocation_obj, t_gene_loc
):
    from plantclusterfinder_py.cluster import read_in_conversion_file

    gene_location = genelocation_obj
    conversion = read_in_conversion_file(
        "./csubellipsoidea/gtpf_CsubellipsoideaC_169_227_v2.0.annotation_info.txt.txt"
    )
    expected_columns = [
        "Gene Name",
        "Gene Start (bp)",
        "Gene End (bp)",
        "Chromosome Name",
        "Strand",
        "Transcript Name",  # added
        "Protein Name",  # added
    ]
    gene_location.remove_protein_transcript_info(conversion)
    assert 7 == len(gene_location.columns)
    assert expected_columns == gene_location.columns
    assert 9_629 == len(gene_location)

    test_gene_location = t_gene_loc
    bad_conversion = read_in_conversion_file("./test/data/bad_conversion.tsv")
    with pytest.raises(Exception) as exception_info:
        test_gene_location.remove_protein_transcript_info(bad_conversion)
        assert "more than one" in str(exception_info.value)


def test_consolidate_gene_info(genelocation_obj, t_gene_loc):
    gene_location = genelocation_obj
    gene_info = gene_location.consolidate_gene_info()
    assert 45 == len(gene_info)
    test_gene_location = t_gene_loc
    test_gene_info = test_gene_location.consolidate_gene_info()
    assert ["A"] == test_gene_info


def test_sort_genes(genelocation_obj, t_gene_loc):
    import polars as pl

    gene_location = genelocation_obj
    gene_location.sort_genes()

    assert 9629 == len(gene_location)

    test_gene_location = t_gene_loc
    test_gene_location.sort_genes()
    sorted_genes = test_gene_location.data.select(pl.col("Gene Name"))
    expected_gene_order = pl.DataFrame(
        {"Gene Name": ["GENE1", "GENE2", "GENE3", "GENE4"]}
    )
    unexpected_gene_order = pl.DataFrame(
        {"Gene Name": ["GENE2", "GENE1", "GENE4", "GENE3"]}
    )
    assert sorted_genes.equals(expected_gene_order)
    assert not sorted_genes.equals(unexpected_gene_order)


def test_map_genelocation_to_conversion_file(genelocation_obj):
    # GeneLocation Test, possibly redundant (mapping takes place in remove_protein_transcript_info)
    from plantclusterfinder_py.cluster import read_in_conversion_file

    conversion = read_in_conversion_file(
        Path(
            "./csubellipsoidea/gtpf_CsubellipsoideaC_169_227_v2.0.annotation_info.txt.txt"
        )
    )
    expected_columns = [
        "Gene Name",
        "Gene Start (bp)",
        "Gene End (bp)",
        "Chromosome Name",
        "Strand",
        "Transcript Name",  # added
        "Protein Name",  # added
    ]
    gene_location = genelocation_obj
    genelocation_obj.remove_protein_transcript_info(conversion)
    BiomG_G = gene_location.map_to_conversion_file(conversion.MapG)
    assert expected_columns == gene_location.data.columns
    assert 9629 == len(BiomG_G)


def test_remove_non_protein_fasta_gene_locations(genelocation_obj):
    from plantclusterfinder_py.cluster import read_in_conversion_file

    conversion = read_in_conversion_file(
        Path(
            "./csubellipsoidea/gtpf_CsubellipsoideaC_169_227_v2.0.annotation_info.txt.txt"
        )
    )
    gene_location = genelocation_obj
    gene_location.remove_protein_transcript_info(conversion)
    BiomG_G = gene_location.map_to_conversion_file(conversion.MapG)
    ok_prot = list(BiomG_G.keys())
    random.shuffle(ok_prot)
    protein_ids = ok_prot[0:100]
    cleaned_gene_location = gene_location.remove_non_protein_fasta_gene_locations(
        conversion, BiomG_G, protein_ids
    )
    unclean_gene_location = gene_location.remove_non_protein_fasta_gene_locations(
        conversion, BiomG_G, []
    )
    assert len(gene_location) == len(unclean_gene_location)
    assert 100 == len(cleaned_gene_location)  # TEST INCOMPLETE


# Part 9
def test_find_intergenic_regions(genelocation_obj, t_gene_loc):
    import polars as pl

    # GeneLocation test
    gene_location = genelocation_obj
    intergenic_regions = gene_location.find_intergenic_regions()
    #    gene_location.data.select_seq(3).unique().to_list()
    # )
    assert 9413 == len(intergenic_regions)  # TEST INCOMPLETE

    intergenic_regions_2 = t_gene_loc.find_intergenic_regions()
    t_gene_result_expected = pl.DataFrame(
        {
            "Chromosome Name": ["A", "A"],
            "Start": [101, 401],
            "End": [149, 499],
            "Gene Before": ["GENE1", "GENE3"],
            "Gene After": ["GENE2", "GENE4"],
        }
    )
    assert intergenic_regions_2.equals(t_gene_result_expected)


# Part 14
def test_annotate_enzyme_info(genelocation_obj):
    pgdb_dir = Path("./csubellipsoidea/pgdb/csubellipsoideacyc/1.0/data")
    from plantclusterfinder_py.cluster import (
        map_pgdb_reactions_to_genes,
        map_PGDBs_to_conversion_file,
        read_in_conversion_file,
        read_in_pgdb_flat,
    )

    conversion = read_in_conversion_file(
        Path(
            "./csubellipsoidea/gtpf_CsubellipsoideaC_169_227_v2.0.annotation_info.txt.txt"
        )
    )
    genes_all = read_in_pgdb_flat(
        pgdb_dir / "genes.dat",
        ("UNIQUE-ID", "ACCESSION-1"),
        unique_keys=("UNIQUE-ID", "ACCESSION-1"),
    )
    gene_info_from_genes_dat = ["ACCESSION-1"]  # TODO
    ids_to_map = "T"
    pgdb_conversions_dict = map_PGDBs_to_conversion_file(
        conversion, genes_all, gene_info_from_genes_dat, ids_to_map
    )
    gene_location = genelocation_obj
    gene_location.remove_protein_transcript_info(conversion)
    gene2rxn = map_pgdb_reactions_to_genes(pgdb_dir)
    genes_with_reactions = gene_location.annotate_enzyme_info(
        pgdb_conversions_dict, gene2rxn
    )
    expected_columns = [
        "Gene Name",
        "Gene Start (bp)",
        "Gene End (bp)",
        "Chromosome Name",
        "Strand",
        "Transcript Name",
        "Protein Name",
        "pgdb_gene_id",  # added
        "has_reaction",  # added
        "pgdb_reaction_id",  # added
        "EC-NUMBER",  # added
    ]
    assert expected_columns == gene_location.data.columns
    assert 9629 == len(genes_with_reactions)
    assert 2053 == len(list(filter(lambda x: x, genes_with_reactions.values())))


def test_introduce_hard_breaks(genelocation_obj):
    from plantclusterfinder_py.cluster import annotate_sequencing_gaps

    sequencing_gaps_file = Path(
        "./csubellipsoidea/Gapanalysis_expected_results/CsubellipsoideaC_169_227_v2.0.hardmasked.fa"
    )
    gene_location = genelocation_obj
    uniq_chr = list(
        gene_location.data.select("Chromosome Name").unique().to_dict().values()
    )[0].to_list()
    sequencing_gaps = annotate_sequencing_gaps(
        sequencing_gaps_file, gene_location.find_intergenic_regions(), uniq_chr
    )
    gene_location.introduce_hard_breaks(
        sequencing_gaps, vMaxSeqGapSize=int(4e4), vPercentileForMaxInterGeneDist=80
    )
    assert 45 == gene_location.data.select(pl.col("Chromosome Name").unique()).height


def test_insert_genes(genelocation_obj):
    print(genelocation_obj.data.head(10))
    test_genes = [
        (f"TEST{i + 1:02}", 50 + 1000 * i, 50 + 1000 * i, "test_scaffold", "-1")
        for i in range(500)
    ]
    gene_location = genelocation_obj.model_copy()
    orig_gene_count = len(gene_location)
    test_gene_df = pl.DataFrame(
        test_genes, schema=gene_location.data.schema, orient="row"
    )
    print(test_gene_df.head(10))
    gene_location.insert_genes(test_gene_df)
    assert orig_gene_count + 500 == len(gene_location)


def test_compute_local_duplication(genelocation_obj):
    from plantclusterfinder_py.cluster import (
        read_in_conversion_file,
        read_MCL_clustering_file,
    )

    conversion = read_in_conversion_file(
        "./csubellipsoidea/gtpf_CsubellipsoideaC_169_227_v2.0.annotation_info.txt.txt"
    )

    mcl_cluster_data = read_MCL_clustering_file(
        "./csubellipsoidea/MCL_clustering_expected_results/dump.CsubellipsoideaC_169_227_v2.0.protein.pcf13.mci.I20",
        "./csubellipsoidea/CsubellipsoideaC_169_227_v2.0.protein.pcf13.fa",
        conversion,
        False,  # do not remove proteins
    )
    gene_location = genelocation_obj
    gene_location.remove_protein_transcript_info(conversion)
    local_duplication = gene_location.compute_local_duplication(
        mcl_cluster_data, {}, 10000, 5
    )
    expected_columns = [
        "Gene Name",
        "Gene Start (bp)",
        "Gene End (bp)",
        "Chromosome Name",
        "Strand",
        "Transcript Name",
        "Protein Name",
        "neighbor_genes",
        "mcl_genes",
        "mcl_clusters",
        "is_locally_duplicated",
    ]
    print(local_duplication.columns)
    assert expected_columns == local_duplication.columns
    assert 848 == local_duplication.filter(pl.col("is_locally_duplicated")).height
