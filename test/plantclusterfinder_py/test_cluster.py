import logging
from pathlib import Path
from tempfile import NamedTemporaryFile

import polars as pl
import pytest

# globals
pgdb_dir = Path("./csubellipsoidea/pgdb/csubellipsoideacyc/1.0/data")
test_input_dir = Path("./Inputs")
_log = logging.getLogger(__name__)

TOP_PERCENT = 5
GENE_COUNT = 9_629
PGDB_GENE_COUNT = 2053
# Part 1


def test_read_in_conversion_file():
    """The conversion file defines the genes, proteins, and reactions
    that many downstream processes rely on."""

    from plantclusterfinder_py.cluster import read_in_conversion_file

    conversion = read_in_conversion_file(
        Path(
            "./csubellipsoidea/gtpf_CsubellipsoideaC_169_227_v2.0.annotation_info.txt.txt"
        )
    )
    assert len(conversion.MapG) == GENE_COUNT
    assert len(conversion.MapT) == GENE_COUNT
    assert len(conversion.MapP) == GENE_COUNT


# Part 2
def test_read_in_pgdb_flat():
    """Used by map_pgdb_reactions_to_genes"""

    from plantclusterfinder_py.cluster import read_in_pgdb_flat

    gene_df = read_in_pgdb_flat(
        pgdb_dir / "genes.dat", ("UNIQUE-ID", "PRODUCT", "ACCESSION-1")
    )
    assert isinstance(gene_df, pl.DataFrame)
    assert len(gene_df) == PGDB_GENE_COUNT
    g4s_3028 = gene_df.filter(pl.col("PRODUCT") == ["G4S-3028-MONOMER"])
    assert len(g4s_3028) == 1


def test_match_pgdb_flat():
    """Used by map_pgdb_reactions_to_genes"""
    from plantclusterfinder_py.cluster import match_pgdb_flat, read_in_pgdb_flat

    flat1_t = ["A", "B", "C", "D"]
    flat2_t = ["A", "D", "B", "C"]
    flat1 = pl.DataFrame([{"UNIQUE-ID": i, "P": j} for i, j in enumerate(flat1_t)])
    flat2 = pl.DataFrame([{"UNIQUE-ID": j, "Q": i} for i, j in enumerate(flat2_t)])
    matched = match_pgdb_flat(flat1, flat2, ("P", "UNIQUE-ID"))
    assert 4 == len(matched)

    genes_flat = integration_flat1 = read_in_pgdb_flat(
        pgdb_dir / "genes.dat", ("UNIQUE-ID", "PRODUCT", "ACCESSION-1", "ACCESSION-2")
    )
    proteins_flat = integration_flat2 = read_in_pgdb_flat(
        pgdb_dir / "proteins.dat", ("UNIQUE-ID", "CATALYZES", "COMPONENTS")
    )
    reaction_flat = read_in_pgdb_flat(
        pgdb_dir / "reactions.dat", ("UNIQUE-ID", "ENZYMATIC-REACTION")
    )
    enzrxns_all = read_in_pgdb_flat(
        pgdb_dir / "enzrxns.dat", ("UNIQUE-ID", "REACTIONS", "ENZYME")
    )
    integration_matched = match_pgdb_flat(
        integration_flat1, integration_flat2, ("PRODUCT", "UNIQUE-ID")
    )
    assert PGDB_GENE_COUNT == len(integration_matched)  # TEST INCOMPLETE
    assert PGDB_GENE_COUNT == len(genes_flat)
    assert 2_054 == len(proteins_flat)
    assert 3_015 == len(enzrxns_all)
    assert len(reaction_flat) > 0


def test_map_pgdb_reactions_to_genes():
    from plantclusterfinder_py.cluster import map_pgdb_reactions_to_genes

    mapping_G2R = map_pgdb_reactions_to_genes(pgdb_dir)
    assert 2195 == len(mapping_G2R.ReactionsAll)
    assert PGDB_GENE_COUNT == len(mapping_G2R.GenesAll)
    assert PGDB_GENE_COUNT == len(mapping_G2R.Gene2Rxn)


# Part 3
def test_get_metabolic_domains():
    from plantclusterfinder_py.cluster import get_metabolic_domains

    metabolic_domains = get_metabolic_domains(
        test_input_dir / "ReactionMetabolicDomainClassification.txt",
        "Specialized Metabolism",
    )
    assert 14_409 == len(metabolic_domains.RxnMetaDom)
    assert 9_682 == len(metabolic_domains.Rxn2MetaDom)


def test_matchPGDBs_to_MDs():
    from plantclusterfinder_py.cluster import (
        get_metabolic_domains,
        map_pgdb_reactions_to_genes,
        match_PGDBs_to_MDs,
    )

    metabolic_domains = get_metabolic_domains(
        test_input_dir / "ReactionMetabolicDomainClassification.txt",
        "Specialized Metabolism",
    )
    mapping_G2R = map_pgdb_reactions_to_genes(pgdb_dir)
    gene2domain = match_PGDBs_to_MDs(mapping_G2R, metabolic_domains)
    assert 2_145 == len(gene2domain)


# part 4
def test_label_reactions_with_sig_and_tail():
    from plantclusterfinder_py.cluster import (
        label_reactions_with_sig_and_tail,
        read_in_pgdb_flat,
    )

    test_signature_enzymes_tsv = (
        test_input_dir / "scaffold-tailoring-reactions-05082016.tab"
    )
    reactions_dat = read_in_pgdb_flat(
        pgdb_dir / "reactions.dat", ("UNIQUE-ID", "ENZYMATIC-REACTION", "EC-NUMBER")
    )

    labeled_reaction = label_reactions_with_sig_and_tail(
        reactions_dat, test_signature_enzymes_tsv
    )
    assert 9 == len(labeled_reaction.STClasses)
    assert 3 == len(labeled_reaction.SigOrTail)
    assert 394 == len(labeled_reaction.Rxn_STClasses)
    assert 394 == len(labeled_reaction.Rxn_SigOrTail)


# part 5
def test_map_PGDBs_to_conversion_file():
    from plantclusterfinder_py.cluster import (
        map_pgdb_gene_ids_to_conversion_file,
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
    ids_to_map = "GTP"
    pgdb_conversions_df = map_pgdb_gene_ids_to_conversion_file(
        conversion, genes_all, gene_info_from_genes_dat, ids_to_map
    )
    assert PGDB_GENE_COUNT == len(pgdb_conversions_df)
    g4s_3028_row = pgdb_conversions_df.filter(pl.col("pgdb_gene_id") == "G4S-3028")
    assert "46517" == g4s_3028_row["Protein Name"][0]
    assert "46517" == g4s_3028_row["Transcript Name"][0]


# part 6
def test_read_MCL_clustering_File():
    from plantclusterfinder_py.cluster import (
        read_in_conversion_file,
        read_MCL_clustering_file,
    )

    conversion = read_in_conversion_file(
        Path(
            "./csubellipsoidea/gtpf_CsubellipsoideaC_169_227_v2.0.annotation_info.txt.txt"
        )
    )
    mcl_cluster_data = read_MCL_clustering_file(
        Path(
            "./csubellipsoidea/MCL_clustering_expected_results/dump.CsubellipsoideaC_169_227_v2.0.protein.pcf13.mci.I20"
        ),
        Path("./csubellipsoidea/CsubellipsoideaC_169_227_v2.0.protein.pcf13.fa"),
        conversion,
        False,  # do not remove proteins
    )
    assert GENE_COUNT - 2 == len(mcl_cluster_data.Gene2MCLClusters)
    assert isinstance(mcl_cluster_data.Protein_IDs, list)
    assert GENE_COUNT - 2 == len(mcl_cluster_data.Protein_IDs)


def test_get_MCL_clustering():
    from plantclusterfinder_py.cluster import (
        get_MCL_clustering,
        read_in_conversion_file,
        read_MCL_clustering_file,
    )

    conversion = read_in_conversion_file(
        Path(
            "./csubellipsoidea/gtpf_CsubellipsoideaC_169_227_v2.0.annotation_info.txt.txt"
        )
    )
    mcl_cluster_diamond = get_MCL_clustering(
        Path("./csubellipsoidea/CsubellipsoideaC_169_227_v2.0.protein.pcf13.fa"),
        conversion,
        True,  # overvrite MCL cluster file
        False,  # remove protein locations
        rebuild_blast_db=True,
        rerun_blastp=True,
        use_diamond=True,
    )
    mcl_cluster_blast = get_MCL_clustering(
        Path("./csubellipsoidea/CsubellipsoideaC_169_227_v2.0.protein.pcf13.fa"),
        conversion,
        True,  # overvrite MCL cluster file
        False,  # remove protein locations
        rebuild_blast_db=True,
        rerun_blastp=True,
        use_diamond=False,
    )
    expected_mcl_cluster_data = read_MCL_clustering_file(
        Path(
            "./csubellipsoidea/MCL_clustering_expected_results/dump.CsubellipsoideaC_169_227_v2.0.protein.pcf13.mci.I20"
        ),
        Path("./csubellipsoidea/CsubellipsoideaC_169_227_v2.0.protein.pcf13.fa"),
        conversion,
        False,  # do not remove proteins
    )
    assert GENE_COUNT - 1 == len(mcl_cluster_diamond.Gene2MCLClusters)
    assert (
        len(expected_mcl_cluster_data.Gene2MCLClusters)
        == len(mcl_cluster_diamond.Gene2MCLClusters) - 1
    )
    assert GENE_COUNT - 2 == len(mcl_cluster_blast.Gene2MCLClusters)
    assert len(expected_mcl_cluster_data.Gene2MCLClusters) == len(
        mcl_cluster_blast.Gene2MCLClusters
    )


# Part 7 through Part 10: No Cluster Functions


# part 11
def test_analyze_plant_cluster_gap_file():
    """Used by get_sequencing_gaps"""
    from plantclusterfinder_py.cluster import analyze_plant_cluster_gap_file

    test_gap_output = Path(
        "./csubellipsoidea/Gapanalysis_expected_results/CsubellipsoideaC_169_227_v2.0.hardmasked.fa_GAPOutput"
    )
    plant_cluster_analysis = analyze_plant_cluster_gap_file(
        test_gap_output,
        False,
        file_output=False,
    )
    assert plant_cluster_analysis is not None
    (filename, gap_output_count) = plant_cluster_analysis
    assert 337 == len(gap_output_count)  # TEST INCOMPLETE
    assert str(test_gap_output) + "_count.txt" == filename  # TEST INCOMPLETE


def test_get_sequencing_gaps():
    """Used by annotate_sequencing_gaps"""
    # this is better tested in get_positions_of_gap
    assert 1 == 1  # TEST INCOMPLETE


def test_annotate_sequencing_gaps(genelocation_obj):
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
    assert 1166 == sequencing_gaps.height  # TEST INCOMPLETE
    assert 6 == len(sequencing_gaps.columns)


# part 12
# hard breaks moved to gene location test


# part 13
def test_introduce_hypo_genes(genelocation_obj):
    from plantclusterfinder_py.cluster import introduce_hypo_genes

    gene_location = genelocation_obj.model_copy()

    orig_gene_count = len(gene_location)
    introduce_hypo_genes(
        gene_location,
        gene_location.find_intergenic_regions(),
        hypothetical_gene_amount=1,
    )
    non_hypo_count = gene_location.data.filter(
        pl.col("Gene Name").str.starts_with("Hypo").not_()
    ).height
    assert orig_gene_count * 2 - 216 == len(gene_location)  # TEST INCOMPLETE
    assert orig_gene_count == non_hypo_count
    assert orig_gene_count != 0


# Part 14
# in gene location test


# Part 15
# in gene location test


# Part 16
def test_compute_background():
    """Used by perform_clustering"""
    from plantclusterfinder_py.cluster import compute_background
    from plantclusterfinder_py.genelocation import GeneLoc

    gene_location = GeneLoc(
        file_path=Path("./test/data/background_computation_test.tsv")
    )
    has_reaction = pl.read_csv(
        Path("./test/data/background_computation_has_reaction.csv")
    )
    gene_location._data = gene_location.data.join(has_reaction, on="Gene Name")

    step_sizes = [1 + x for x in range(4) for _ in range(100)]
    enzyme_counts = [3 for _ in range(400)]
    window_sizes = [3 + x for _ in range(100) for x in range(4)]

    cluster = pl.DataFrame(
        (
            {
                "step_size": x,
                "enzyme_count": y,
                "non_enzyme_gene_count": z - y,
                "test_window_size": z,
            }
            for x, y, z in zip(step_sizes, enzyme_counts, window_sizes)
        )
    )
    cluster.write_csv("test_cluster.csv")
    expected_background = pl.DataFrame(
        [
            {"window_size": x, "enzyme_count": y, "tested": True}
            for x in range(1, 7)
            for y in range(x + 1)
        ]
    ).filter(  # filtering out background cases not seen in data
        (pl.col("window_size") <= 4) | (pl.col("window_size") != pl.col("enzyme_count"))
    )
    background = compute_background(cluster, gene_location, pl.DataFrame())
    assert expected_background.height == background.height
    for window in (
        background.group_by("window_size")
        .agg(
            pl.sum("winzyme_occurrence"),
            pl.sum("winzyme_frequency"),
            pl.last("window_occurrence"),
        )
        .iter_rows()
    ):
        (size, occurrence, frequency, total_occurrence) = window
        assert occurrence == total_occurrence
        assert occurrence == 186 - (size - 1) * 6
    return


def test_perform_clustering(genelocation_transformed_obj):
    from plantclusterfinder_py.cluster import (
        perform_clustering,
        read_in_conversion_file,
        read_MCL_clustering_file,
    )

    conversion = read_in_conversion_file(
        Path(
            "./csubellipsoidea/gtpf_CsubellipsoideaC_169_227_v2.0.annotation_info.txt.txt"
        )
    )
    gene_location = genelocation_transformed_obj
    mcl_cluster_data = read_MCL_clustering_file(
        Path(
            "./csubellipsoidea/MCL_clustering_expected_results/dump.CsubellipsoideaC_169_227_v2.0.protein.pcf13.mci.I20"
        ),
        Path("./csubellipsoidea/CsubellipsoideaC_169_227_v2.0.protein.pcf13.fa"),
        conversion,
        False,  # do not remove proteins
    )

    _log.debug(f"Gene Location before clustering: {gene_location}")
    clusters = perform_clustering(
        gene_location,
        mcl_cluster_data,
        3,  # enzymes min for clusters
        2,  # criterion
        5,  # min step size
        20,  # max step size
    )
    clusters.perform_cutoff(TOP_PERCENT)
    print(clusters)
    assert 3 == clusters.FinalStepSize
    assert 5 == clusters.MaxStepSize
    # assert 86 == clusters.ClusterMeta.height # This was the old test value
    # assert 482 == clusters.ClusterBackground.height
    assert 1192 == clusters.ClusterMeta.height
    assert 830 == clusters.ClusterBackground.height


def test_evaluate_criterion():
    from collections import namedtuple

    from plantclusterfinder_py.cluster import evaluate_criterion

    # named tuple to simplify creation of clusterdf
    cluster_tuple = namedtuple(
        "cluster_tuple",
        [
            "chromosome",
            "step_size",
            "start",
            "end",
            "enzyme_count",
            "non_enzyme_gene_count",  # includes hypotheticals
            "hypothetical_gene_count",
        ],
    )

    # Evaluate 1:
    #   calculated as enzymes / (non enzymes) > 0.5
    #   AKA: keep searching as long as there are more enzymes than 2x non-enzymes in clusters.
    #   Continue search case: More enzymes than non-enzymes.
    case_one_continue_clusters = pl.DataFrame(
        [
            cluster_tuple(
                "TEST",
                3,
                x * 50,
                (x + 1) * 50 + 12,
                6,  # 6 enzymes per cluster
                12,  # 12 non-enzyme
                0,  # 0 hypothetical gene
            )
            for x in range(10)
        ]
    )
    case_one_continue = evaluate_criterion(1, case_one_continue_clusters, 3, 24, 3)
    #   Stop search case: Enzymes are fewer than or equal to the number of non-enzymes.
    case_one_stop_clusters = pl.DataFrame(
        [
            cluster_tuple(
                "TEST",
                3,
                x * 50,
                (x + 1) * 50 + 12,
                6,  # 6 enzymes per cluster
                13,  # 13 non-enzymes
                0,  # of which 0 are hypothetical
            )
            for x in range(10)
        ]
    )
    case_one_stop = evaluate_criterion(1, case_one_stop_clusters, 3, 24, 3)
    assert (2, True) == case_one_continue
    assert (2, False) == case_one_stop
    # Evaluate 2: continue search until step_size is 3
    empty_clusters = pl.DataFrame([cluster_tuple("TEST", 0, 0, 0, 0, 0, 0)])
    #    Continue search: step size less than 3
    case_two_continue = evaluate_criterion(2, empty_clusters, 2, 24, 3)
    #    Stop search: step size 3 or greater
    case_two_stop = evaluate_criterion(2, empty_clusters, 3, 24, 3)
    assert (1, True) == case_two_continue
    assert (3, False) == case_two_stop
    # Evaluate 3:
    #    Continue Search: enzymes > non-enzymes + 2*hypos
    case_three_continue_clusters = pl.DataFrame(
        [
            cluster_tuple(
                "TEST",
                3,
                x * 50,
                (x + 1) * 50 + 12,
                12,  # 12 enzymes per cluster
                9,  # 9 non-enzymes
                3,  # of which 3 are hypothetical
            )
            for x in range(10)
        ]
    )
    case_three_continue = evaluate_criterion(3, case_three_continue_clusters, 3, 24, 3)
    #    Stop Search: enzymes <= non-enzymes + 2*hypos
    case_three_stop_clusters = pl.DataFrame(
        [
            cluster_tuple(
                "TEST",
                3,
                x * 50,
                (x + 1) * 50 + 12,
                12,  # 12 enzymes per cluster
                10,  # 10 non-enzyme
                3,  # of which 3 are hypothetical
            )
            for x in range(10)
        ]
    )
    case_three_stop = evaluate_criterion(3, case_three_stop_clusters, 3, 24, 3)
    assert (2, True) == case_three_continue
    assert (2, False) == case_three_stop
    # Evaluate 4:
    #    Continue case: there are more clusters for step size 3 than 2
    case_four_continue_clusters = pl.DataFrame(
        [
            cluster_tuple(
                "TEST",
                3,
                x * 50,
                (x + 1) * 50 + 12,
                12,  # 12 enzymes per cluster
                6,  # 1 non-enzyme
                3,  # 1 hypothetical gene
            )
            for x in range(10)
        ]
        + [
            cluster_tuple(
                "TEST",
                2,
                x * 50,
                (x + 1) * 50 + 12,
                12,  # 12 enzymes per cluster
                6,  # 1 non-enzyme
                3,  # 1 hypothetical gene
            )
            for x in range(9)
        ]
    )
    case_four_continue = evaluate_criterion(4, case_four_continue_clusters, 3, 24, 3)
    case_four_stop_clusters = pl.DataFrame(
        [
            cluster_tuple(
                "TEST",
                3,
                x * 50,
                (x + 1) * 50 + 12,
                12,  # 12 enzymes per cluster
                6,  # 1 non-enzyme
                3,  # 1 hypothetical gene
            )
            for x in range(9)
        ]
        + [
            cluster_tuple(
                "TEST",
                2,
                x * 50,
                (x + 1) * 50 + 12,
                12,  # 12 enzymes per cluster
                6,  # 1 non-enzyme
                3,  # 1 hypothetical gene
            )
            for x in range(10)
        ]
    )
    case_four_stop = evaluate_criterion(4, case_four_stop_clusters, 3, 24, 3)
    assert (2, True) == case_four_continue
    assert (2, False) == case_four_stop
    # Test invalid criterion:
    with pytest.raises(Exception) as exception_info:
        evaluate_criterion(5, pl.DataFrame(), 1, 1, 1)
        assert "valid options are" in str(exception_info.value)


# Part 17
# function moved to cluster_datatypes
# test is here
def test_perform_cutoff():
    from plantclusterfinder_py.cluster import compute_background
    from plantclusterfinder_py.cluster_datatypes import Clusters
    from plantclusterfinder_py.genelocation import GeneLoc

    gene_location = GeneLoc(
        file_path=Path("./test/data/background_computation_test.tsv")
    )
    has_reaction = pl.read_csv(
        Path("./test/data/background_computation_has_reaction.csv")
    )
    gene_location._data = gene_location.data.join(has_reaction, on="Gene Name")

    step_sizes = [1 + x for x in range(4) for _ in range(100)]
    enzyme_counts = [3 for _ in range(400)]
    window_sizes = [3 + x for _ in range(100) for x in range(4)]

    cluster_df = pl.DataFrame(
        (
            {
                "step_size": x,
                "enzyme_count": y,
                "non_enzyme_gene_count": z - y,
                "test_window_size": z,
            }
            for x, y, z in zip(step_sizes, enzyme_counts, window_sizes)
        )
    )
    background = compute_background(cluster_df, gene_location, pl.DataFrame())

    cluster = Clusters(
        ClusterMeta=cluster_df,
        ClusterBackground=background,
        FinalStepSize=6,
        MaxStepSize=5,
    )
    top_clusters = cluster.perform_cutoff(TOP_PERCENT)
    assert len(top_clusters) == len(cluster.ClusterMeta)
    assert (
        abs(
            100
            * len(cluster.ClusterMeta.filter("top_percent_clusters"))
            / len(top_clusters)
        )
        - TOP_PERCENT
        < TOP_PERCENT**0.5  # about TOP_PERCENT of clusters should be returned
    )


# Part 18
def test_write_gene_output(genelocation_transformed_obj):
    from plantclusterfinder_py.cluster import (
        get_metabolic_domains,
        label_reactions_with_sig_and_tail,
        map_pgdb_reactions_to_genes,
        match_PGDBs_to_MDs,
        perform_clustering,
        read_in_conversion_file,
        read_in_pgdb_flat,
        read_MCL_clustering_file,
        write_gene_output,
    )

    top_percent_clusters = 5
    list_of_MD_to_annotate = [
        "Amines and Polyamines Metabolism",
        "Amino Acids Metabolism",
        "Carbohydrates Metabolism",
        "Cofactors Metabolism",
        "Detoxification Metabolism",
        "Energy Metabolism",
        "Fatty Acids and Lipids Metabolism",
        "Hormones Metabolism",
        "Inorganic Nutrients Metabolism",
        "Nitrogen-Containing Compounds",
        "Nucleotides Metabolism",
        "Phenylpropanoid Derivatives",
        "Polyketides",
        "Primary-Specialized Interface Metabolism",
        "Redox Metabolism",
        "Specialized Metabolism",
        "Sugar Derivatives",
        "Terpenoids",
    ]
    metabolc_domain_reactions_filename = Path(
        "./Inputs/ReactionMetabolicDomainClassification.txt"
    )

    conversion = read_in_conversion_file(
        Path(
            "./csubellipsoidea/gtpf_CsubellipsoideaC_169_227_v2.0.annotation_info.txt.txt"
        )
    )
    gene_location = genelocation_transformed_obj
    # needs mcl clusters
    mcl_cluster_data = read_MCL_clustering_file(
        Path(
            "./csubellipsoidea/MCL_clustering_expected_results/dump.CsubellipsoideaC_169_227_v2.0.protein.pcf13.mci.I20"
        ),
        Path("./csubellipsoidea/CsubellipsoideaC_169_227_v2.0.protein.pcf13.fa"),
        conversion,
        False,  # do not remove proteins
    )

    clusters = perform_clustering(
        gene_location,
        mcl_cluster_data,
        3,  # enzymes min for clusters
        2,  # criterion
        5,  # min step size
        20,  # max step size
    )
    clusters.perform_cutoff(TOP_PERCENT)
    genes_to_reactions = map_pgdb_reactions_to_genes(pgdb_dir)
    pgdb_to_metabolic_domain_map = match_PGDBs_to_MDs(
        genes_to_reactions,
        get_metabolic_domains(
            metabolc_domain_reactions_filename, ",".join(list_of_MD_to_annotate)
        ),
    )
    test_signature_enzymes_tsv = (
        test_input_dir / "scaffold-tailoring-reactions-05082016.tab"
    )
    reactions_dat = read_in_pgdb_flat(
        pgdb_dir / "reactions.dat", ("UNIQUE-ID", "ENZYMATIC-REACTION", "EC-NUMBER")
    )

    reaction_labels = label_reactions_with_sig_and_tail(
        reactions_dat, test_signature_enzymes_tsv
    )
    # test "old" output
    # test "simple" output (default)
    simple_gene_output = NamedTemporaryFile()
    write_gene_output(
        Path(simple_gene_output.name),
        gene_location,
        clusters,
        top_percent_clusters,
        list_of_MD_to_annotate,
        pgdb_to_metabolic_domain_map,
        reaction_labels,
        "simple",
    )
    simple_gene_df = pl.read_csv(simple_gene_output.name, separator="\t")
    assert GENE_COUNT == simple_gene_df.height
    # test "verbose" output


# Part 19
def test_write_cluster_output(genelocation_transformed_obj):
    from plantclusterfinder_py.cluster import (
        perform_clustering,
        read_in_conversion_file,
        read_MCL_clustering_file,
        write_cluster_output,
    )

    conversion = read_in_conversion_file(
        Path(
            "./csubellipsoidea/gtpf_CsubellipsoideaC_169_227_v2.0.annotation_info.txt.txt"
        )
    )
    gene_location = genelocation_transformed_obj
    # needs mcl clusters
    mcl_cluster_data = read_MCL_clustering_file(
        Path(
            "./csubellipsoidea/MCL_clustering_expected_results/dump.CsubellipsoideaC_169_227_v2.0.protein.pcf13.mci.I20"
        ),
        Path("./csubellipsoidea/CsubellipsoideaC_169_227_v2.0.protein.pcf13.fa"),
        conversion,
        False,  # do not remove proteins
    )
    # using the defaults from main.py
    clusters = perform_clustering(
        gene_location,
        mcl_cluster_data,
        3,  # minimum enzymes for clusters
        2,  # criterion
        5,  # min step size
        20,  # max step size
    )
    clusters.perform_cutoff(TOP_PERCENT)

    # test "old" output
    # test "simple" output (default)
    simple_cluster_file = NamedTemporaryFile()
    write_cluster_output(
        Path(simple_cluster_file.name), gene_location, clusters, "simple"
    )
    # verify file was written
    simple_cluster_df = pl.read_csv(simple_cluster_file.name, separator="\t")
    assert 109 == simple_cluster_df.height
    # test "verbose" output
