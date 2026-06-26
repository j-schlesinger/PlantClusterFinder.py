"""Tests for PlantClusterFinder cluster.py functions."""

import logging
import time
from pathlib import Path

import polars as pl
import pytest

from tests.plantclusterfinder_py.assets import (
    GAPANALYSIS_EXPECTED_RESULTS,
    GENE_COUNT,
    GENOME_FASTA_PATH,
    MCL_EXPECTED_RESULTS_PATH,
    PGDB_DIR,
    PGDB_GENE_COUNT,
    PROTEIN_FASTA_PATH,
)

# globals
test_input_dir = Path("./Inputs")
_log = logging.getLogger(__name__)

TOP_PERCENT = 5
# Part 1


# Part 2 -- in test_pgdb_utils

# Part 3 match PGDB flat in test_pgdb_utils


@pytest.mark.benchmark
def test_map_pgdb_reactions_to_genes(conversion_obj):
    from plantclusterfinder_py.cluster import map_pgdb_reactions_to_genes

    mapping_genes_to_reactions = map_pgdb_reactions_to_genes(PGDB_DIR, conversion_obj)
    assert len(mapping_genes_to_reactions.ReactionsAll) == 2195
    assert len(mapping_genes_to_reactions.GenesAll) == PGDB_GENE_COUNT
    assert len(mapping_genes_to_reactions.Gene2Rxn) == PGDB_GENE_COUNT


# Part 3
def test_get_metabolic_domains():
    from plantclusterfinder_py.cluster import get_metabolic_domains

    metabolic_domains = get_metabolic_domains(
        test_input_dir / "ReactionMetabolicDomainClassification.txt",
        ["Specialized Metabolism"],
    )
    # number of rows with "Specialized Metabolism"
    assert len(metabolic_domains.RxnMetaDom) == 3_113
    assert len(metabolic_domains.Rxn2MetaDom) == 3_113


@pytest.mark.benchmark
def test_match_pgdb_to_md(conversion_obj):
    from plantclusterfinder_py.cluster import (
        get_metabolic_domains,
        map_pgdb_reactions_to_genes,
        match_pgdbs_to_mds,
    )

    metabolic_domains = get_metabolic_domains(
        test_input_dir / "ReactionMetabolicDomainClassification.txt",
        ["Specialized Metabolism"],
    )
    mapping_gene_to_reactions = map_pgdb_reactions_to_genes(PGDB_DIR, conversion_obj)
    gene2domain = match_pgdbs_to_mds(mapping_gene_to_reactions, metabolic_domains)
    assert len(gene2domain) == 261


# part 4
@pytest.mark.benchmark
def test_label_reactions_with_sig_and_tail():
    from plantclusterfinder_py.cluster import label_reactions_with_sig_and_tail
    from plantclusterfinder_py.pgdb_utils import read_in_pgdb_flat

    test_signature_enzymes_tsv = test_input_dir / "scaffold-tailoring-reactions-05082016.tab"
    reactions_dat = read_in_pgdb_flat(
        PGDB_DIR / "reactions.dat", ("UNIQUE-ID", "ENZYMATIC-REACTION", "EC-NUMBER")
    )

    labeled_reaction = label_reactions_with_sig_and_tail(reactions_dat, test_signature_enzymes_tsv)
    assert len(labeled_reaction.STClasses) == 9
    assert len(labeled_reaction.SigOrTail) == 3
    assert len(labeled_reaction.Rxn_STClasses) == 394
    assert len(labeled_reaction.Rxn_SigOrTail) == 394


# part 5
def test_map_pgdb_gene_ids_to_conversion_file(conversion_obj):
    from plantclusterfinder_py.cluster import map_pgdb_gene_ids_to_conversion_file
    from plantclusterfinder_py.pgdb_utils import read_in_pgdb_flat

    conversion = conversion_obj
    genes_all = read_in_pgdb_flat(
        PGDB_DIR / "genes.dat",
        ("UNIQUE-ID", "ACCESSION-1", "PRODUCT", "COMMON-NAME"),
        unique_keys=("UNIQUE-ID",),
    )
    gene_info_from_genes_dat = ["ACCESSION-1", "COMMON-NAME"]
    ids_to_map = "GTP"
    pgdb_conversions_df = map_pgdb_gene_ids_to_conversion_file(
        conversion, genes_all, gene_info_from_genes_dat, ids_to_map
    )
    assert len(pgdb_conversions_df) == PGDB_GENE_COUNT
    g4s_3028_row = pgdb_conversions_df.filter(pl.col("pgdb_gene_id") == "G4S-3028")
    print(g4s_3028_row)
    assert g4s_3028_row["Protein Name"][0] == "46517"
    assert g4s_3028_row["Transcript Name"][0] == "46517"


# part 6
def test_read_mcl_clustering_file(conversion_obj):
    from plantclusterfinder_py.cluster import read_mcl_clustering_file

    conversion = conversion_obj
    mcl_cluster_data = read_mcl_clustering_file(
        MCL_EXPECTED_RESULTS_PATH,
        PROTEIN_FASTA_PATH,
        conversion,
        False,  # do not remove proteins
    )
    assert len(mcl_cluster_data.Gene2MCLClusters) == GENE_COUNT - 2
    assert isinstance(mcl_cluster_data.Protein_IDs, list)
    assert len(mcl_cluster_data.Protein_IDs) == GENE_COUNT - 2
    assert mcl_cluster_data.MCL_Clusters.schema == {
        "ID": pl.String,
        "Transcript_IDs": pl.List(pl.String),
        "Gene_Names": pl.List(pl.String),
        "Cluster_Size": pl.Int64,
    }


def test_get_mcl_clustering(conversion_obj):
    from plantclusterfinder_py.cluster import (
        get_mcl_clustering,
        read_mcl_clustering_file,
    )

    conversion = conversion_obj
    mcl_cluster_diamond = get_mcl_clustering(
        PROTEIN_FASTA_PATH,
        conversion,
        True,  # overwrite MCL cluster file
        False,  # remove protein locations
        rebuild_blast_db=True,
        rerun_blastp=True,
        use_diamond=True,
    )
    mcl_cluster_blast = get_mcl_clustering(
        PROTEIN_FASTA_PATH,
        conversion,
        True,  # overwrite MCL cluster file
        False,  # remove protein locations
        rebuild_blast_db=True,
        rerun_blastp=True,
        use_diamond=False,
    )
    expected_mcl_cluster_data = read_mcl_clustering_file(
        MCL_EXPECTED_RESULTS_PATH,
        PROTEIN_FASTA_PATH,
        conversion,
        False,  # do not remove proteins
    )
    assert len(mcl_cluster_diamond.Gene2MCLClusters) == GENE_COUNT - 1
    # Note: Diamond's algorithm yields slightly different results.
    assert (
        len(expected_mcl_cluster_data.Gene2MCLClusters)
        == len(mcl_cluster_diamond.Gene2MCLClusters) - 1
    )
    assert len(mcl_cluster_blast.Gene2MCLClusters) == GENE_COUNT - 2
    assert len(expected_mcl_cluster_data.Gene2MCLClusters) == len(mcl_cluster_blast.Gene2MCLClusters)


# Part 7 through Part 10: No Cluster Functions


# part 11
def test_analyze_plant_cluster_gap_file():
    """Used by get_sequencing_gaps."""
    from plantclusterfinder_py.cluster import analyze_plant_cluster_gap_file

    test_gap_output = GAPANALYSIS_EXPECTED_RESULTS
    plant_cluster_analysis = analyze_plant_cluster_gap_file(
        test_gap_output,
        False,
        file_output=False,
    )
    assert isinstance(plant_cluster_analysis, tuple)
    (filename, gap_output_count) = plant_cluster_analysis
    assert len(gap_output_count) == 337  # TEST INCOMPLETE
    assert str(test_gap_output) + "_count.txt" == filename  # TEST INCOMPLETE


@pytest.mark.benchmark
def test_annotate_sequencing_gaps(genelocation_obj):
    """Tests annotate_sequencing_gaps and get_sequencing_gaps."""
    from plantclusterfinder_py.cluster import annotate_sequencing_gaps

    sequencing_gaps_file = GAPANALYSIS_EXPECTED_RESULTS
    sequencing_gaps_fasta = sequencing_gaps_file.with_suffix(
        sequencing_gaps_file.suffix.split("_")[0]
    )
    gene_location = genelocation_obj
    uniq_chr = list(gene_location.data.select("Chromosome Name").unique().to_dict().values())[
        0
    ].to_list()
    sequencing_gaps = annotate_sequencing_gaps(
        sequencing_gaps_fasta, gene_location.find_intergenic_regions(), uniq_chr
    )
    assert sequencing_gaps.height == 1166
    assert len(sequencing_gaps.columns) == 6
    sequencing_gaps_test = annotate_sequencing_gaps(
        GENOME_FASTA_PATH, gene_location.find_intergenic_regions(), uniq_chr
    )
    assert sequencing_gaps_test.height == sequencing_gaps.height
    sequencing_gaps_eval = sequencing_gaps.join(
        sequencing_gaps_test, on=(pl.col("scaffold"), pl.col("start index"))
    ).filter(pl.col("seq_gaps_length") == pl.col("seq_gaps_length_right"))
    assert sequencing_gaps_eval.height == sequencing_gaps.height


# part 12
# hard breaks moved to gene location test


# part 13
@pytest.mark.benchmark
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
    assert len(gene_location) == orig_gene_count * 2 - 216  # TEST INCOMPLETE
    assert orig_gene_count == non_hypo_count
    assert orig_gene_count != 0


# Part 14
# in gene location test


# Part 15
# in gene location test


# Part 16
@pytest.mark.benchmark
def test_compute_background():
    """Used by perform_clustering."""
    from plantclusterfinder_py.cluster import compute_background
    from plantclusterfinder_py.genelocation import GeneLoc

    gene_location = GeneLoc(file_path=Path("./tests/data/background_computation_test.tsv"))
    has_reaction = pl.read_csv(Path("./tests/data/background_computation_has_reaction.csv"))
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
            for x, y, z in zip(step_sizes, enzyme_counts, window_sizes, strict=True)
        )
    )
    # cluster.write_csv(tmp_path / "test_cluster.csv")
    expected_background = pl.DataFrame(
        [
            {
                "window_size": x,
                "enzyme_count": y,
                "window_occurrence": 192 - 6 * x,
                "tested": True,
            }
            for x in range(1, 7)
            for y in range(x + 1)
        ]
    ).filter(  # filtering out background cases not seen in data
        (pl.col("window_size") <= 4) | (pl.col("window_size") != pl.col("enzyme_count"))
    )
    print(expected_background)
    background = compute_background(cluster, gene_location, pl.DataFrame())

    assert expected_background.height == background.height
    ex_grouped = expected_background.group_by("window_size").agg(
        pl.last("window_occurrence"),
    )
    for window in (
        background.group_by("window_size")
        .agg(
            pl.sum("winzyme_occurrence"),
            pl.sum("winzyme_frequency"),
            pl.min("window_occurrence"),
        )
        .iter_rows()
    ):
        (size, winzyme_occurrence_sum, _, total_window_occurrence) = (
            window  # not testing winzyme frequency
        )
        print(size, total_window_occurrence)
        assert winzyme_occurrence_sum == total_window_occurrence
        assert total_window_occurrence == 192 - size * 6
        assert (
            total_window_occurrence
            == ex_grouped.filter(pl.col("window_size") == size).select(pl.col("window_occurrence"))[
                "window_occurrence"
            ][0]
        )


@pytest.mark.benchmark
def test_perform_clustering(genelocation_transformed_obj, expected_mcl_cluster_obj):
    from plantclusterfinder_py.cluster import perform_clustering

    gene_location = genelocation_transformed_obj
    mcl_cluster_data = expected_mcl_cluster_obj

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
    assert clusters.FinalStepSize == 3
    assert clusters.MaxStepSize == 5
    assert clusters.ClusterMeta.height == 1699
    assert clusters.ClusterBackground.height == 848


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
    assert case_one_continue == (2, True)
    assert case_one_stop == (2, False)
    # Evaluate 2: continue search until step_size is 3
    empty_clusters = pl.DataFrame([cluster_tuple("TEST", 0, 0, 0, 0, 0, 0)])
    #    Continue search: step size less than 3
    case_two_continue = evaluate_criterion(2, empty_clusters, 2, 24, 3)
    #    Stop search: step size 3 or greater
    case_two_stop = evaluate_criterion(2, empty_clusters, 3, 24, 3)
    assert case_two_continue == (1, True)
    assert case_two_stop == (3, False)
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
    assert case_three_continue == (2, True)
    assert case_three_stop == (2, False)
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
    assert case_four_continue == (2, True)
    assert case_four_stop == (2, False)
    # Test invalid criterion:
    with pytest.raises(ValueError, match="valid options are") as _:
        evaluate_criterion(5, pl.DataFrame(), 1, 1, 1)


# Part 17
# function moved to cluster_datatypes
# test is here
def test_perform_cutoff():
    from plantclusterfinder_py.cluster import compute_background
    from plantclusterfinder_py.cluster_datatypes import Clusters
    from plantclusterfinder_py.genelocation import GeneLoc

    gene_location = GeneLoc(file_path=Path("./tests/data/background_computation_test.tsv"))
    has_reaction = pl.read_csv(Path("./tests/data/background_computation_has_reaction.csv"))
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
            for x, y, z in zip(step_sizes, enzyme_counts, window_sizes, strict=True)
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
        abs(100 * len(cluster.ClusterMeta.filter("top_percent_clusters")) / len(top_clusters))
        - TOP_PERCENT
        < TOP_PERCENT**0.5  # about TOP_PERCENT of clusters should be returned
    )


# Part 18
def test_write_gene_output(
    genelocation_transformed_obj, expected_mcl_cluster_obj, conversion_obj, tmp_path
):
    from plantclusterfinder_py.cluster import (
        get_metabolic_domains,
        label_reactions_with_sig_and_tail,
        map_pgdb_reactions_to_genes,
        match_pgdbs_to_mds,
        perform_clustering,
        write_gene_output,
    )
    from plantclusterfinder_py.pgdb_utils import read_in_pgdb_flat

    gene_location = genelocation_transformed_obj
    mcl_cluster_data = expected_mcl_cluster_obj
    write_gene_start_time = int(time.time())

    top_percent_clusters = 5
    list_of_md_to_annotate = [
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
    metabolc_domain_reactions_filename = Path("./Inputs/ReactionMetabolicDomainClassification.txt")

    clusters = perform_clustering(
        gene_location,
        mcl_cluster_data,
        3,  # enzymes min for clusters
        2,  # criterion
        5,  # min step size
        20,  # max step size
    )
    clusters.perform_cutoff(TOP_PERCENT)
    genes_to_reactions = map_pgdb_reactions_to_genes(PGDB_DIR, conversion_obj)
    pgdb_to_metabolic_domain_map = match_pgdbs_to_mds(
        genes_to_reactions,
        get_metabolic_domains(metabolc_domain_reactions_filename, list_of_md_to_annotate),
    )
    test_signature_enzymes_tsv = test_input_dir / "scaffold-tailoring-reactions-05082016.tab"
    reactions_dat = read_in_pgdb_flat(
        PGDB_DIR / "reactions.dat", ("UNIQUE-ID", "ENZYMATIC-REACTION", "EC-NUMBER")
    )

    reaction_labels = label_reactions_with_sig_and_tail(reactions_dat, test_signature_enzymes_tsv)
    # test "old" output
    # test "simple" output (default)
    simple_gene_output = tmp_path / f"tmp_gene_output_test_write.{write_gene_start_time}.tsv"
    write_gene_output(
        Path(simple_gene_output.name),
        gene_location,
        clusters,
        top_percent_clusters,
        list_of_md_to_annotate,
        pgdb_to_metabolic_domain_map,
        reaction_labels,
        "simple",
    )
    simple_gene_df = pl.read_csv(simple_gene_output.name, separator="\t")
    assert simple_gene_df.height == GENE_COUNT
    # test "verbose" output


# Part 19
def test_write_cluster_output(genelocation_transformed_obj, tmp_path):
    from plantclusterfinder_py.cluster import (
        perform_clustering,
        read_mcl_clustering_file,
        write_cluster_output,
    )
    from plantclusterfinder_py.conversion import Conversion

    conversion = Conversion.read_in_conversion_file(
        Path("./csubellipsoidea/gtpf_CsubellipsoideaC_169_227_v2.0.annotation_info.txt.txt")
    )
    write_cluster_start_time = int(time.time())
    gene_location = genelocation_transformed_obj
    # needs mcl clusters
    mcl_cluster_data = read_mcl_clustering_file(
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
    simple_cluster_file = tmp_path / f"tmp_cluster_output_test_write.{write_cluster_start_time}.tsv"
    write_cluster_output(Path(simple_cluster_file.name), gene_location, clusters, "simple")
    # verify file was written
    simple_cluster_df = pl.read_csv(simple_cluster_file.name, separator="\t")
    assert simple_cluster_df.height == 258
    # test "verbose" output
    verbose_cluster_file = (
        tmp_path / f"tmp_cluster_output_test_write_verbose.{write_cluster_start_time}.tsv"
    )
    write_cluster_output(Path(verbose_cluster_file.name), gene_location, clusters, "verbose")
    simple_cluster_df = pl.read_csv(verbose_cluster_file.name, separator="\t")
    assert simple_cluster_df.height == 477
