from pathlib import Path

import polars as pl

# Part 1

# globals
pgdb_dir = Path("./csubellipsoidea/pgdb/csubellipsoideacyc/1.0/data")
test_input_dir = Path("./Inputs")


def test_read_in_conversion_file():
    """The conversion file defines the genes, proteins, and reactions
    that many downstream processes rely on."""

    from plantclusterfinder_py.cluster import read_in_conversion_file

    conversion = read_in_conversion_file(
        Path(
            "./csubellipsoidea/gtpf_CsubellipsoideaC_169_227_v2.0.annotation_info.txt.txt"
        )
    )
    assert len(conversion.MapG) == 9629
    assert len(conversion.MapT) == 9629
    assert len(conversion.MapP) == 9629


# Part 2
def test_read_in_pgdb_flat():
    """Used by map_pgdb_reactions_to_genes"""

    from plantclusterfinder_py.cluster import read_in_pgdb_flat

    gene_df = read_in_pgdb_flat(
        pgdb_dir / "genes.dat", ("UNIQUE-ID", "PRODUCT", "ACCESSION-1")
    )
    assert isinstance(gene_df, pl.DataFrame)
    assert len(gene_df) == 2053
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
    assert 2_053 == len(integration_matched)  # TEST INCOMPLETE
    assert 2_053 == len(genes_flat)
    assert 2_054 == len(proteins_flat)
    assert 3_015 == len(enzrxns_all)
    assert len(reaction_flat) > 0


def test_map_pgdb_reactions_to_genes():
    from plantclusterfinder_py.cluster import map_pgdb_reactions_to_genes

    mapping_G2R = map_pgdb_reactions_to_genes(pgdb_dir)
    assert 2195 == len(mapping_G2R.ReactionsAll)
    assert 2053 == len(mapping_G2R.GenesAll)
    assert 2053 == len(mapping_G2R.Gene2Rxn)


# Part 3
def test_get_metabolic_domains():
    from plantclusterfinder_py.cluster import get_metabolic_domains

    metabolic_domains = get_metabolic_domains(
        test_input_dir / "ReactionMetabolicDomainClassification.txt",
        "Specialized Metabolism",
    )
    assert 11_435 == len(metabolic_domains.RxnMetaDom)


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
    from plantclusterfinder_py.cluster import label_reactions_with_sig_and_tail

    reactions_dat = pgdb_dir / "reactions.dat"
    test_signature_enzymes_tsv = (
        test_input_dir / "scaffold-tailoring-reactions-05082016.tab"
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
        map_PGDBs_to_conversion_file,
        read_in_conversion_file,
        read_in_pgdb_flat,
    )

    # from plantclusterfinder_py.cluster_datatypes import Conversion

    conversion = read_in_conversion_file(
        "./csubellipsoidea/gtpf_CsubellipsoideaC_169_227_v2.0.annotation_info.txt.txt"
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
    assert 2_053 == len(pgdb_conversions_dict)  # TEST INCOMPLETE
    assert ["46517"] == pgdb_conversions_dict["G4S-3028"]


# part 6
def test_read_MCL_clustering_File():
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
    assert 9627 == len(mcl_cluster_data.Gene2MCLClusters)  # TEST INCOMPLETE
    assert 9_627 == len(mcl_cluster_data.Protein_IDs)  # TEST INCOMPLETE


def test_get_MCL_clustering():
    from plantclusterfinder_py.cluster import (
        get_MCL_clustering,
        read_in_conversion_file,
    )

    conversion = read_in_conversion_file(
        "./csubellipsoidea/gtpf_CsubellipsoideaC_169_227_v2.0.annotation_info.txt.txt"
    )
    mcl_cluster = get_MCL_clustering(
        "./csubellipsoidea/CsubellipsoideaC_169_227_v2.0.protein.pcf13.fa",
        conversion,
        True,  # overvrite MCL cluster file
        False,  # remove protein locations
        vParallelize=1,
    )
    assert 9_627 == len(mcl_cluster.Gene2MCLClusters)


# part 7


# Part 10: No Cluster Functions


# part 11
def test_analyze_plant_cluster_gap_file():
    """Used by get_sequencing_gaps"""
    from plantclusterfinder_py.cluster import analyze_plant_cluster_gap_file

    test_gap_output = Path(
        "./csubellipsoidea/Gapanalysis_expected_results/CsubellipsoideaC_169_227_v2.0.hardmasked.fa_GAPOutput"
    )
    (filename, gap_output_count) = analyze_plant_cluster_gap_file(
        test_gap_output,
        False,
        file_output=False,
    )
    assert 337 == len(gap_output_count)  # TEST INCOMPLETE
    assert str(test_gap_output) + "_count.txt" == filename  # TEST INCOMPLETE


def test_get_sequencing_gaps():
    """Used by annotate_sequencing_gaps"""
    # this is better tested in get_positions_of_gap
    assert 1 == 1  # TEST INCOMPLETE


def test_annotate_sequencing_gaps():
    # from plantclusterfinder_py.cluster import annotate_sequencing_gaps
    assert 1 == 0  # TEST INCOMPLETE


# part 12


def test_introduce_hard_breaks():
    assert 1 == 0  # TEST INCOMPLETE


# part 13
def test_introduce_hypo_genes():
    assert 1 == 0  # TEST INCOMPLETE


# Part 14
# gene location test


# Part 15
def test_compute_local_duplication():
    assert 1 == 0  # TEST INCOMPLETE


# Part 16
def test_compute_background():
    """Used by perform_clustering"""
    assert 1 == 0  # TEST INCOMPLETE


def test_perform_clustering():
    assert 1 == 0  # TEST INCOMPLETE


# Part 17
def test_perform_cutoff():
    assert 1 == 0  # TEST INCOMPLETE


# Part 18
def test_write_gene_output():
    assert 1 == 0  # TEST INCOMPLETE


# Part 19
def test_write_cluster_output():
    assert 1 == 0  # TEST INCOMPLETE
