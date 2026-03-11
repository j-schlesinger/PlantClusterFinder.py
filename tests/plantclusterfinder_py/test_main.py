import time


def test_run_main(tmp_path):
    from plantclusterfinder_py.main import main

    test_start = int(time.time())
    gene_output_temp_file = tmp_path / f"gene_output_test.{test_start}.tsv"
    cluster_output_temp_file = tmp_path / f"cluster_output_test.{test_start}.tsv"

    args = [
        "--GenesDat_GeneInfo",
        "ACCESSION-1",
        "-md",
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
        "-dnaf",
        "./csubellipsoidea/CsubellipsoideaC_169_227_v2.0.hardmasked.fa",
        "-gff3",
        "./csubellipsoidea/CsubellipsoideaC_169_227_v2.0.gene.gff3",
        "-pgdb",
        "./csubellipsoidea/pgdb/csubellipsoideacyc/1.0/data/",
        "-sitf",
        "./Inputs/scaffold-tailoring-reactions-05082016.tab",
        "-rmdf",
        "./Inputs/ReactionMetabolicDomainClassification.txt",
        "-gout",
        gene_output_temp_file.name,
        "-cout",
        cluster_output_temp_file.name,
        "--SeqGapSizesChromBreak",
        "10000",
        "--PGDBIdsToMap",
        "GTP",
        "--OutputType",
        "simple",
        "--Verbose",
        "1",
        "--MaxInterGeneDistByMedian",
        "120",
        "--OutputGeneLocationParquet",
        "./Outputs/csubellipsoidea_gene_loc.parquet",
    ]
    main(args)
