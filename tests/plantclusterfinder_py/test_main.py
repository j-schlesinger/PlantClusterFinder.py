import hashlib
import time
from pathlib import Path

import polars as pl
import pytest


@pytest.mark.benchmark
def test_run_main(tmp_path):
    from plantclusterfinder_py.main import main

    test_start = int(time.time())
    gene_output_temp_file = tmp_path / f"gene_output_test.{test_start}.tsv"
    cluster_output_temp_file = tmp_path / f"cluster_output_test.{test_start}.tsv"
    parquet_output_temp_file = tmp_path / f"parquet_output_test.{test_start}.tsv"

    args = [
        "-r",
        "PyTest_Execution",
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
        parquet_output_temp_file.name,
    ]
    main(args)

    gene_location = pl.read_parquet(parquet_output_temp_file.name)
    assert len(gene_location) == 9629
    genes_file = pl.read_csv(gene_output_temp_file.name, separator="\t")
    assert len(genes_file) == 9629
    clusters_file = pl.read_csv(cluster_output_temp_file.name, separator="\t")
    assert len(clusters_file) == 123


def test_download_input(tmp_path):
    from plantclusterfinder_py.main import download_input_files

    assert download_input_files(tmp_path)
    with (tmp_path / "scaffold-tailoring-reactions-05082016.tab").open("rb") as f:
        print(tmp_path)
        test_digest_s = hashlib.sha256(f.read()).digest()
    with (tmp_path / "ReactionMetabolicDomainClassification.txt").open("rb") as f:
        test_digest_r = hashlib.sha256(f.read()).digest()
    with (Path("./Inputs") / "scaffold-tailoring-reactions-05082016.tab").open("rb") as f:
        expect_digest_s = hashlib.sha256(f.read()).digest()
    with (Path("./Inputs") / "ReactionMetabolicDomainClassification.txt").open("rb") as f:
        expect_digest_r = hashlib.sha256(f.read()).digest()
    assert test_digest_s == expect_digest_s
    assert test_digest_r == expect_digest_r
