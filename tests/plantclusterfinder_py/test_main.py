import hashlib
import logging
import time
from pathlib import Path
from unittest.mock import patch

import polars as pl
import pytest

log = logging.getLogger()


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
    assert len(clusters_file) == 322


def test_download_input(tmp_path):
    """Download the static Input files."""
    from plantclusterfinder_py.main import download_input_files

    with patch("requests.get") as mock_get:
        mock_get.return_value.status_code = 200
        mock_get.return_value.content = b"TEST_CONTENT"
        assert download_input_files(tmp_path)
    """
    with (tmp_path / "scaffold-tailoring-reactions-05082016.tab").open("rb") as f:
        test_digest_s = hashlib.sha256(f.read()).digest()
    with (tmp_path / "ReactionMetabolicDomainClassification.txt").open("rb") as f:
        test_digest_r = hashlib.sha256(f.read()).digest()
    with (Path("./Inputs") / "scaffold-tailoring-reactions-05082016.tab").open("rb") as f:
        expect_digest_s = hashlib.sha256(f.read()).digest()
    with (Path("./Inputs") / "ReactionMetabolicDomainClassification.txt").open("rb") as f:
        expect_digest_r = hashlib.sha256(f.read()).digest()
    assert test_digest_s == expect_digest_s
    assert test_digest_r == expect_digest_r
    """


def test_download_csub(tmp_path):
    """Download the Csubellipsoidea data."""
    import subprocess

    from plantclusterfinder_py.main import download_csub_files

    temp_tgz = tmp_path / "csubellipsoidea_old.tar.gz"
    subprocess.run(["/usr/bin/env", "tar", "-czvf", temp_tgz, Path("./csubellipsoidea/")])

    with (
        patch("requests.get") as mock_get,
        temp_tgz.open("rb") as csub_data,
    ):
        mock_get.return_value.status_code = 200
        mock_get.return_value.content = csub_data.read()
        assert download_csub_files(tmp_path)
    """
    test_digests = {}
    expected_digests = {}
    for file in (tmp_path / "csubellipsoidea").iterdir():
        if not file.is_file():
            continue
        with file.open("rb") as f:
            test_digests[file.name] = hashlib.sha256(f.read()).digest()
    for file in (Path("./csubellipsoidea")).iterdir():
      m if not file.is_file():
            continue
        with file.open("rb") as f:
            expected_digests[file.name] = hashlib.sha256(f.read()).digest()
    for filename in test_digests:
        if filename not in expected_digests:
            log.warning("Filename %s not in expected digests", filename)
            print(list(expected_digests.keys())[:5])
        elif test_digests[filename] != expected_digests[filename]:
            print(f"Hashes do not match for file {filename}")
            with (tmp_path / filename).open("rb") as file:
                print(file.read())
            assert test_digests[filename] == expected_digests[filename]
    """
