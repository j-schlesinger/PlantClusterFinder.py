from itertools import zip_longest
from pathlib import Path

import polars as pl

from tests.plantclusterfinder_py.assets import (
    GAPANALYSIS_EXPECTED_RESULTS,
    GENOME_FASTA_PATH,
)


def test_char_run_search():
    from plantclusterfinder_py.get_positions_of_gap import char_run_search

    assert list(char_run_search("A", "BANANA")) == list(zip([2, 4, 6], [1, 1, 1], strict=True))
    assert (
        len(
            list(
                char_run_search(
                    "aoeiu",
                    "Now is the time for all good men to come to the aid of their community.",
                )
            )
        )
        == 20
    )
    assert 18 in [
        x
        for x, _ in char_run_search(
            " \n", "Haikus can be fun\nbut sometimes they don't make sense\nrefigerator"
        )
    ]
    assert (
        list(
            char_run_search("259", "You have either chosen, or been chosen to relocate to City 17.")
        )
        == []
    )
    test_start_end = list(char_run_search("A", "AAAAAAABAAAAAABBBBBBBABAAA"))
    print(test_start_end)
    assert test_start_end == [(1, 7), (9, 6), (22, 1), (24, 3)]


def test_read_fasta():
    from plantclusterfinder_py.get_positions_of_gap import read_fasta

    input_path = Path("tests/data/test_fasta.fa")
    with input_path.open() as fa_fh:
        fasta_gen = read_fasta(fa_fh)
        seq = dict(fasta_gen)
    assert ">Test Sequence 1" in seq
    assert ">Test Sequence 2" in seq
    assert len(seq[">Test Sequence 1"]) == 143
    assert len(seq[">Test Sequence 2"]) == 80 * 16


def test_get_positions_of_gap(tmp_path):
    from plantclusterfinder_py.get_positions_of_gap import get_positions_of_gap

    test_d = {}
    output_path = tmp_path / "groupbyseq.csv"  # output
    input_path = Path("tests/data/test_fasta.fa")
    get_positions_of_gap(input_path, output_path)
    with output_path.open() as test_output:
        for line in test_output:
            nuc = line[5]
            test_d[nuc] = test_d.get(nuc, 0) + 1
    assert test_d["N"] == 18
    assert test_d["A"] == 230

    output_path2 = tmp_path / "groupbyseq2.csv"  # output
    get_positions_of_gap(input_path, output_path2, parallelize=10)
    with output_path.open() as test_output, output_path2.open() as test_parallel:
        for test_line, para_line in zip_longest(test_output, test_parallel):
            assert test_line == para_line
    output_path.unlink()
    output_path2.unlink()


def test_get_positions_of_gap_against_expected(tmp_path):
    from plantclusterfinder_py.get_positions_of_gap import get_positions_of_gap

    output_path = tmp_path / "groupbyseq.csv"  # output
    input_path = GENOME_FASTA_PATH
    expected_output = GAPANALYSIS_EXPECTED_RESULTS
    get_positions_of_gap(input_path, output_path)
    test_d = {}
    with output_path.open() as test_output:
        for line in test_output:
            (
                _,  # chr
                nuc,
                _,  # pos
                _,  # len
            ) = line.split(" ")
            test_d[nuc] = test_d.get(nuc, 0) + 1
    ex_d = {}
    with expected_output.open() as ex_output, output_path.open() as test_out:
        for ex_line, test_line in zip_longest(ex_output, test_out):
            assert test_line == ex_line
            (
                _,  # chr
                nuc,
                _,  # pos
                _,  # len
            ) = ex_line.split(" ")
            ex_d[nuc] = ex_d.get(nuc, 0) + 1
    for nuc, value in ex_d.items():
        assert test_d[nuc] == value
    output_path.unlink()


def test_positions_of_gap_main(tmp_path):
    from plantclusterfinder_py.get_positions_of_gap import main

    test_output = tmp_path / "test.csv"
    args = [
        "-i",
        "csubellipsoidea/CsubellipsoideaC_169_227_v2.0.hardmasked.fa",
        "-o",
        str(test_output),
        "-p",
        "10",
    ]
    main(args)
    gap_positions_df = pl.read_csv(
        test_output,
        schema={
            "chromosome": pl.String,
            "nucleotide": pl.String,
            "start": pl.UInt64,
            "length": pl.UInt64,
        },
        separator=" ",
        has_header=False,
    )
    assert len(gap_positions_df) == 8828
