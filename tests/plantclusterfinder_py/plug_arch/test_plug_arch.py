import pytest

from tests.plantclusterfinder_py.assets import GENOME_FASTA_PATH


def test_plug_arch(tmp_path):
    from plantclusterfinder_py.plugin_arch.plugin_main import main

    main([str(GENOME_FASTA_PATH), tmp_path])
