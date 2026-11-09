from pathlib import Path

import pytest

from tests.plantclusterfinder.assets import GENOME_FASTA_PATH


def test_plug_arch(tmp_path, request):
    from plantclusterfinder.plugin_arch.plugin_main import test_plugins

    original_root = request.config.rootpath
    cluster_filename = original_root / "tests/data/example_cluster.tsv"
    plugin_outputs = test_plugins(tmp_path, GENOME_FASTA_PATH, cluster_filename)
    for path in plugin_outputs.values():
        assert path.is_file()
