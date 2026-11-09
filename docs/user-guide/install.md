# Installation

[TOC]

## Conda

Difficulty: Easy
Customization: None

```bash
# mamba or pixi are also ok
# FIXME:
conda env create bgc  bioconda::plantclusterfinder_py
```

!!! note
    If install is successful, the command `plantclusterfinder` will be added to your $PATH.

Run `plantclusterfinder --help` to test if it ran, and see the help message. The first time you run
PlantClusterFinder, it will download and install necessary data files to your $USER directory.

## Python Installation

Difficulty: Medium
Customization: Minimal

!!! important
    To run PlantClusterFinder.py, you need the following prerequisite programs installed:
    - [Diamond](https://github.com/bbuchfink/diamond/wiki/2.-Installation)
    - [MCL](https://github.com/micans/mcl#installation-mcl-versions-and-linux-packages)
    You may need to use `apt install diamond-aligner mcl` or `yum install diamond mcl` depending on your
    package manager.

```bash
# using pip
pip install plantclusterfinder_py
# OR using uv
uv pip install plantclusterfinder_py
```

## Build from Source

Difficulty: High
Customization: High

You will need to install the prerequisite programs from the [[#Python Installation|Python Installation Section]].
In addition, you will need [UV](https://docs.astral.sh/uv) installed.

Then clone the repository, and sync UV:

```bash
git clone https://github.com/j-schlesinger/PlantClusterFinder.py.git
cd PlantClusterFinder.py
uv sync
```

Now you can test your installation was successful using:

```bash
uv run plantclusterfinder --help
```

## Testing Your Installation

If your installation was successful, the output of the above should be the help message, but you will
want to confirm the output of your installation matches with the test output before proceeding.

```bash
plantclusterfinder --test-install
```

This will run PlantClusterFinder on a known input and compare the results to known good results.
If this process fails, there may be an issue with your installation.
Check the [Troubleshooting](troubleshooting.md) section to see if your issue is known.
