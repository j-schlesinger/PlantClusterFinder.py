# About PlantClusterFinder.py

PlantClusterFinder.py performs genome-wide biosynthetic gene cluster prediction and annotation.
It is a complete rewrite of the MATLAB tool [PlantClusterFinder](https://github.com/carnegie/PlantClusterFinder) in Python with paralellism and plugin support.

## News

### PlantClusterFinder.py Released
**PlantClusterFinder.py** is being released as the second major version of PlantClusterFinder.
It accepts similar inputs to, and produces compatible outputs with PlantClusterFinder, but the execution arguments are not compatible.
However, any researcher using PlantClusterFinder is recommended to switch to using this version for future research, as all future updates will only apply to the Python version.

New Features in PlantClusterFinder.py include:
- Easy install using Conda or Pip
- Support for plugins to automate or extend analysis
- Significant performance improvements over PlantClusterFinder
- Parquet export of raw gene metadata
