# Built-in Plugins
Plugins have been built-in to PlantClusterFinder.py to facilitate research, and provide an example for further development of plugins.
## Biosynthetic Gene Cluster Statistics Aggregator Plugin
This plugin parses the cluster and gene output files and generates an additional file summarizing key statistics about the results.
## Graph Biosynthetic Gene Clusters
This plugin plots the location of discovered BGCs on the chromosome.
By default, this plugin group all BGCs together, but additional configuration may group BGCs differently.
## Automated Filtering
This plugin requires a set of "positive" BGCs to guarantee are included.
It then adjusts the filtering parameters to include all positive BGCs while still being as restrictive as possible.
