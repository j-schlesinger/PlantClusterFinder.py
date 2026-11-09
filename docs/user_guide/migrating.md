# Migrating from PlantClusterFinder v1.0

For users of the previous PlantClusterFinder, many parameters may seem similar but there are some major differences to be aware of.
First off, this is a *major* version change, and so any scripts you have are not guaranteed to work on this version of PlantClusterFinder.
Second off, generation of the GLOF, GTPF, and protein FASTA files are now optional.
Third, execution is parallel (CPU count - 1) by default, so use the `-para` parameter to modify it if needed.
Lastly, the convention of setting parameters to `-1` to disable them should be considered as deprecated and will be removed in a future version.
