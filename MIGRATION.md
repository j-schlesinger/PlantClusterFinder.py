# Migrating from PlantClusterFinder 1.3
## Major Differences
The following are the major differences in using PlantClusterFinder.py (2.0) vs. PlantClusterFinder 1.3.
1. The command `plantclusterfinder` is added to $PATH on install (via uv, conda, or pip). Once installing it, you should be able to run it from anywhere in your environment.
2. GTPF, GLOF, and Protein FASTA are no longer required, instead provide a GFF3 annotation file to generate them automatically. Instead of specifying `-glof`, `-gtpf`, and `-psf`, you now specify `-gff3`.
  1. Depending on your GFF3 formatting, you may need to provide `-gff3_gene` (Gene Name parameter),
  2. `-gff3_prot` (Protein Name parameter), and.
  3. `-gff3_regex` (regular expression for protein name).
3. Boolean arguments no longer require 1 or 0 to specify if it is enabled.
  1. `--ComputedMCLClusterFile`
  2. `--OverwriteSeqGapInfo`
  3. `--OverwriteMCLClustering`
  4. `--UnmaskedDNA`
  5. `--RemoveNonProtLocations`
  6. `--InsertHypos`
4. Diamond BLAST is preferred over NCBI-BLAST. (No change to user experience)
5. No longer reliant on AWK, and no intermediate files generated.
6. No option to stop and restart process partway through execution.
7. Can now optionally skip BLAST processing steps (if output already exists)
  1. `--NoRebuildBlastDB`
  2. `--NoRerunBlastP`
8. Parallel processing is enabled by default (automatically selects CPU count - 1)
## Backward Compatibility
In an effort to maintain backward compatibility as best as possible, the following options are available:
- `-glof` `-gtpf` and `-psf` are supported if **all three** arguments are provided
- NCBI-BLAST can be enabled with the `--use-ncbi-blast` argument
- `get_positions_of_gap.awk` replaced with the `get_positions_of_gap` command, added to $PATH
