# PlantClusterFinder.py CLI Arguments
```ascii
usage: PlantClusterFinder [-h] --RunName RUNNAME -gout FILENAME -cout FILENAME
                          [--RunInfoFile FILENAME]
                          [--LogFile FILENAME]
                          --PathwayGenomeDatabasePath FILENAME
                          --GeneAnnoatationGFF3 FILENAME
                          --DNAFastaFilename FILENAME
                          --MetabolicDomainsToAnnotate DOMAIN [DOMAIN ...]
                          [--MaxSeqGapSize GAP]
                          [--MaxInterGeneDist DISTANCE]
                          [--MaxInterGeneDistByMedian DISTANCE]
                          [--PercentileForMaxInterGeneDist DISTANCE]
                          [--SeqGapSizesChromBreak SIZE [SIZE ...]]
                          [--EnzMinForClusters MINIMUM]
                          [--KbRangeLD KB_RANGE]
                          [--GeneRangeLD GENES]
                          [--TopPercentClusters PERCENT]
                          [--MinStepSize STEP]
                          [--MaxStepSize STEP]
                          [--Criterion {1,2,3,4}]
                          [--HypoGenePercentile PERCENTILE]
                          [--InsertHypos]
                          [--HypoAmount AMOUNT]
                          [--GFF3GeneAttribute ATTRIBUTE]
                          [--GFF3ProteinAttribute ATTRIBUTE]
                          [--GFF3ProteinRegex REGEX]
                          [--ProteinSequenceFasta FILENAME]
                          [--GeneLocationTSV FILENAME]
                          [--GeneTranscriptProteinFilename FILENAME]
                          [--GenesDat_GeneInfo COLUMN_NAME [COLUMN_NAME ...]]
                          [--ComputedMCLClusterFile]
                          [--OverwriteSeqGapInfo]
                          [--OverwriteMCLClustering]
                          [--Verbose VERBOSITY]
                          [-para PARALLELIZATION]
                          [--UnmaskedDNA]
                          [--PGDBIdsToMap PGDBIDSTOMAP]
                          [--RemoveNonProtLocations]
                          [--OutputType TYPE]
                          [--OutputGeneLocationParquet FILENAME]
                          [--use_ncbi_blast]
                          [--ReactionMetabolicDomainFilename FILENAME]
                          [--SignatureTailoringFilename FILENAME]
Find Metabolic Gene Clusters in Plant Genomes
```

## Mandatory Arguments
All of these arguments must be set for every run.

- `RunName`: Name of run of PlantClusterFinder. Helpful for identifying execution across different organisms/samples.
- `PathwayGenomeDatabasePath`: Full path of PGDB flat file folder, where the flat files of the species pgdb is stored.
- `GeneAnnoatationGFF3`: Gene Annotation GFF3 Filename.
- `DNAFastaFilename`: Fasta file of the HARDMASKED genome nucleotide sequences.
- `MetabolicDomainsToAnnotate`: A list of Metabolic domains that should be analyzed. For example {'Amines and Polyamines Metabolism'; 'Amino Acids '}
- `gout`: Path and file name to the gene-Outputfile, make sure you have access to the file and folder. This is the first results output file that will be generated.
- `cout`: Path and file name to the cluster-Outputfile, make sure you have access to the file and folder. This is the second Results output file that will be generated.
## Input Parsing Arguments
These arguments specifically deal with the parsing and processing of input.
### GFF3 Parsing Arguments
When parsing the GFF3 file, certain attributes may need to be overridden to match the data in PGDB.
If the GFF3 file structure is unclear, you may need to view it directly to determine how to set these arguments.

- `GFF3GeneAttribute`: Attribute in GFF3 file to extract for Gene Name Typically "Name", but possibly "locus_tag". (Default: Name)
- `GFF3ProteinAttribute`: Attribute in GFF3 file to extract for transcript ID. Typically "Name", but possibly "original_protein_id". (Default: Name)
- `GFF3ProteinRegex`: Regex to find the protein ID shared with PGDB data. If the column contains jnl|GCVI|AT1G5560.1 and PGDB protein ID is AT1G5560.1, you can pass the regex "[^|]+$" to select everything after the last '|' (Default: ".*")

### Compatability Options
By default, PlantClusterFinder.py uses Diamond for blastp-compatible alignment, and has default metabolic domains and signature/tailoring files it pulls from.
However, if you need to override those defaults, these are the options to use.

- `use_ncbi_blast`: Use NCBI Blast rather than Diamond for BlastP (Default False)
- `ReactionMetabolicDomainFilename`: TSV file of the metabolic domains of reactions. Check if file is outdated (needs to contain all reactions that are present in the PGDB that have genes to be annotated). Defaults to $HOME/.plantclusterfinder/ReactionMetabolicDomainClassification.txt
- `SignatureTailoringFilename`: Signature and tailoring reactions TSV file describing the classes of enzymes that should be classified as signature or tailoring enzymes to identify clusters containing such.
#### GFF3 Parsing Override
If you have the protein sequence fasta, a gene location TSV, and the Gene Transcript Protein file, it is possible to override the preprocessing step.
!!! important
    All three arguments must be set to override GFF3 parsing.

- `psf`: [DEPRECATED] Protein sequence file of the genome. Can only contain Protein IDs as header e.g. >protein1-ID This is auto-generated from Genome FASTA and GFF3 file.
- `GeneLocationTSV`: [DEPRECATED] Gene Location (position) TSV file with gene ID, start bp, end bp, chromosome / scaffold, strand (encoded as 1 or -1). This file can be generated with a biomart or out of a gff3 file.
- `GeneTranscriptProteinFilename`: [DEPRECATED] TSV file with a header describing gene, transcript, and protein identifier and then for all the identfiers the listings (gene-ID, transcript-ID, protein-ID). This is used to map information regarding transcripts (currently not used) and proteins (used for MCL clustering, potentially also used in PGDBs) to the genes (e.g. the gene location information).

## Biosynthetic Gene Cluster Search Parameters
The following arguments configure the parameters of searching for, and predicting biosynthetic gene clusters.
The defaults are suitable for most applications, but they may need to be adjusted depending on the genome.

- `MaxSeqGapSize`: (*number*) Maximal masked nucleotide region (in bp) that a cluster is allowed to bridge (default is 100000).
- `MaxInterGeneDist`: (*number*) Maximal intergenic distance (in bp) that still can be crossed by a gene cluster. (default inactive).
- `MaxInterGeneDistByMedian`: (*number*) Maximal intergenic distance that still can be crossed by a gene cluster, here defined by this number times the median of the gene sizes. (default is 50. Set it to -1 if it should be inactive).
- `PercentileForMaxInterGeneDist`: (*number*) Calculates the maximal intergenic distance that still can be crossed by a gene cluster. This number determines the percentile to choose from all the distances. If this is 99.9, then a cluster is allowed to bridge if the largest intergenic distance in the cluster is below 99.9% of the background size of all intergenic distances of the genome (default is -1,thus inactive).
- `SeqGapSizesChromBreak`: (*number*) The masked nucleotides sequences of these specific lengths will prevent cluster formation (Default is an empty list: []).
- `EnzMinForClusters`: (*number*) Defines the minimum number of enyzmes that need to be present in a cluster (default is 3).
- `KbRangeLD`: (*number*) Genes that are separated by more than this size (in bp) are not considered as local duplicates. (default is 100000)
- `GeneRangeLD`: (*number*) Number of intervening genes that are allowed to separate two local duplicated genes. Two genes that are separated by more than this many intervening genes are not called local duplicates. (default is 10)
- `TopPercentClusters`: (*number*) Clusters are labeled if they are within this top X% of enzyme dense clusters compared to the background. (default is 5)
- `MinStepSize`: (*number*) Minimal intervening gene size that should be computed (default is 5).
- `MaxStepSize`: (*number*) Maximal intervening gene size that should be computed (default is 20).
- `Criterion`: (*number*) Criterion after which stepsize (intervening gene size) is chosen
    1. as published in Schlapfer et al 2017 (PMID:28228535)'' sum(enzymes in clusters) - sum(non-enzymes in clusters) > 0.
    2. **(default)** stepsize 3 (because most organisms take this one)
    3. sum(enzymes in clusters) - sum(non-enzymes in clusters, without hypo genes) -2*sum(hypotetical genes)
    4. sum(clusters in stepsize n) > sum(clusters in stepsize n+1).

### Parameters Regarding the Insertion of Hypothetical Genes:
In the case of sequence gaps, hypothetical genes may be introduced to act as a buffer for BGC prediction.
These parameters should be configured based heavily on the nature of the genome being studied.

- `HypoGenePercentile`: (*number*) Defines the length of sequencing gap that should be populateded by two hypothetical genes. If this is 10, then two hypothetical genes are introduced as soon as the sequencing gap exceeds the (100-10) lower percentile of the background size of all genes of the genome (default is 5).
- `InsertHypos`: Defines if the genome should be populated by hypothetical genes in sequencing gaps that are not physical breaks. (Default False)
- `HypoAmount`: (*number*) Defines by how many hypothetical genes a sequencing gap that is not a physical break should be populated by artificial genes to punish such sequencing gaps from becoming clusters. (Default 0)

## Other Optional Parameters:
- `OutputType`: Defines the format of the Outputfiles: 'old' is using the old format to report all clusters of all stepsizes, 'verbose' defines that top percent clusters of all stepsizes are reported. 'simple' just reports the top percent clusters of the chosen stepsize (default simple).
- `OutputGeneLocationParquet`: The file location of the parquet file to store the GeneLocation data. The file is only written if the argument is set.
- `RunInfoFile`: Store information (arguments, settings, timing) for a run. (Default: STDOUT)
- `LogFile`: Store log output to a file. (Default: STDERR)
  - Note: the logs of external programs (MCL, Diamond) are always sent to STDERR
- `GenesDat_GeneInfo`: Attributes in PGDB genes.dat where gene ID information should be searched for. (Default is "ACCESSION-1, COMMON-NAME")
- `ComputedMCLClusterFile`: Indicates that the input file is not a protein fasta file, but a precomputed MCL clustering file (precomputing can be useful for speed). This file can only contain Gene Names, tabs, and newline characters.
- `OverwriteSeqGapInfo`: Force masked nucleotide analysis for finding sequencing gaps (Default False)
- `OverwriteMCLClustering`: Forces (if set) that the MCL clustering is redone even if the result files are already present (Default False).
- `Verbose`: Verbosity Level
    - **0 (default)** only logs warnings and errors
    - **1** logs info-level information
    - **2** logs debug information
- `para`: How many cpus can be used by the algorithm (used for MCL clustering and Sequencing gap search). The default is 0 (auto-scale).
- `UnmaskedDNA`: Defines if PlantClusterFinder should continue, despite a non-masked genome sequence file being used. (Default False)
- `PGDBIdsToMap`: Can contain the Letters G and/or T and/or P. In that case (the pgdb is mapped to Gene (G) and/or Transcript (T) and/or Protein (P) Ids of the Gene conversion file. The default value is 'G'.
- `RemoveNonProtLocation`: Defines if the protein-fasta file should be used to clean out the gene position file from sequences that do not show up on the protein fasta file, and thus do not have had a chance to be predicted to be an enzyme. (Default False)
## Test Arguments
!!! important
    This argument overrides all others.
    Do not try to use it in conjunction with other arguments.
Use the argument `--test-install` to test an execution against Csubellipsoidea data.
