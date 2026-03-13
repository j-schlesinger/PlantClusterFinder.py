#!/usr/bin/env bash
mkdir -p Outputs
rm ./Outputs/csubellipsoidea_Clust_v2_0.txt 2>/dev/null
rm ./Outputs/csubellipsoidea_Gene_v2_0.txt 2>/dev/null
#   '-psf'  './csubellipsoidea/CsubellipsoideaC_169_227_v2.0.protein.pcf13.fa'\
#   '-glof' './csubellipsoidea/glof_CsubellipsoideaC_169_227_v2.0.gene.gff3.txt'\
#   '-gtpf' './csubellipsoidea/gtpf_CsubellipsoideaC_169_227_v2.0.annotation_info.txt.txt'\
uv run python -m plantclusterfinder_py\
   '-r' 'Csubellipsoidea_Test_Run'\
   '--GenesDat_GeneInfo' 'ACCESSION-1'\
   '-md'\
   "Amines and Polyamines Metabolism" "Amino Acids Metabolism"\
   "Carbohydrates Metabolism" "Cofactors Metabolism"\
   "Detoxification Metabolism" "Energy Metabolism"\
   "Fatty Acids and Lipids Metabolism" "Hormones Metabolism"\
   "Inorganic Nutrients Metabolism" "Nitrogen-Containing Compounds"\
   "Nucleotides Metabolism" "Phenylpropanoid Derivatives"\
   "Polyketides" "Primary-Specialized Interface Metabolism"\
   "Redox Metabolism" "Specialized Metabolism" "Sugar Derivatives"\
   "Terpenoids"\
   '-dnaf' './csubellipsoidea/CsubellipsoideaC_169_227_v2.0.hardmasked.fa'\
   '-gff3' './csubellipsoidea/CsubellipsoideaC_169_227_v2.0.gene.gff3'\
   '-pgdb' './csubellipsoidea/pgdb/csubellipsoideacyc/1.0/data/'\
   '-sitf' './Inputs/scaffold-tailoring-reactions-05082016.tab'\
   '-rmdf' './Inputs/ReactionMetabolicDomainClassification.txt'\
   '-gout' './Outputs/csubellipsoidea_Gene_v2_0.txt'\
   '-cout' './Outputs/csubellipsoidea_Clust_v2_0.txt'\
   '--SeqGapSizesChromBreak' 10000 '--PGDBIdsToMap' 'GTP'\
   '--OutputType' 'simple'\
   --Verbose 2\
   "--OutputGeneLocationParquet" "./Outputs/csubellipsoidea_gene_loc.parquet"\
   --NoRerunBlastP\
   --NoRebuildBlastDB
# the last 2 options should be removed for *real* testing
