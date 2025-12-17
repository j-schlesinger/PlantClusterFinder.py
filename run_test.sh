#!/usr/bin/env bash
rm ./csubellipsoidea/csubellipsoidea_Clust_v1_3.txt
rm ./csubellipsoidea/csubellipsoidea_Gene_v1_3.txt
rye run python -m plantclusterfinder_py -pgdb './csubellipsoidea/pgdb/csubellipsoideacyc/1.0/data/'\
   '--GenesDat_GeneInfo' 'ACCESSION-1'\
   '-rmdf' './Inputs/ReactionMetabolicDomainClassification.txt'\
   '-md'\
   "Amines and Polyamines Metabolism, Amino Acids Metabolism, Carbohydrates Metabolism, Cofactors Metabolism, Detoxification Metabolism, Energy Metabolism, Fatty Acids and Lipids Metabolism, Hormones Metabolism, Inorganic Nutrients Metabolism, Nitrogen-Containing Compounds, Nucleotides Metabolism, Phenylpropanoid Derivatives, Polyketides, Primary-Specialized Interface Metabolism, Redox Metabolism, Specialized Metabolism, Sugar Derivatives, Terpenoids"\
   '-psf' './csubellipsoidea/CsubellipsoideaC_169_227_v2.0.protein.pcf13.fa'\
   '-gtpf' './csubellipsoidea/gtpf_CsubellipsoideaC_169_227_v2.0.annotation_info.txt.txt'\
   '-glof' './csubellipsoidea/glof_CsubellipsoideaC_169_227_v2.0.gene.gff3.txt'\
   '-dnaf' './csubellipsoidea/CsubellipsoideaC_169_227_v2.0.hardmasked.fa'\
   '-sitf' './Inputs/scaffold-tailoring-reactions-05082016.tab'\
   '-gout' './csubellipsoidea/csubellipsoidea_Gene_v1_3.txt'\
   '-cout' './csubellipsoidea/csubellipsoidea_Clust_v1_3.txt'\
   '--SeqGapSizesChromBreak' 10000 '--PGDBIdsToMap' 'GTP'\
   --Verbose 2
