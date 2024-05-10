## janus

Janus has been implemented in both Golang and C, and the source code (and documentation) is available at https://git.sr.ht/~hms/janus and https://git.sr.ht/~hms/hringhorni

**files**

This folder contains five subdirectories containing janus input/output files. The output files are imported by the primary R script ([constrained_species_tree.R](../constrained_species_tree.R)) for post-processing. The sequence data files that these runs used are located and described [here](../AHE_REASSEMBLY/sequences).

**/files/exons_MFP_MERGE_MRL3_constraint_G_UE_UL_M3**

Contains input/output files of the final analysis run for exon data. Inside this directory the *.call.txt file contains the specific command used to run the janus software.

**/files/introns_MFP_MERGE_MRL3_constraint_G_UE_UL_M3**

Contains input/output files of the final analysis run for intron data. Inside this directory the *.call.txt file contains the specific command used to run the janus software.

**/files/utrs_MFP_MERGE_MRL3_constraint_G_UE_UL_M3**

Contains input/output files of the final analysis run for UTR data. Inside this directory the *.call.txt file contains the specific command used to run the janus software.

**/files/ALL_MFP_MERGE_MRL3_constraint_G_UE_UL_M3**

Contains input/output files of an analysis run that considered the entire concatenated datset. Discussion of these results were excluded from the manuscript because we decided to focus on data type. Inside this directory the *.call.txt file contains the specific command used to run the janus software.

**/files/mtdnas**

Contains input/output files for three separate analyses (in three sub-directories) of mtDNA sequences. The 'concat_all' directory contains analyses of all mtDNA sequences concatenated together. The 'concat_proteins' directory contains analyses of all mtDNA protein sequences. Lastly, the 'concat_rRNA' directory contains analyses of mtDNA rRNA sequences. Inside each sub-directory the *.call.txt file contains the specific command used to run the janus software.

