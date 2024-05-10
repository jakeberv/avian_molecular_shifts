
## Files related to assembly of off-target mtDNA. 

### Description of mtDNA reassembly

Using the [reassembled contigs](../AHE_REASSEMBLY), we ran [Mitofinder 1.4](https://github.com/RemiAllio/MitoFinder). For reference mitogenomes, we used complete mitogenomes available in GenBank [Table S3](../Supplementary%20Table%203.xlsx). When available, we used a reference from the same order (though for Passeriformes, we used different references for oscines [Passeri] and suboscines [Tyranni]); in a few cases, it was necessary to use a reference from a closely related order [Table S3](../Supplementary%20Table%203.xlsx). We then extracted the 13 protein-coding genes and 2 rRNAs from the mitofinder output. In some cases, limited mitochondrial data were recovered (Table S3). In those cases, we searched GenBank for the same or a phylogenetically equivalent species that could be substituted. When no suitable alternative was available from GenBank, we also used mitogenomes assembled from the raw data collected as part of other studies (referenced in the manuscript). To increase data coverage in five cases, we generated chimeric sequences using available GenBank data from multiple individuals of the same species [Table S3](../Supplementary%20Table%203.xlsx). These sequence data are used as input into [Janus](../janus)

### Processed data files

**unaligned.zip**

This zip file contains unaligned sequences generated with MitoFinder

**alignments.zip**

Zip file containing sequences aligned with MACSE/FSA (described in the manuscript).

**concat_all_paup.fas**

All mtDNA loci in a concatenated fasta file

**concat_proteins_no_nd6.fas**

mtDNA proteins in a concatenated fasta file

**concat_rRNAs.fas**

mtDNA rRNAs in a concatenated fasta file

### Data pre-processing

**mtDNA-processing.R**

R script that assembles the unaligned sequences into fasta files (calls MITOGENOME-Berv-Aug-2021.xlsx, below)

**MITOGENOME-Berv-Aug-2021.xlsx**

Assembly spreadsheet (only called by internal R scripts and not intended for reader use).

### Data assembly