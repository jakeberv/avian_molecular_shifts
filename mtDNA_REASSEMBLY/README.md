
## Assembly of off-target mtDNA from Prum et al. 2015

---

### Processed data files

**unaligned.zip**

* This zip file contains unaligned sequences generated with MitoFinder

**alignments.zip**

* Zip file containing sequences aligned with MACSE/FSA (described in the manuscript).

**concat_all_paup.fas**

* All mtDNA loci in a concatenated fasta file

**concat_proteins_no_nd6.fas**

* mtDNA proteins in a concatenated fasta file

**concat_rRNAs.fas**

* mtDNA rRNAs in a concatenated fasta file

---

### Data pre-processing

**mtDNA-processing.R**

* R script that assembles the unaligned sequences into fasta files (calls MITOGENOME-Berv-Aug-2021.xlsx, below)

**MITOGENOME-Berv-Aug-2021.xlsx**

* Assembly spreadsheet (only called by internal R scripts and not intended for reader use).

---

### Summary of mtDNA reassembly (from methods section in our manuscript)

* Using the [reassembled contigs](../AHE_REASSEMBLY/bioinformatics_pipeline), we ran [Mitofinder 1.4](https://github.com/RemiAllio/MitoFinder). For reference mitogenomes, we used complete mitogenomes available in GenBank [Table S3](../Supplementary%20Table%203.xlsx). When available, we used a reference from the same order (though for Passeriformes, we used different references for oscines [Passeri] and suboscines [Tyranni]); in a few cases, it was necessary to use a reference from a closely related order [Table S3](../Supplementary%20Table%203.xlsx). We then extracted the 13 protein-coding genes and 2 rRNAs from the mitofinder output. In some cases, limited mitochondrial data were recovered (Table S3). In those cases, we searched GenBank for the same or a phylogenetically equivalent species that could be substituted. When no suitable alternative was available from GenBank, we also used mitogenomes assembled from the raw data collected as part of other studies (referenced in the manuscript and Table S3). To increase data coverage in five cases, we generated chimeric sequences using available GenBank data from multiple individuals of the same species [Table S3](../Supplementary%20Table%203.xlsx). Final alignments were used as input into [Janus](../janus)

**reannotate.pl**

* This is a perl script developed by Ed Braun (ebraun68@ufl.edu) and Rebecca Kimball (rkimball@ufl.edu) for running MitoFinder as summarized above and in our manuscript. Instructions on how to use it are included as commented annotations.

* This script generates a tab-delimited output file which is converted to fasta format using awk. Our preliminary alignment used to vet the sequences was produced with a simple mafft command:

`mafft --adjustdirection xxx.fasta > xxx.aln.fasta`

From there we concatenated the alignments of individual regions using a tool here
`https://github.com/ebraun68/RYcode/blob/main/simple_concat.pl`

* After preliminary data checking and assembly, we re-aligned using the more sophisticated approaches described in the manuscript (e.g., with the FSA/MACSE alignment software).

* As noted above, due to data quality variation, the final data set we used is a combination of off-target mitochondrial data from Prum et al. 2015 (most of the data), as well as data from other studies. The specific Trinity assemblies (derived from Prum et al. 2015) that contributed to our final mtDNA data set are noted in Supplementary Table 3 (noted as specific *.fasta files in the leftmost column, [available here](../AHE_REASSEMBLY/)).