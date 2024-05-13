pre-print available: https://doi.org/10.1101/2022.10.21.513146

## Jacob S. Berv<sup>1,2,3,*</sup>, Sonal Singhal<sup>4</sup>, Daniel J. Field<sup>5,6</sup>, Nathanael Walker-Hale<sup>7</sup>, Sean W. McHugh<sup>8</sup>, J. Ryan Shipley<sup>9</sup>, Eliot T. Miller<sup>10</sup>, Rebecca T. Kimball<sup>11</sup>, Edward L. Braun<sup>11</sup>, Alex Dornburg<sup>12</sup>, C. Tomomi Parins-Fukuchi<sup>13</sup>, Richard O. Prum<sup>14,15</sup>, Benjamin M. Winger<sup>1,3</sup>, Matt Friedman<sup>2, 16</sup>, Stephen A. Smith<sup>1</sup>

Corresponding Author: *Jacob S. Berv, jberv@umich.edu, jacob.berv@gmail.com

Author affiliations:

1.	Department of Ecology and Evolutionary Biology, 1105 North University Avenue, Biological Sciences Building, University of Michigan, Ann Arbor, Michigan, 48109-1085, USA
2.	University of Michigan Museum of Paleontology, 1105 North University Avenue, Biological Sciences Building, University of Michigan, Ann Arbor, Michigan, 48109-1085, USA
3.	University of Michigan Museum of Zoology, 1105 North University Avenue, Biological Sciences Building, University of Michigan, Ann Arbor, Michigan, 48109-1085, USA
4.	Department of Biology, California State University, Dominguez Hills, Carson, California 90747, USA
5.	Department of Earth Sciences, Downing Street, University of Cambridge, Cambridge CB2 3EQ, UK
6.	Museum of Zoology, University of Cambridge, Downing Street, Cambridge CB2 3EJ, UK
7.	Department of Plant Sciences, University of Cambridge, Downing Street, Cambridge, CB2 3EA, UK
8.	Department of Evolution, Ecology, and Population Biology, Washington University in St Louis, St Louis, Missouri, USA
9.	Department of Forest Dynamics Swiss Federal Institute for Forest, Snow, and Landscape Research WSL Zürcherstrasse 111 8903 Birmensdorf, Switzerland
10.	Macaulay Library, Cornell Lab of Ornithology, Ithaca, New York, 14850, USA
11.	Department of Biology, University of Florida, Gainesville, Florida 32611, USA
12.	Department of Bioinformatics and Genomics, University of North Carolina at Charlotte, Charlotte, North Carolina, USA
13.	Department of Ecology and Evolutionary Biology, University of Toronto, Toronto, Ontario, Canada, M5S 3B2
14.	Department of Ecology and Evolutionary Biology, Yale University, New Haven, Connecticut, 06520, USA
15.	Peabody Museum of Natural History, Yale University, New Haven, Connecticut, 06520, USA
16.	Department of Earth and Environmental Sciences, 1100 North University Avenue, University of Michigan, Ann Arbor, Michigan, 48109-1085, USA


## Code and datasets

---

**constrained_species_tree.R**

This is a large, complex R script containing the code for the primary analyses in our manuscript. It is separated into 17 sub-sections delimited using curly-bracket notation {}. Each section is commented and annotated for reproducibility. This script is mostly for processing and analyzing the *output files* from *Janus* and other software applied herein.

Here we provide a brief overview of the various sub-sections of the script for reader reference.

### System report and package references

Here we provide a overview of the R package versions and system information used when the analyses were performed.

### Section 1 - Setting up 1

This section sets up the appropriate file paths for importing data and performs initial data pre-processing.

### Section 2 - Setting up 2

This section performs additional data pre-processing, including reading the output data from Janus, reading in particular tree output files, and performing some pre-processing related to quantifying GC content.

### Section 3 - Generating new datasets for logistic regression

This section sets up input data objects for phylogenetic logistic regression, as described in the manuscript.

### Section 4 - Running phylogenetic logistic regression 

This section runs various logistic regression models, summarizes the model fits, and generates plots that are provided in the supplementary material (Figure S1).

### Section 5 - Get data for GC content and codon usage

This section imports the sequence datasets and estimates various relevant parameters, including aspects of GC content across various slices of the dataset and other nucleotide statistics. The output is an intermetdiate data frame summarizing the results for downstream processing.

### Section 6 - Estimating phylogenetic signal

This section generates estimates of phylogenetic signal across various nucleotide statistics. Largely exploratory.

### Section 7 - Codon usage calculations

This section includes analyses of patterns of codon usage, as described in the manuscript. We estimate different metrics of codon usage and then compare them across groups identified by Janus, using phylogenetic anove (see manuscript text).

### Section 8 - Life history trait data setup

This section sets up the life history character datasets and does initial pre-processing (eg. data imputation described in the manuscript).

### Section 9 - Setting up input objects for mvMORPH/OUwie analyses

Setting up character and tree data sets in a format that can be read by the R packages for fitting models of character evolution.

### Section 10 - OUwie model fitting

This section fits various models (eg, BM, OU) to life history trait data, as described in the manuscript. These analyses contribute to the graphics depicted in Figure 2 of the manuscript.

### Section 11 - Setting up plots of paramter estimates

This section sets up the data objects for generating the radar plots provided in the supplemental text (Figure S7a-d). These plots depict the parameter estimates from Janus (e.g., estimates of equilibrium base frequencies).

### Section 12 - Setting up BMR datasets

This section does pre-processing on the metabolic scaling data sets. The analyses of these data are described in [another script](./BMR) in the BMR directory.

### Section 13 - Setting up datasets for the random forest classifier

This section does does some processing on the cleaned up data sets, to pass to the random forest classifier (tidymodels) described in [another script](./RandomForest_var_imp.R)

### Section 14 - Analyses with l1ou

This section sets up and performs analyses using the [l1ou R package](https://github.com/khabbazian/l1ou), which fits a pseudo-multivariate OU model with shifting optimum. See the supplementary text in the manuscript for methods details.

### Section 15 - Processing l1ou results

This section takes the models fit with l1ou and generates output products. These results contribute to Figure 1 and Figure S4a-b.

### Section 16 - Additional BMR models

This section tests some additional exploratory models of BMR allometry. Most of these are not discussed in the manuscript and are provided as exploratory for interested readers.

### Section 17 - Miscellaneous

This section includes code related to diagnostic plots, as well as some additional experimental analyses relative to multivariate model fitting. These are mostly not discussed in the manuscript and are provided as exploratory for interested readers.

--- 

Most lines in the primary script where analyses are actually executed are commented out so that processed RDS objects ([provided here](./RDS), see below) can be loaded instead. RDS (R Data Serialization) files are a common format for saving R objects, and they allow R users to preserve the state of an object between R sessions. Thus, the script is provided in a format that should allow anyone to download the GitHub repository and go through every sub-section without needing to re-fit or re-process any data locally.

**Functions_consensus_genetrees.R**

This file contains many custom or new functions called by the primary analysis script. This functions file is sourced by the primary analysis script when executed in order.

**RandomForest_var_imp.R**

This script contains the code for RandomForest supervised classification using the tidymodels framework. This code is adapted from the Physalia course examples for machine learning in R.

**shift_model_sims.R**

This script contains the code for assessing the statistical performance of Janus; including assessing false positive/negative rates for substitution model shifts. It include various R functions defining a comprehensive pipeline that, in combination with IQ-tree, can be used to assess the performance of Janus on arbitrary topologies. See annotations within this script to get a sense of how to run it for your own topologies or contact the authors for additional detail.

**Supplementary Table 3.xlsx**

Supplementary Table 3, as noted in the manuscript. This file describes assembly details for assembly and extraction of mtDNA genome data from the original AHE target capture data.

**LHT_DATA**

LHT_reference.xlsx contains a database of life-history traits analyzed in the present study, which are collated from a varietyof sources (as described in the manuscript). This table contains only raw data, and does not include estimated values (See primary R script and RDS files).

**RDS**

This folder contains RDS files (R data objects). These files correspond to intermediate or final R data objects generated in the course of analysis and which are called/loaded by the [primary analysis script](/constrained_species_tree.R)

**AHE_REASSEMBLY**

This folder contains files related to assembly of the AHE nuclear DNA dataset. See the detailed README in this directory for additional detail.

**mtDNA_REASSEMBLY**

This folder contains files related to assembly of the mtDNA dataset. See the detailed README in this directory for additional detail.

**BMR**

This folder contains scripts and files related to analysis of the BMR dataset. See the detailed readme in this directory for additional detail.

**Janus**

This folder contains the relevant results from janus (e.g.,input/output files). The source tree file is in ‘timetree_all_taxa_OW_2019.nextree.tre’ available in the supplementary materials from Kimball et al. (2019) here: https://www.mdpi.com/1424-2818/11/7/109/s1, labeled as ‘MRL_3backbone.’ It is also provided in our GitHub repository [here](/trees/MRL_3backbone.tre) for reference. Janus has been implemented in both Golang and C, and the source code (and documentation) is available at https://git.sr.ht/~hms/janus and https://git.sr.ht/~hms/hringhorni.
