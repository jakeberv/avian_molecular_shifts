Release v1.0.0 has been archived at Zenodo

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.12513738.svg)](https://doi.org/10.5281/zenodo.12513738)

# Supplemental Data Repository

## Genome and life-history evolution link bird diversification to the end-Cretaceous mass extinction

*in press* at Science Advances

The pre-print for our research is available at [this link](https://doi.org/10.1101/2022.10.21.513146).

---

## Authors

### Lead and Corresponding author

-   **Jacob S. Berv**<sup>1,2,3,</sup>
    -   Email: [jberv\@umich.edu](mailto:jberv@umich.edu)

### Co-authors

-   Sonal Singhal<sup>4</sup>
-   Daniel J. Field<sup>5,6</sup>
-   Nathanael Walker-Hale<sup>7</sup>
-   Sean W. McHugh<sup>8</sup>
-   J. Ryan Shipley<sup>9</sup>
-   Eliot T. Miller<sup>10</sup>
-   Rebecca T. Kimball<sup>11</sup>
-   Edward L. Braun<sup>11</sup>
-   Alex Dornburg<sup>12</sup>
-   C. Tomomi Parins-Fukuchi<sup>13</sup>
-   Richard O. Prum<sup>14,15</sup>
-   Benjamin M. Winger<sup>1,3</sup>
-   Matt Friedman<sup>2,16</sup>
-   Stephen A. Smith<sup>1</sup>

### Author Affiliations

1.  Department of Ecology and Evolutionary Biology, University of Michigan, Ann Arbor, MI, USA
2.  University of Michigan Museum of Paleontology, Ann Arbor, MI, USA
3.  University of Michigan Museum of Zoology, Ann Arbor, MI, USA
4.  Department of Biology, California State University, Dominguez Hills, Carson, CA, USA
5.  Department of Earth Sciences, University of Cambridge, Cambridge, UK
6.  Museum of Zoology, University of Cambridge, Cambridge, UK
7.  Department of Plant Sciences, University of Cambridge, Cambridge, UK
8.  Department of Evolution, Ecology, and Population Biology, Washington University in St Louis, St Louis, MO, USA
9.  Department of Forest Dynamics, Swiss Federal Institute for Forest, Snow, and Landscape Research WSL, Birmensdorf, Switzerland
10. Macaulay Library, Cornell Lab of Ornithology, Ithaca, NY, USA
11. Department of Biology, University of Florida, Gainesville, FL, USA
12. Department of Bioinformatics and Genomics, University of North Carolina at Charlotte, Charlotte, NC, USA
13. Department of Ecology and Evolutionary Biology, University of Toronto, Toronto, Ontario, Canada
14. Department of Ecology and Evolutionary Biology, Yale University, New Haven, CT, USA
15. Peabody Museum of Natural History, Yale University, New Haven, CT, USA
16. Department of Earth and Environmental Sciences, University of Michigan, Ann Arbor, MI, USA

---

## Code and Datasets

### Primary Analysis Script Overview

[**constrained_species_tree.R**](./constrained_species_tree.R)

This complex R script is central to our manuscript's analysis. It comprises 17 detailed sub-sections, each delimited by curly brackets `{}` for clear demarcation. Annotations and comments within each section ensure the script is comprehensible and reproducible. The primary focus of this script is to process and analyze output files from *Janus* and other utilized software.

#### Sub-section Descriptions

-   **System report and package references**: Documents the R package versions and system configurations used during the analysis.

-   **Section 1 - Setting up 1**: Establishes file paths and conducts initial data pre-processing.

-   **Section 2 - Setting up 2**: Continues data pre-processing, includes reading output from Janus, tree output files, and pre-processing for GC content quantification.

-   **Section 3 - Generating new datasets for logistic regression**: Prepares input data objects for phylogenetic logistic regression as outlined in the manuscript.

-   **Section 4 - Running phylogenetic logistic regression**: Executes logistic regression models, summarizes the fits, and produces plots for the supplementary materials (Figure S1).

-   **Section 5 - Get data for GC content and codon usage**: Imports sequence datasets and estimates parameters such as GC content and other nucleotide statistics, summarizing results in an intermediate data frame.

-   **Section 6 - Estimating phylogenetic signal**: Generates estimates of phylogenetic signal for various nucleotide statistics in an exploratory fashion.

-   **Section 7 - Codon usage calculations**: Analyzes codon usage patterns, estimating different metrics and comparing them across groups identified by Janus using phylogenetic ANOVA.

-   **Section 8 - Life history trait data setup**: Prepares life history trait datasets and performs initial preprocessing, including data imputation as described in the manuscript.

-   **Section 9 - Setting up input objects for mvMORPH/OUwie analyses**: Prepares datasets for analysis with R packages that model character evolution.

-   **Section 10 - OUwie model fitting**: Fits various models (e.g., BM, OU) to life history trait data, contributing to the visuals in Figure 2 of the manuscript.

-   **Section 11 - Setting up plots of parameter estimates**: Arranges data objects for generating radar plots (Figure S7a-d) that depict parameter estimates from Janus.

-   **Section 12 - Setting up BMR datasets**: Processes data for metabolic rate analysis, detailed further in [another script](./BMR) located in the BMR directory.

-   **Section 13 - Setting up datasets for the random forest classifier**: Processes cleaned data sets for use with the random forest classifier, detailed in [another script](./RandomForest_var_imp.R).

-   **Section 14 - Analyses with l1ou**: Utilizes the [l1ou R package](https://github.com/khabbazian/l1ou) to perform analyses with a pseudo-multivariate OU model with shifting optima, as detailed in the supplementary manuscript text.

-   **Section 15 - Processing l1ou results**: Processes output from l1ou models, contributing to Figure 1 and Figure S4a-b.

-   **Section 16 - Additional BMR models**: Explores additional, exploratory models of basal metabolic rate (BMR) allometry, not extensively discussed in the manuscript.

-   **Section 17 - Miscellaneous**: Includes code for diagnostic plots and additional experimental analyses concerning multivariate model fitting. These sections are provided as exploratory content for interested readers.

*Note: To facilitate seamless reproducibility, most execution lines in the script are commented out, allowing the script to load processed RDS objects directly from [here](./RDS). RDS files are a standard format for storing R objects, enabling the preservation of data states across sessions.*

---

### Additional Resources and Directories

This section outlines the other crucial files and directories included in our repository that support our main analysis or provide further insights and data.

#### Essential Scripts and Files

-   [**Functions_consensus_genetrees.R**](./Functions_consensus_genetrees.R): Contains custom and new functions utilized by the primary analysis script. This file is sourced in sequence during the script's execution.

-   [**RandomForest_var_imp.R**](./RandomForest_var_imp.R): Implements RandomForest supervised classification using the tidymodels framework. This implementation is adapted from machine learning examples provided by the Physalia course.

-   [**shift_model_sims.R**](./shift_model_sims.R): Features code for evaluating the statistical performance of Janus, specifically addressing false positive and negative rates for substitution model shifts. The script includes a variety of R functions that establish a pipeline with IQ-tree for assessing Janus's performance across different topologies. Detailed annotations within the script provide additional insights.

-   [**Supplementary Table 3.xlsx**](./Supplementary%20Table%203.xlsx): Contains assembly details for the mtDNA genome data extracted from the original AHE target capture data, as referenced in our manuscript.

#### Important Directories

-   [**LHT_DATA**](./LHT_DATA): This directory houses 'LHT_reference.xlsx', a database of life-history traits examined in our study. These traits are compiled from various sources, with the database containing only raw data, devoid of estimated values (refer to the primary R script and RDS files for more details).

-   [**RDS**](./RDS): Includes RDS files, which are R data objects that represent either intermediate or final data states generated during our analyses. These objects are crucial for the operations carried out by the [primary analysis script](/constrained_species_tree.R).

-   [**AHE_REASSEMBLY**](./AHE_REASSEMBLY): Contains files pertinent to the assembly of the AHE nuclear DNA dataset. A detailed README within this directory provides extensive information about the processes and methodologies employed.

-   [**mtDNA_REASSEMBLY**](./mtDNA_REASSEMBLY): Hosts files necessary for the assembly of the mtDNA dataset. Detailed documentation in this directory elaborates on the specific assembly steps and procedures.

-   [**BMR**](./BMR): This directory includes scripts and files related to the analysis of the Basal Metabolic Rate (BMR) dataset. A comprehensive readme file in this directory offers detailed insights into the analysis techniques and results.

-   [**Janus**](./janus): Contains all relevant results from Janus, including input and output files. The directory also includes the ‘timetree_all_taxa_OW_2019.nextree.tre’ from the supplementary materials of Kimball et al. (2019), accessible [here](https://www.mdpi.com/1424-2818/11/7/109/s1) and in our repository [here](/trees/MRL_3backbone.tre). Janus is developed in both Golang and C, with its source code and documentation available at <https://git.sr.ht/~hms/janus> and <https://git.sr.ht/~hms/hringhorni>.
