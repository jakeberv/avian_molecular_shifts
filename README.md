pre-print available: https://doi.org/10.1101/2022.10.21.513146

## Jacob S. Berv<sup>1,2,3,*</sup>, Sonal Singhal<sup>4</sup>, Daniel J. Field<sup>5,6</sup>, Nathanael Walker-Hale<sup>7</sup>, Sean W. McHugh<sup>8</sup>, J. Ryan Shipley<sup>9</sup>, Eliot T. Miller<sup>10</sup>, Rebecca T. Kimball<sup>11</sup>, Edward L. Braun<sup>11</sup>, Alex Dornburg<sup>12</sup>, C. Tomomi Parins-Fukuchi<sup>13</sup>, Richard O. Prum<sup>14,15</sup>, Benjamin M. Winger<sup>1,3</sup>, Matt Friedman<sup>2, 16</sup>, Stephen A. Smith<sup>1</sup>

Corresponding Author: *Jacob S. Berv, jberv@umich.edu, jacob.berv@gmail.com

Author affiliations:

1.	Department of Ecology and Evolutionary Biology, 1105 North University Avenue, Biological Sciences Building, University of Michigan, Ann Arbor, Michigan, 48109-1085, USA

2.	University of Michigan Museum of Paleontology, 1105 North University Avenue, Biological Sciences Building, University of Michigan, Ann Arbor, Michigan, 48109-1085, USA

3.	University of Michigan Museum of Zoology, 1105 North University Avenue, Biological Sciences Building, University of Michigan, Ann Arbor, Michigan, 48109-1085, USA

4.	Department of Biology, California State University, Dominguez Hills, Carson, California 90747, USA

5.	Department of Earth Sciences, Downing Street, University of Cambridge, Cambridge CB2 3EQ, UK

6.	Museum of Zoology, Downing Street, University of Cambridge, Cambridge CB2 3EJ, UK

7.	Department of Plant Sciences, Downing Street, University of Cambridge, Cambridge, CB2 3EA, UK

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

**constrained_species_tree.R**

This is a large R script containing the code for the primary 
series of analyses. It is divided into 16 sub-sections delimited 
using curly-bracket notation {}. Each section is commented and 
annotated for reproducibility. Many lines where analyses are 
executed are commented out so that processed RDS objects (provided 
here, see below) can be loaded. This script focuses on processing
and analyzing the *output* from *janus* and other software applied
herein. The Go code for *janus* can be found here: https://git.sr.ht/~hms/janus

**Functions_consensus_genetrees.R**

This file contains many custom or new functions called by the 
primary analysis script.

**RandomForest_var_imp.R**

This script contains the code for RandomForest supervised 
classification using the tidymodels framework. This code is 
adapted from the Physalia course examples for machine learning 
in R.

**shift_model_sims.R**
This script contains the code for assessing the janus
false positive/false negative performance. We include
functions defining a generic pipeline that can be used to
assess model performance on arbitrary topologies.

**Supplementary Table 3.xlsx**

Supplementary Table 3, as noted in the manuscript. This file
describes assembly details for assembly and extraction of 
mtDNA genome data from the original AHE target capture data.

**LHT_DATA**

LHT_reference.xlsx contains a database of life-history traits
analyzed in the present study, which are collated from a variety 
of sources (as described in the manuscript). This table 
contains only raw data, and does not include estimated values
(See R code and RDS files).

**RDS**

This folder contains RDS files (R data objects). These files 
correspond to intermediate or final R data objects 
generated in the course of analysis and which are called/loaded 
by the primary analysis script "constrained_species_tree.R".

**AHE_REASSEMBLY**

This folder contains files related to assembly of the AHE dataset.

**mtDNA_REASSEMBLY**

This folder contains files related to assembly of the mtDNA dataset.

**BMR**

This folder contains scripts and files related to 
analysis of the BMR dataset.

**janus**

This folder contains the relevant results from janus 
(e.g.,input/output files). The source tree file is in 
‘timetree_all_taxa_OW_2019.nextree.tre’ available in 
the supplementary materials from Kimball et al. (2019) 
here: https://www.mdpi.com/1424-2818/11/7/109/s1, 
labeled as ‘MRL_3backbone.’

