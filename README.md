# Molecular early burst associated with the diversification of birds at the K–Pg boundary
Author details: Jacob S. Berv (1,2,*), Sonal Singhal (3), Daniel J.
Field (4,5), Nathanael Walker-Hale (6), Sean W. McHugh (7), J. Ryan
Shipley (8), Eliot T. Miller (9), Rebecca T. Kimball (10), Edward L.
Braun (10), Alex Dornburg (11), C. Tomomi Parins-Fukuchi (12), Richard
O. Prum (13,14), Matt Friedman (2, 15), Stephen A. Smith (1)

Corresponding Author: *Jacob S. Berv, jberv@umich.edu,
jacob.berv@gmail.com

Author affiliations:

1.	Department of Ecology and Evolutionary Biology, 1105 North
University Avenue, Biological Sciences Building, University of Michigan,
Ann Arbor, Michigan, 48109-1085, USA

2.	University of Michigan Museum of Paleontology, 1105 North
University Avenue, Biological Sciences Building, University of Michigan,
Ann Arbor, Michigan, 48109-1085, USA

3.	Department of Biology, California State University, Dominguez
Hills, Carson, California 90747, USA

4.	Department of Earth Sciences, Downing Street, University of
Cambridge, Cambridge CB2 3EQ, UK

5.	Museum of Zoology, Downing Street, University of Cambridge,
Cambridge CB2 3EJ, UK

6.	Department of Plant Sciences, Downing Street, University of
Cambridge, Cambridge, CB2 3EA, UK

7.	Department of Evolution, Ecology, and Population Biology,
Washington University in St Louis, St Louis, Missouri, USA

8.	Department of Forest Dynamics Swiss Federal Institute for Forest,
Snow, and Landscape Research WSL Zürcherstrasse 111 8903 Birmensdorf,
Switzerland

9.	Macaulay Library, Cornell Lab of Ornithology, Ithaca, New York,
14850, USA

10.	Department of Biology, University of Florida, Gainesville,
Florida 32611, USA

11.	Department of Bioinformatics and Genomics, University of North
Carolina at Charlotte, Charlotte, North Carolina, USA

12.	Department of Ecology and Evolutionary Biology, University of
Toronto, Toronto, Ontario, Canada, M5S 3B2

13.	Department of Ecology and Evolutionary Biology, Yale University,
New Haven, Connecticut, 06520, USA

14.	Peabody Museum of Natural History, Yale University, New Haven,
Connecticut, 06520, USA

15.	Department of Earth and Environmental Sciences, 1100 North
University Avenue, University of Michigan, Ann Arbor, Michigan,
48109-1085, USA

## Last updated 3 November 2022

## Primary files

**constrained_species_tree.R**

This is a large R script containing the code for the primary 
series of analyses. It is divided into 16 sub-sections delimited 
using curly-bracket notation {}. Each section is commented and 
annotated for reproducibility. Many lines where analyses are 
executed are commented out so that processed RDS objects (provided 
here, see below) can be loaded.

**Functions_consensus_genetrees.R**

This file contains many custom or new functions called by the 
primary analysis script.

**RandomForest_var_imp.R**

This script contains the code for RandomForest supervised 
classification using the tidymodels framework. This code is 
adapted from the Physalia course examples for machine learning 
in R.

**.RDS files**

These files correspond to intermediate or final R data objects 
generated in the course of analysis and which are called/loaded 
by the primary analysis script "constrained_species_tree.R".

**mtDNA_supp_table.xlsx**

Supplementary Table 3, as noted in the manuscript. This file
describes assembly details for assembly and extraction of 
mtDNA genome data from the original AHE target capture data.
