---
output:
  pdf_document: default
  html_document: default
---

pre-print available: <https://doi.org/10.1101/2022.10.21.513146>

## Jacob S. Berv<sup>1,2,3,\*</sup>, Sonal Singhal<sup>4</sup>, Daniel J. Field<sup>5,6</sup>, Nathanael Walker-Hale<sup>7</sup>, Sean W. McHugh<sup>8</sup>, J. Ryan Shipley<sup>9</sup>, Eliot T. Miller<sup>10</sup>, Rebecca T. Kimball<sup>11</sup>, Edward L. Braun<sup>11</sup>, Alex Dornburg<sup>12</sup>, C. Tomomi Parins-Fukuchi<sup>13</sup>, Richard O. Prum<sup>14,15</sup>, Benjamin M. Winger<sup>1,3</sup>, Matt Friedman<sup>2, 16</sup>, Stephen A. Smith<sup>1</sup>

Corresponding Author: \*Jacob S. Berv, [jberv\@umich.edu](mailto:jberv@umich.edu){.email}, [jacob.berv\@gmail.com](mailto:jacob.berv@gmail.com){.email} Relevant co-author: Sean Mchugh [sean.mchugh4\@gmail.com](mailto:sean.mchugh4@gmail.com){.email}

Author affiliations:

1.  Department of Ecology and Evolutionary Biology, 1105 North University Avenue, Biological Sciences Building, University of Michigan, Ann Arbor, Michigan, 48109-1085, USA

2.  University of Michigan Museum of Paleontology, 1105 North University Avenue, Biological Sciences Building, University of Michigan, Ann Arbor, Michigan, 48109-1085, USA

3.  University of Michigan Museum of Zoology, 1105 North University Avenue, Biological Sciences Building, University of Michigan, Ann Arbor, Michigan, 48109-1085, USA

4.  Department of Biology, California State University, Dominguez Hills, Carson, California 90747, USA

5.  Department of Earth Sciences, Downing Street, University of Cambridge, Cambridge CB2 3EQ, UK

6.  Museum of Zoology, Downing Street, University of Cambridge, Cambridge CB2 3EJ, UK

7.  Department of Plant Sciences, Downing Street, University of Cambridge, Cambridge, CB2 3EA, UK

8.  Department of Evolution, Ecology, and Population Biology, Washington University in St Louis, St Louis, Missouri, USA

9.  Department of Forest Dynamics Swiss Federal Institute for Forest, Snow, and Landscape Research WSL Zürcherstrasse 111 8903 Birmensdorf, Switzerland

10. Macaulay Library, Cornell Lab of Ornithology, Ithaca, New York, 14850, USA

11. Department of Biology, University of Florida, Gainesville, Florida 32611, USA

12. Department of Bioinformatics and Genomics, University of North Carolina at Charlotte, Charlotte, North Carolina, USA

13. Department of Ecology and Evolutionary Biology, University of Toronto, Toronto, Ontario, Canada, M5S 3B2

14. Department of Ecology and Evolutionary Biology, Yale University, New Haven, Connecticut, 06520, USA

15. Peabody Museum of Natural History, Yale University, New Haven, Connecticut, 06520, USA

16. Department of Earth and Environmental Sciences, 1100 North University Avenue, University of Michigan, Ann Arbor, Michigan, 48109-1085, USA

## Metabolic allometry

---
Bayou allometry shifting regimes model pipeline written by: Sean McHugh, Josef Uyeda, Jacob Berv

**full_restrict_equal_10Mx3_ME01_L8.R**

This script estimates the evolution of an allometric slope and intercept 
as a multiregime OU process using the R package bayou (Uyeda and Harmon 
2014; Uyeda et al. 2017), with shifts  either provided a priori or 
estimated via reversible jump MCMC. Additional details describing how
this script works are provided in the supplementary methods section
of our manuscript. After performing the analysis, our script includes code 
to generate a version of Figure 3 in our manuscript.

In the top level directory, we include several files that are accessed 
by the provided R script:

*nomiss.RDS* 
'R Data Serialization' formated file containing the metabolic rate dataset (described in the manuscript)

*simmap.janus.nuc.alldata.RDS*
'R Data Serialization' formated file containing the tree object (described in the manuscript)

*fixed.zip*
Zip file containing bayou output files for each tested configuration of
'fixed' shifts -- eg shifts that are defined on the basis of molecular
shifts identified in each distinct data type.

*global1.zip, global2.zip, global3.zip*
These zip files contain the bayou output files for analyses
(repeated three times) of an 'unconstrained' global shift search

*saved_objects* directory 
This directory contains the mcmc chain data generated in the 
bayou model search described in the provided R script. These files
are provided so that an interseted reader can load these files into
our script and be able to reproduce our Figure 3.
---

References Josef C. Uyeda, Luke J. Harmon, A Novel Bayesian Method for Inferring and Interpreting the Dynamics of Adaptive Landscapes from Phylogenetic Comparative Data, Systematic Biology, Volume 63, Issue 6, November 2014, Pages 902–918, <https://doi.org/10.1093/sysbio/syu057>

The Evolution of Energetic Scaling across the Vertebrate Tree of Life Josef C. Uyeda, Matthew W. Pennell, Eliot T. Miller, Rafael Maia, and Craig R. McClain The American Naturalist 2017 190:2, 185-199
