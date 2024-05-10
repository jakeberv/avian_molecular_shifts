## Metabolic allometry

---
### Bayou allometry shifting regimes model pipeline written by: Sean McHugh, Josef Uyeda, Jacob Berv

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

**nomiss.RDS**

'R Data Serialization' formatted file containing the metabolic rate dataset (described in the manuscript)

**simmap.janus.nuc.alldata.RDS**

'R Data Serialization' formatted file containing the tree object (described in the manuscript)

**fixed.zip**

Zip file containing bayou output files for each tested configuration of
'fixed' shifts -- eg shifts that are defined on the basis of molecular
shifts identified in each distinct data type.

**global1.zip, global2.zip, global3.zip**

These zip files contain the bayou output files for analyses
(repeated three times) of an 'unconstrained' global shift search

**saved_objects**

This directory contains the mcmc chain data generated in the 
bayou model search described in the provided R script. These files
are provided so that an interseted reader can load these files into
our script and be able to reproduce our Figure 3.

---

### References

Josef C. Uyeda, Luke J. Harmon, A Novel Bayesian Method for Inferring and Interpreting the Dynamics of Adaptive Landscapes from Phylogenetic Comparative Data, Systematic Biology, Volume 63, Issue 6, November 2014, Pages 902–918, <https://doi.org/10.1093/sysbio/syu057>

The Evolution of Energetic Scaling across the Vertebrate Tree of Life Josef C. Uyeda, Matthew W. Pennell, Eliot T. Miller, Rafael Maia, and Craig R. McClain The American Naturalist 2017 190:2, 185-199
