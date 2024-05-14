## Primary sequence datasets

Here we provide files containing our processed sequence data sets. Details of data set processing are described in the supplementary methods and include sequence trimming and alignment steps.

[**unaligned**](/unaligned)

* This directory contains two files. 'locus-filtering.R' is an R script used for data filtering after re-assembly. min50bp_min10p.zip contains filtered assemblies such that 1) any locus shorter than 50bp is removed, 2) any locus with fewer than 10% of the taxa is removed.

[**aligned**](/aligned)

* This directory contains two subdirectories. 'sample2haplo' contains the full aligned datasets including both reconstructed alleles from each taxon. 'sample1haplo' contains the sub-sampled datasets passed to analyses with Janus (described in our manuscript).

[**concat**](/concat)

* This directory contains concatenated sequence data files