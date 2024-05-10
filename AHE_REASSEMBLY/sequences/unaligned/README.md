**locus-filtering.R** is an R script used for data filtering after re-assembly (see supplementary methods):

Before alignment, we applied a series of sequential filtering steps to remove remaining short or low-quality fragments. We removed (1) leading and trailing N characters from each fragment, and resulting sequences that were zero length (see locus-filtering.R script), (2) fragments with > 40% N characters, (3) fragments that were < 50 bp long, and (4) whole loci that lacked coverage for at least 10% of the taxa in the dataset.

**min50bp_min10p.zip** contains filtered assemblies such that 1) any locus shorter than 50bp is removed, 2) any locus with fewer than 10% of the taxa is removed.