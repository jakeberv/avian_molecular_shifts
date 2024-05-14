
## R objects

**RDS folder**

* RDS (R Data Serialization) files are a common format for saving R objects, and they allow R users to preserve the state of an object between R sessions. For reproducibility, we provide RDS objects corresponding to intermediate or final R data objects generated in the course of analysis. For example, various model fits, transformed data frames or data manipulations, or other outputs from many analyses performed in our analysis pipeline [constrained_species_tree.R](../constrained_species_tree.R)). These files are required to reproduce our analyses or figures, and will be loaded automatically (where necessary) by the primary script if a user downloads the entire GitHub repository.