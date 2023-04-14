# Last updated 13 April 2023
# Master R code for executing analyses for
# Molecular early burst associated with the diversification of birds at the K–Pg boundary

#load packages to get started
require(phytools)
require(ape)
require(treeio)
require(RColorBrewer)
require(ggtree)
library(doSNOW)
library(doParallel)
library(parallel)
require(geiger)
require(viridis)
require(coRdon)
require(seqinr)
require(Rphylopars)
require(pegas)
require(mvMORPH)
require(english)
require(car)
require(phylolm)
require(words2number)
require(OUwie)
require(GISTools)
require(Quartet)
require(readxl)
require(pbmcapply)
require(caper)
require(Biobase)
require(rr2)
require(dplyr)
require(ggridges)
require(ggplot2)
require(ggtree)
require(pheatmap)
require(ggpubr)
require(ggradar)
require(fmsb)
require(TeachingDemos)
require(ratematrix)
require(mapplots)
require(autoimage)

#read in function definitions
source("Functions_consensus_genetrees.R")

#system report and package references
{
#require(report)
#session <- sessionInfo()
#r <- report(session)
#saveRDS(r, file='./RDS/session_report.RDS')

##Note, treePL must be installed in your path for congruification to work
  
#report_system()
#Analyses were conducted using the R Statistical language (version 4.2.2; R Core Team, 2022) on macOS Ventura 13.2.1

#report_packages(include_R = FALSE)
#- geiger (version 2.0.10; Alfaro M et al., 2009)
#- iterators (version 1.0.14; Analytics R, Weston S, 2022)
#- magrittr (version 2.0.3; Bache S, Wickham H, 2022)
#- OUwie (version 2.10; Beaulieu JM, O'Meara B, 2022)
#- maps (version 3.4.1; Becker OScbRA et al., 2022)
#- ggradar (version 0.2; Bion R, 2022)
#- maptools (version 1.1.6; Bivand R, Lewin-Koh N, 2022)
#- rgeos (version 0.6.1; Bivand R, Rundel C, 2022)
#- GISTools (version 0.7.4; Brunsdon C, Chen H, 2014)
#- ratematrix (version 1.2.4; Caetano D, Harmon L, 2022)
#- seqinr (version 4.2.23; Charif D, Lobry J, 2007)
#- mvMORPH (version 1.1.6; Clavel J et al., 2015)
#- doParallel (version 1.0.17; Corporation M, Weston S, 2022)
#- doSNOW (version 1.0.20; Corporation M, Weston S, 2022)
#- coRdon (version 1.16.0; Elek A et al., 2022)
#- english (version 1.2.6; Fox J et al., 2021)
#- car (version 3.1.1; Fox J, Weisberg S, 2019)
#- carData (version 3.0.5; Fox J et al., 2022)
#- autoimage (version 2.2.3; French JP, 2017)
#- viridis (version 0.6.2; Garnier et al., 2021)
#- viridisLite (version 0.4.1; Garnier et al., 2022)
#- mvtnorm (version 1.1.3; Genz A et al., 2021)
#- mapplots (version 1.5.1; Gerritsen H, 2018)
#- Rphylopars (version 0.3.9; Goolsby E et al., 2022)
#- phylolm (version 2.6.2; Ho LST, Ane C, 2014)
#- Biobase (version 2.58.0; Huber W et al., 2015)
#- BiocGenerics (version 0.44.0; Huber et al., 2015)
#- rr2 (version 1.1.0; Ives AR, 2018)
#- nloptr (version 2.0.3; Johnson SG, ?)
#- ggpubr (version 0.5.0; Kassambara A, 2022)
#- subplex (version 1.8; King AA, Rowan T, 2022)
#- pheatmap (version 1.0.12; Kolde R, 2019)
#- pbmcapply (version 1.5.1; Kuang K et al., 2022)
#- report (version 0.5.7; Makowski D et al., 2023)
#- foreach (version 1.5.2; Microsoft, Weston S, 2022)
#- fmsb (version 0.7.5; Nakazawa M, 2023)
#- RColorBrewer (version 1.1.3; Neuwirth E, 2022)
#- caper (version 1.0.1; Orme D et al., 2018)
#- pegas (version 1.1; Paradis E, 2010)
#- ape (version 5.7; Paradis E, Schliep K, 2019)
#- sp (version 1.6.0; Pebesma EJ, Bivand RS, 2005)
#- phytools (version 1.2.0; Revell LJ, 2012)
#- corpcor (version 1.6.10; Schafer J et al., 2021)
#- Quartet (version 1.2.5; Smith MR, 2019)
#- TreeTools (version 1.9.0; Smith MR, 2019)
#- TeachingDemos (version 2.12; Snow G, 2020)
#- snow (version 0.4.4; Tierney L et al., 2021)
#- MASS (version 7.3.58.1; Venables WN, Ripley BD, 2002)
#- words2number (version 0.2.0; Vera Oteo J, Marwick B, 2018)
#- ggplot2 (version 3.4.1; Wickham H, 2016)
#- readxl (version 1.4.1; Wickham H, Bryan J, 2022)
#- dplyr (version 1.1.0; Wickham H et al., 2023)
#- ggridges (version 0.5.4; Wilke C, 2022)
#- treeio (version 1.22.0; Yu G, 2022)
#- ggtree (version 3.6.2; Yu G, 2022)

#cite_packages()
#- Alfaro M, Santini F, Brock C, Alamillo H, Dornburg A, Rabosky D, Carnevale G, Harmon L (2009). “Nine exceptional radiations plus high turnover explain species diversity in jawed vertebrates.” _Proceedings of the National Academy of Sciences of the United States of America_, *106*, 13410-13414. Eastman J, Alfaro M, Joyce P, Hipp A, Harmon L (2011). “A novel comparative method for identifying shifts in the rate of character evolution on trees.” _Evolution_, *65*, 3578-3589. Slater G, Harmon L, Wegmann D, Joyce P, Revell L, Alfaro M (2012). “Fitting models of continuous trait evolution to incompletely sampled comparative data using approximate Bayesian computation.” _Evolution_, *66*, 752-762. Harmon L, Weir J, Brock C, Glor R, Challenger W (2008). “GEIGER: investigating evolutionary radiations.” _Bioinformatics_, *24*, 129-131. Pennell M, Eastman J, Slater G, Brown J, Uyeda J, Fitzjohn R, Alfaro M, Harmon L (2014). “geiger v2.0: an expanded suite of methods for fitting macroevolutionary models to phylogenetic trees.” _Bioinformatics_, *30*, 2216-2218.
#- Analytics R, Weston S (2022). _iterators: Provides Iterator Construct_. R package version 1.0.14, <https://CRAN.R-project.org/package=iterators>.
#- Bache S, Wickham H (2022). _magrittr: A Forward-Pipe Operator for R_. R package version 2.0.3, <https://CRAN.R-project.org/package=magrittr>.
#- Beaulieu JM, O'Meara B (2022). _OUwie: Analysis of Evolutionary Rates in an OU Framework_. R package version 2.10, <https://github.com/thej022214/OUwie>.
#- Becker OScbRA, Minka ARWRvbRBEbTP, Deckmyn. A (2022). _maps: Draw Geographical Maps_. R package version 3.4.1, <https://CRAN.R-project.org/package=maps>.
#- Bion R (2022). _ggradar: Create radar charts using ggplot2_. R package version 0.2.
#- Bivand R, Lewin-Koh N (2022). _maptools: Tools for Handling Spatial Objects_. R package version 1.1-6, <https://CRAN.R-project.org/package=maptools>.
#- Bivand R, Rundel C (2022). _rgeos: Interface to Geometry Engine - Open Source ('GEOS')_. R package version 0.6-1, <https://CRAN.R-project.org/package=rgeos>.
#- Brunsdon C, Chen H (2014). _GISTools: Some further GIS capabilities for R_. R package version 0.7-4, <https://CRAN.R-project.org/package=GISTools>.
#- Caetano D, Harmon L (2022). _ratematrix: Bayesian Estimation of the Evolutionary Rate Matrix_. R package version 1.2.4, <https://CRAN.R-project.org/package=ratematrix>.
#- Charif D, Lobry J (2007). “SeqinR 1.0-2: a contributed package to the R project for statistical computing devoted to biological sequences retrieval and analysis.” In Bastolla U, Porto M, Roman H, Vendruscolo M (eds.), _Structural approaches to sequence evolution: Molecules, networks, populations_, series Biological and Medical Physics, Biomedical Engineering, 207-232. Springer Verlag, New York. ISBN : 978-3-540-35305-8.
#- Clavel J, Escarguel G, Merceron G (2015). “mvMORPH: an R package for fitting multivariate evolutionary models to morphometric data.” _Methods in Ecology and Evolution_, *6*, 1311-1319.
#- Corporation M, Weston S (2022). _doParallel: Foreach Parallel Adaptor for the 'parallel' Package_. R package version 1.0.17, <https://CRAN.R-project.org/package=doParallel>.
#- Corporation M, Weston S (2022). _doSNOW: Foreach Parallel Adaptor for the 'snow' Package_. R package version 1.0.20, <https://CRAN.R-project.org/package=doSNOW>.
#- Elek A, Kuzman M, Vlahovicek K (2022). _coRdon: Codon Usage Analysis and Prediction of Gene Expressivity_. R package version 1.16.0, <https://github.com/BioinfoHR/coRdon>.
#- Fox J, Venables B, Damico A, Salverda AP (2021). _english: Translate Integers into English_. R package version 1.2-6, <https://CRAN.R-project.org/package=english>.
#- Fox J, Weisberg S (2019). _An R Companion to Applied Regression_, Third edition. Sage, Thousand Oaks CA. <https://socialsciences.mcmaster.ca/jfox/Books/Companion/>.
#- Fox J, Weisberg S, Price B (2022). _carData: Companion to Applied Regression Data Sets_. R package version 3.0-5, <https://CRAN.R-project.org/package=carData>.
#- French JP (2017). “autoimage: Multiple Heat Maps for Projected Coordinates.” _The R Journal_, *9*, 284-297. <https://journal.r-project.org/archive/2017/RJ-2017-025/RJ-2017-025.pdf>.
#- Garnier, Simon, Ross, Noam, Rudis, Robert, Camargo, Pedro A, Sciaini, Marco, Scherer, Cédric (2021). _viridis - Colorblind-Friendly Color Maps for R_. doi:10.5281/zenodo.4679424 <https://doi.org/10.5281/zenodo.4679424>, R package version 0.6.2, <https://sjmgarnier.github.io/viridis/>.
#- Garnier, Simon, Ross, Noam, Rudis, Robert, Camargo, Pedro A, Sciaini, Marco, Scherer, Cédric (2022). _viridis - Colorblind-Friendly Color Maps for R_. doi:10.5281/zenodo.4679424 <https://doi.org/10.5281/zenodo.4679424>, R package version 0.4.1, <https://sjmgarnier.github.io/viridis/>.
#- Genz A, Bretz F, Miwa T, Mi X, Leisch F, Scheipl F, Hothorn T (2021). _mvtnorm: Multivariate Normal and t Distributions_. R package version 1.1-3, <https://CRAN.R-project.org/package=mvtnorm>. Genz A, Bretz F (2009). _Computation of Multivariate Normal and t Probabilities_, series Lecture Notes in Statistics. Springer-Verlag, Heidelberg. ISBN 978-3-642-01688-2.
#- Gerritsen H (2018). _mapplots: Data Visualisation on Maps_. R package version 1.5.1, <https://CRAN.R-project.org/package=mapplots>.
#- Goolsby E, Bruggeman J, Ane C (2022). _Rphylopars: Phylogenetic Comparative Tools for Missing Data and Within-Species Variation_. R package version 0.3.9, <https://CRAN.R-project.org/package=Rphylopars>.
#- Ho LST, Ane C (2014). “A linear-time algorithm for Gaussian and non-Gaussian trait evolution models.” _Systematic Biology_, *63*, 397-408.
#- Huber W, Carey VJ, Gentleman R, Anders S, Carlson M, Carvalho BS, Bravo HC, Davis S, Gatto L, Girke T, Gottardo R, Hahne F, Hansen KD, Irizarry RA, Lawrence M, Love MI, MacDonald J, Obenchain V, Ole's AK, Pag`es H, Reyes A, Shannon P, Smyth GK, Tenenbaum D, Waldron L, Morgan M (2015). “Orchestrating high-throughput genomic analysis with Bioconductor.” _Nature Methods_, *12*(2), 115-121. <http://www.nature.com/nmeth/journal/v12/n2/full/nmeth.3252.html>.
#- Huber, W., Carey, J. V, Gentleman, R., Anders, S., Carlson, M., Carvalho, S. B, Bravo, C. H, Davis, S., Gatto, L., Girke, T., Gottardo, R., Hahne, F., Hansen, D. K, Irizarry, A. R, Lawrence, M., Love, I. M, MacDonald, J., Obenchain, V., Ole's, K. A, Pag`es, H., Reyes, A., Shannon, P., Smyth, K. G, Tenenbaum, D., Waldron, L., Morgan, M. (2015). “Orchestrating high-throughput genomic analysis with Bioconductor.” _Nature Methods_, *12*(2), 115-121. <http://www.nature.com/nmeth/journal/v12/n2/full/nmeth.3252.html>.
#- Ives AR (2018). “R^2s for Correlated Data: Phylogenetic Models, LMMs, and GLMMs.” _Systematic Biology_, syy060. <https://doi.org/10.1093/sysbio/syy060>. Ives AR, Li D (2018). “rr2: An R package to calculate R^2s for regression models.” _The Journal of Open Source Software_, *3*(30), 1028. <https://doi.org/10.21105/joss.01028>.
#- Johnson SG (?). “The NLopt nonlinear-optimization package.” _?_, *?*(?), ?
#- Kassambara A (2022). _ggpubr: 'ggplot2' Based Publication Ready Plots_. R package version 0.5.0, <https://CRAN.R-project.org/package=ggpubr>.
#- King AA, Rowan T (2022). _subplex: Unconstrained Optimization using the Subplex Algorithm_. R package version 1.8, <https://CRAN.R-project.org/package=subplex>.
#- Kolde R (2019). _pheatmap: Pretty Heatmaps_. R package version 1.0.12, <https://CRAN.R-project.org/package=pheatmap>.
#- Kuang K, Kong Q, Napolitano F (2022). _pbmcapply: Tracking the Progress of Mc*pply with Progress Bar_. R package version 1.5.1, <https://CRAN.R-project.org/package=pbmcapply>.
#- Makowski D, Lüdecke D, Patil I, Thériault R, Ben-Shachar M, Wiernik B (2023). “Automated Results Reporting as a Practical Tool to Improve Reproducibility and Methodological Best Practices Adoption.” _CRAN_. <https://easystats.github.io/report/>.
#- Microsoft, Weston S (2022). _foreach: Provides Foreach Looping Construct_. R package version 1.5.2, <https://CRAN.R-project.org/package=foreach>.
#- Nakazawa M (2023). _fmsb: Functions for Medical Statistics Book with some Demographic Data_. R package version 0.7.5, <https://CRAN.R-project.org/package=fmsb>.
#- Neuwirth E (2022). _RColorBrewer: ColorBrewer Palettes_. R package version 1.1-3, <https://CRAN.R-project.org/package=RColorBrewer>.
#- Orme D, Freckleton R, Thomas G, Petzoldt T, Fritz S, Isaac N, Pearse W (2018). _caper: Comparative Analyses of Phylogenetics and Evolution in R_. R package version 1.0.1, <https://CRAN.R-project.org/package=caper>.
#- Paradis E (2010). “pegas: an R package for population genetics with an integrated-modular approach.” _Bioinformatics_, *26*, 419-420.
#- Paradis E, Schliep K (2019). “ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R.” _Bioinformatics_, *35*, 526-528. doi:10.1093/bioinformatics/bty633 <https://doi.org/10.1093/bioinformatics/bty633>.
#- Pebesma EJ, Bivand RS (2005). “Classes and methods for spatial data in R.” _R News_, *5*(2), 9-13. <https://CRAN.R-project.org/doc/Rnews/>. Bivand RS, Pebesma E, Gomez-Rubio V (2013). _Applied spatial data analysis with R, Second edition_. Springer, NY. <https://asdar-book.org/>.
#- R Core Team (2022). _R: A Language and Environment for Statistical Computing_. R Foundation for Statistical Computing, Vienna, Austria. <https://www.R-project.org/>.
#- Revell LJ (2012). “phytools: An R package for phylogenetic comparative biology (and other things).” _Methods in Ecology and Evolution_, *3*, 217-223.
#- Schafer J, Opgen-Rhein R, Zuber V, Ahdesmaki M, Silva APD, Strimmer. K (2021). _corpcor: Efficient Estimation of Covariance and (Partial) Correlation_. R package version 1.6.10, <https://CRAN.R-project.org/package=corpcor>.
#- Smith MR (2019). _Quartet: comparison of phylogenetic trees using quartet and split measures_. doi:10.5281/zenodo.2536318 <https://doi.org/10.5281/zenodo.2536318>, R package version 1.2.5. Sand A, Holt MK, Johansen J, Brodal GS, Mailund T, Pedersen CNS (2014). “tqDist: a library for computing the quartet and triplet distances between binary or general trees.” _Bioinformatics_, *30*(14), 2079-2080. doi:10.1093/bioinformatics/btu157 <https://doi.org/10.1093/bioinformatics/btu157>. Smith MR (2019). “Bayesian and parsimony approaches reconstruct informative trees from simulated morphological datasets.” _Biology Letters_, *15*(2), 20180632. doi:10.1098/rsbl.2018.0632 <https://doi.org/10.1098/rsbl.2018.0632>.
#- Smith MR (2019). _TreeTools: create, modify and analyse phylogenetic trees_. Comprehensive R Archive Network. doi:10.5281/zenodo.3522725 <https://doi.org/10.5281/zenodo.3522725>, R package version 1.9.0.
#- Snow G (2020). _TeachingDemos: Demonstrations for Teaching and Learning_. R package version 2.12, <https://CRAN.R-project.org/package=TeachingDemos>.
#- Tierney L, Rossini AJ, Li N, Sevcikova H (2021). _snow: Simple Network of Workstations_. R package version 0.4-4, <https://CRAN.R-project.org/package=snow>.
#- Venables WN, Ripley BD (2002). _Modern Applied Statistics with S_, Fourth edition. Springer, New York. ISBN 0-387-95457-0, <https://www.stats.ox.ac.uk/pub/MASS4/>.
#- Vera Oteo J, Marwick B (2018). _words2number: Convert number words to numeric digits_. R package version 0.2.0, <https://github.com/benmarwick/words2number>.
#- Wickham H (2016). _ggplot2: Elegant Graphics for Data Analysis_. Springer-Verlag New York. ISBN 978-3-319-24277-4, <https://ggplot2.tidyverse.org>.
#- Wickham H, Bryan J (2022). _readxl: Read Excel Files_. R package version 1.4.1, <https://CRAN.R-project.org/package=readxl>.
#- Wickham H, François R, Henry L, Müller K, Vaughan D (2023). _dplyr: A Grammar of Data Manipulation_. R package version 1.1.0, <https://CRAN.R-project.org/package=dplyr>.
#- Wilke C (2022). _ggridges: Ridgeline Plots in 'ggplot2'_. R package version 0.5.4, <https://CRAN.R-project.org/package=ggridges>.
#- Yu G (2022). _Data Integration, Manipulation and Visualization of Phylogenetic Treess_, 1st edition edition. Chapman and Hall/CRC. <https://www.amazon.com/Integration-Manipulation-Visualization-Phylogenetic-Computational-ebook/dp/B0B5NLZR1Z/>. Wang L, Lam TT, Xu S, Dai Z, Zhou L, Feng T, Guo P, Dunn CW, Jones BR, Bradley T, Zhu H, Guan Y, Jiang Y, Yu G (2020). “treeio: an R package for phylogenetic tree input and output with richly annotated and associated data.” _Molecular Biology and Evolution_, *37*, 599-603. doi:10.1093/molbev/msz240 <https://doi.org/10.1093/molbev/msz240>.
#- Yu G (2022). _Data Integration, Manipulation and Visualization of Phylogenetic Treess_, 1st edition edition. Chapman and Hall/CRC. <https://www.amazon.com/Integration-Manipulation-Visualization-Phylogenetic-Computational-ebook/dp/B0B5NLZR1Z/>. Yu G (2020). “Using ggtree to Visualize Data on Tree-Like Structures.” _Current Protocols in Bioinformatics_, *69*(1), e96. doi:10.1002/cpbi.96 <https://doi.org/10.1002/cpbi.96>, <https://currentprotocols.onlinelibrary.wiley.com/doi/abs/10.1002/cpbi.96>. Yu G, Lam TT, Zhu H, Guan Y (2018). “Two methods for mapping and visualizing associated data on phylogeny using ggtree.” _Molecular Biology and Evolution_, *35*, 3041-3043. doi:10.1093/molbev/msy194 <https://doi.org/10.1093/molbev/msy194>, <https://academic.oup.com/mbe/article/35/12/3041/5142656>. Yu G, Smith D, Zhu H, Guan Y, Lam TT (2017). “ggtree: an R package for visualization and annotation of phylogenetic trees with their covariates and other associated data.” _Methods in Ecology and Evolution_, *8*, 28-36. doi:10.1111/2041-210X.12628 <https://doi.org/10.1111/2041-210X.12628>, <http://onlinelibrary.wiley.com/doi/10.1111/2041-210X.12628/abstract>.
  
}

#Section 1, setting up
#this code uses bracket notation {} to delineate discrete sections
#each section is labeled with a comment, 
#section contents are described within each section
#code is formatted to read .RDS files representing 
#intermediate data objects
{

#generate file paths
{
##run analyses##
#import time scaled reference from kimball for congruification
#this is generated by running kimball_time_parser.R
#local file path (Jake's computer)
#time_scale<-read.tree(file="/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/unmasked/min2x/qual20_cov2_haplo_bestonly/initial_test_filters/min50bp_min10p_aligned/ALIGNED/phased/sample1haplo/phyhetnucbf/backbone_constraint/timetree.prune.nameswap.tre")
#git repo file path
time_scale<-read.tree(file="./trees/timetree.prune.nameswap.tre")

#important file paths

#file path for janus directory reflecting analysis of the entire dataset
#consensus.all.path<-"/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/unmasked/min2x/qual20_cov2_haplo_bestonly/initial_test_filters/min50bp_min10p_aligned/ALIGNED/phased/sample1haplo/phyhetnucbf/backbone_constraint/janus/8952e31d/gamma/redo_m3/ALL_MFP_MERGE_MRL3_constraint_G_UE_UL_M3/ALL_MFP_MERGE_MRL3_constraint.rooted.treefile.gophy.results.tre"
consensus.all.path<-"./janus/files/ALL_MFP_MERGE_MRL3_constraint_G_UE_UL_M3/ALL_MFP_MERGE_MRL3_constraint.rooted.treefile.gophy.results.tre"

#file path for janus directory containing exploratory (not noted in text) analysis of individual gene trees
loci.standard.path<-"/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/unmasked/min2x/qual20_cov2_haplo_bestonly/initial_test_filters/min50bp_min10p_aligned/ALIGNED/phased/sample1haplo/phyhetnucbf/unconstrained_new/janus/0cba26ae/loci_standard/tmp_results/gophy"

#file path for janus directory containing analysis of the exon dataset
#consensus.exons.path<-"/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/unmasked/min2x/qual20_cov2_haplo_bestonly/initial_test_filters/min50bp_min10p_aligned/ALIGNED/phased/sample1haplo/phyhetnucbf/backbone_constraint/janus/8952e31d/gamma/redo_m3/exons_MFP_MERGE_MRL3_constraint_G_UE_UL_M3/exons_MFP_MERGE_MRL3_constraint.rooted.treefile.gophy.results.tre"
consensus.exons.path<-"./janus/files/exons_MFP_MERGE_MRL3_constraint_G_UE_UL_M3/exons_MFP_MERGE_MRL3_constraint.rooted.treefile.gophy.results.tre"

#file path for janus directory containing analysis of the exon dataset
#consensus.introns.path<-"/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/unmasked/min2x/qual20_cov2_haplo_bestonly/initial_test_filters/min50bp_min10p_aligned/ALIGNED/phased/sample1haplo/phyhetnucbf/backbone_constraint/janus/8952e31d/gamma/redo_m3/introns_MFP_MERGE_MRL3_constraint_G_UE_UL_M3/introns_MFP_MERGE_MRL3_constraint.rooted.treefile.gophy.results.tre"
consensus.introns.path<-"./janus/files/introns_MFP_MERGE_MRL3_constraint_G_UE_UL_M3/introns_MFP_MERGE_MRL3_constraint.rooted.treefile.gophy.results.tre"

#file path for janus directory containing analysis of the utr dataset
#consensus.utrs.path<-"/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/unmasked/min2x/qual20_cov2_haplo_bestonly/initial_test_filters/min50bp_min10p_aligned/ALIGNED/phased/sample1haplo/phyhetnucbf/backbone_constraint/janus/8952e31d/gamma/redo_m3/utrs_MFP_MERGE_MRL3_constraint_G_UE_UL_M3/utrs_MFP_MERGE_MRL3_constraint.rooted.treefile.gophy.results.tre"
consensus.utrs.path<-"./janus/files/utrs_MFP_MERGE_MRL3_constraint_G_UE_UL_M3/utrs_MFP_MERGE_MRL3_constraint.rooted.treefile.gophy.results.tre"

#file path for janus directory containing analysis of the mtdna (all) dataset
#consensus.mtdna.all.path<-"/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/mtdnas/janus/concat_all/mtDNA_all_MRL3_constraint.janus.tre.gophy.results.tre"
consensus.mtdna.all.path<-"./janus/files/mtdnas/concat_all/mtDNA_all_MRL3_constraint.janus.tre.gophy.results.tre"

#file path for janus directory containing analysis of the mtdna (protein only) dataset
#consensus.mtdna.proteins.path<-"/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/mtdnas/janus/concat_proteins/mtDNA_proteins_MRL3_constraint.janus.treefile.gophy.results.tre"
consensus.mtdna.proteins.path<-"./janus/files/mtdnas/concat_proteins/mtDNA_proteins_MRL3_constraint.janus.treefile.gophy.results.tre"

#file path for janus directory containing analysis of the mtdna (protein only) dataset
#consensus.mtdna.rRNAs.path<-"/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/mtdnas/janus/concat_rRNA/mtDNA_rRNAs_MRL3_constraint.janus.treefile.gophy.results.tre"
consensus.mtdna.rRNAs.path<-"./janus/files/mtdnas/concat_rRNA/mtDNA_rRNAs_MRL3_constraint.janus.treefile.gophy.results.tre"

#file path for tree estimated with both alleles
#allele.tree.path<-"/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/unmasked/min2x/qual20_cov2_haplo_bestonly/initial_test_filters/min50bp_min10p_aligned/ALIGNED/phased/sample2haplo/tree_building/ALL_MFP_2haplo/ALL_MFP_MERGE_2haplo.treefile"
allele.tree.path<-"./trees/ALL_MFP_2haplo/ALL_MFP_MERGE_2haplo.treefile"

#directories for gene trees with low support edges collapsed
#loci.standard.50.path<-"/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/unmasked/min2x/qual20_cov2_haplo_bestonly/initial_test_filters/min50bp_min10p_aligned/ALIGNED/phased/sample1haplo/tree_building/loci_standard/collapsed/50"
#loci.standard.75.path<-"/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/unmasked/min2x/qual20_cov2_haplo_bestonly/initial_test_filters/min50bp_min10p_aligned/ALIGNED/phased/sample1haplo/tree_building/loci_standard/collapsed/75"
#loci.standard.90.path<-"/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/unmasked/min2x/qual20_cov2_haplo_bestonly/initial_test_filters/min50bp_min10p_aligned/ALIGNED/phased/sample1haplo/tree_building/loci_standard/collapsed/90"
#loci.standard.95.path<-"/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/unmasked/min2x/qual20_cov2_haplo_bestonly/initial_test_filters/min50bp_min10p_aligned/ALIGNED/phased/sample1haplo/tree_building/loci_standard/collapsed/95"

loci.standard.50.path<-"./trees/loci_standard/collapsed/50"
loci.standard.75.path<-"./trees/loci_standard/collapsed/75"
loci.standard.90.path<-"./trees/loci_standard/collapsed/90"
loci.standard.95.path<-"./trees/loci_standard/collapsed/95"

}

#read and process data
{
##### ALL GENE TREES

#read in the consensus tree from BF+G analysis
#v8952e31d
consensus.all<-read.beast.fixed(consensus.all.path, ladderize=T)
#genematch.all<-match_generator(refpath=consensus.all.path, targetpath=loci.standard.path)
#saveRDS(genematch.all, file="./RDS/genematch.all.RDS")
genematch.all<-readRDS(file="./RDS/genematch.all.RDS")

#concordance calculation on gene trees with low support edges collapsed
#all_concordance.95<-gene_concordance_perc(refpath=consensus.all.path, targetpath=loci.standard.95.path)
#saveRDS(all_concordance.95, file="./RDS/all_concordance.95.RDS")
all_concordance.95<-readRDS(file="./RDS/all_concordance.95.RDS")

##### EXONS ONLY

#read in the consensus tree from BF+G analysis
#v8952e31d 
consensus.exons<-read.beast.fixed(consensus.exons.path, ladderize=T)
#consensus.exons.bl<-read.beast.fixed("~/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/unmasked/min2x/qual20_cov2_haplo_bestonly/initial_test_filters/min50bp_min10p_aligned/ALIGNED/phased/sample1haplo/phyhetnucbf/backbone_constraint/janus/5d81a027/branch_lengths/exons_MFP_MERGE_MRL3_constraint_B_RM_UE_UL_M4/exons_MFP_MERGE_MRL3_constraint.rooted.treefile.gophy.results.tre", ladderize=T)
#genematch.exons<-match_generator(refpath=consensus.exons.path, targetpath = loci.standard.path, kind="exon")
#saveRDS(genematch.exons, file="./RDS/genematch.exons.RDS")
genematch.exons<-readRDS(file="./RDS/genematch.exons.RDS")

#genematch.exons.allref<-match_generator(refpath=consensus.all.path, targetpath = loci.standard.path, kind="exon")
#saveRDS(genematch.exons.allref, file="./RDS/genematch.exons.allref.RDS")
genematch.exons.allref<-readRDS(file="./RDS/genematch.exons.allref.RDS")

#concordance data on gene trees with low support edges collapsed
#exons_concordance.allref<-gene_concordance_perc(refpath=consensus.all.path, targetpath = loci.standard.95.path, kind="exon")
#saveRDS(exons_concordance.allref, file="./RDS/exons_concordance.allref.RDS")
exons_concordance.allref<-readRDS(file="./RDS/exons_concordance.allref.RDS")

##### EXONS with MFP+MERGE codons

#read in the consensus tree from BFRM analysis
#consensus.exons<-read.beast.fixed("~/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/unmasked/min2x/qual20_cov2_haplo_bestonly/initial_test_filters/min50bp_min10p_aligned/ALIGNED/phased/sample1haplo/phyhetnucbf/backbone_constraint/janus/bfrm_2_Dec_2020/regular/exons_MFP_MERGE_MRL3_constraint_RM_UE_UL_M4/exons_MFP_MERGE_MRL3_constraint.rooted.treefile.gophy.results.tre", ladderize=T)
#consensus.exons.bl<-read.beast.fixed("~/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/unmasked/min2x/qual20_cov2_haplo_bestonly/initial_test_filters/min50bp_min10p_aligned/ALIGNED/phased/sample1haplo/phyhetnucbf/backbone_constraint/janus/bfrm_2_Dec_2020/branch_lengths/exons_MFP_MERGE_MRL3_constraint_B_RM_UE_UL_M4/exons_MFP_MERGE_MRL3_constraint.rooted.treefile.gophy.results.tre", ladderize=T)
#genematch.exons.MFP<-match_generator(refpath="~/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/unmasked/min2x/qual20_cov2_haplo_bestonly/initial_test_filters/min50bp_min10p_aligned/ALIGNED/phased/sample1haplo/phyhetnucbf/backbone_constraint/janus/5d81a027/regular/exons_MFP_MERGE_MRL3_constraint_RM_UE_UL_M4/exons_MFP_MERGE_MRL3_constraint.rooted.treefile.gophy.results.tre", targetpath = "/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/unmasked/min2x/qual20_cov2_haplo_bestonly/initial_test_filters/min50bp_min10p_aligned/ALIGNED/phased/sample1haplo/phyhetnucbf/unconstrained/exons_partitioned/results/gophytrees/gophy", kind="exon")
#concordance data
#exons_MFP_concordance<-gene_concordance_perc(refpath="~/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/unmasked/min2x/qual20_cov2_haplo_bestonly/initial_test_filters/min50bp_min10p_aligned/ALIGNED/phased/sample1haplo/phyhetnucbf/backbone_constraint/janus/5d81a027/regular/exons_MFP_MERGE_MRL3_constraint_RM_UE_UL_M4/exons_MFP_MERGE_MRL3_constraint.rooted.treefile.gophy.results.tre", targetpath = "/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/unmasked/min2x/qual20_cov2_haplo_bestonly/initial_test_filters/min50bp_min10p_aligned/ALIGNED/phased/sample1haplo/phyhetnucbf/unconstrained/exons_partitioned/results/gophytrees/gophy", kind="exon")


#trying with exons with failed partitions removed
#genematch.exons.MFP.sym<-match_generator(refpath="~/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/unmasked/min2x/qual20_cov2_haplo_bestonly/initial_test_filters/min50bp_min10p_aligned/ALIGNED/phased/sample1haplo/phyhetnucbf/backbone_constraint/janus/5d81a027/regular/exons_MFP_MERGE_MRL3_constraint_RM_UE_UL_M4/exons_MFP_MERGE_MRL3_constraint.rooted.treefile.gophy.results.tre", targetpath = "/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/unmasked/min2x/qual20_cov2_haplo_bestonly/initial_test_filters/min50bp_min10p_aligned/ALIGNED/phased/sample1haplo/phyhetnucbf/unconstrained/exons_partitioned_symtest/results/gophytrees/gophy", kind="exon")
#concordance data
#exons_MFP_concordance<-gene_concordance_perc(refpath="~/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/unmasked/min2x/qual20_cov2_haplo_bestonly/initial_test_filters/min50bp_min10p_aligned/ALIGNED/phased/sample1haplo/phyhetnucbf/backbone_constraint/janus/5d81a027/regular/exons_MFP_MERGE_MRL3_constraint_RM_UE_UL_M4/exons_MFP_MERGE_MRL3_constraint.rooted.treefile.gophy.results.tre", targetpath = "/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/unmasked/min2x/qual20_cov2_haplo_bestonly/initial_test_filters/min50bp_min10p_aligned/ALIGNED/phased/sample1haplo/phyhetnucbf/unconstrained/exons_partitioned_symtest/results/gophytrees/gophy", kind="exon")



##### INTRONS ONLY

#read in the consensus tree from BF+G analysis
#v8952e31d
consensus.introns<-read.beast.fixed(consensus.introns.path, ladderize=T)
#consensus.introns.bl<-read.beast.fixed("~/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/unmasked/min2x/qual20_cov2_haplo_bestonly/initial_test_filters/min50bp_min10p_aligned/ALIGNED/phased/sample1haplo/phyhetnucbf/backbone_constraint/janus/5d81a027/branch_lengths/introns_MFP_MERGE_MRL3_constraint_B_RM_UE_UL_M4/introns_MFP_MERGE_MRL3_constraint.rooted.treefile.gophy.results.tre", ladderize=T)
#genematch.introns<-match_generator(refpath = consensus.introns.path, targetpath = loci.standard.path, kind = "intron|Unknown")
#saveRDS(genematch.introns, file="./RDS/genematch.introns.RDS")
genematch.introns<-readRDS(file="./RDS/genematch.introns.RDS")

#genematch.introns.allref<-match_generator(refpath = consensus.all.path, targetpath = loci.standard.path, kind = "intron|Unknown")
#aveRDS(genematch.introns.allref, file="./RDS/genematch.introns.allref.RDS")
genematch.introns.allref<-readRDS(file="./RDS/genematch.introns.allref.RDS")

#concordance calcs
#introns_concordance.allref<-gene_concordance_perc(refpath = consensus.all.path, targetpath = loci.standard.95.path, kind = "intron|Unknown")
#saveRDS(introns_concordance.allref, file="./RDS/introns_concordance.allref.RDS")
introns_concordance.allref<-readRDS(file="./RDS/introns_concordance.allref.RDS")

##### UTRs ONLY

#read in the consensus tree from BF+G analysis
#v8952e31d
consensus.utrs<-read.beast.fixed(consensus.utrs.path, ladderize=T)
#consensus.utrs.bl<-read.beast.fixed("~/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/unmasked/min2x/qual20_cov2_haplo_bestonly/initial_test_filters/min50bp_min10p_aligned/ALIGNED/phased/sample1haplo/phyhetnucbf/backbone_constraint/janus/5d81a027/branch_lengths/utrs_MFP_MERGE_MRL3_constraint_B_RM_UE_UL_M4_stderror/utrs_MFP_MERGE_MRL3_constraint.rooted.treefile.gophy.results.tre", ladderize=T)
#genematch.utrs<-match_generator(refpath = consensus.utrs.path, targetpath = loci.standard.path, kind = "utr")
#saveRDS(genematch.utrs, file="./RDS/genematch.utrs.RDS")
genematch.utrs<-readRDS(file="./RDS/genematch.utrs.RDS")

#genematch.utrs.allref<-match_generator(refpath = consensus.utrs.path, targetpath = loci.standard.path, kind = "utr")
#saveRDS(genematch.utrs.allref, file="./RDS/genematch.utrs.allref.RDS")
genematch.utrs.allref<-readRDS(file="./RDS/genematch.utrs.allref.RDS")

#concordance calcs
#utrs_concordance.allref<-gene_concordance_perc(refpath = consensus.all.path, targetpath = loci.standard.95.path, kind = "utr")
#saveRDS(utrs_concordance.allref, file="./RDS/utrs_concordance.allref.RDS")
utrs_concordance.allref<-readRDS(file="./RDS/utrs_concordance.allref.RDS")

#reading in the mtDNA datasets

#mtDNAs ONLY
#mtDNA.data<-read_excel(path="/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/mtdnas/MITOGENOME-Berv-Aug-2021.xlsx", sheet = "jake_mod")
mtDNA.data<-read_excel(path="./mtDNA_REASSEMBLY/MITOGENOME-Berv-Aug-2021.xlsx", sheet = "jake_mod")
mtDNA.data<-data.frame(treeID=mtDNA.data$treeID, mtdnaNames = mtDNA.data$`NAMES FOR EXPORT`)
mtDNA.data<-mtDNA.data[complete.cases(mtDNA.data),]

#read in the consensus tree from BF+G analysis
#ran with build 0cba26ae
consensus.mtdnas.all<-read.beast.fixed(consensus.mtdna.all.path, ladderize=T)
#need to translate the names to those which are phylogenetically equivalent 
#in the nuclear dataset
nameswap<-merge(data.frame(mtdnaNames=consensus.mtdnas.all@phylo$tip.label), mtDNA.data, by='mtdnaNames', sort=F)
consensus.mtdnas.all@phylo$tip.label<-nameswap$treeID

#mtdna protein only analysis
consensus.mtdnas.proteins<-read.beast.fixed(consensus.mtdna.proteins.path, ladderize=T)
#need to translate the names to those which are phylogenetically equivalent 
#in the nuclear dataset
nameswap<-merge(data.frame(mtdnaNames=consensus.mtdnas.proteins@phylo$tip.label), mtDNA.data, by='mtdnaNames', sort=F)
consensus.mtdnas.proteins@phylo$tip.label<-nameswap$treeID


#now the rRNAs 
consensus.mtdnas.rRNAs<-read.beast.fixed(consensus.mtdna.rRNAs.path, ladderize = T)
#need to translate the names to those which are phylogenetically equivalent 
#in the nuclear dataset
nameswap<-merge(data.frame(mtdnaNames=consensus.mtdnas.rRNAs@phylo$tip.label), mtDNA.data, by='mtdnaNames', sort=F)
consensus.mtdnas.rRNAs@phylo$tip.label<-nameswap$treeID

#presently not worrying about concordance on the mtDNA because 
#there are very few shifts, and therefore we would not expect any significant signal
}

#generate condordance summaries for each node relative to "all data" analysis
{
con_counts<-data.frame(exons=exons_concordance.allref[,2]*length(genematch.exons.allref[1,]), introns=introns_concordance.allref[,2]*length(genematch.introns.allref[1,]), utrs=utrs_concordance.allref[,2]*length(genematch.utrs.allref[1,]))
con_prop<-con_counts/length(genematch.all[1,])
#there seems to be one exon for which aves is not monophyletic, not sure how that can happen?

#calclate the discordance
con_prop$discordance <- 1- (con_prop$exons+con_prop$introns+con_prop$utrs)
con_prop$discordance.exons <- exons_concordance.allref$discordance
con_prop$discordance.introns <- introns_concordance.allref$discordance
con_prop$discordance.utrs <- utrs_concordance.allref$discordance


#visualize patterns of phylogenomic discordance across loci
#pairs(con_prop)

#pairs(con_prop)

#setting up the dataset for logistic regression
con_prop_logit<-con_prop
con_prop_logit$node.all<-seq(from=length(consensus.all@phylo$tip.label)+1, length.out=consensus.all@phylo$Nnode)

#make a time tree version of consensus.all
consensus.all.timetree<-make_timetree(ref = time_scale, tar = consensus.all@phylo)
#plot(consensus.all.timetree$phy, cex=0.001)
#nodelabels()
}

#collect the shift nodes from each dataset
{
shifts.allnuc<-as.data.frame(consensus.all@data)[,c(3,5)][seq(from=length(consensus.all@phylo$tip.label)+1, length.out=consensus.all@phylo$Nnode),]
shifts.allnuc$node.all<-match_phylo_nodes(consensus.all@phylo, consensus.all.timetree$phy)[,2]
shifts.allnuc<-shifts.allnuc[order(shifts.allnuc$node.all),]

shifts.exons<-as.data.frame(consensus.exons@data)[,c(3,5)][seq(from=length(consensus.exons@phylo$tip.label)+1, length.out=consensus.exons@phylo$Nnode),]
shifts.exons$node.all<-match_phylo_nodes(consensus.exons@phylo, consensus.all.timetree$phy)[,2]
shifts.exons<-shifts.exons[order(shifts.exons$node.all),]

shifts.introns<-as.data.frame(consensus.introns@data)[,c(3,5)][seq(from=length(consensus.introns@phylo$tip.label)+1, length.out=consensus.introns@phylo$Nnode),]
shifts.introns$node.all<-match_phylo_nodes(consensus.introns@phylo, consensus.all.timetree$phy)[,2]
shifts.introns<-shifts.introns[order(shifts.introns$node.all),]

shifts.utrs<-as.data.frame(consensus.utrs@data)[,c(3,5)][seq(from=length(consensus.utrs@phylo$tip.label)+1, length.out=consensus.utrs@phylo$Nnode),]
shifts.utrs$node.all<-match_phylo_nodes(consensus.utrs@phylo, consensus.all.timetree$phy)[,2]
shifts.utrs<-shifts.utrs[order(shifts.utrs$node.all),]


#mtdna
shifts.mtdnas.all<-as.data.frame(consensus.mtdnas.all@data)[,c(3,5)][seq(from=length(consensus.mtdnas.all@phylo$tip.label)+1, length.out=consensus.mtdnas.all@phylo$Nnode),]
shifts.mtdnas.all$node.all<-match_phylo_nodes(consensus.mtdnas.all@phylo, consensus.all.timetree$phy)[,2]
shifts.mtdnas.all<-shifts.mtdnas.all[order(shifts.mtdnas.all$node.all),]

shifts.mtdnas.proteins<-as.data.frame(consensus.mtdnas.proteins@data)[,c(3,5)][seq(from=length(consensus.mtdnas.proteins@phylo$tip.label)+1, length.out=consensus.mtdnas.proteins@phylo$Nnode),]
shifts.mtdnas.proteins$node.all<-match_phylo_nodes(consensus.mtdnas.proteins@phylo, consensus.all.timetree$phy)[,2]
shifts.mtdnas.proteins<-shifts.mtdnas.proteins[order(shifts.mtdnas.proteins$node.all),]

shifts.mtdnas.rRNAs<-as.data.frame(consensus.mtdnas.rRNAs@data)[,c(3,5)][seq(from=length(consensus.mtdnas.rRNAs@phylo$tip.label)+1, length.out=consensus.mtdnas.rRNAs@phylo$Nnode),]
shifts.mtdnas.rRNAs$node.all<-match_phylo_nodes(consensus.mtdnas.rRNAs@phylo, consensus.all.timetree$phy)[,2]
shifts.mtdnas.rRNAs<-shifts.mtdnas.rRNAs[order(shifts.mtdnas.rRNAs$node.all),]
#add back in missing data rows for those braches which are missing in the rRNA dataset
shifts.mtdnas.rRNAs<-merge(shifts.mtdnas.rRNAs, data.frame(node.all=shifts.mtdnas.proteins$node.all), by='node.all', all.y=T)
}


#checking uncex and unloc
as.data.frame(consensus.exons@data[!is.na(consensus.exons@data$uncex),])
as.data.frame(consensus.introns@data[!is.na(consensus.introns@data$uncex),])
as.data.frame(consensus.utrs@data[!is.na(consensus.utrs@data$uncex),])
as.data.frame(consensus.all@data[!is.na(consensus.all@data$uncex),])
as.data.frame(consensus.mtdnas.all@data[!is.na(consensus.mtdnas.all@data$uncex),])
as.data.frame(consensus.mtdnas.proteins@data[!is.na(consensus.mtdnas.proteins@data$uncex),])
as.data.frame(consensus.mtdnas.rRNAs@data[!is.na(consensus.mtdnas.rRNAs@data$uncex),])
#all cases have uncex and unloc at 100% model weight


#concordance processing
{
con_prop_logit<-cbind(con_prop_logit, 
                      uncex.allnudatc=shifts.allnuc$uncex,
                      uncex.exons=shifts.exons$uncex, 
                      uncex.introns=shifts.introns$uncex, 
                      uncex.utrs=shifts.utrs$uncex, 
                      uncex.mtdnas.all=shifts.mtdnas.all$uncex, 
                      uncex.mtdnas.proteins=shifts.mtdnas.proteins$uncex, 
                      uncex.mtdnas.rrnas=shifts.mtdnas.rRNAs$uncex)


#which nodes are considered in the janus estimation?
con_prop_logit$considered.nodes<-as.numeric(node_indices_N(tree=consensus.all.timetree$phy, min=4) >=4)


#create merged uncex column with merged uncex exons+introns+utrs
con_prop_logit$uncex.merged<-dplyr::coalesce(con_prop_logit$uncex.exons, con_prop_logit$uncex.introns, con_prop_logit$uncex.utrs)
con_prop_logit$uncex.merged.mtdnas<-dplyr::coalesce(con_prop_logit$uncex.merged, con_prop_logit$uncex.mtdnas.all, con_prop_logit$uncex.mtdnas.proteins, con_prop_logit$uncex.mtdnas.rrnas)
#con_prop_logit$uncex.merged.all<-dplyr::coalesce(con_prop_logit$uncex.allnucdat, con_prop_logit$uncex.merged, con_prop_logit$uncex.merged.mtdnas)

#cbind(con_prop_logit$node.all, con_prop_logit$uncex.merged, con_prop_logit$uncex.merged.mtdnas)
#there is one extra shift added by the mtDNA -- the shift on Neognathae

#con_prop_logit$node[con_prop_logit$uncex==1]
con_prop_logit[is.na(con_prop_logit)] <- 0
}

#notes
{
  # 
  # #plot concordance % on model shift tree
  # pdf(file="concordance_test.pdf", height=11, width=8.5)
  # #plot(compute.brlen(consensus.all@phylo, 1), cex=0.5, no.margin=T, edge.width=0.75)
  # 
  # plot_janus(treepath=consensus.all.path, ladderize=T, xlim=c(0, 0.2), width=1.5)
  # 
  # par(lwd = 0.2)
  # #nodelabels(pie=as.matrix(exons_concordance[,2:3]), piecol=c("black", "white"), cex=0.2, adj=c(0.725,0.5))
  # #nodelabels(pie=as.matrix(introns_concordance[,2:3]), piecol=c("black", "white"), cex=0.2, adj=c(0.5,0.5))
  # #nodelabels(pie=as.matrix(utrs_concordance[,2:3]), piecol=c("black", "white"), cex=0.2, adj=c(0.275,0.5))
  # #nodelabels(pie=as.matrix(exons_concordance[,2:3]), piecol=c("black", "white"), cex=0.2, adj=c(0.50,0.5))
  # #nodelabels(pie=as.matrix(introns_concordance[,2:3]), piecol=c("black", "white"), cex=0.2, adj=c(0.5,0.5))
  # #nodelabels(pie=as.matrix(utrs_concordance[,2:3]), piecol=c("black", "white"), cex=0.2, adj=c(0.5,0.5))
  # 
  # #labelplot(shifts=genematch.all, minscale=0.65)
  # nodelabels(pie=as.matrix(con_prop), piecol=c(viridis_pal()(3), "white"), cex=0.3)
  # 
  # dev.off()
  # 
  # 
  # #plot_janus(treepath=consensus.all.path, ladderize=T, xlim=c(0, 0.2), width=1.5)
  # #this looks OK
  # plot(consensus.all.timetree$phy, cex=0.0001)
  # nodelabels(pie=as.matrix(con_prop_logit$uncex.merged), piecol=c(viridis_pal()(3), "white"), cex=0.3)
  # nodelabels(pie=as.matrix(con_prop_logit$uncex.merged.mtdnas), piecol=c(viridis_pal()(3), "white"), cex=0.3)


########################################

  # 
  # #plot big comparison
  # 
  # setwd("~/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/unmasked/min2x/qual20_cov2_haplo_bestonly/initial_test_filters/min50bp_min10p_aligned/ALIGNED/phased/sample1haplo/phyhetnucbf/")
  # 
  # pdf(file="test.pdf", width=24, height=18)
  # 
  # par(mfrow=c(3,5))
  # #regular, estimate base freq shift
  # 
  # par(mar = c(1.5,0.25,1,1))
  # #plot an empty plot
  # plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  # #plot the label
  # text(x = 0.5, y = 0.5, "base freq shift",cex = 3.0, col = "black")
  # 
  # plot_janus(treepath="~/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/unmasked/min2x/qual20_cov2_haplo_bestonly/initial_test_filters/min50bp_min10p_aligned/ALIGNED/phased/sample1haplo/phyhetnucbf/backbone_constraint/janus/a65a5d46/bf/ALL_MFP_MERGE_MRL3/ALL_MFP_MERGE_MRL3_constraint.rooted.treefile.gophy.results.tre", ladderize=T, xlim=c(0, 0.2), width=1.5)
  # 
  # plot_janus(treepath="~/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/unmasked/min2x/qual20_cov2_haplo_bestonly/initial_test_filters/min50bp_min10p_aligned/ALIGNED/phased/sample1haplo/phyhetnucbf/backbone_constraint/janus/a65a5d46/bf/exons_MFP_MERGE_MRL3/exons_MFP_MERGE_MRL3_constraint.rooted.treefile.gophy.results.tre", ladderize=T, xlim=c(0, 0.15), width=1.5)
  # 
  # plot_janus(treepath="~/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/unmasked/min2x/qual20_cov2_haplo_bestonly/initial_test_filters/min50bp_min10p_aligned/ALIGNED/phased/sample1haplo/phyhetnucbf/backbone_constraint/janus/a65a5d46/bf/introns_MFP_MERGE_MRL3/introns_MFP_MERGE_MRL3_constraint.rooted.treefile.gophy.results.tre", ladderize=T, xlim=c(0, 0.3), width=1.5)
  # 
  # plot_janus(treepath="~/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/unmasked/min2x/qual20_cov2_haplo_bestonly/initial_test_filters/min50bp_min10p_aligned/ALIGNED/phased/sample1haplo/phyhetnucbf/backbone_constraint/janus/a65a5d46/bf/utrs_MFP_MERGE_MRL3/utrs_MFP_MERGE_MRL3_constraint.rooted.treefile.gophy.results.tre", ladderize=T, xlim=c(0, 0.3), width=1.5)
  # 
  # 
  # #regular, estimate base freq shift + rm shift
  # #par(mfrow=c(1,5))
  # #par(mar = c(1,0.25,1,0.25))
  # 
  # #plot an empty plot
  # plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  # #plot the label
  # text(x = 0.5, y = 0.5, "base freq +\n rate matrix shift",cex = 3.0, col = "black")
  # 
  # plot_janus(treepath="~/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/unmasked/min2x/qual20_cov2_haplo_bestonly/initial_test_filters/min50bp_min10p_aligned/ALIGNED/phased/sample1haplo/phyhetnucbf/backbone_constraint/janus/5d81a027/regular/ALL_MFP_MERGE_MRL3_constraint_RM_UE_UL_M4/ALL_MFP_MERGE_MRL3_constraint.rooted.treefile.gophy.results.tre", xlim=c(0,0.2), ladderize=T, width=1.5)
  # #title(main=paste("model shifts on", length(genematch.all[1,]), "loci"))
  # labelplot(shifts=genematch.all, minscale=0.65)
  # #labelplot.simple(shifts=genematch.all, minscale=0.1, matches.summary=matches.summary.all)
  # 
  # plot_janus(treepath="~/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/unmasked/min2x/qual20_cov2_haplo_bestonly/initial_test_filters/min50bp_min10p_aligned/ALIGNED/phased/sample1haplo/phyhetnucbf/backbone_constraint/janus/5d81a027/regular/exons_MFP_MERGE_MRL3_constraint_RM_UE_UL_M4/exons_MFP_MERGE_MRL3_constraint.rooted.treefile.gophy.results.tre", xlim=c(0,0.15), ladderize=T, width=1.5)
  # #title(main=paste("model shifts on", length(genematch.exons[1,]), "exons"))
  # #labelplot(shifts=genematch.exons.MFP, minscale=0.65, "red")
  # labelplot(shifts=genematch.exons, minscale=0.65, "black")
  # 
  # #labelplot.simple(shifts=genematch.exons, minscale=0.3, matches.summary=matches.summary.exons)
  # 
  # plot_janus(treepath="~/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/unmasked/min2x/qual20_cov2_haplo_bestonly/initial_test_filters/min50bp_min10p_aligned/ALIGNED/phased/sample1haplo/phyhetnucbf/backbone_constraint/janus/5d81a027/regular/introns_MFP_MERGE_MRL3_constraint_RM_UE_UL_M4/introns_MFP_MERGE_MRL3_constraint.rooted.treefile.gophy.results.tre", xlim=c(0,0.3), ladderize=T, width=1.5)
  # #title(main=paste("model shifts on", length(genematch.introns[1,]), "introns"))
  # labelplot(shifts=genematch.introns, minscale=0.65)
  # 
  # plot_janus(treepath="~/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/unmasked/min2x/qual20_cov2_haplo_bestonly/initial_test_filters/min50bp_min10p_aligned/ALIGNED/phased/sample1haplo/phyhetnucbf/backbone_constraint/janus/5d81a027/regular/utrs_MFP_MERGE_MRL3_constraint_RM_UE_UL_M4_stderror/utrs_MFP_MERGE_MRL3_constraint.rooted.treefile.gophy.results.tre", xlim=c(0,0.3), ladderize=T, width=1.5)
  # #title(main=paste("model shifts on", length(genematch.utrs[1,]), "UTRs"))
  # labelplot(shifts=genematch.utrs, minscale=0.65)
  # 
  # 
  # #branch lengths
  # #par(mar = c(0,0,0,0))
  # #par(mar = c(2,0.25,1,0.25))
  # #plot an empty plot
  # plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  # #plot the label
  # text(x = 0.5, y = 0.5, "base freq +\n rate matrix shift +\n re-estimate branch lengths",cex = 3.0, col = "black")
  # 
  # compare.phylo(t1=(consensus.all@phylo), t2=(consensus.all.bl@phylo), limit=c(0,0.2), width=1.5)
  # #plot_janus(treepath="~/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/unmasked/min2x/qual20_cov2_haplo_bestonly/initial_test_filters/min50bp_min10p_aligned/ALIGNED/phased/sample1haplo/phyhetnucbf/backbone_constraint/janus/bfrm_2_Dec_2020/branch_lengths/ALL_MFP_MERGE_MRL3_constraint_RM_UE_UL_M4/ALL_MFP_MERGE_MRL3_constraint.rooted.treefile.gophy.results.tre", xlim=c(0,0.2), ladderize=T, width=1.5)
  # 
  # compare.phylo(t1=(consensus.exons@phylo), t2=(consensus.exons.bl@phylo), limit=c(0,0.15), width=1.5)
  # #plot_janus(treepath="~/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/unmasked/min2x/qual20_cov2_haplo_bestonly/initial_test_filters/min50bp_min10p_aligned/ALIGNED/phased/sample1haplo/phyhetnucbf/backbone_constraint/janus/bfrm_2_Dec_2020/branch_lengths/exons_MFP_MERGE_MRL3_constraint_B_RM_UE_UL_M4/exons_MFP_MERGE_MRL3_constraint.rooted.treefile.gophy.results.tre", xlim=c(0,0.15), ladderize=T, width=1.5)
  # 
  # compare.phylo(t1=(consensus.introns@phylo), t2=(consensus.introns.bl@phylo), limit=c(0,0.3), width=1.5)
  # #plot_janus(treepath="~/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/unmasked/min2x/qual20_cov2_haplo_bestonly/initial_test_filters/min50bp_min10p_aligned/ALIGNED/phased/sample1haplo/phyhetnucbf/backbone_constraint/janus/bfrm_2_Dec_2020/branch_lengths/introns_MFP_MERGE_MRL3_constraint_B_RM_UE_UL_M4/introns_MFP_MERGE_MRL3_constraint.rooted.treefile.gophy.results.tre", xlim=c(0,0.3), ladderize=T, width=1.5)
  # 
  # compare.phylo(t1=(consensus.utrs@phylo), t2=(consensus.utrs.bl@phylo), limit=c(0,0.3), width=1.5)
  # #plot_janus(treepath="~/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/unmasked/min2x/qual20_cov2_haplo_bestonly/initial_test_filters/min50bp_min10p_aligned/ALIGNED/phased/sample1haplo/phyhetnucbf/backbone_constraint/janus/bfrm_2_Dec_2020/branch_lengths/utrs_MFP_MERGE_MRL3_constraint_B_RM_UE_UL_M4_stderror/utrs_MFP_MERGE_MRL3_constraint.rooted.treefile.gophy.results.tre", xlim=c(0,0.3), ladderize=T, width=1.5)
  # 
  # 
  # dev.off()

}

}

#intermediate plotting
#generating time scaled versions of the trees, 
#then painting molecular regimes
{
par(mfrow=c(1,6))
par(mar = c(1.5,0.25,1,1))
#plot_janus_time_rect(reference = time_scale, target = consensus.all)
#title(main="all nuclear data")
#labelplot(shifts=genematch.all, minscale=0.65, make.transparent("black", 0.75), size=0.75)

plot_janus_time_rect(reference = time_scale, target = consensus.exons)
title(main="exons")
#labelplot(shifts=genematch.exons, minscale=0.65, make.transparent("black", 0.75), size=0.75)

plot_janus_time_rect(reference = time_scale, target = consensus.introns)
title(main="introns")
#labelplot(shifts=genematch.introns, minscale=0.65, make.transparent("black", 0.75), size=0.75)

plot_janus_time_rect(reference = time_scale, target = consensus.utrs)
title(main="UTRs")
#labelplot(shifts=genematch.utrs, minscale=0.65, make.transparent("black", 0.75), size=0.75)

plot_janus_time_rect(reference = time_scale, target = consensus.mtdnas.all)
title(main="all mtDNA")

plot_janus_time_rect(reference = time_scale, target = consensus.mtdnas.proteins)
title(main="mtDNA proteins")

plot_janus_time_rect(reference = time_scale, target = consensus.mtdnas.rRNAs)
title(main="mtDNA rRNAs")
}

#Section 2
#more set up, data formatting etc
{
### new section
#generating estimates of GC content, codon usage bias, and nucleotide diversity
{
#getting tip labels and associated node numbers
#read in the data for the "all data" analysis
consensus.all.data<-model2tipdata(output_path=consensus.all.path)
colnames(consensus.all.data)[1]<-"all_nuc_models"

#read in the data for the exon only analysis
consensus.exons.data<-model2tipdata(output_path=consensus.exons.path)
colnames(consensus.exons.data)[1]<-"exon_models"

#read in the data for the intron only analysis
consensus.introns.data<-model2tipdata(output_path = consensus.introns.path)
colnames(consensus.introns.data)[1]<-"intron_models"

#read in the data for the utr only analysis
consensus.utrs.data<-model2tipdata(output_path = consensus.utrs.path)
colnames(consensus.utrs.data)[1]<-"utr_models"


#read in the data for the "all" mtdna analysis
consensus.mtdnas.all.data<-model2tipdata(output_path = consensus.mtdna.all.path)
colnames(consensus.mtdnas.all.data)[1]<-"mtdna_all_models"
#need to translate the names to those which are phylogenetically equivalent 
#in the nuclear dataset
nameswap<-merge(data.frame(mtdnaNames=consensus.mtdnas.all.data$label), mtDNA.data, by='mtdnaNames', sort=F)
consensus.mtdnas.all.data$label<-nameswap$treeID

#read in the data for the "protein" mtdna analysis
consensus.mtdnas.proteins.data<-model2tipdata(output_path = consensus.mtdna.proteins.path)
colnames(consensus.mtdnas.proteins.data)[1]<-"mtdna_proteins_models"
#need to translate the names to those which are phylogenetically equivalent 
#in the nuclear dataset
nameswap<-merge(data.frame(mtdnaNames=consensus.mtdnas.proteins.data$label), mtDNA.data, by='mtdnaNames', sort=F)
consensus.mtdnas.proteins.data$label<-nameswap$treeID

#read in the data for the "rRNA" mtdna analysis
consensus.mtdnas.rRNAs.data<-model2tipdata(output_path = consensus.mtdna.rRNAs.path)
colnames(consensus.mtdnas.rRNAs.data)[1]<-"mtdna_rRNAs_models"
#need to translate the names to those which are phylogenetically equivalent 
#in the nuclear dataset
nameswap<-merge(data.frame(mtdnaNames=consensus.mtdnas.rRNAs.data$label), mtDNA.data, by='mtdnaNames', sort=F)
consensus.mtdnas.rRNAs.data$label<-nameswap$treeID

#merge(x= consensus.mtdnas.rRNAs.data,y= ,  by.x = "node", by.y="node.all") 
}

#estimating diversity from total branch length between alleles
{
two_haplo<-readrooter(treepath = allele.tree.path, outgroup = c("Caiman_croccodilus_1", "Caiman_croccodilus_2", "Crocodylus_porosus_1", "Crocodylus_porosus_2"), drop=T)
two_haplo.dist<-allele_dist(two_haplo)
}
#AUGUST 1## the merging below may not be valid if the node numbers across the consensus trees are not identical

#create new data matrix representing the model IDs and 
#start building up the tip dataset
#add the allele distances
{
tip.data <- merge(x=consensus.all.data[,c(1,6)], y=two_haplo.dist, by="label")

#add in the exon model shifts
tip.data <- merge(consensus.exons.data[,c(1,6)], y = tip.data, by="label")

#add in the intron model shifts
tip.data <- merge(consensus.introns.data[,c(1,6)], y = tip.data, by="label")

#add in the utr model shifts
tip.data <- merge(consensus.utrs.data[,c(1,6)], y = tip.data, by="label")


#add in the mtdna "all" model shifts
tip.data <- merge(consensus.mtdnas.all.data[,c(1,6)], y = tip.data, by="label")

#add in the mtdna "protein" model shifts
tip.data <- merge(consensus.mtdnas.proteins.data[,c(1,6)], y = tip.data, by="label")


#add in the mtdna "rRNA" model shifts
tip.data <- merge(consensus.mtdnas.rRNAs.data[,c(1,6)], y = tip.data, by="label", all.y=T)

#fill in the missing models for the rRNAs assuming phylogenetic equivalence+signal
#it is identical to the protein model set, so just overwrite it
tip.data$mtdna_rRNAs_models<-tip.data$mtdna_proteins_models
#[is.na(tip.data$mtdna_rRNAs_models)]<-c(0,1,1,1,1,1,1,1,1)
}

#convert model shifts to characters for downstream anova
{
tip.data$all_nuc_models<-as.character(consensus.all.data$all_nuc_models)
tip.data$exon_models<-as.character(tip.data$exon_models)
tip.data$intron_models<-as.character(tip.data$intron_models)
tip.data$utr_models<-as.character(tip.data$utr_models)
tip.data$mtdna_all_models<-as.character(tip.data$mtdna_all_models)
tip.data$mtdna_proteins_models<-as.character(tip.data$mtdna_proteins_models)
tip.data$mtdna_rRNAs_models<-as.character(tip.data$mtdna_rRNAs_models)

#add rownames
rownames(tip.data)<-tip.data$label

#reorder to match tree order
tip.data <- tip.data[consensus.all.timetree$phy$tip.label,]
}
}

#Section 3
########################################################
### generating new datasets for logistic regression ####
########################################################
{
#try phylogenetic logistic regression with each of the model shift data types
#is phylogemic discordance related to the identification of a model shift?
{
logisticreg_tree<-consensus.all.timetree$phy
#logisticreg_tree$root.edge<-0.1
#nodelabels()
#plot(logisticreg_tree, cex=0.1, no.margin=T)

#graft in a zero length branch to each node (now set to 1 Ma)
graftLength<-1
logisticreg_tree<-grafter(logisticreg_tree, edge.length = graftLength, position = 0.0)

#fixing the graft
logisticreg_tree<-multi2di(logisticreg_tree, random=F)
logisticreg_tree$edge.length[logisticreg_tree$edge.length==0]<-graftLength
#now we have a tree with grafted tips like fossils

logisticreg.newdat<-data.frame(node=logisticreg_tree$tip.label, label=logisticreg_tree$tip.label)

#add the discordance states
logisticreg.newdat<-merge(x=con_prop_logit, y=logisticreg.newdat, by.x="node.all", by.y="node", all.y=T)
#is.unsorted(con_prop_logit$node.all)

#add row labels
rownames(logisticreg.newdat)<-logisticreg.newdat$label

#reorder to match tip label order
logisticreg.newdat<-logisticreg.newdat[logisticreg_tree$tip.label,]

#replace NAs with zeros
logisticreg.newdat[is.na(logisticreg.newdat)]<-0


#add in node ages
ages<-paleotree::dateNodes(consensus.all.timetree$phy, labelDates = F)
ages<-c(setNames(ages[1:198], consensus.all.timetree$phy$tip.label), ages[199:395])
logisticreg.newdat$ages<-ages[logisticreg.newdat$label]
#create alternative version with tips set to 0.0001 years
logisticreg.newdat$ages.nozero<-logisticreg.newdat$ages+0.0001
#logisticreg.newdat$ages.nozero[logisticreg.newdat$ages.nozero==0]<-0.0001

#same for discordance
logisticreg.newdat$discordance.nozeros<-logisticreg.newdat$discordance+0.0001
#logisticreg.newdat$discordance.nozeros[logisticreg.newdat$discordance.nozeros==0]<-0.0001

#is a node within 10 Ma of the K-Pg
logisticreg.newdat$ages.thresh10<-logisticreg.newdat$ages
logisticreg.newdat$ages.thresh10[logisticreg.newdat$ages >= 56 & logisticreg.newdat$ages <= 76] <- "yes"
logisticreg.newdat$ages.thresh10[logisticreg.newdat$ages <= 56 | logisticreg.newdat$ages >= 76] <- "no"

#is a node within 5 Ma of the K-Pg
logisticreg.newdat$ages.thresh5<-logisticreg.newdat$ages
logisticreg.newdat$ages.thresh5[logisticreg.newdat$ages >= 61 & logisticreg.newdat$ages <= 71] <- "yes"
logisticreg.newdat$ages.thresh5[logisticreg.newdat$ages <= 61 | logisticreg.newdat$ages >= 71] <- "no"

#time since KPg
logisticreg.newdat$ages.age_since<-abs(logisticreg.newdat$ages-66)

#time since Kpg, but stem ages
logisticreg.newdat$stem.ages<-abs(stem_age(consensus.all.timetree$phy, all=T)[logisticreg.newdat$label])
logisticreg.newdat$stem.ages_since<-abs(stem_age(consensus.all.timetree$phy, all=T)[logisticreg.newdat$label]-66)

#is a stem node within 10 Ma of the K-Pg
logisticreg.newdat$stem.ages.thresh10<-logisticreg.newdat$stem.ages
logisticreg.newdat$stem.ages.thresh10[logisticreg.newdat$stem.ages >= 56 & logisticreg.newdat$stem.ages <= 76] <- "yes"
logisticreg.newdat$stem.ages.thresh10[logisticreg.newdat$stem.ages <= 56 | logisticreg.newdat$stem.ages >= 76] <- "no"

#is a node within 5 Ma of the K-Pg
logisticreg.newdat$stem.ages.thresh5<-logisticreg.newdat$stem.ages
logisticreg.newdat$stem.ages.thresh5[logisticreg.newdat$stem.ages >= 61 & logisticreg.newdat$stem.ages <= 71] <- "yes"
logisticreg.newdat$stem.ages.thresh5[logisticreg.newdat$stem.ages <= 61 | logisticreg.newdat$stem.ages >= 71] <- "no"


#create alternative versions which exclude the contemporary terminals
logisticreg_tree.alt<-drop.tip(logisticreg_tree, tip=consensus.all.timetree$phy$tip.label)
logisticreg.newdat.alt<-logisticreg.newdat[!rownames(logisticreg.newdat) %in% consensus.all.timetree$phy$tip.label,]

#save.image("~/jsb439@cornell.edu/Code/avian_molecular_shifts/workspace_Feb_3_2022.RData")
}

########################################################
#### end section datasets for logistic regression #####
########################################################

#notes
{

# 
# test<-phylolm(log(ages.age_since) ~ as.factor(uncex.merged), phy = (logisticreg_tree), data = logisticreg.newdat)
# plot(log(ages.age_since) ~ as.factor(uncex.merged), data= logisticreg.newdat)
# 
# test2<-phylolm(log(ages.age_since) ~ 1, phy = (logisticreg_tree), data = logisticreg.newdat)
# plot(log(ages.age_since) ~ 1, data= logisticreg.newdat)
# 
# test   #AIC 787.6
# test2  #AIC 799.7
# 10 AIC units improvement


#pageltest<-fitPagel(logisticreg_tree, x = setNames(logisticreg.newdat$ages.thresh10, logisticreg.newdat$label) , y= setNames(as.character(logisticreg.newdat$uncex.merged), logisticreg.newdat$label))
#pageltest2<-fitPagel(logisticreg_tree, x = setNames(logisticreg.newdat$ages.thresh5, logisticreg.newdat$label) , y= setNames(as.character(logisticreg.newdat$uncex.merged), logisticreg.newdat$label))

# tmp<-(logisticreg.newdat.alt[,c("label", "ages.thresh10","uncex.merged" )])
# tmp$uncex.merged<-as.character(tmp$uncex.merged)
# tmp2<-corDISC(phy = logisticreg_tree.alt, data = tmp, ntraits=2, model="ARD", node.states="marginal", diagn=FALSE)
# tmp2$states


# phylo.logisticreg.fit<- phyloglm(uncex.merged~(log(ages.age_since))+log(discordance), data=logisticreg.newdat.alt, phy=(logisticreg_tree.alt), boot=0, "logistic_IG10", btol=1000, log.alpha.bound=6)
# summary(phylo.logisticreg.fit)
#library(future); plan(multiprocess)

#using age since (MRCA)
# 
# #log transformations
# phylo.logisticreg.fit.ages<- phyloglm(uncex.merged~(log(ages.age_since)), data=logisticreg.newdat, phy=(logisticreg_tree), method="logistic_IG10", boot=0)
# summary(phylo.logisticreg.fit.ages)
# 
# phylo.logisticreg.fit.ages<- phyloglm(uncex.merged~(log(ages.age_since))+log(discordance.nozeros), data=logisticreg.newdat, phy=(logisticreg_tree), method="logistic_IG10", boot=0)
# summary(phylo.logisticreg.fit.ages)
# 
# #arcsin sqrt transformation for discordance
# phylo.logisticreg.fit.ages<- phyloglm(uncex.merged~(log(ages.age_since))+asinTransform(discordance), data=logisticreg.newdat, phy=(logisticreg_tree), method="logistic_IG10")
# summary(phylo.logisticreg.fit.ages)
# 
# #using age since (stem ages)
# 
# #log transformations
# 
# phylo.logisticreg.fit.stem.ages<- phyloglm(uncex.merged~(log(stem.ages_since)), data=logisticreg.newdat, phy=(logisticreg_tree.alt), method="logistic_IG10", btol=16)
# summary(phylo.logisticreg.fit.stem.ages)
# 
# 
# #models with discordance
# phylo.logisticreg.fit.stem.ages<- phyloglm(uncex.merged~(log(stem.ages_since))+log(discordance.nozeros), data=logisticreg.newdat, phy=(logisticreg_tree), method="logistic_IG10", btol=16)
# summary(phylo.logisticreg.fit.stem.ages)
# 
# #arcsing sqrt transformation for discordance
# phylo.logisticreg.fit.stem.ages<- phyloglm(uncex.merged~(log(stem.ages_since))+asinTransform(discordance), data=logisticreg.newdat, phy=(logisticreg_tree), method="logistic_IG10", btol=16)
# summary(phylo.logisticreg.fit.stem.ages)
}
}

#Section 4
########################################################
###########running logistic regression #################
########################################################
{
#individual models for each data type
library(future); plan(multisession, workers=10)
#library(future); plan(sequential)
#source("/Users/cotinga/Downloads/phylolm-master\ copy/R/phyloglm.R")
#this custom version doesn't work now for some reason-- using boot mean values

#run logistic regression
{
#phylo.logisticreg.fit.stem.ages.allnucdata <- phyloglm(uncex.allnudatc~(log(stem.ages_since)), data=logisticreg.newdat, phy=(logisticreg_tree.alt), boot = 1000)
phylo.logisticreg.fit.stem.ages.exons <- phyloglm(uncex.exons~(log(stem.ages_since)), data=logisticreg.newdat, phy=(logisticreg_tree.alt), boot = 1000, start.alpha=0.1, log.alpha.bound=4)
phylo.logisticreg.fit.stem.ages.introns <- phyloglm(uncex.introns~(log(stem.ages_since)), data=logisticreg.newdat, phy=(logisticreg_tree.alt), boot=1000, start.alpha=0.15, log.alpha.bound=4)#, log.alpha.bound=4, boot = 1000)
phylo.logisticreg.fit.stem.ages.utrs<- phyloglm(uncex.utrs~(log(stem.ages_since)), data=logisticreg.newdat, phy=(logisticreg_tree.alt), boot=1000, start.alpha=0.1)#, start.alpha=0.2, log.alpha.bound=3, boot = 1000)
phylo.logisticreg.fit.stem.ages.mtdnas.all<- phyloglm(uncex.mtdnas.all~(log(stem.ages_since)), data=logisticreg.newdat, phy=(logisticreg_tree.alt), boot = 1000)
phylo.logisticreg.fit.stem.ages.mtdnas.proteins<- phyloglm(uncex.mtdnas.proteins~(log(stem.ages_since)), data=logisticreg.newdat, phy=(logisticreg_tree.alt), boot=1000)#, start.alpha=0.2, boot = 1000)
phylo.logisticreg.fit.stem.ages.mtdnas.rrnas<- phyloglm(uncex.mtdnas.rrnas~(log(stem.ages_since)), data=logisticreg.newdat, phy=(logisticreg_tree.alt), boot=1000)#, start.alpha=0.2, boot = 1000)
phylo.logisticreg.fit.stem.ages.merged<- phyloglm(uncex.merged~(log(stem.ages_since)), data=logisticreg.newdat, phy=(logisticreg_tree.alt), boot = 1000)
phylo.logisticreg.fit.stem.ages.merged.mtdnas<- phyloglm(uncex.merged.mtdnas~(log(stem.ages_since)), data=logisticreg.newdat, phy=(logisticreg_tree.alt), boot = 1000)

#individual models for each data type + discordance
#phylo.logisticreg.fit.stem.ages.allnucdata.discordance <- phyloglm(uncex.allnudatc~(log(stem.ages_since))+(asinTransform(discordance)), data=logisticreg.newdat, phy=(logisticreg_tree.alt), start.alpha=0.001, log.alpha.bound = 10, boot = 100)
phylo.logisticreg.fit.stem.ages.exons.discordance  <- phyloglm(uncex.exons~(log(stem.ages_since))+(asinTransform(discordance.exons)), data=logisticreg.newdat, phy=(logisticreg_tree.alt), boot=1000, start.alpha=0.1, log.alpha.bound = 4.5)
#phylo.logisticreg.fit.stem.ages.introns.discordance  <- phyloglm(uncex.introns~(log(stem.ages_since))+(asinTransform(discordance.introns)), data=logisticreg.newdat, phy=(logisticreg_tree.alt), boot=1000, start.alpha=1, log.alpha.bound = 10)
phylo.logisticreg.fit.stem.ages.introns.discordance  <- phyloglm(uncex.introns~(log(stem.ages_since))+(asinTransform(discordance.introns)), data=logisticreg.newdat, phy=(logisticreg_tree.alt), boot=1000)
phylo.logisticreg.fit.stem.ages.utrs.discordance <- phyloglm(uncex.utrs~(log(stem.ages_since))+(asinTransform(discordance.utrs)), data=logisticreg.newdat, phy=(logisticreg_tree.alt), boot = 1000)#, start.alpha=0.2, log.alpha.bound = 5)
phylo.logisticreg.fit.stem.ages.merged.discordance <- phyloglm(uncex.merged~(log(stem.ages_since))+(asinTransform(discordance)), data=logisticreg.newdat, phy=(logisticreg_tree.alt), boot=1000)#, start.alpha=0.1, log.alpha.bound = 4.5)


#min(fitted.values(phylo.logisticreg.fit.stem.ages.exons.discordance))
#max(fitted.values(phylo.logisticreg.fit.stem.ages.exons.discordance))


# 
# phylo.logisticreg.fit.stem.ages.merged.mtdnas.discordance <- phyloglm(uncex.merged.mtdnas~scale(log(stem.ages_since))*scale(asinTransform(discordance)), data=logisticreg.newdat, phy=(logisticreg_tree.alt), boot = 100)
# summary(phylo.logisticreg.fit.stem.ages.merged.mtdnas.discordance)
# 
# test <- glm(uncex.merged.mtdnas~(log(stem.ages_since))+asinTransform(discordance), data=logisticreg.newdat.alt, family="binomial", correlation)
# summary(test)
# LogisticDx:::dx(test)
# plot_logistic_curve(test)
# 
# plot(test$residuals~test$fitted.values)
# length(test$fitted.values)
# length(logisticreg.newdat.alt$uncex.merged.mtdnas[-197])
# 
# 
}

# #run logistic regression (assuming fixed alpha) -- commented out
# {
#   #phylo.logisticreg.fit.stem.ages.allnucdata <- phyloglm(uncex.allnudatc~(log(stem.ages_since)), data=logisticreg.newdat, phy=(logisticreg_tree.alt), boot = 1000)
#   phylo.logisticreg.fit.stem.ages.exons <- phyloglm(uncex.exons~(log(stem.ages_since)), data=logisticreg.newdat, phy=(logisticreg_tree.alt), boot = 1000, log.alpha.bound=0)#start.alpha=0.1)
#   phylo.logisticreg.fit.stem.ages.introns <- phyloglm(uncex.introns~(log(stem.ages_since)), data=logisticreg.newdat, phy=(logisticreg_tree.alt), boot=1000, log.alpha.bound=0)#start.alpha=0.15, log.alpha.bound=4)#, log.alpha.bound=4, boot = 1000)
#   phylo.logisticreg.fit.stem.ages.utrs<- phyloglm(uncex.utrs~(log(stem.ages_since)), data=logisticreg.newdat, phy=(logisticreg_tree.alt), boot=1000, log.alpha.bound=0)#, start.alpha=0.1)#, start.alpha=0.2, log.alpha.bound=3, boot = 1000)
#   phylo.logisticreg.fit.stem.ages.mtdnas.all<- phyloglm(uncex.mtdnas.all~(log(stem.ages_since)), data=logisticreg.newdat, phy=(logisticreg_tree.alt), boot = 1000, log.alpha.bound=0)
#   phylo.logisticreg.fit.stem.ages.mtdnas.proteins<- phyloglm(uncex.mtdnas.proteins~(log(stem.ages_since)), data=logisticreg.newdat, phy=(logisticreg_tree.alt), boot=1000, log.alpha.bound=0)#, start.alpha=0.2, boot = 1000)
#   phylo.logisticreg.fit.stem.ages.mtdnas.rrnas<- phyloglm(uncex.mtdnas.rrnas~(log(stem.ages_since)), data=logisticreg.newdat, phy=(logisticreg_tree.alt), boot=1000, log.alpha.bound=0)#, start.alpha=0.2, boot = 1000)
#   phylo.logisticreg.fit.stem.ages.merged<- phyloglm(uncex.merged~(log(stem.ages_since)), data=logisticreg.newdat, phy=(logisticreg_tree.alt), boot = 1000, log.alpha.bound=0)
#   phylo.logisticreg.fit.stem.ages.merged.mtdnas<- phyloglm(uncex.merged.mtdnas~(log(stem.ages_since)), data=logisticreg.newdat, phy=(logisticreg_tree.alt), boot = 1000, log.alpha.bound=0)
#   
#   #individual models for each data type + discordance
#   #phylo.logisticreg.fit.stem.ages.allnucdata.discordance <- phyloglm(uncex.allnudatc~(log(stem.ages_since))+(asinTransform(discordance)), data=logisticreg.newdat, phy=(logisticreg_tree.alt), start.alpha=0.001, log.alpha.bound = 10, boot = 100)
#   phylo.logisticreg.fit.stem.ages.exons.discordance  <- phyloglm(uncex.exons~(log(stem.ages_since))+(asinTransform(discordance.exons)), data=logisticreg.newdat, phy=(logisticreg_tree.alt), boot=1000, log.alpha.bound=0)#, start.alpha=0.1, log.alpha.bound = 4.5)
#   #phylo.logisticreg.fit.stem.ages.introns.discordance  <- phyloglm(uncex.introns~(log(stem.ages_since))+(asinTransform(discordance.introns)), data=logisticreg.newdat, phy=(logisticreg_tree.alt), boot=1000, start.alpha=1, log.alpha.bound = 10)
#   phylo.logisticreg.fit.stem.ages.introns.discordance  <- phyloglm(uncex.introns~(log(stem.ages_since))+(asinTransform(discordance.introns)), data=logisticreg.newdat, phy=(logisticreg_tree.alt), log.alpha.bound=0)#, start.alpha=0.1, boot=1000, log.alpha.bound = 4.5)
#   phylo.logisticreg.fit.stem.ages.utrs.discordance <- phyloglm(uncex.utrs~(log(stem.ages_since))+(asinTransform(discordance.utrs)), data=logisticreg.newdat, phy=(logisticreg_tree.alt), boot = 1000, log.alpha.bound=0)#, start.alpha=0.2, log.alpha.bound = 5)
#   phylo.logisticreg.fit.stem.ages.merged.discordance <- phyloglm(uncex.merged~(log(stem.ages_since))+(asinTransform(discordance)), data=logisticreg.newdat, phy=(logisticreg_tree.alt), boot=1000, log.alpha.bound=0)# start.alpha=0.1, log.alpha.bound = 4.5)
#   
#   #min(fitted.values(phylo.logisticreg.fit.stem.ages.exons.discordance))
#   #max(fitted.values(phylo.logisticreg.fit.stem.ages.exons.discordance))
#   
#   # 
#   # phylo.logisticreg.fit.stem.ages.merged.mtdnas.discordance <- phyloglm(uncex.merged.mtdnas~scale(log(stem.ages_since))*scale(asinTransform(discordance)), data=logisticreg.newdat, phy=(logisticreg_tree.alt), boot = 100)
#   # summary(phylo.logisticreg.fit.stem.ages.merged.mtdnas.discordance)
#   # 
#   # test <- glm(uncex.merged.mtdnas~(log(stem.ages_since))+asinTransform(discordance), data=logisticreg.newdat.alt, family="binomial", correlation)
#   # summary(test)
#   # LogisticDx:::dx(test)
#   # plot_logistic_curve(test)
#   # 
#   # plot(test$residuals~test$fitted.values)
#   # length(test$fitted.values)
#   # length(logisticreg.newdat.alt$uncex.merged.mtdnas[-197])
#   # 
#   # 
# }
  
#load some function from this mrhelmus github repo
#calculate AICc
{
devtools:::source_url("https://raw.githubusercontent.com/mrhelmus/phylogeny_manipulation/master/AIC_func.r")
AICc.phyloglm(phylo.logisticreg.fit.stem.ages.merged.discordance) #69.36564
AICc.phyloglm(phylo.logisticreg.fit.stem.ages.merged) #67.90789
#dAIC < 1
aictab(list(phylo.logisticreg.fit.stem.ages.merged.discordance,phylo.logisticreg.fit.stem.ages.merged ))


AICc.phyloglm(phylo.logisticreg.fit.stem.ages.exons.discordance) #61.27305
AICc.phyloglm(phylo.logisticreg.fit.stem.ages.exons) #61.15363
#dAIC ~ NA
aictab(list(phylo.logisticreg.fit.stem.ages.exons.discordance,phylo.logisticreg.fit.stem.ages.exons ))


AICc.phyloglm(phylo.logisticreg.fit.stem.ages.introns.discordance) #44.31081
AICc.phyloglm(phylo.logisticreg.fit.stem.ages.introns) #39.28601
#dAIC = 5.51
aictab(list(phylo.logisticreg.fit.stem.ages.introns.discordance,phylo.logisticreg.fit.stem.ages.introns ))


AICc.phyloglm(phylo.logisticreg.fit.stem.ages.utrs.discordance) #31.15809
AICc.phyloglm(phylo.logisticreg.fit.stem.ages.utrs) #29.99632
#dAIC < 1
#in all cases, models excluding discordance have lower AICc
aictab(list(phylo.logisticreg.fit.stem.ages.utrs.discordance,phylo.logisticreg.fit.stem.ages.utrs ))


}



#what is the average VIF for predictors across
#models that exclude mtDNA
{
mean(vif(phylo.logisticreg.fit.stem.ages.exons.discordance)[1],
vif(phylo.logisticreg.fit.stem.ages.introns.discordance)[1],
vif(phylo.logisticreg.fit.stem.ages.utrs.discordance)[1],
vif(phylo.logisticreg.fit.stem.ages.merged.discordance)[1])
}

#what is the relationship between time and discordance?
model.time.discordance<-phylolm(log(stem.ages_since)~asinTransform(discordance), data=logisticreg.newdat, phy=(logisticreg_tree.alt), boot = 100)
summary(model.time.discordance)
#not significant after phylogeny

model.time.discordance.lm<-lm(log(stem.ages_since)~asinTransform(discordance), data=logisticreg.newdat)
summary(model.time.discordance.lm)
#significant before phylogeny

#plot(model.time.discordance)
#plot(log(stem.ages_since)~asinTransform(discordance), data=logisticreg.newdat)
#abline(model.time.discordance)

#notes -- trying rr2
{

#individual models for each data type + discordance + considered nodes
    # phylo.logisticreg.fit.stem.ages.allnucdata<- phyloglm(uncex.allnudatc~(log(stem.ages_since))+asinTransform(discordance)+as.factor(considered.nodes), data=logisticreg.newdat, phy=(logisticreg_tree.alt), start.alpha=0.1, btol=10)
    # phylo.logisticreg.fit.stem.ages.exons <- phyloglm(uncex.exons~(log(stem.ages_since))+asinTransform(discordance)+as.factor(considered.nodes), data=logisticreg.newdat, phy=(logisticreg_tree.alt))
    # phylo.logisticreg.fit.stem.ages.introns <- phyloglm(uncex.introns~(log(stem.ages_since))+asinTransform(discordance)+as.factor(considered.nodes), data=logisticreg.newdat, phy=(logisticreg_tree.alt))
    # phylo.logisticreg.fit.stem.ages.utrs<- phyloglm(uncex.utrs~(log(stem.ages_since))+asinTransform(discordance)+as.factor(considered.nodes), data=logisticreg.newdat, phy=(logisticreg_tree.alt))
    # phylo.logisticreg.fit.stem.ages.mtdnas.all<- phyloglm(uncex.mtdnas.all~(log(stem.ages_since))+asinTransform(discordance)+as.factor(considered.nodes), data=logisticreg.newdat, phy=(logisticreg_tree.alt))
    # phylo.logisticreg.fit.stem.ages.mtdnas.proteins<- phyloglm(uncex.mtdnas.proteins~(log(stem.ages_since))+asinTransform(discordance)+as.factor(considered.nodes), data=logisticreg.newdat, phy=(logisticreg_tree.alt))
    # phylo.logisticreg.fit.stem.ages.mtdnas.rrnas<- phyloglm(uncex.mtdnas.rrnas~(log(stem.ages_since))+asinTransform(discordance)+as.factor(considered.nodes), data=logisticreg.newdat, phy=(logisticreg_tree.alt))
    # phylo.logisticreg.fit.stem.ages.merged<- phyloglm(uncex.merged~(log(stem.ages_since))+asinTransform(discordance)+as.factor(considered.nodes), data=logisticreg.newdat, phy=(logisticreg_tree.alt))
    # phylo.logisticreg.fit.stem.ages.merged.mtdnas<- phyloglm(uncex.merged.mtdnas~(log(stem.ages_since))+asinTransform(discordance)+as.factor(considered.nodes), data=logisticreg.newdat, phy=(logisticreg_tree.alt))


# 
# phylo.logisticreg.fit.null <- phyloglm(uncex.merged~1, data=logisticreg.newdat, phy=(logisticreg_tree.alt), method="logistic_IG10")
# #phylo.logisticreg.fit.null <- phyloglm(uncex.merged~1, data=logisticreg.newdat, phy=(logisticreg_tree), btol=10, start.alpha=0.1, log.alpha.bound=5)
# summary(phylo.logisticreg.fit.null)
# 
# z.f.plog2 <- phylo.logisticreg.fit.stem.ages
# z.x.plog2 <- phylo.logisticreg.fit.null
# 
# # R2.resid and R2.pred do not apply for phyloglm
# require(rr2)
# R2(z.f.plog2, z.x.plog2)

#revisit this 
#https://github.com/cran/rr2
}
}

#plotting for supp figur 1

#picking values for nodevalues
#setting up params for supp figure 1

#setting up alternative values for X2 in the logistic regressions with discordance
{
X2s.exons <- (asinTransform(logisticreg.newdat[logisticreg_tree.alt$tip.label,]$discordance.exons))
X2s.introns <- (asinTransform(logisticreg.newdat[logisticreg_tree.alt$tip.label,]$discordance.introns))
X2s.utrs <- (asinTransform(logisticreg.newdat[logisticreg_tree.alt$tip.label,]$discordance.utrs))
X2s.merged <- (asinTransform(logisticreg.newdat[logisticreg_tree.alt$tip.label,]$discordance))

X2_l.exon<- mean(X2s.exons) - sd(X2s.exons)#*2
X2_l.intron<- mean(X2s.introns) - sd(X2s.introns)#*2
X2_l.utr<- mean(X2s.utrs) - sd(X2s.utrs)#*2
X2_l.merged<- mean(X2s.merged) - sd(X2s.merged)#*2
mean(c(sin(X2_l.exon)^2, sin(X2_l.intron)^2, sin(X2_l.utr)^2, sin(X2_l.merged)^2))

X2_m.exon <- mean(X2s.exons)
X2_m.intron <- mean(X2s.introns)
X2_m.utr <- mean(X2s.utrs)
X2_m.merged <- mean(X2s.merged)
mean(c(sin(X2_m.exon)^2, sin(X2_m.intron)^2, sin(X2_m.utr)^2, sin(X2_m.merged)^2))

X2_h.exon<- mean(X2s.exons) + sd(X2s.exons)#*2
X2_h.intron<- mean(X2s.introns) + sd(X2s.introns)#*2
X2_h.utr<- mean(X2s.utrs) + sd(X2s.utrs)#*2
X2_h.merged<- mean(X2s.merged) + sd(X2s.merged)#*2
mean(c(sin(X2_h.exon)^2, sin(X2_h.intron)^2, sin(X2_h.utr)^2, sin(X2_h.merged)^2))

}

#these parameters control the confidence interval and the spline smoothing param
int=0.7
smooth=0.95

#set up base plot with points
baseplot<-function(title="Model Shift Binomial pGLM", gradient=F){
  # create a color ramp function that goes from white to black
  color_ramp <- colorRampPalette(c("gray100","gray90","grey", "black"), bias=0.5)
  
  # normalize the values in logisticreg.newdat.alt$discordance to lie between 0 and 1
  normalized_discordance <- asinTransform(logisticreg.newdat.alt$discordance) / asinTransform(max(logisticreg.newdat.alt$discordance))
  
  normalized_discordance <- round(normalized_discordance, 2)*100
  names(normalized_discordance) <- rownames(logisticreg.newdat.alt)
  
  # apply the color ramp to the normalized values to get a corresponding color for each point
  point_colors <- color_ramp(100)
  
  # generate the scatterplot
  plot(x=log(logisticreg.newdat$stem.ages_since), y=jitter(logisticreg.newdat$uncex.merged.mtdnas,factor=0,amount=0.02), 
       xlim=c(-3,4.2), ylim=c(-0.05,1.05),
       xlab="log(time) to T=66 Ma", ylab="probability", 
       main=title, pch=21, col=make.transparent("gray", 0.0), 
       bty='n', cex=0.000001, cex.main=1.0)
  
  # # add points to the scatterplot with the colors determined by the discordance values
  # points(x=log(logisticreg.newdat$stem.ages_since), y=jitter(logisticreg.newdat$uncex.merged.mtdnas,factor=0,amount=0.075), 
  #        pch=21, cex=1.0, bg=scales::alpha(point_colors[normalized_discordance], alpha=0.8), lwd=0.5)
  # 
  if (gradient==T){
    points(x=log(logisticreg.newdat$stem.ages_since), y=jitter(logisticreg.newdat$uncex.merged.mtdnas,factor=0,amount=0.05), 
           pch=21, cex=0.9, bg=scales::alpha(point_colors[normalized_discordance], alpha=0.8), lwd=0.5)
    {# First, create a vector of values that covers the full range of the data
      values <- seq(min(normalized_discordance), max(normalized_discordance), length.out=100)
      
      # Next, create a matrix of colors corresponding to the values in the 'values' vector
      colors <- color_ramp(length(values))
      
      # legend.scale(zlim=range(values), col = colors, horizontal = F)
      
      # Now, create the legend
      par(lend=2)
      legend("topright", pch="",
             legend=values, text.col = scales::alpha("white",0),
             col = rev(colors), 
             title="Color scale",
             lty=1, lwd=4, bty='n', box.lty = 1, box.lwd=2,
             cex=0.05, text.width=0.00001, seg.len=8
      )
      par(lend=1)
    }
    
  }
  
  if (gradient==F){
    points(x=log(logisticreg.newdat$stem.ages_since), y=jitter(logisticreg.newdat$uncex.merged.mtdnas,factor=0,amount=0.05), 
           pch=21, cex=0.9, bg=scales::alpha('white', alpha=0.8), lwd=0.5)
  }

  
}

pdf(file="model_shift_pGLM.pdf", height=9, width = 7.5)
{
par(mfrow=c(2,2))

###set up plotting for shift vs time ###
{
    # plot(x=log(logisticreg.newdat$stem.ages_since),y=jitter(logisticreg.newdat$uncex.merged.mtdnas,factor=0,amount=0.02), xlim=c(-3,4.2), ylim=c(-0.05,1.05),
    #      xlab="log(time) to T=66 Ma", ylab="probability", main="Model Shift Binomial pGLM", pch=21, col=make.transparent("gray", 0.0), bty='n', cex=0.000001, cex.main=1.0)
    # #axis(side=1, at = c(-3,-2,-1,0,1,2,3, 4.189655), labels = round(exp(c(-3,-2,-1,0,1,2,3, 4.189655)), 1))
    # points(lwd=0.5, x=log(logisticreg.newdat$stem.ages_since),y=jitter(logisticreg.newdat$uncex.merged.mtdnas,factor=0,amount=0.075), pch=21, bg=rgb(colorRamp(c("black", "white"))(asinTransform(logisticreg.newdat.alt$discordance)/max(asinTransform(logisticreg.newdat.alt$discordance)))/255, alpha = 0.8), cex=1.0)#make.transparent("gray", 0.3)
    # 
    baseplot()
    #plot the exon data
    {
      # #add the bootstrap lines for exons
      # for(i in 1:(phylo.logisticreg.fit.stem.ages.exons$boot/10)){
      #   cc <- (phylo.logisticreg.fit.stem.ages.exons$bootstrap[i,])
      #   curve(plogis(cc[1]+cc[2]*x),col=make.transparent(palette.colors(palette = "Okabe-Ito")[2], 0.25),add=TRUE, lwd=0.25, xlim=c(-0.8726974, 4.2))
      #   curve(plogis(cc[1]+cc[2]*x),col=make.transparent(palette.colors(palette = "Okabe-Ito")[2], 0.25),add=TRUE, lwd=0.25, xlim=c(-3.5, -0.8726974), lty=1)
      # }
      # 
      exonplot<-phylolm_bootstrap_CI(input = phylo.logisticreg.fit.stem.ages.exons, lims = c(-3.5, 4.2), breaks=1000, interval=int)
      phylolm_plot_CIs(input=exonplot, type='hdi', smooth=smooth, color=palette.colors(palette = "Okabe-Ito")[2], alpha=0.1, outline=F)
      
      #paste(coef(phylo.logisticreg.fit.stem.ages.exons), names(coef(phylo.logisticreg.fit.stem.ages.exons)), sep = ' * ', collapse = ' + ')
      
      #cc <- coef(phylo.logisticreg.fit.stem.ages.exons)
      cc <- (phylo.logisticreg.fit.stem.ages.exons$bootmean)
      curve(plogis(cc[1]+cc[2]*x),col=make.transparent(palette.colors(palette = "Okabe-Ito")[2], 0.95),add=TRUE, lwd=2, xlim=c(-0.8726974, 4.2))
      curve(plogis(cc[1]+cc[2]*x),col=make.transparent(palette.colors(palette = "Okabe-Ito")[2], 0.95),add=TRUE, lwd=2, xlim=c(-3.5, -0.8), lty=2)
      
    }
    
    max_y <- max(curvemod(plogis(cc[1] + cc[2] * x), xlim=c(-0.8726974, 4.2))$y)
    lines(x=c(-0.8726974, 4.2), y=c(max_y,max_y), col = palette.colors(palette = "Okabe-Ito")[2], lwd=0.5)
    #abline(h = max_y, col = palette.colors(palette = "Okabe-Ito")[2], xlim=c(-0.8726974, 4.2))
    text(x = 4.05, y = max_y, label = bquote(.(round(max_y,2))), pos = 4, cex=0.65)
    
    #plot the intron data
    {
      # for(i in 1:(phylo.logisticreg.fit.stem.ages.introns$boot/10)){
      #   cc <- (phylo.logisticreg.fit.stem.ages.introns$bootstrap[i,])
      #   curve(plogis(cc[1]+cc[2]*x),col=make.transparent(palette.colors(palette = "Okabe-Ito")[3], 0.25),add=TRUE, lwd=0.25, xlim=c(-0.8726974, 4.2))
      #   curve(plogis(cc[1]+cc[2]*x),col=make.transparent(palette.colors(palette = "Okabe-Ito")[3], 0.25),add=TRUE, lwd=0.25, xlim=c(-3.5, -0.8726974), lty=1)
      # }
      
      intronplot<-phylolm_bootstrap_CI(input = phylo.logisticreg.fit.stem.ages.introns, lims = c(-3.5, 4.2), breaks=1000, interval=int)
      phylolm_plot_CIs(input=intronplot, type='hdi', smooth=smooth, color=palette.colors(palette = "Okabe-Ito")[3], alpha=0.1, outline=F)
      
      
      #cc <- coef(phylo.logisticreg.fit.stem.ages.introns)
      cc <- (phylo.logisticreg.fit.stem.ages.introns$bootmean)
      curve(plogis(cc[1]+cc[2]*x),col=make.transparent(palette.colors(palette = "Okabe-Ito")[3], 0.95),add=TRUE, lwd=2, xlim=c(-0.8726974, 4.2))
      curve(plogis(cc[1]+cc[2]*x),col=make.transparent(palette.colors(palette = "Okabe-Ito")[3], 0.95),add=TRUE, lwd=2, xlim=c(-3.5, -0.8), lty=2)
      
    }
    
    max_y <- max(curvemod(plogis(cc[1] + cc[2] * x), xlim=c(-0.8726974, 4.2))$y)
    lines(x=c(-0.8726974, 4.2), y=c(max_y,max_y), col = palette.colors(palette = "Okabe-Ito")[3], lwd=0.5)
    #abline(h = max_y, col = palette.colors(palette = "Okabe-Ito")[3], xlim=c(-0.8726974, 4.2))
    text(x = 4.05, y = max_y, label = bquote(.(round(max_y,2))), pos = 4, cex=0.65)
    
    #plot the utr data
    {
      # for(i in 1:(phylo.logisticreg.fit.stem.ages.utrs$boot/10)){
      #   cc <- (phylo.logisticreg.fit.stem.ages.utrs$bootstrap[i,])
      #   curve(plogis(cc[1]+cc[2]*x),col=make.transparent(palette.colors(palette = "Okabe-Ito")[4], 0.25),add=TRUE, lwd=0.25, xlim=c(-0.8726974, 4.2))
      #   curve(plogis(cc[1]+cc[2]*x),col=make.transparent(palette.colors(palette = "Okabe-Ito")[4], 0.25),add=TRUE, lwd=0.25, xlim=c(-3.5, -0.8726974), lty=1)
      # }
      
      utrplot<-phylolm_bootstrap_CI(input = phylo.logisticreg.fit.stem.ages.utrs, lims = c(-3.5, 4.2), breaks=1000, interval=int)
      phylolm_plot_CIs(input=utrplot, type='hdi', smooth=smooth, color=palette.colors(palette = "Okabe-Ito")[4], alpha=0.1, outline=F)
      
      #cc <- coef(phylo.logisticreg.fit.stem.ages.utrs)
      cc <- (phylo.logisticreg.fit.stem.ages.utrs$bootmean)
      curve(plogis(cc[1]+cc[2]*x),col=make.transparent(palette.colors(palette = "Okabe-Ito")[4], 0.95),add=TRUE, lwd=2, xlim=c(-0.8726974, 4.2))
      curve(plogis(cc[1]+cc[2]*x),col=make.transparent(palette.colors(palette = "Okabe-Ito")[4], 0.95),add=TRUE, lwd=2, xlim=c(-3.5, -0.8), lty=2)
      
    }
    
    max_y <- max(curvemod(plogis(cc[1] + cc[2] * x), xlim=c(-0.8726974, 4.2))$y)
    lines(x=c(-0.8726974, 4.2), y=c(max_y,max_y), col = palette.colors(palette = "Okabe-Ito")[4], lwd=0.5)
    #abline(h = max_y, col = palette.colors(palette = "Okabe-Ito")[4], xlim=c(-0.8726974, 4.2))
    text(x = 4.05, y = max_y, label = bquote(.(round(max_y,2))), pos = 4, cex=0.65)
    
    #plot the mtdna data (all)
    {
      # for(i in 1:(phylo.logisticreg.fit.stem.ages.mtdnas.all$boot/10)){
      #   cc <- (phylo.logisticreg.fit.stem.ages.mtdnas.all$bootstrap[i,])
      #   curve(plogis(cc[1]+cc[2]*x),col=make.transparent(palette.colors(palette = "Okabe-Ito")[6], 0.25),add=TRUE, lwd=0.25, xlim=c(-0.8726974, 4.2))
      #   curve(plogis(cc[1]+cc[2]*x),col=make.transparent(palette.colors(palette = "Okabe-Ito")[6], 0.25),add=TRUE, lwd=0.25, xlim=c(-3.5, -0.8726974), lty=1)
      # }
      
      mt_all_plot<-phylolm_bootstrap_CI(input = phylo.logisticreg.fit.stem.ages.mtdnas.all, lims = c(-3.5, 4.2), breaks=1000, interval=int)
      phylolm_plot_CIs(input=mt_all_plot, type='hdi', smooth=smooth, color=palette.colors(palette = "Okabe-Ito")[6], alpha=0.1, outline=F)
      
      #cc <- coef(phylo.logisticreg.fit.stem.ages.mtdnas.all)
      cc <- (phylo.logisticreg.fit.stem.ages.mtdnas.all$bootmean)
      curve(plogis(cc[1]+cc[2]*x),col=make.transparent(palette.colors(palette = "Okabe-Ito")[6], 0.95),add=TRUE, lwd=2, xlim=c(-0.8726974, 4.2))
      curve(plogis(cc[1]+cc[2]*x),col=make.transparent(palette.colors(palette = "Okabe-Ito")[6], 0.95),add=TRUE, lwd=2, xlim=c(-3.5, -0.8726974), lty=2)
      
    }
    
    #plot the mtdna data (proteins)
    {
      # for(i in 1:(phylo.logisticreg.fit.stem.ages.mtdnas.proteins$boot/10)){
      #   cc <- (phylo.logisticreg.fit.stem.ages.mtdnas.proteins$bootstrap[i,])
      #   curve(plogis(cc[1]+cc[2]*x),col=make.transparent(palette.colors(palette = "Okabe-Ito")[7], 0.25),add=TRUE, lwd=0.25, xlim=c(-0.8726974, 4.2))
      #   curve(plogis(cc[1]+cc[2]*x),col=make.transparent(palette.colors(palette = "Okabe-Ito")[7], 0.25),add=TRUE, lwd=0.25, xlim=c(-3.5, -0.8726974), lty=1)
      # }
      
      mt_prot_plot<-phylolm_bootstrap_CI(input = phylo.logisticreg.fit.stem.ages.mtdnas.proteins, lims = c(-3.5, 4.2), breaks=1000, interval=int)
      phylolm_plot_CIs(input=mt_prot_plot, type='hdi', smooth=smooth, color=palette.colors(palette = "Okabe-Ito")[7], alpha=0.1, outline=F)
      
      #cc <- coef(phylo.logisticreg.fit.stem.ages.mtdnas.proteins)
      cc <- (phylo.logisticreg.fit.stem.ages.mtdnas.proteins$bootmean)
      curve(plogis(cc[1]+cc[2]*x),col=make.transparent(palette.colors(palette = "Okabe-Ito")[7], 0.95),add=TRUE, lwd=2, xlim=c(-0.8726974, 4.2))
      curve(plogis(cc[1]+cc[2]*x),col=make.transparent(palette.colors(palette = "Okabe-Ito")[7], 0.95),add=TRUE, lwd=2, xlim=c(-3.5, -0.8726974), lty=2)
      
    }
    
    #plot the merged signal + mtdna signal
    {
      # for(i in 1:(phylo.logisticreg.fit.stem.ages.merged.mtdnas$boot/10)){
      #   cc <- (phylo.logisticreg.fit.stem.ages.merged.mtdnas$bootstrap[i,])
      #   curve(plogis(cc[1]+cc[2]*x),col=make.transparent("gray", 0.25),add=TRUE, lwd=0.25, xlim=c(-0.8726974, 4.2))
      #   curve(plogis(cc[1]+cc[2]*x),col=make.transparent("gray", 0.25),add=TRUE, lwd=0.25, xlim=c(-3.5, -0.8726974), lty=1)
      # }
      
      mt_merged_plot<-phylolm_bootstrap_CI(input = phylo.logisticreg.fit.stem.ages.merged.mtdnas, lims = c(-3.5, 4.2), breaks=1000, interval=int)
      phylolm_plot_CIs(input=mt_merged_plot, type='hdi', smooth=smooth, color="gray", alpha=0.1, outline=F)
      
      #cc <- coef(phylo.logisticreg.fit.stem.ages.merged.mtdnas)
      cc <- phylo.logisticreg.fit.stem.ages.merged.mtdnas$bootmean
      curve(plogis(cc[1]+cc[2]*x),col=make.transparent("gray", 0.95),add=TRUE, lwd=2, xlim=c(-0.8726974, 4.2))
      curve(plogis(cc[1]+cc[2]*x),col=make.transparent("gray", 0.95),add=TRUE, lwd=2, xlim=c(-3.5, -0.8726974), lty=2)
      
    }
    
    #plot the merged data signal
    {
      # for(i in 1:(phylo.logisticreg.fit.stem.ages.merged$boot/10)){
      #   cc <- (phylo.logisticreg.fit.stem.ages.merged$bootstrap[i,])
      #   curve(plogis(cc[1]+cc[2]*x),col=make.transparent("black", 0.25),add=TRUE, lwd=0.25, xlim=c(-0.8726974, 4.2))
      #   curve(plogis(cc[1]+cc[2]*x),col=make.transparent("black", 0.25),add=TRUE, lwd=0.25, xlim=c(-3.5, -0.8726974), lty=1)
      # }
      
      merged_plot<-phylolm_bootstrap_CI(input = phylo.logisticreg.fit.stem.ages.merged, lims = c(-3.5, 4.2), breaks=1000, interval=int)
      phylolm_plot_CIs(input=merged_plot, type='hdi', smooth=smooth, color="black", alpha=0.1, outline=F)
      
      #cc <- coef(phylo.logisticreg.fit.stem.ages.merged)
      cc <- (phylo.logisticreg.fit.stem.ages.merged$bootmean)
      curve(plogis(cc[1]+cc[2]*x),col=make.transparent("black", 0.95),add=TRUE, lwd=2, xlim=c(-0.8726974, 4.2))
      curve(plogis(cc[1]+cc[2]*x),col=make.transparent("black", 0.95),add=TRUE, lwd=2, xlim=c(-3.5, -0.8726974), lty=2)
      
    }
    
    abline(v=2.302585, col="black", lty=3)
    max_y <- max(curvemod(plogis(cc[1] + cc[2] * x), xlim=c(-0.8726974, 4.2))$y)
    lines(x=c(-0.8726974, 4.2), y=c(max_y,max_y), col = 'black', lwd=0.5)
    #abline(h = max_y, col = "black", xlim=c(-0.8726974, 4.2))
    text(x = 4.05, y = max_y, label = bquote(.(round(max_y,2))), pos = 4, cex=0.65)
    
}

  
### set up plot for shift vs time + discordance ### (low discordance)
{
    # plot(x=log(logisticreg.newdat$stem.ages_since),y=jitter(logisticreg.newdat$uncex.merged,factor=0,amount=0.02), xlim=c(-3,4.2), ylim=c(-0.05,1.05),
    #      xlab="log(time) to T=66 Ma", ylab="probability", main="Model Shift Binomial pGLM (low discordance)", pch=21, col=make.transparent("gray", 0.0), bty='n', cex=0.000001, cex.main=1.0)
    # #axis(side=1, at = c(-3,-2,-1,0,1,2,3, 4.189655), labels = round(exp(c(-3,-2,-1,0,1,2,3, 4.189655)), 1))
    # points(lwd=0.5, x=log(logisticreg.newdat$stem.ages_since),y=jitter(logisticreg.newdat$uncex.merged,factor=0,amount=0.075), pch=21, bg=rgb(colorRamp(c("black", "white"))(asinTransform(logisticreg.newdat.alt$discordance)/max(asinTransform(logisticreg.newdat.alt$discordance)))/255, alpha=0.8), cex=1.0)#make.transparent("gray", 0.3)
    # 
    baseplot(title="Model Shift Binomial pGLM (low discordance)", gradient=T)
    #plot the data for exons
    X2<-X2_l.exon
    {
      # #add the bootstrap lines for exons
      # for(i in 1:(phylo.logisticreg.fit.stem.ages.exons.discordance$boot/10)){
      #   cc <- (phylo.logisticreg.fit.stem.ages.exons.discordance$bootstrap[i,])
      #   curve(plogis(cc[1]+cc[2]*x+cc[3]*X2),col=make.transparent(palette.colors(palette = "Okabe-Ito")[2], 0.25),add=TRUE, lwd=0.25, xlim=c(-0.8726974, 4.2))
      #   curve(plogis(cc[1]+cc[2]*x+cc[3]*X2),col=make.transparent(palette.colors(palette = "Okabe-Ito")[2], 0.25),add=TRUE, lwd=0.25, xlim=c(-3.5, -0.8726974), lty=1)
      # }
      
      exonplot<-phylolm_bootstrap_CI(input = phylo.logisticreg.fit.stem.ages.exons.discordance, lims = c(-3.5, 4.2), breaks=1000, interval=int, X2)
      phylolm_plot_CIs(input=exonplot, type='hdi', smooth=smooth, color=palette.colors(palette = "Okabe-Ito")[2], alpha=0.1, outline=F)
      
      #cc <- coef(phylo.logisticreg.fit.stem.ages.exons.discordance)
      cc <- (phylo.logisticreg.fit.stem.ages.exons.discordance$bootmean)
      curve(plogis(cc[1]+cc[2]*x+cc[3]*X2),col=make.transparent(palette.colors(palette = "Okabe-Ito")[2], 0.95),add=TRUE, lwd=2, xlim=c(-0.8726974, 4.2), type='l')
      curve(plogis(cc[1]+cc[2]*x+cc[3]*X2),col=make.transparent(palette.colors(palette = "Okabe-Ito")[2], 0.95),add=TRUE, lwd=2, xlim=c(-3.5, -0.8726974), lty=2, type='l')
    }
    
    max_y <- max(curvemod(plogis(cc[1]+cc[2]*x+cc[3]*X2), xlim=c(-0.8726974, 4.2))$y)
    lines(x=c(-0.8726974, 4.2), y=c(max_y,max_y), col = palette.colors(palette = "Okabe-Ito")[2], lwd=0.5)
    #abline(h = max_y, col = "black", lwd=0.5)
    text(x = 4.05, y = max_y, label = bquote(.(round(max_y,2))), pos = 4, cex=0.65)
    
    #plot data for introns
    X2<-X2_l.intron
    {
      # #add the bootstrap lines for introns
      # for(i in 1:(phylo.logisticreg.fit.stem.ages.introns.discordance$boot/10)){
      #   cc <- (phylo.logisticreg.fit.stem.ages.introns.discordance$bootstrap[i,])
      #   curve(plogis(cc[1]+cc[2]*x+cc[3]*x),col=make.transparent(palette.colors(palette = "Okabe-Ito")[3], 0.25),add=TRUE, lwd=0.25, xlim=c(-0.8726974, 4.2))
      #   curve(plogis(cc[1]+cc[2]*x+cc[3]*x),col=make.transparent(palette.colors(palette = "Okabe-Ito")[3], 0.25),add=TRUE, lwd=0.25, xlim=c(-3.5, -0.8726974), lty=1)
      # }
      
      intronplot<-phylolm_bootstrap_CI(input = phylo.logisticreg.fit.stem.ages.introns.discordance, lims = c(-3.5, 4.2), breaks=1000, interval=int, X2)
      phylolm_plot_CIs(input=intronplot, type='hdi', smooth=smooth, color=palette.colors(palette = "Okabe-Ito")[3], alpha=0.1, outline=F)
      
      #cc <- coef(phylo.logisticreg.fit.stem.ages.introns.discordance)
      cc <- (phylo.logisticreg.fit.stem.ages.introns.discordance$bootmean)
      curve(plogis(cc[1]+cc[2]*x+cc[3]*X2),col=make.transparent(palette.colors(palette = "Okabe-Ito")[3], 0.95),add=TRUE, lwd=2, xlim=c(-0.8726974, 4.2))
      curve(plogis(cc[1]+cc[2]*x+cc[3]*X2),col=make.transparent(palette.colors(palette = "Okabe-Ito")[3], 0.95),add=TRUE, lwd=2, xlim=c(-3.5, -0.8726974), lty=2)
      
    }
    
    max_y <- max(curvemod(plogis(cc[1]+cc[2]*x+cc[3]*X2), xlim=c(-0.8726974, 4.2))$y)
    lines(x=c(-0.8726974, 4.2), y=c(max_y,max_y), col = palette.colors(palette = "Okabe-Ito")[3], lwd=0.5)
    #abline(h = max_y, col = "black", lwd=0.5)
    text(x = 4.05, y = max_y, label = bquote(.(round(max_y,2))), pos = 4, cex=0.65)
    
    #plot data for utrs
    X2<-X2_l.utr
    {
      # #add the bootstrap lines for utrs
      # for(i in 1:(phylo.logisticreg.fit.stem.ages.utrs.discordance$boot/10)){
      #   cc <- (phylo.logisticreg.fit.stem.ages.utrs.discordance$bootstrap[i,])
      #   curve(plogis(cc[1]+cc[2]*x+cc[3]*x),col=make.transparent(palette.colors(palette = "Okabe-Ito")[4], 0.25),add=TRUE, lwd=0.25, xlim=c(-0.8726974, 4.2))
      #   curve(plogis(cc[1]+cc[2]*x+cc[3]*x),col=make.transparent(palette.colors(palette = "Okabe-Ito")[4], 0.25),add=TRUE, lwd=0.25, xlim=c(-3.5, -0.8726974), lty=1)
      # }
      
      utrplot<-phylolm_bootstrap_CI(input = phylo.logisticreg.fit.stem.ages.utrs.discordance, lims = c(-3.5, 4.2), breaks=1000, interval=int, X2)
      phylolm_plot_CIs(input=utrplot, type='hdi', smooth=smooth, color=palette.colors(palette = "Okabe-Ito")[4], alpha=0.1, outline=F)
      
      #cc <- coef(phylo.logisticreg.fit.stem.ages.utrs.discordance)
      cc <- (phylo.logisticreg.fit.stem.ages.utrs.discordance$bootmean)
      curve(plogis(cc[1]+cc[2]*x+cc[3]*X2),col=make.transparent(palette.colors(palette = "Okabe-Ito")[4], 0.95),add=TRUE, lwd=2, xlim=c(-0.8726974, 4.2))
      curve(plogis(cc[1]+cc[2]*x+cc[3]*X2),col=make.transparent(palette.colors(palette = "Okabe-Ito")[4], 0.95),add=TRUE, lwd=2, xlim=c(-3.5, -0.8726974), lty=2)
      
    }
    
    max_y <- max(curvemod(plogis(cc[1]+cc[2]*x+cc[3]*X2), xlim=c(-0.8726974, 4.2))$y)
    lines(x=c(-0.8726974, 4.2), y=c(max_y,max_y), col = palette.colors(palette = "Okabe-Ito")[4], lwd=0.5)
    #abline(h = max_y, col = "black", lwd=0.5)
    text(x = 4.05, y = max_y, label = bquote(.(round(max_y,2))), pos = 4, cex=0.65)
    
    #plot the data for merged
    X2<-X2_l.merged
    {
      # #add the bootstrap lines for merged
      # for(i in 1:(phylo.logisticreg.fit.stem.ages.merged.discordance$boot/10)){
      #   cc <- (phylo.logisticreg.fit.stem.ages.merged.discordance$bootstrap[i,])
      #   curve(plogis(cc[1]+cc[2]*x+cc[3]*x),col=make.transparent("black", 0.25),add=TRUE, lwd=0.25, xlim=c(-0.8726974, 4.2))
      #   curve(plogis(cc[1]+cc[2]*x+cc[3]*x),col=make.transparent("black", 0.25),add=TRUE, lwd=0.25, xlim=c(-3.5, -0.8726974), lty=1)
      # }
      #cc <- coef(phylo.logisticreg.fit.stem.ages.merged.discordance)
      
      mergedplot<-phylolm_bootstrap_CI(input = phylo.logisticreg.fit.stem.ages.merged.discordance, lims = c(-3.5, 4.2), breaks=1000, interval=int, X2)
      phylolm_plot_CIs(input=mergedplot, type='hdi', smooth=smooth, color='black', alpha=0.1, outline=F)
      
      cc <- (phylo.logisticreg.fit.stem.ages.merged.discordance$bootmean)
      curve(plogis(cc[1]+cc[2]*x+cc[3]*X2),col=make.transparent("black", 0.95),add=TRUE, lwd=2, xlim=c(-0.8726974, 4.2))
      curve(plogis(cc[1]+cc[2]*x+cc[3]*X2),col=make.transparent("black", 0.95),add=TRUE, lwd=2, xlim=c(-3.5, -0.8726974), lty=2)
      
    }
    
    abline(v=2.302585, col="black", lty=3)
    max_y <- max(curvemod(plogis(cc[1]+cc[2]*x+cc[3]*X2), xlim=c(-0.8726974, 4.2))$y)
    lines(x=c(-0.8726974, 4.2), y=c(max_y,max_y), col = 'black', lwd=0.5)
    #abline(h = max_y, col = "black", lwd=0.5)
    text(x = 4.05, y = max_y, label = bquote(.(round(max_y,2))), pos = 4, cex=0.65)
    
    ### add more points to curve function plotting --- try to figure out why the max isn't matching up
    
    
}
  
### set up plot for shift vs time + discordance ### (mean discordance)
{
  # plot(x=log(logisticreg.newdat$stem.ages_since),y=jitter(logisticreg.newdat$uncex.merged,factor=0,amount=0.02), xlim=c(-3,4.2), ylim=c(-0.05,1.05),
  #      xlab="log(time) to T=66 Ma", ylab="probability", main="Model Shift Binomial pGLM (mean discordance)", pch=21, col=make.transparent("gray", 0.0), bty='n', cex=0.000001, cex.main=1.0)
  # #axis(side=1, at = c(-3,-2,-1,0,1,2,3, 4.189655), labels = round(exp(c(-3,-2,-1,0,1,2,3, 4.189655)), 1))
  # points(lwd=0.5, x=log(logisticreg.newdat$stem.ages_since),y=jitter(logisticreg.newdat$uncex.merged,factor=0,amount=0.075), pch=21, bg=rgb(colorRamp(c("black", "white"))(asinTransform(logisticreg.newdat.alt$discordance)/max(asinTransform(logisticreg.newdat.alt$discordance)))/255, alpha=0.8), cex=1.0)#make.transparent("gray", 0.3)
  # 
  baseplot(title="Model Shift Binomial pGLM (mean discordance)", gradient=T)
  
  X2<-X2_m.exon
  #plot the data for exons
  {
    # #add the bootstrap lines for exons
    # for(i in 1:(phylo.logisticreg.fit.stem.ages.exons.discordance$boot/10)){
    #   cc <- (phylo.logisticreg.fit.stem.ages.exons.discordance$bootstrap[i,])
    #   curve(plogis(cc[1]+cc[2]*x+cc[3]*X2),col=make.transparent(palette.colors(palette = "Okabe-Ito")[2], 0.25),add=TRUE, lwd=0.25, xlim=c(-0.8726974, 4.2))
    #   curve(plogis(cc[1]+cc[2]*x+cc[3]*X2),col=make.transparent(palette.colors(palette = "Okabe-Ito")[2], 0.25),add=TRUE, lwd=0.25, xlim=c(-3.5, -0.8726974), lty=1)
    # }
    
    exonplot<-phylolm_bootstrap_CI(input = phylo.logisticreg.fit.stem.ages.exons.discordance, lims = c(-3.5, 4.2), breaks=1000, interval=int, X2)
    phylolm_plot_CIs(input=exonplot, type='hdi', smooth=smooth, color=palette.colors(palette = "Okabe-Ito")[2], alpha=0.1, outline=F)
    
    #cc <- coef(phylo.logisticreg.fit.stem.ages.exons.discordance)
    cc <- (phylo.logisticreg.fit.stem.ages.exons.discordance$bootmean)
    curve(plogis(cc[1]+cc[2]*x+cc[3]*X2),col=make.transparent(palette.colors(palette = "Okabe-Ito")[2], 0.95),add=TRUE, lwd=2, xlim=c(-0.8726974, 4.2), type='l')
    curve(plogis(cc[1]+cc[2]*x+cc[3]*X2),col=make.transparent(palette.colors(palette = "Okabe-Ito")[2], 0.95),add=TRUE, lwd=2, xlim=c(-3.5, -0.8726974), lty=2, type='l')
  }
  
  max_y <- max(curvemod(plogis(cc[1]+cc[2]*x+cc[3]*X2), xlim=c(-0.8726974, 4.2))$y)
  lines(x=c(-0.8726974, 4.2), y=c(max_y,max_y), col = palette.colors(palette = "Okabe-Ito")[2], lwd=0.5)
  #abline(h = max_y, col = "black", lwd=0.5)
  text(x = 4.05, y = max_y, label = bquote(.(round(max_y,2))), pos = 4, cex=0.65)
  
  #plot data for introns
  X2<-X2_m.intron
  {
    # #add the bootstrap lines for introns
    # for(i in 1:(phylo.logisticreg.fit.stem.ages.introns.discordance$boot/10)){
    #   cc <- (phylo.logisticreg.fit.stem.ages.introns.discordance$bootstrap[i,])
    #   curve(plogis(cc[1]+cc[2]*x+cc[3]*x),col=make.transparent(palette.colors(palette = "Okabe-Ito")[3], 0.25),add=TRUE, lwd=0.25, xlim=c(-0.8726974, 4.2))
    #   curve(plogis(cc[1]+cc[2]*x+cc[3]*x),col=make.transparent(palette.colors(palette = "Okabe-Ito")[3], 0.25),add=TRUE, lwd=0.25, xlim=c(-3.5, -0.8726974), lty=1)
    # }
    
    intronplot<-phylolm_bootstrap_CI(input = phylo.logisticreg.fit.stem.ages.introns.discordance, lims = c(-3.5, 4.2), breaks=1000, interval=int, X2)
    phylolm_plot_CIs(input=intronplot, type='hdi', smooth=smooth, color=palette.colors(palette = "Okabe-Ito")[3], alpha=0.1, outline=F)
    
    #cc <- coef(phylo.logisticreg.fit.stem.ages.introns.discordance)
    cc <- (phylo.logisticreg.fit.stem.ages.introns.discordance$bootmean)
    curve(plogis(cc[1]+cc[2]*x+cc[3]*X2),col=make.transparent(palette.colors(palette = "Okabe-Ito")[3], 0.95),add=TRUE, lwd=2, xlim=c(-0.8726974, 4.2))
    curve(plogis(cc[1]+cc[2]*x+cc[3]*X2),col=make.transparent(palette.colors(palette = "Okabe-Ito")[3], 0.95),add=TRUE, lwd=2, xlim=c(-3.5, -0.8726974), lty=2)
    
  }
  
  max_y <- max(curvemod(plogis(cc[1]+cc[2]*x+cc[3]*X2), xlim=c(-0.8726974, 4.2))$y)
  lines(x=c(-0.8726974, 4.2), y=c(max_y,max_y), col = palette.colors(palette = "Okabe-Ito")[3], lwd=0.5)
  #abline(h = max_y, col = "black", lwd=0.5)
  text(x = 4.05, y = max_y, label = bquote(.(round(max_y,2))), pos = 4, cex=0.65)
  
  #plot data for utrs
  X2<-X2_m.utr
  {
    # #add the bootstrap lines for utrs
    # for(i in 1:(phylo.logisticreg.fit.stem.ages.utrs.discordance$boot/10)){
    #   cc <- (phylo.logisticreg.fit.stem.ages.utrs.discordance$bootstrap[i,])
    #   curve(plogis(cc[1]+cc[2]*x+cc[3]*x),col=make.transparent(palette.colors(palette = "Okabe-Ito")[4], 0.25),add=TRUE, lwd=0.25, xlim=c(-0.8726974, 4.2))
    #   curve(plogis(cc[1]+cc[2]*x+cc[3]*x),col=make.transparent(palette.colors(palette = "Okabe-Ito")[4], 0.25),add=TRUE, lwd=0.25, xlim=c(-3.5, -0.8726974), lty=1)
    # }
    
    utrplot<-phylolm_bootstrap_CI(input = phylo.logisticreg.fit.stem.ages.utrs.discordance, lims = c(-3.5, 4.2), breaks=1000, interval=int, X2)
    phylolm_plot_CIs(input=utrplot, type='hdi', smooth=smooth, color=palette.colors(palette = "Okabe-Ito")[4], alpha=0.1, outline=F)
    
    #cc <- coef(phylo.logisticreg.fit.stem.ages.utrs.discordance)
    cc <- (phylo.logisticreg.fit.stem.ages.utrs.discordance$bootmean)
    curve(plogis(cc[1]+cc[2]*x+cc[3]*X2),col=make.transparent(palette.colors(palette = "Okabe-Ito")[4], 0.95),add=TRUE, lwd=2, xlim=c(-0.8726974, 4.2))
    curve(plogis(cc[1]+cc[2]*x+cc[3]*X2),col=make.transparent(palette.colors(palette = "Okabe-Ito")[4], 0.95),add=TRUE, lwd=2, xlim=c(-3.5, -0.8726974), lty=2)
    
  }
  
  max_y <- max(curvemod(plogis(cc[1]+cc[2]*x+cc[3]*X2), xlim=c(-0.8726974, 4.2))$y)
  lines(x=c(-0.8726974, 4.2), y=c(max_y,max_y), col = palette.colors(palette = "Okabe-Ito")[4], lwd=0.5)
  #abline(h = max_y, col = "black", lwd=0.5)
  text(x = 4.05, y = max_y, label = bquote(.(round(max_y,2))), pos = 4, cex=0.65)
  
  #plot the data for merged
  X2<-X2_m.merged
  {
    # #add the bootstrap lines for merged
    # for(i in 1:(phylo.logisticreg.fit.stem.ages.merged.discordance$boot/10)){
    #   cc <- (phylo.logisticreg.fit.stem.ages.merged.discordance$bootstrap[i,])
    #   curve(plogis(cc[1]+cc[2]*x+cc[3]*x),col=make.transparent("black", 0.25),add=TRUE, lwd=0.25, xlim=c(-0.8726974, 4.2))
    #   curve(plogis(cc[1]+cc[2]*x+cc[3]*x),col=make.transparent("black", 0.25),add=TRUE, lwd=0.25, xlim=c(-3.5, -0.8726974), lty=1)
    # }
    #cc <- coef(phylo.logisticreg.fit.stem.ages.merged.discordance)
    
    mergedplot<-phylolm_bootstrap_CI(input = phylo.logisticreg.fit.stem.ages.merged.discordance, lims = c(-3.5, 4.2), breaks=1000, interval=int, X2)
    phylolm_plot_CIs(input=mergedplot, type='hdi', smooth=smooth, color='black', alpha=0.1, outline=F)
    
    cc <- (phylo.logisticreg.fit.stem.ages.merged.discordance$bootmean)
    curve(plogis(cc[1]+cc[2]*x+cc[3]*X2),col=make.transparent("black", 0.95),add=TRUE, lwd=2, xlim=c(-0.8726974, 4.2))
    curve(plogis(cc[1]+cc[2]*x+cc[3]*X2),col=make.transparent("black", 0.95),add=TRUE, lwd=2, xlim=c(-3.5, -0.8726974), lty=2)
    
  }
  
  abline(v=2.302585, col="black", lty=3)
  max_y <- max(curvemod(plogis(cc[1]+cc[2]*x+cc[3]*X2), xlim=c(-0.8726974, 4.2))$y)
  lines(x=c(-0.8726974, 4.2), y=c(max_y,max_y), col = 'black', lwd=0.5)
  #abline(h = max_y, col = "black", lwd=0.5)
  text(x = 4.05, y = max_y, label = bquote(.(round(max_y,2))), pos = 4, cex=0.65)
  
}
  
### set up plot for shift vs time + discordance ### (high discordance)
{
  # plot(x=log(logisticreg.newdat$stem.ages_since),y=jitter(logisticreg.newdat$uncex.merged,factor=0,amount=0.02), xlim=c(-3,4.2), ylim=c(-0.05,1.05),
  #      xlab="log(time) to T=66 Ma", ylab="probability", main="Model Shift Binomial pGLM (high discordance)", pch=21, col=make.transparent("gray", 0.0), bty='n', cex=0.000001, cex.main=1.0)
  # #axis(side=1, at = c(-3,-2,-1,0,1,2,3, 4.189655), labels = round(exp(c(-3,-2,-1,0,1,2,3, 4.189655)), 1))
  # points(lwd=0.5, x=log(logisticreg.newdat$stem.ages_since),y=jitter(logisticreg.newdat$uncex.merged,factor=0,amount=0.075), pch=21, bg=rgb(colorRamp(c("black", "white"))(asinTransform(logisticreg.newdat.alt$discordance)/max(asinTransform(logisticreg.newdat.alt$discordance)))/255, alpha=0.8), cex=1.0)#make.transparent("gray", 0.3)
  # 
  baseplot(title="Model Shift Binomial pGLM (high discordance)", gradient=T)
  
  X2<-X2_h.exon
  #plot the data for exons
  {
    # #add the bootstrap lines for exons
    # for(i in 1:(phylo.logisticreg.fit.stem.ages.exons.discordance$boot/10)){
    #   cc <- (phylo.logisticreg.fit.stem.ages.exons.discordance$bootstrap[i,])
    #   curve(plogis(cc[1]+cc[2]*x+cc[3]*X2),col=make.transparent(palette.colors(palette = "Okabe-Ito")[2], 0.25),add=TRUE, lwd=0.25, xlim=c(-0.8726974, 4.2))
    #   curve(plogis(cc[1]+cc[2]*x+cc[3]*X2),col=make.transparent(palette.colors(palette = "Okabe-Ito")[2], 0.25),add=TRUE, lwd=0.25, xlim=c(-3.5, -0.8726974), lty=1)
    # }
    
    exonplot<-phylolm_bootstrap_CI(input = phylo.logisticreg.fit.stem.ages.exons.discordance, lims = c(-3.5, 4.2), breaks=1000, interval=int, X2)
    phylolm_plot_CIs(input=exonplot, type='hdi', smooth=smooth, color=palette.colors(palette = "Okabe-Ito")[2], alpha=0.1, outline=F)
    
    #cc <- coef(phylo.logisticreg.fit.stem.ages.exons.discordance)
    cc <- (phylo.logisticreg.fit.stem.ages.exons.discordance$bootmean)
    curve(plogis(cc[1]+cc[2]*x+cc[3]*X2),col=make.transparent(palette.colors(palette = "Okabe-Ito")[2], 0.95),add=TRUE, lwd=2, xlim=c(-0.8726974, 4.2), type='l')
    curve(plogis(cc[1]+cc[2]*x+cc[3]*X2),col=make.transparent(palette.colors(palette = "Okabe-Ito")[2], 0.95),add=TRUE, lwd=2, xlim=c(-3.5, -0.8726974), lty=2, type='l')
  }
  
  max_y <- max(curvemod(plogis(cc[1]+cc[2]*x+cc[3]*X2), xlim=c(-0.8726974, 4.2))$y)
  lines(x=c(-0.8726974, 4.2), y=c(max_y,max_y), col = palette.colors(palette = "Okabe-Ito")[2], lwd=0.5)
  #abline(h = max_y, col = "black", lwd=0.5)
  text(x = 4.05, y = max_y, label = bquote(.(round(max_y,2))), pos = 4, cex=0.65)
  
  #plot data for introns
  X2<-X2_h.intron
  {
    # #add the bootstrap lines for introns
    # for(i in 1:(phylo.logisticreg.fit.stem.ages.introns.discordance$boot/10)){
    #   cc <- (phylo.logisticreg.fit.stem.ages.introns.discordance$bootstrap[i,])
    #   curve(plogis(cc[1]+cc[2]*x+cc[3]*x),col=make.transparent(palette.colors(palette = "Okabe-Ito")[3], 0.25),add=TRUE, lwd=0.25, xlim=c(-0.8726974, 4.2))
    #   curve(plogis(cc[1]+cc[2]*x+cc[3]*x),col=make.transparent(palette.colors(palette = "Okabe-Ito")[3], 0.25),add=TRUE, lwd=0.25, xlim=c(-3.5, -0.8726974), lty=1)
    # }
    
    intronplot<-phylolm_bootstrap_CI(input = phylo.logisticreg.fit.stem.ages.introns.discordance, lims = c(-3.5, 4.2), breaks=1000, interval=int, X2)
    phylolm_plot_CIs(input=intronplot, type='hdi', smooth=smooth, color=palette.colors(palette = "Okabe-Ito")[3], alpha=0.1, outline=F)
    
    #cc <- coef(phylo.logisticreg.fit.stem.ages.introns.discordance)
    cc <- (phylo.logisticreg.fit.stem.ages.introns.discordance$bootmean)
    curve(plogis(cc[1]+cc[2]*x+cc[3]*X2),col=make.transparent(palette.colors(palette = "Okabe-Ito")[3], 0.95),add=TRUE, lwd=2, xlim=c(-0.8726974, 4.2))
    curve(plogis(cc[1]+cc[2]*x+cc[3]*X2),col=make.transparent(palette.colors(palette = "Okabe-Ito")[3], 0.95),add=TRUE, lwd=2, xlim=c(-3.5, -0.8726974), lty=2)
    
  }
  
  max_y <- max(curvemod(plogis(cc[1]+cc[2]*x+cc[3]*X2), xlim=c(-0.8726974, 4.2))$y)
  lines(x=c(-0.8726974, 4.2), y=c(max_y,max_y), col = palette.colors(palette = "Okabe-Ito")[3], lwd=0.5)
  #abline(h = max_y, col = "black", lwd=0.5)
  text(x = 4.05, y = max_y, label = bquote(.(round(max_y,2))), pos = 4, cex=0.65)
  
  #plot data for utrs
  X2<-X2_h.utr
  {
    # #add the bootstrap lines for utrs
    # for(i in 1:(phylo.logisticreg.fit.stem.ages.utrs.discordance$boot/10)){
    #   cc <- (phylo.logisticreg.fit.stem.ages.utrs.discordance$bootstrap[i,])
    #   curve(plogis(cc[1]+cc[2]*x+cc[3]*x),col=make.transparent(palette.colors(palette = "Okabe-Ito")[4], 0.25),add=TRUE, lwd=0.25, xlim=c(-0.8726974, 4.2))
    #   curve(plogis(cc[1]+cc[2]*x+cc[3]*x),col=make.transparent(palette.colors(palette = "Okabe-Ito")[4], 0.25),add=TRUE, lwd=0.25, xlim=c(-3.5, -0.8726974), lty=1)
    # }
    
    utrplot<-phylolm_bootstrap_CI(input = phylo.logisticreg.fit.stem.ages.utrs.discordance, lims = c(-3.5, 4.2), breaks=1000, interval=int, X2)
    phylolm_plot_CIs(input=utrplot, type='hdi', smooth=smooth, color=palette.colors(palette = "Okabe-Ito")[4], alpha=0.1, outline=F)
    
    #cc <- coef(phylo.logisticreg.fit.stem.ages.utrs.discordance)
    cc <- (phylo.logisticreg.fit.stem.ages.utrs.discordance$bootmean)
    curve(plogis(cc[1]+cc[2]*x+cc[3]*X2),col=make.transparent(palette.colors(palette = "Okabe-Ito")[4], 0.95),add=TRUE, lwd=2, xlim=c(-0.8726974, 4.2))
    curve(plogis(cc[1]+cc[2]*x+cc[3]*X2),col=make.transparent(palette.colors(palette = "Okabe-Ito")[4], 0.95),add=TRUE, lwd=2, xlim=c(-3.5, -0.8726974), lty=2)
    
  }
  
  max_y <- max(curvemod(plogis(cc[1]+cc[2]*x+cc[3]*X2), xlim=c(-0.8726974, 4.2))$y)
  lines(x=c(-0.8726974, 4.2), y=c(max_y,max_y), col = palette.colors(palette = "Okabe-Ito")[4], lwd=0.5)
  #abline(h = max_y, col = "black", lwd=0.5)
  text(x = 4.05, y = max_y, label = bquote(.(round(max_y,2))), pos = 4, cex=0.65)
  
  #plot the data for merged
  X2<-X2_h.merged
  {
    # #add the bootstrap lines for merged
    # for(i in 1:(phylo.logisticreg.fit.stem.ages.merged.discordance$boot/10)){
    #   cc <- (phylo.logisticreg.fit.stem.ages.merged.discordance$bootstrap[i,])
    #   curve(plogis(cc[1]+cc[2]*x+cc[3]*x),col=make.transparent("black", 0.25),add=TRUE, lwd=0.25, xlim=c(-0.8726974, 4.2))
    #   curve(plogis(cc[1]+cc[2]*x+cc[3]*x),col=make.transparent("black", 0.25),add=TRUE, lwd=0.25, xlim=c(-3.5, -0.8726974), lty=1)
    # }
    #cc <- coef(phylo.logisticreg.fit.stem.ages.merged.discordance)
    
    mergedplot<-phylolm_bootstrap_CI(input = phylo.logisticreg.fit.stem.ages.merged.discordance, lims = c(-3.5, 4.2), breaks=1000, interval=int, X2)
    phylolm_plot_CIs(input=mergedplot, type='hdi', smooth=smooth, color='black', alpha=0.1, outline=F)
    
    cc <- (phylo.logisticreg.fit.stem.ages.merged.discordance$bootmean)
    curve(plogis(cc[1]+cc[2]*x+cc[3]*X2),col=make.transparent("black", 0.95),add=TRUE, lwd=2, xlim=c(-0.8726974, 4.2))
    curve(plogis(cc[1]+cc[2]*x+cc[3]*X2),col=make.transparent("black", 0.95),add=TRUE, lwd=2, xlim=c(-3.5, -0.8726974), lty=2)
    
  }
  
  abline(v=2.302585, col="black", lty=3)
  max_y <- max(curvemod(plogis(cc[1]+cc[2]*x+cc[3]*X2), xlim=c(-0.8726974, 4.2))$y)
  lines(x=c(-0.8726974, 4.2), y=c(max_y,max_y), col = 'black', lwd=0.5)
  #abline(h = max_y, col = "black", lwd=0.5)
  text(x = 4.05, y = max_y, label = bquote(.(round(max_y,2))), pos = 4, cex=0.65)
  
}

#add a bit of distance 

### end plot ####
  
# 
# ### set up plot for shift vs time + discordance ### Bootmodes version
# {
#   plot(x=log(logisticreg.newdat$stem.ages_since),y=jitter(logisticreg.newdat$uncex.merged,factor=0,amount=0.02), xlim=c(-3,4.2), ylim=c(-0.05,1.05),
#        xlab="log(time) to T=66 Ma", ylab="probability", main="Model Shift Binomial pGLM", pch=21, col=make.transparent("gray", 0.0), bty='n', cex=0.000001, cex.main=1.0)
#   #axis(side=1, at = c(-3,-2,-1,0,1,2,3, 4.189655), labels = round(exp(c(-3,-2,-1,0,1,2,3, 4.189655)), 1))
#   points(x=log(logisticreg.newdat$stem.ages_since),y=jitter(logisticreg.newdat$uncex.merged,factor=0,amount=0.075), pch=21, bg=rgb(colorRamp(c("red", "white"))(asinTransform(logisticreg.newdat.alt$discordance)/max(asinTransform(logisticreg.newdat.alt$discordance)))/255), cex=1.0)#make.transparent("gray", 0.3)
#   
#   #plot the data for exons
#   {
#     #add the bootstrap lines for exons
#     for(i in 1:(phylo.logisticreg.fit.stem.ages.exons.discordance$boot/10)){
#       cc <- (phylo.logisticreg.fit.stem.ages.exons.discordance$bootstrap[i,])
#       curve(plogis(cc[1]+cc[2]*x+cc[3]*x),col=make.transparent(palette.colors(palette = "Okabe-Ito")[2], 0.25),add=TRUE, lwd=0.25, xlim=c(-0.8726974, 4.2))
#       curve(plogis(cc[1]+cc[2]*x+cc[3]*x),col=make.transparent(palette.colors(palette = "Okabe-Ito")[2], 0.25),add=TRUE, lwd=0.25, xlim=c(-3.5, -0.8726974), lty=1)
#     }
#     
#     #cc <- coef(phylo.logisticreg.fit.stem.ages.exons.discordance)
#     cc <- (phylo.logisticreg.fit.stem.ages.exons.discordance$bootmode)
#     curve(plogis(cc[1]+cc[2]*x+cc[3]*x),col=make.transparent(palette.colors(palette = "Okabe-Ito")[2], 0.95),add=TRUE, lwd=2, xlim=c(-0.8726974, 4.2))
#     curve(plogis(cc[1]+cc[2]*x+cc[3]*x),col=make.transparent(palette.colors(palette = "Okabe-Ito")[2], 0.95),add=TRUE, lwd=2, xlim=c(-3, -0.8726974), lty=2)
#   }
#   
#   #plot data for introns
#   {
#     #add the bootstrap lines for introns
#     for(i in 1:(phylo.logisticreg.fit.stem.ages.introns.discordance$boot/10)){
#       cc <- (phylo.logisticreg.fit.stem.ages.introns.discordance$bootstrap[i,])
#       curve(plogis(cc[1]+cc[2]*x+cc[3]*x),col=make.transparent(palette.colors(palette = "Okabe-Ito")[3], 0.25),add=TRUE, lwd=0.25, xlim=c(-0.8726974, 4.2))
#       curve(plogis(cc[1]+cc[2]*x+cc[3]*x),col=make.transparent(palette.colors(palette = "Okabe-Ito")[3], 0.25),add=TRUE, lwd=0.25, xlim=c(-3.5, -0.8726974), lty=1)
#     }
#     #cc <- coef(phylo.logisticreg.fit.stem.ages.introns.discordance)
#     cc <- (phylo.logisticreg.fit.stem.ages.introns.discordance$bootmode)
#     curve(plogis(cc[1]+cc[2]*x+cc[3]*x),col=make.transparent(palette.colors(palette = "Okabe-Ito")[3], 0.95),add=TRUE, lwd=2, xlim=c(-0.8726974, 4.2))
#     curve(plogis(cc[1]+cc[2]*x+cc[3]*x),col=make.transparent(palette.colors(palette = "Okabe-Ito")[3], 0.95),add=TRUE, lwd=2, xlim=c(-3, -0.8726974), lty=2)
#     
#   }
#   
#   #plot data for utrs
#   {
#     #add the bootstrap lines for utrs
#     for(i in 1:(phylo.logisticreg.fit.stem.ages.utrs.discordance$boot/10)){
#       cc <- (phylo.logisticreg.fit.stem.ages.utrs.discordance$bootstrap[i,])
#       curve(plogis(cc[1]+cc[2]*x+cc[3]*x),col=make.transparent(palette.colors(palette = "Okabe-Ito")[4], 0.25),add=TRUE, lwd=0.25, xlim=c(-0.8726974, 4.2))
#       curve(plogis(cc[1]+cc[2]*x+cc[3]*x),col=make.transparent(palette.colors(palette = "Okabe-Ito")[4], 0.25),add=TRUE, lwd=0.25, xlim=c(-3.5, -0.8726974), lty=1)
#     }
#     #cc <- coef(phylo.logisticreg.fit.stem.ages.utrs.discordance)
#     cc <- (phylo.logisticreg.fit.stem.ages.utrs.discordance$bootmode)
#     curve(plogis(cc[1]+cc[2]*x+cc[3]*x),col=make.transparent(palette.colors(palette = "Okabe-Ito")[4], 0.95),add=TRUE, lwd=2, xlim=c(-0.8726974, 4.2))
#     curve(plogis(cc[1]+cc[2]*x+cc[3]*x),col=make.transparent(palette.colors(palette = "Okabe-Ito")[4], 0.95),add=TRUE, lwd=2, xlim=c(-3, -0.8726974), lty=2)
#     
#   }
#   
#   #plot the data for merged
#   {
#     #add the bootstrap lines for merged
#     for(i in 1:(phylo.logisticreg.fit.stem.ages.merged.discordance$boot/10)){
#       cc <- (phylo.logisticreg.fit.stem.ages.merged.discordance$bootstrap[i,])
#       curve(plogis(cc[1]+cc[2]*x+cc[3]*x),col=make.transparent("black", 0.25),add=TRUE, lwd=0.25, xlim=c(-0.8726974, 4.2))
#       curve(plogis(cc[1]+cc[2]*x+cc[3]*x),col=make.transparent("black", 0.25),add=TRUE, lwd=0.25, xlim=c(-3.5, -0.8726974), lty=1)
#     }
#     #cc <- coef(phylo.logisticreg.fit.stem.ages.merged.discordance)
#     cc <- (phylo.logisticreg.fit.stem.ages.merged.discordance$bootmode)
#     curve(plogis(cc[1]+cc[2]*x+cc[3]*x),col=make.transparent("black", 0.95),add=TRUE, lwd=4, xlim=c(-0.8726974, 4.2))
#     curve(plogis(cc[1]+cc[2]*x+cc[3]*x),col=make.transparent("black", 0.95),add=TRUE, lwd=4, xlim=c(-3, -0.8726974), lty=2)
#     
#   }
#   
#   
#   abline(v=-0.8726974, col="black", lwd=1)
#   #abline(v=1.609438, col="black", lty=2)
#   abline(v=2.302585, col="black", lty=2)
#   abline(h=0.9, col="black")
# }
# #plotting modal values for introns
# ### end plot ####
# 
# #trying separation plot
# {
# par(mfrow=c(2,1))
# separationplot(pred=phylo.logisticreg.fit.stem.ages.merged.mtdnas$fitted.values, actual=logisticreg.newdat.alt$uncex.merged.mtdnas[-197], type="rect", line=T, show.expected = T, shuffle=T, newplot=F, BW=T)
# separationplot(pred=phylo.logisticreg.fit.stem.ages.merged.mtdnas.discordance$fitted.values, actual=logisticreg.newdat.alt$uncex.merged.mtdnas[-197], type="rect", line=T, show.expected = T, shuffle=T, newplot=F, BW=T)
# }


#notes - logistic regression
{
# 
# #see if we can make a prediction from the phyloglm
# #predict(phylo.logisticreg.fit, data.frame(discordance=0.9), type="response")
# #this doesn't work, no predict method for phyloglm
# 
# #use wolfram alpha to find the intercept
# 
# p<- ? 
# x <- (log(p/(1-p)) - coef(phylo.logisticreg.fit)[1]) / coef(phylo.logisticreg.fit)[2]
# x
# 
# -1=(log(p/(1-p)) - -1.12931 ) / -1.414319
# 
# 0.570774

}
}
dev.off()



#Section 5 
###############d###########################
#get data for gc content and codon usage #
##########################################
{
#consensus.all.timetree$phy$tip.label
{
#read in the exons for coRdon analysis
#exons<-readSet(file="~/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/unmasked/min2x/qual20_cov2_haplo_bestonly/initial_test_filters/min50bp_min10p_aligned/ALIGNED/phased/sample1haplo/CODING/exon_concat.fas")
exons<-readSet(file="./AHE_REASSEMBLY/sequences/concat/exon_concat.fas")

#reorder to match tree
exons<-exons[consensus.all.timetree$phy$tip.label,]
#convert to character vector
exonchar<-as.character(exons)
#this may need to be resorted
#exon_bin<-read.dna("/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/unmasked/min2x/qual20_cov2_haplo_bestonly/initial_test_filters/min50bp_min10p_aligned/ALIGNED/phased/sample2haplo/CODING/exons_concat.fasta", format="fasta")
#as.DNAbin(exons)

#reading in individual exons (codon usage looping)
#exons_separate<-exonload(path='/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/unmasked/min2x/qual20_cov2_haplo_bestonly/initial_test_filters/min50bp_min10p_aligned/ALIGNED/phased/sample1haplo/CODING/exons/NT_fasta', treeNames=consensus.all.timetree$phy$tip.label)
#uncomment the next line to load from github path
#exons_separate<-exonload(path='./AHE_REASSEMBLY/sequences/aligned/sample1haplo/CODING/exons/NT_fasta', treeNames=consensus.all.timetree$phy$tip.label)

#saveRDS(exons_separate, file="./RDS/exons_separate.RDS")
exons_separate<-readRDS(file="./RDS/exons_separate.RDS")

#read in the introns for GC content calculation
#introns<-readSet(file="~/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/unmasked/min2x/qual20_cov2_haplo_bestonly/initial_test_filters/min50bp_min10p_aligned/ALIGNED/phased/sample1haplo/NONCODING/TRIMAL_GT_0.05/intron_concat.fas")
introns<-readSet(file="./AHE_REASSEMBLY/sequences/concat/intron_concat.fas")

#reorder to match tree
introns<-introns[consensus.all.timetree$phy$tip.label,]
#convert to character vector
intronchar<-as.character(introns)
#this may need to be resorted
#intron_bin<- as.DNAbin(introns)
#read.dna("/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/unmasked/min2x/qual20_cov2_haplo_bestonly/initial_test_filters/min50bp_min10p_aligned/ALIGNED/phased/sample2haplo/NONCODING/intron_concat.fasta", format="fasta")

#read in the UTRs for GC content calculation
#utrs<-readSet(file="~/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/unmasked/min2x/qual20_cov2_haplo_bestonly/initial_test_filters/min50bp_min10p_aligned/ALIGNED/phased/sample1haplo/NONCODING/TRIMAL_GT_0.05/utr_concat.fas")
utrs<-readSet(file="./AHE_REASSEMBLY/sequences/concat/utr_concat.fas")

#reorder to match tree
utrs<-utrs[consensus.all.timetree$phy$tip.label,]
#convert to character vector
utrchar<-as.character(utrs)
#this may need to be resorted
#utr_bin <- as.DNAbin(utrs)
#<-read.dna("/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/unmasked/min2x/qual20_cov2_haplo_bestonly/initial_test_filters/min50bp_min10p_aligned/ALIGNED/phased/sample2haplo/NONCODING/utr_concat.fasta", format="fasta")

#read in all data alignmenet for GC content calculation
#file is too big for github, so you must find the zipped file in ./AHE_REASSEMBLY/sequences/concat/ and unzip it, then change the path
alldata<-readSet(file="~/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/unmasked/min2x/qual20_cov2_haplo_bestonly/initial_test_filters/min50bp_min10p_aligned/ALIGNED/phased/sample1haplo/concat_all_0.05_haplo1.fasta")

#reorder to match tree
alldata<-alldata[consensus.all.timetree$phy$tip.label,]
#convert to character vector
alldatachar<-as.character(alldata)
#this may need to be resorted
#all_bin <- as.DNAbin(alldata)
#<-read.dna("/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/unmasked/min2x/qual20_cov2_haplo_bestonly/initial_test_filters/min50bp_min10p_aligned/ALIGNED/phased/sample2haplo/concat_all_0.05.fasta", format = "fasta")


###mtDNAs###
#need a new translation table 
mtDNA.data.seq<-rbind(mtDNA.data, c("Caiman_croccodilus_1", "Caiman_crocodilus"), c("Crocodylus_porosus_1", "Crocodylus_porosus"))

#read in the mtdna dataset for coRdon analysis
#mtdnas<-readSet(file="/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/mtdnas/concat/concat_all/concat_all_paup.fas")
mtdnas<-readSet(file="./mtDNA_REASSEMBLY/concat_all_paup.fas")

#swap names
mtdnas<-translation_nameswap_set(object=mtdnas, names = mtDNA.data.seq$mtdnaNames, swap = mtDNA.data.seq$treeID)
#reorder to match the nuclear dataset
mtdnas<-mtdnas[consensus.all.timetree$phy$tip.label,]
#convert to character vector
mtdnaschar<-as.character(mtdnas)
#swap names
#mtdnaschar<-translation_nameswap(object=mtdnaschar, names = mtDNA.data.seq$mtdnaNames, swap = mtDNA.data.seq$treeID)

#read in the proteins for coRdon analysis
#mtdnaproteins<-readSet(file="/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/mtdnas/concat/concat_proteins_no_nd6/concat_proteins_no_nd6.fas")
mtdnaproteins<-readSet(file="./mtDNA_REASSEMBLY/concat_proteins_no_nd6.fas")

#swap names
mtdnaproteins<-translation_nameswap_set(object=mtdnaproteins, names = mtDNA.data.seq$mtdnaNames, swap = mtDNA.data.seq$treeID)
#reorder to match the nuclear dataset
mtdnaproteins<-mtdnaproteins[consensus.all.timetree$phy$tip.label,]
#convert to character vector
mtdnaproteinschar<-as.character(mtdnaproteins)

#read in the rRNAs for coRdon analysis
#mtdnas_rRNAs<-readSet(file="/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/mtdnas/concat/concat_rRNAs/concat_rRNAs.fas")
mtdnas_rRNAs<-readSet(file="./mtDNA_REASSEMBLY/concat_rRNAs.fas")

#swap names
mtdnas_rRNAs<-translation_nameswap_set(object=mtdnas_rRNAs, names = mtDNA.data.seq$mtdnaNames, swap = mtDNA.data.seq$treeID)
#filter out crocodiles
mtdnas_rRNAs<-mtdnas_rRNAs[names(mtdnas_rRNAs)[names(mtdnas_rRNAs) %in% names(alldatachar)],]
#resort to match tree order
#this is going to be a problem because there are empty rows, figure out later
mtdnas_rRNAs<-mtdnas_rRNAs[consensus.all.timetree$phy$tip.label[consensus.all.timetree$phy$tip.label %in% names(mtdnas_rRNAs)],]
mtdnas_rRNAschar<-as.character(mtdnas_rRNAs)
}

#calculate nucleotide statistics
{
#generate table of codon counts for nuclear data
nuc.codon_table<-codonTable(exons)
#calcualte various codon usage statistics
nuc.milc<-MILC(nuc.codon_table)
nuc.b<-B(nuc.codon_table)
nuc.enc<- ENC(nuc.codon_table)
nuc.encprime<-ENCprime(nuc.codon_table)
nuc.mcb<-MCB(nuc.codon_table)
nuc.scuo<- SCUO(nuc.codon_table)
}
  
#calcuate nucleotide statistics (per exon)
{
  #generate table of codon counts for nuclear data
  #nuc.codon_table.separate<-pbmclapply(exons_separate, codonTable, mc.cores=3)
  #saveRDS(nuc.codon_table.separate, file="./RDS/nuc.codon_table.separate.RDS")
  nuc.codon_table.separate<-readRDS("./RDS/nuc.codon_table.separate.RDS")
  
  #calculate various codon usage statistics
  ENCprime.alt<-function(input){
    result<-log(ENCprime(input)[,1])
    result<-as.data.frame(setNames(result,input@ID))
  }
  SCUO.alt<-function(input){
    result<-log(SCUO(input))
    result<-as.data.frame(setNames(result, input@ID))
  }
  #calcualte values
  {
  #nuc.encprime.separate<-pbmclapply(nuc.codon_table.separate, ENCprime.alt, mc.cores=3)
  #saveRDS(nuc.encprime.separate, file="./RDS/nuc.encprime.separate.RDS")
  nuc.encprime.separate<-readRDS(file="./RDS/nuc.encprime.separate.RDS")
  
  for(i in 1:length(nuc.encprime.separate)){
    names(nuc.encprime.separate[[i]])<-names(nuc.encprime.separate[i])
    nuc.encprime.separate[[i]]$names<-rownames(nuc.encprime.separate[[i]])
  }
  
  #nuc.scuo.separate<- pbmclapply(nuc.codon_table.separate, SCUO.alt, mc.cores=3)
  #saveRDS(nuc.scuo.separate, file="./RDS/nuc.scuo.separate.RDS")
  nuc.scuo.separate<-readRDS(file="./RDS/nuc.scuo.separate.RDS")
  
  for(i in 1:length(nuc.scuo.separate)){
    names(nuc.scuo.separate[[i]])<-names(nuc.scuo.separate[i])
    nuc.scuo.separate[[i]]$names<-rownames(nuc.scuo.separate[[i]])
  }
  }

  nuc.encprime.separate <- Reduce(function(x, y) merge(x, y, by="names", all=TRUE), nuc.encprime.separate)
  rownames(nuc.encprime.separate)<-nuc.encprime.separate$names
  nuc.encprime.separate$names<-NULL
  nuc.encprime.separate<-nuc.encprime.separate[consensus.all.timetree$phy$tip.label,]
  
  nuc.scuo.separate <- Reduce(function(x, y) merge(x, y, by="names", all=TRUE), nuc.scuo.separate)
  rownames(nuc.scuo.separate)<-nuc.scuo.separate$names
  nuc.scuo.separate$names<-NULL
  nuc.scuo.separate<-nuc.scuo.separate[consensus.all.timetree$phy$tip.label,]
  
}

#generate table of codon counts for mtDNA
#generate statistics for mtDNA
{mtdna.codon_table<-codonTable(mtdnaproteins)
#calcualte various codon usage statistics
{mtdna.milc<-MILC(mtdna.codon_table)
mtdna.b<-B(mtdna.codon_table)
mtdna.enc<- ENC(mtdna.codon_table)
mtdna.encprime<-ENCprime(mtdna.codon_table)
mtdna.mcb<-MCB(mtdna.codon_table)
mtdna.scuo<- SCUO(mtdna.codon_table)}
}

#combine these into a data frame
{mol_stats<-data.frame(nuc.milc=nuc.milc[,1], nuc.b=nuc.b[,1], nuc.enc=nuc.enc, nuc.encprime=nuc.encprime[,1], nuc.mcb=nuc.mcb[,1], nuc.scuo=nuc.scuo, 
                      mtdna.milc=mtdna.milc[,1], mtdna.b=mtdna.b[,1], mtdna.enc=mtdna.enc, mtdna.encprime=mtdna.encprime[,1], mtdna.mcb=mtdna.mcb[,1], mtdna.scuo=mtdna.scuo)
rownames(mol_stats)<-names(alldatachar)}

#calculate GC content on the exon data
{#nuc.gc1<-pbmclapply(strsplit(exonchar, split=""), GC1, exact=T, mc.cores = 4)
#saveRDS(nuc.gc1, file="./RDS/nuc.gc1.RDS")
nuc.gc1<-readRDS(file="./RDS/nuc.gc1.RDS")

#nuc.gc2<-pbmclapply(strsplit(exonchar, split=""), GC2, exact=T, mc.cores = 4)
#saveRDS(nuc.gc2, file="./RDS/nuc.gc2.RDS")
nuc.gc2<-readRDS(file="./RDS/nuc.gc2.RDS")

#nuc.gc3<-pbmclapply(strsplit(exonchar, split=""), GC3, exact=T, mc.cores = 4)
#saveRDS(nuc.gc3, file="./RDS/nuc.gc3.RDS")
nuc.gc3<-readRDS(file="./RDS/nuc.gc3.RDS")

#nuc.exonGC<-pbmclapply((exonchar), function(dat){GC(s2c(dat), exact=T)}, mc.cores = 4)
#saveRDS(nuc.exonGC, file="./RDS/nuc.exonGC.RDS")
nuc.exonGC<-readRDS(file="./RDS/nuc.exonGC.RDS")
}

#calculate intron GC
{#nuc.intronGC<-pbmclapply((intronchar), function(dat){GC(s2c(dat), exact=T)}, mc.cores = 4)
#saveRDS(nuc.intronGC, file="./RDS/nuc.intronGC.RDS")
nuc.intronGC<-readRDS(file="./RDS/nuc.intronGC.RDS")
#hist(unlist(intronGC))
}

#calculate utr GC
{#nuc.utrGC<-pbmclapply((utrchar), function(dat){GC(s2c(dat), exact=T)}, mc.cores = 4)
#saveRDS(nuc.utrGC, file="./RDS/nuc.utrGC.RDS")
nuc.utrGC<-readRDS(file="./RDS/nuc.utrGC.RDS")
#hist(unlist(utrGC))
}

#calcuate GC on the "whole" dataset
{#nuc.allGC<-pbmclapply((alldatachar), function(dat){GC(s2c(dat), exact=T)}, mc.cores=4)
#saveRDS(nuc.allGC, file="./RDS/nuc.allGC.RDS")
nuc.allGC<-readRDS(file="./RDS/nuc.allGC.RDS")
}

#calculate GC content on the mtdna protein data
{mtdna.protein.gc1<-pbmclapply(strsplit(mtdnaproteinschar, split=""), GC1, exact=T, mc.cores = 3)
mtdna.protein.gc2<-pbmclapply(strsplit(mtdnaproteinschar, split=""), GC2, exact=T, mc.cores = 3)
mtdna.protein.gc3<-pbmclapply(strsplit(mtdnaproteinschar, split=""), GC3, exact=T, mc.cores = 3)
mtdna.protein.gc<-pbmclapply((mtdnaproteinschar), function(dat){GC(s2c(dat), exact=T)}, mc.cores = 3)
}

#calculate rRNA
{mtdna.rrnaGC<-pbmclapply((mtdnas_rRNAschar), function(dat){GC(s2c(dat), exact=T)}, mc.cores = 3)
}

#calcuate GC on the "whole" dataset
{mtdna.allGC<-mclapply((mtdnaschar), function(dat){GC(s2c(dat), exact=T)}, mc.cores=3)
}

#calculate A T G C composition for each data type
{
nuc.exonACGT <- base_comps(exons)
nuc.intronACGT <- base_comps(introns)
nuc.utrACGT <- base_comps(utrs)
mtdna.ACGT <- base_comps(mtdnas)
}

# tmp<-as.data.frame(base_comps(exons))
# abs(mean(tmp$a+tmp$t)-mean(tmp$g+tmp$c)) #SW mean difference is about 3%
# abs(mean(tmp$g+tmp$t)-mean(tmp$a+tmp$c)) #KM mean difference is about 6%
# abs(mean(tmp$g+tmp$a)-mean(tmp$c+tmp$t)) #RY mean difference is about 4%

#add gc content to master data frame
{
mol_stats$nuc.gc1<-unlist(nuc.gc1)
mol_stats$nuc.gc2<-unlist(nuc.gc2)
mol_stats$nuc.gc3<-unlist(nuc.gc3)
mol_stats$nuc.allGC<-unlist(nuc.allGC)
}

#add the intron and utr GC to the dataset
{
mol_stats$nuc.exonGC<-unlist(nuc.exonGC)
mol_stats$nuc.intronGC<-unlist(nuc.intronGC)
mol_stats$nuc.utrGC<-unlist(nuc.utrGC)
}

#add gc contents from the mtDNA
{
mol_stats$mtdna.protein.gc1<-unlist(mtdna.protein.gc1)
mol_stats$mtdna.protein.gc2<-unlist(mtdna.protein.gc2)
mol_stats$mtdna.protein.gc3<-unlist(mtdna.protein.gc3)
mol_stats$mtdna.protein.gc<-unlist(mtdna.protein.gc)

#mol_stats$scuo.log<-log(mol_stats$scuo)
#mol_stats$label<-rownames(mol_stats)
#mol_stats$encprime.log<-log(mol_stats$encprime)
}

#check names for intron and utr gc
#rownames(mol_stats) == names(unlist(nuc.gc)) #exongc, OK
#rownames(mol_stats) == names(unlist(nuc.intronGC))#, OK
#rownames(mol_stats) == names(unlist(nuc.utrGC))#, OK

#generate merged data frame with the model shifts, and nucleotide stats
{
tip.data.molstats<-merge(x=tip.data, y=mol_stats, by="row.names", sort=F)
rownames(tip.data.molstats)<-tip.data.molstats$Row.names
tip.data.molstats$Row.names<-NULL
}

#notes
{
#namecheck
#rownames(tip.data.molstats) == consensus.all.timetree$phy$tip.label


#generate merged data frame with the model shifts, and more nucleotide stats
#consensus.all.data.molstats<-merge(x=consensus.all.data.molstats, y=div_calcs, by="label")
#rownames(consensus.all.data.molstats)<-consensus.all.data.molstats$label


#reorder to match tree order
#consensus.all.data.molstats <- consensus.all.data.molstats[consensus.all.timetree$phy$tip.label,]


# 
# #checking some plots
# plot(tip.data.molstats$nuc.gc3 ~ tip.data.molstats$nuc.scuo)
# #plot(log(tip.data.molstats$nuc_div) ~ consensus.all.data.molstats$scuo)
# #plot(log(consensus.all.data.molstats$nuc_div) ~ consensus.all.data.molstats$gc3)
# plot(tip.data.molstats$nuc.intronGC ~ tip.data.molstats$nuc.gc3)
# plot(tip.data.molstats$nuc.utrGC ~ tip.data.molstats$nuc.gc3)


  
}
}

#Section 6
##########################################
### phylogenetic signal calculations #####
##########################################
{
  #blomberg's K
  #trying with timetree
  
  phylosig(tree=consensus.all.timetree$phy, x=log(setNames(tip.data.molstats$nuc_div, rownames(tip.data.molstats))),  method="K", test=T, nsim=10000)
  #Phylogenetic signal K : 0.00467956
  #P-value (based on 10000 randomizations) : 0.5021
  
  #trying with phylogram
  phylosig(tree=consensus.all@phylo, x=log(setNames(tip.data.molstats$nuc_div, rownames(tip.data.molstats))),  method="K", test=T, nsim=10000)
  #Phylogenetic signal K : 0.423024
  #P-value (based on 10000 randomizations) : 0.0104
  
  #trying with phylogram and picante version of K estimator
  picante::phylosignal(x = log(setNames(tip.data.molstats$nuc_div, rownames(tip.data.molstats))), phy=consensus.all@phylo, reps=10000)
  #K = 0.4230242
  #p = 9.999e-05
  
  #trying lambda estimator
  phylosig(tree=consensus.all.timetree$phy, x=log(setNames(tip.data.molstats$nuc_div, rownames(tip.data.molstats))),  method="lambda", test=T, nsim=10000)
  #lambda =  0.234789
  #p=0.121315
  fitContinuous(consensus.all.timetree$phy, dat = log(setNames(tip.data.molstats$nuc_div, rownames(tip.data.molstats))), model="OU", bounds=list(alpha = c(min = exp(-500), max = exp(10))))
  
  phylosig(tree=consensus.all@phylo, x=log(setNames(tip.data.molstats$nuc_div, rownames(tip.data.molstats))),  method="lambda", test=T, nsim=10000)
  #lambda = 0.873677
  #0.00516263
  #options(mc.cores=8)
  #fitContinuous(consensus.all@phylo, dat = log(setNames(tip.data.molstats$nuc_div, rownames(tip.data.molstats))), model="OU", bounds=list(alpha = c(min = exp(-500), max = exp(3))))
  
  #trying phylogenetic anova to see if patterns of nucleotide diversity
  #reflect anything about the model shifts
  
  #trying with the time tree
  phylANOVA(tree = consensus.all.timetree$phy, x = setNames(tip.data.molstats$exon_models, rownames(tip.data.molstats)), y = setNames(log(tip.data.molstats$nuc_div), rownames(tip.data.molstats)), nsim=10000)
  #no evidence of association between model regimes and patterns of nucleotide diversity
  
  # ANOVA table: Phylogenetic ANOVA
  # 
  # Response: y
  # Sum Sq  Mean Sq  F value Pr(>F)
  # x          3.959784 0.565683 0.934263 0.9874
  # Residual 115.042425 0.605486    
  
  #trying with the phylogram
  phylANOVA(tree = consensus.all@phylo, x = setNames(tip.data.molstats$exon_models, rownames(tip.data.molstats)), y = setNames(log(tip.data.molstats$nuc_div), rownames(tip.data.molstats)), nsim=10000)
  #no evidence of association between model regimes and patterns of nucleotide diversity
  
  # ANOVA table: Phylogenetic ANOVA
  # 
  # Response: y
  # Sum Sq  Mean Sq  F value Pr(>F)
  # x          3.959784 0.565683 0.934263 0.9965
  # Residual 115.042425 0.605486   
  
  {
    #is there phylogenetic signal in codon usage and gc content
    
    # #trying on the timetree, blom K
    # #testing gc
    # phylosig(tree=consensus.all.timetree$phy, x=log(setNames(tip.data.molstats$nuc.allGC, rownames(tip.data.molstats))),  method="K", test=T, nsim=10000)
    # #Phylogenetic signal K : 0.0359857
    # #P-value (based on 10000 randomizations) : 0.0014
    # 
    # #trying with phylogram, blom K
    # phylosig(tree=consensus.all@phylo, x=log(setNames(tip.data.molstats$nuc.allGC, rownames(tip.data.molstats))),  method="K", test=T, nsim=10000)
    # #Phylogenetic signal K : 2.60728
    # #P-value (based on 10000 randomizations) : 1e-04
    # 
    # 
    # #testing gc, pagel's lambda
    # phylosig(tree=consensus.all.timetree$phy, x=log(setNames(tip.data.molstats$nuc.allGC, rownames(tip.data.molstats))),  method="lambda", test=T, nsim=10000)
    # #Phylogenetic signal lambda : 0.884006
    # #logL(lambda) : 633.197
    # #LR(lambda=0) : 137.694
    # #P-value (based on LR test) : 8.50317e-32
    # 
    # #trying with phylogram, pagel's lambda
    # phylosig(tree=consensus.all@phylo, x=log(setNames(tip.data.molstats$nuc.allGC, rownames(tip.data.molstats))),  method="lambda", test=T, nsim=10000)
    # 
    # 
    # #testing gc3 on exons, K
    # #trying on the timetree
    # phylosig(tree=consensus.all.timetree$phy, x=log(setNames(tip.data.molstats$nuc.gc3, rownames(tip.data.molstats))),  method="K", test=T, nsim=10000)
    # #Phylogenetic signal K : 0.0499484
    # #P-value (based on 10000 randomizations) : 2e-04
    # 
    # #trying with phylogram, K
    # phylosig(tree=consensus.all@phylo, x=log(setNames(tip.data.molstats$nuc.gc3, rownames(tip.data.molstats))),  method="K", test=T, nsim=10000)
    # #Phylogenetic signal K : 4.13493
    # #P-value (based on 10000 randomizations) : 1e-04
    # 
    # #testing gc3, pagel's lambda
    # phylosig(tree=consensus.all.timetree$phy, x=log(setNames(tip.data.molstats$nuc.gc3, rownames(tip.data.molstats))),  method="lambda", test=T, nsim=10000)
    # #Phylogenetic signal lambda : 0.947243
    # #logL(lambda) : 488.819
    # #LR(lambda=0) : 196.912
    # #P-value (based on LR test) : 9.85646e-45
    # 
    # #trying with phylogram
    # phylosig(tree=consensus.all@phylo, x=log(setNames(tip.data.molstats$nuc.gc3, rownames(tip.data.molstats))),  method="lambda", test=T, nsim=10000)
    # #Phylogenetic signal lambda : 0.999934
    # #logL(lambda) : 482.937
    # #LR(lambda=0) : 195.021
    # #P-value (based on LR test) : 2.54912e-44
    # 
    # 
    # #trying on the timetree
    # #testing intron gc
    # phylosig(tree=consensus.all.timetree$phy, x=log(setNames(tip.data.molstats$nuc.intronGC, rownames(tip.data.molstats))),  method="K", test=T, nsim=10000)
    # #Phylogenetic signal K : 0.0454585
    # #P-value (based on 10000 randomizations) : 8e-04
    # 
    # #trying with phylogram
    # phylosig(tree=consensus.all@phylo, x=log(setNames(tip.data.molstats$nuc.intronGC, rownames(tip.data.molstats))),  method="K", test=T, nsim=10000)
    # #Phylogenetic signal K : 1.61121
    # #P-value (based on 10000 randomizations) : 1e-04
    # 
    # #testing timetree, lambda
    # phylosig(tree=consensus.all.timetree$phy, x=log(setNames(tip.data.molstats$nuc.intronGC, rownames(tip.data.molstats))),  method="lambda", test=T, nsim=10000)
    # #Phylogenetic signal lambda : 0.784473
    # #logL(lambda) : 582.959
    # #LR(lambda=0) : 72.1018
    # #P-value (based on LR test) : 2.04375e-17
    # 
    # #trying with phylogram, lambda
    # phylosig(tree=consensus.all@phylo, x=log(setNames(tip.data.molstats$nuc.intronGC, rownames(tip.data.molstats))),  method="lambda", test=T, nsim=10000)
    # #Phylogenetic signal lambda : 0.932228
    # #logL(lambda) : 572.926
    # #LR(lambda=0) : 84.4774
    # #P-value (based on LR test) : 3.88623e-20
    # 
    # 
    # #trying on the timetree
    # #testing utr gc, K
    # phylosig(tree=consensus.all.timetree$phy, x=log(setNames(tip.data.molstats$nuc.utrGC, rownames(tip.data.molstats))),  method="K", test=T, nsim=10000)
    # #Phylogenetic signal K : 0.00373311
    # #P-value (based on 10000 randomizations) : 0.5515
    # 
    # #trying with phylogram, K
    # phylosig(tree=consensus.all@phylo, x=log(setNames(tip.data.molstats$nuc.utrGC, rownames(tip.data.molstats))),  method="K", test=T, nsim=10000)
    # #Phylogenetic signal K : 0.519275
    # #P-value (based on 10000 randomizations) : 0.0208
    # 
    # #testing timetree lambda
    # phylosig(tree=consensus.all.timetree$phy, x=log(setNames(tip.data.molstats$nuc.utrGC, rownames(tip.data.molstats))),  method="lambda", test=T, nsim=10000)
    # #Phylogenetic signal lambda : 0.572429
    # #logL(lambda) : 494.655
    # #LR(lambda=0) : 29.3251
    # #P-value (based on LR test) : 6.11963e-08
    # 
    # #trying with phylogram, lambda
    # phylosig(tree=consensus.all@phylo, x=log(setNames(tip.data.molstats$nuc.utrGC, rownames(tip.data.molstats))),  method="lambda", test=T, nsim=10000)
    # #Phylogenetic signal lambda : 0.638472
    # #logL(lambda) : 492.262
    # #LR(lambda=0) : 27.9993
    # #P-value (based on LR test) : 1.21358e-07
  }
  
  #testing codon usage
  
  #trying on the timetree
  #testing scuo, kappa
  phylosig(tree=consensus.all.timetree$phy, x=log(setNames(tip.data.molstats$nuc.scuo, rownames(tip.data.molstats))),  method="K", test=T, nsim=10000)
  #Phylogenetic signal K : 0.112858
  #P-value (based on 10000 randomizations) : 1e-04
  
  #trying with phylogram, kappa
  phylosig(tree=consensus.all@phylo, x=log(setNames(tip.data.molstats$nuc.scuo, rownames(tip.data.molstats))),  method="K", test=T, nsim=10000)
  #Phylogenetic signal K : 2.19923
  #P-value (based on 10000 randomizations) : 1e-04
  
  #testing scuo, lambda
  phylosig(tree=consensus.all.timetree$phy, x=log(setNames(tip.data.molstats$nuc.scuo, rownames(tip.data.molstats))),  method="lambda", test=T, nsim=10000)
  #Phylogenetic signal lambda : 0.886638
  #logL(lambda) : 239.448
  #LR(lambda=0) : 195.492
  #P-value (based on LR test) : 2.01192e-44
  
  #trying with phylogram, lambda
  phylosig(tree=consensus.all@phylo, x=log(setNames(tip.data.molstats$nuc.scuo, rownames(tip.data.molstats))),  method="lambda", test=T, nsim=10000)
  #Phylogenetic signal lambda : 0.999934
  #logL(lambda) : 300.429
  #LR(lambda=0) : 293.928
  #P-value (based on LR test) : 6.92803e-66
  
  
  #testing encprime
  #trying on the timetree
  #testing encprime, kappa
  phylosig(tree=consensus.all.timetree$phy, x=log(setNames(tip.data.molstats$nuc.encprime, rownames(tip.data.molstats))),  method="K", test=T, nsim=10000)
  #Phylogenetic signal K : 0.058319
  #P-value (based on 10000 randomizations) : 0.0066
  
  #trying with phylogram, kappa
  phylosig(tree=consensus.all@phylo, x=log(setNames(tip.data.molstats$nuc.encprime, rownames(tip.data.molstats))),  method="K", test=T, nsim=10000)
  #Phylogenetic signal K : 2.27968
  #P-value (based on 10000 randomizations) : 1e-04
  
  
  #testing encprime
  phylosig(tree=consensus.all.timetree$phy, x=log(setNames(tip.data.molstats$nuc.encprime, rownames(tip.data.molstats))),  method="lambda", test=T, nsim=10000)
  #Phylogenetic signal lambda : 0.878424
  #logL(lambda) : 1000.69
  #LR(lambda=0) : 117.902
  #P-value (based on LR test) : 1.82174e-27
  
  #trying with phylogram
  phylosig(tree=consensus.all@phylo, x=log(setNames(tip.data.molstats$nuc.encprime, rownames(tip.data.molstats))),  method="lambda", test=T, nsim=10000)
  #Phylogenetic signal lambda : 0.999934
  #logL(lambda) : 1028.48
  #LR(lambda=0) : 183.119
  #P-value (based on LR test) : 1.01034e-41
  
  #result: extremely high phylogenetic signal for gc3, scuo, encprime
  
  
  #now testing to see if there is any associations between gc3, scuo, encprime and the model shifts
  # #trying with the all data timetree
  # phylANOVA(tree=consensus.all.timetree$phy, x = setNames(tip.data.molstats$all_nuc_models, rownames(tip.data.molstats)), y=log(setNames(tip.data.molstats$nuc.allGC, rownames(tip.data.molstats))), nsim=10000)
  # # ANOVA table: Phylogenetic ANOVA
  # #
  # # Response: y
  # # Sum Sq  Mean Sq F value Pr(>F)
  # # x        0.003038 0.000506 2.70661 0.0184
  # # Residual 0.035732 0.000187
  # #
  # # P-value based on simulation.
  # 
  # 
  # #trying with phylogram
  # phylANOVA(tree=consensus.all@phylo, x = setNames(tip.data.molstats$all_nuc_models, rownames(tip.data.molstats)), y=log(setNames(tip.data.molstats$nuc.allGC, rownames(tip.data.molstats))), nsim=10000)
  # # ANOVA table: Phylogenetic ANOVA
  # #
  # # Response: y
  # # Sum Sq  Mean Sq F value Pr(>F)
  # # x        0.003038 0.000506 2.70661 0.0317
  # # Residual 0.035732 0.000187
  # #
  # # P-value based on simulation.
  # 
  # 
  # #gc3
  # #testing gc3 in exons by model shifts
  # phylANOVA(tree=consensus.all.timetree$phy, x = setNames(tip.data.molstats$exon_models, rownames(tip.data.molstats)), y=log(setNames(tip.data.molstats$nuc.gc3, rownames(tip.data.molstats))), nsim=10000)
  # # ANOVA table: Phylogenetic ANOVA
  # #
  # # Response: y
  # # Sum Sq  Mean Sq   F value Pr(>F)
  # # x        0.118998 0.017000 30.536072  3e-04
  # # Residual 0.105774 0.000557
  # #
  # # P-value based on simulation.
  # 
  # 
  # #trying with phylogram
  # phylANOVA(tree=consensus.exons@phylo,  x = setNames(tip.data.molstats$exon_models, rownames(tip.data.molstats)), y=log(setNames(tip.data.molstats$nuc.gc3, rownames(tip.data.molstats))), nsim=10000)
  # # ANOVA table: Phylogenetic ANOVA
  # #
  # # Response: y
  # # Sum Sq  Mean Sq   F value Pr(>F)
  # # x        0.118998 0.017000 30.536072  0.004
  # # Residual 0.105774 0.000557
  # #
  # # P-value based on simulation.
  # 
  # 
  # 
  # ## testing intron GC vs model shifts
  # #trying timetree
  # phylANOVA(tree=consensus.all.timetree$phy, x = setNames(tip.data.molstats$intron_models, rownames(tip.data.molstats)), y=log(setNames(tip.data.molstats$nuc.intronGC, rownames(tip.data.molstats))), nsim=10000)
  # # ANOVA table: Phylogenetic ANOVA
  # #
  # # Response: y
  # # Sum Sq  Mean Sq   F value Pr(>F)
  # # x        0.015607 0.003121 19.565445 0.0078
  # # Residual 0.030632 0.000160
  # #
  # # P-value based on simulation.
  # #
  # 
  # 
  # #trying phylogram
  # phylANOVA(tree=consensus.introns@phylo, x = setNames(tip.data.molstats$intron_models, rownames(tip.data.molstats)), y=log(setNames(tip.data.molstats$nuc.intronGC, rownames(tip.data.molstats))), nsim=10000)
  # # ANOVA table: Phylogenetic ANOVA
  # #
  # # Response: y
  # # Sum Sq  Mean Sq   F value Pr(>F)
  # # x        0.015607 0.003121 19.565445 0.0541
  # # Residual 0.030632 0.000160
  # #
  # # P-value based on simulation.
  # 
  # 
  # ## testing utr GC vs model shifts
  # phylANOVA(tree=consensus.all.timetree$phy,  x = setNames(tip.data.molstats$utr_models, rownames(tip.data.molstats)), y=log(setNames(tip.data.molstats$nuc.utrGC, rownames(tip.data.molstats))), nsim=10000)
  # # ANOVA table: Phylogenetic ANOVA
  # #
  # # Response: y
  # # Sum Sq  Mean Sq  F value Pr(>F)
  # # x        0.009715 0.003238 7.738865 0.2698
  # # Residual 0.081183 0.000418
  # #
  # # P-value based on simulation.
  # 
  # 
  # ## testing phylogram
  # phylANOVA(tree=consensus.utrs@phylo,  x = setNames(tip.data.molstats$utr_models, rownames(tip.data.molstats)), y=log(setNames(tip.data.molstats$nuc.utrGC, rownames(tip.data.molstats))), nsim=10000)
  # # ANOVA table: Phylogenetic ANOVA
  # #
  # # Response: y
  # # Sum Sq  Mean Sq  F value Pr(>F)
  # # x        0.009715 0.003238 7.738865 0.4772
  # # Residual 0.081183 0.000418
  # #
  # # P-value based on simulation.
  # 
  # 
  # 
  # #testing scuo vs model shifts
  # #testing with the time tree
  # phylANOVA(tree=consensus.all.timetree$phy, x = setNames(tip.data.molstats$exon_models, rownames(tip.data.molstats)), y=log(setNames(tip.data.molstats$nuc.scuo, rownames(tip.data.molstats))), nsim=10000)
  # # ANOVA table: Phylogenetic ANOVA
  # #
  # # Response: y
  # # Sum Sq  Mean Sq   F value Pr(>F)
  # # x        1.091425 0.155918 17.642515 0.0049
  # # Residual 1.679148 0.008838
  # #
  # # P-value based on simulation.
  # 
  # 
  # #trying with phylogram
  # phylANOVA(tree=consensus.exons@phylo, x = setNames(tip.data.molstats$exon_models, rownames(tip.data.molstats)), y=log(setNames(tip.data.molstats$nuc.scuo, rownames(tip.data.molstats))), nsim=10000)
  # # ANOVA table: Phylogenetic ANOVA
  # #
  # # Response: y
  # # Sum Sq  Mean Sq   F value Pr(>F)
  # # x        1.091425 0.155918 17.642515 0.0753
  # # Residual 1.679148 0.008838
  # #
  # # P-value based on simulation.
  
  
  #testing encprime vs model shifts
  #testing with the time tree
  aov.time.encprime<-phylANOVA(tree=consensus.all.timetree$phy, x = setNames(tip.data.molstats$exon_models, rownames(tip.data.molstats)), y=log(setNames(tip.data.molstats$nuc.encprime, rownames(tip.data.molstats))), nsim=10000, p.adj="BH")
  # ANOVA table: Phylogenetic ANOVA
  #
  # Response: y
  # Sum Sq Mean Sq F value Pr(>F)
  # x        0.000422   6e-05 27.2314  3e-04
  # Residual 0.000420   2e-06
  contMap(tree=consensus.all.timetree$phy, x = log(setNames(tip.data.molstats$nuc.encprime, rownames(tip.data.molstats))), ftype="off")
  
  
  #testing encprime vs model shifts
  #testing with the time tree
  aov.time.enc<-phylANOVA(tree=consensus.all.timetree$phy, x = setNames(tip.data.molstats$exon_models, rownames(tip.data.molstats)), y=log(setNames(tip.data.molstats$nuc.enc, rownames(tip.data.molstats))), nsim=10000, p.adj="BH")
  # ANOVA table: Phylogenetic ANOVA
  # 
  # Response: y
  # Sum Sq  Mean Sq   F value Pr(>F)
  # x        0.004567 0.000652 17.232435 0.0054
  # Residual 0.007194 0.000038  
  contMap(tree=consensus.all.timetree$phy, x = log(setNames(tip.data.molstats$nuc.enc, rownames(tip.data.molstats))), ftype="off")
  
  
  #aov.phyl.encprime<-phylANOVA(tree=consensus.exons@phylo,  x = setNames(tip.data.molstats$exon_models, rownames(tip.data.molstats)), y=log(setNames(tip.data.molstats$nuc.encprime, rownames(tip.data.molstats))), nsim=10000, p.adj="BH")
  # ANOVA table: Phylogenetic ANOVA
  # 
  # Response: y
  # Sum Sq Mean Sq F value Pr(>F)
  # x        0.000422   6e-05 27.2314 0.0093
  # Residual 0.000420   2e-06
  
  
  # #trying manova for exons
  # dat<-as.matrix(tip.data.molstats[,c(9:16)])
  # #dat<-as.matrix(consensus.all.data.molstats[,c(14,16)])
  # grp<-tip.data.molstats$exon_models
  # names(grp)<-rownames(tip.data.molstats)
  # grp<-as.factor(grp)
  # 
  # manova.time<-aov.phylo(formula = log(dat) ~ grp, phy = consensus.all.timetree$phy, nsim = 1000, test="Wilks")
  # print(attributes(manova.time)$summary)
  
  # Multivariate Analysis of Variance Table
  #
  # Response: dat
  # Df    Wilks approx-F num-Df den-Df     Pr(>F) Pr(>F) given phy
  # group       7 0.065495   11.658     56  990.8 2.3022e-75         0.000999 ***
  #   Residuals 190
  # ---
  #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  
  
  # #trying with phylogram
  # manova.phylogram<-aov.phylo(formula = log(dat) ~ grp, phy = consensus.exons@phylo, nsim = 1000, test="Wilks")
  # print(attributes(manova.phylogram)$summary)
  
  
  # Multivariate Analysis of Variance Table
  #
  # Response: dat
  # Df    Wilks approx-F num-Df den-Df     Pr(>F) Pr(>F) given phy
  # group       7 0.065495   11.658     56  990.8 2.3022e-75          0.02098 *
  #   Residuals 190
  # ---
  #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  
  
  
}

#Section 7
##########################################
####### patterns of codon usage ##########
##########################################
{
#testing codon usage across exon model shifts
{

#testing scuo vs model shifts
#testing with the time tree
aov.time.scuo<-phylANOVA(tree=consensus.all.timetree$phy, x = setNames(tip.data.molstats$exon_models, rownames(tip.data.molstats)), y=log(setNames(tip.data.molstats$nuc.scuo, rownames(tip.data.molstats))), nsim=10000, posthoc = T, p.adj = "BH")
# ANOVA table: Phylogenetic ANOVA
#
# Response: y
# Sum Sq  Mean Sq   F value Pr(>F)
# x        1.091425 0.155918 17.642515 0.0049
# Residual 1.679148 0.008838
#
# P-value based on simulation.
contMap(tree=consensus.all.timetree$phy, x = log(setNames(tip.data.molstats$nuc.scuo, rownames(tip.data.molstats))), ftype="off")

#testing encprime vs model shifts
#testing with the time tree
#aov.phyl.scuo<-phylANOVA(tree=consensus.all.timetree$phy, x = setNames(tip.data.molstats$exon_models, rownames(tip.data.molstats)), y=log(setNames(tip.data.molstats$nuc.encprime, rownames(tip.data.molstats))), nsim=10000, posthoc = T)
# ANOVA table: Phylogenetic ANOVA
#
# Response: y
# Sum Sq Mean Sq F value Pr(>F)
# x        0.000422   6e-05 27.2314  3e-04
# Residual 0.000420   2e-06
}

#plotting set up for supplemental codon usage
{
{
  colnames(aov.time.encprime$Pt)<- c('Gallo+Columbea+S', "Aequornithes", "Coraciimorphae", "Passeri", "Tinamiformes", "Psittaciformes", "Passerea", "Sister to Tinamiformes")
  rownames(aov.time.encprime$Pt)<- c('Gallo+Columbea+S', "Aequornithes", "Coraciimorphae", "Passeri", "Tinamiformes", "Psittaciformes", "Passerea", "Sister to Tinamiformes")
  Ncp<-pheatmap(aov.time.encprime$Pt, display_numbers = T, 
                color = rev(colorRampPalette(c(rev(colorRampPalette(c("grey75", "grey100"))(10)), '#0079bc'))(1000)), 
                cluster_rows = F, cluster_cols = F, 
                fontsize_number = 9, angle_col=45, fontsize_row=8, fontsize_col=8,
                number_format="%.2f" , main="phylANOVA Nc' BH p-values",
                legend=F, number_color = 'black', show_rownames = F, show_colnames=F)
  
  colnames(aov.time.enc$Pt)<- c('Gallo+Columbea+S', "Aequornithes", "Coraciimorphae", "Passeri", "Tinamiformes", "Psittaciformes", "Passerea", "Sister to Tinamiformes")
  rownames(aov.time.enc$Pt)<- c('Gallo+Columbea+S', "Aequornithes", "Coraciimorphae", "Passeri", "Tinamiformes", "Psittaciformes", "Passerea", "Sister to Tinamiformes")
  Nc<-pheatmap(aov.time.enc$Pt, display_numbers = T, 
               color = rev(colorRampPalette(c(rev(colorRampPalette(c("grey75", "grey100"))(10)), '#0079bc'))(1000)),
               cluster_rows = F, cluster_cols = F, 
               fontsize_number = 9, angle_col=45, fontsize_row=8, fontsize_col=8,
               number_format="%.2f" , main="phylANOVA Nc BH p-values",
               legend=F, number_color = 'black', show_rownames=F, show_colnames=F)
  
  colnames(aov.time.scuo$Pt)<- c('Gallo+Columbea+S', "Aequornithes", "Coraciimorphae", "Passeri", "Tinamiformes", "Psittaciformes", "Passerea", "Sister to Tinamiformes")
  rownames(aov.time.scuo$Pt)<- c('Gallo+Columbea+S', "Aequornithes", "Coraciimorphae", "Passeri", "Tinamiformes", "Psittaciformes", "Passerea", "Sister to Tinamiformes")
  SCUO<-pheatmap(aov.time.scuo$Pt, display_numbers = T, 
                 color = rev(colorRampPalette(c(rev(colorRampPalette(c("grey75", "grey100"))(10)), '#0079bc'))(1000)), 
                 cluster_rows = F, cluster_cols = F, 
                 fontsize_number = 9, angle_col=45, fontsize_row=8, fontsize_col=8,
                 number_format="%.2f" , main='phylANOVA SCUO BH p-values',
                 legend=T, number_color = 'black', show_colnames=F)
}
{
  colnames(aov.time.encprime$T)<- c('Gallo+Columbea+S', "Aequornithes", "Coraciimorphae", "Passeri", "Tinamiformes", "Psittaciformes", "Passerea", "Sister to Tinamiformes")
  rownames(aov.time.encprime$T)<- c('Gallo+Columbea+S', "Aequornithes", "Coraciimorphae", "Passeri", "Tinamiformes", "Psittaciformes", "Passerea", "Sister to Tinamiformes")
  Ncp.t<-pheatmap(aov.time.encprime$T, display_numbers = T, 
                color = rev(colorRampPalette(c('#c43129', 'grey','white', 'grey', '#0079bc'))(1000)), 
                cluster_rows = F, cluster_cols = F, 
                fontsize_number = 9, angle_col=45, fontsize_row=8, fontsize_col=8,
                number_format="%.2f" , main="phylANOVA Nc' T-values",
                legend=F, number_color = 'black', show_rownames = F)
  
  colnames(aov.time.enc$T)<- c('Gallo+Columbea+S', "Aequornithes", "Coraciimorphae", "Passeri", "Tinamiformes", "Psittaciformes", "Passerea", "Sister to Tinamiformes")
  rownames(aov.time.enc$T)<- c('Gallo+Columbea+S', "Aequornithes", "Coraciimorphae", "Passeri", "Tinamiformes", "Psittaciformes", "Passerea", "Sister to Tinamiformes")
  Nc.t<-pheatmap(aov.time.enc$T, display_numbers = T, 
               color = rev(colorRampPalette(c('#c43129', 'grey','white', 'grey', '#0079bc'))(1000)), 
               cluster_rows = F, cluster_cols = F, 
               fontsize_number = 9, angle_col=45, fontsize_row=8, fontsize_col=8,
               number_format="%.2f" , main="phylANOVA Nc T-values",
               legend=F, number_color = 'black', show_rownames=F)
  
  colnames(aov.time.scuo$T)<- c('Gallo+Columbea+S', "Aequornithes", "Coraciimorphae", "Passeri", "Tinamiformes", "Psittaciformes", "Passerea", "Sister to Tinamiformes")
  rownames(aov.time.scuo$T)<- c('Gallo+Columbea+S', "Aequornithes", "Coraciimorphae", "Passeri", "Tinamiformes", "Psittaciformes", "Passerea", "Sister to Tinamiformes")
  SCUO.t<-pheatmap(aov.time.scuo$T, display_numbers = T, 
                 color = rev(colorRampPalette(c('#c43129', 'grey','white', 'grey', '#0079bc'))(1000)), 
                 cluster_rows = F, cluster_cols = F, 
                 fontsize_number = 9, angle_col=45, fontsize_row=8, fontsize_col=8,
                 number_format="%.2f" , main='phylANOVA SCUO T-values',
                 legend=T, number_color = 'black')
  
}
}

#plotting for supplemental figure codon usage
pdf(file="codon_usage.pdf", width=13, height=9)
ggarrange(Nc$gtable, Ncp$gtable, SCUO$gtable,
          Nc.t$gtable, Ncp.t$gtable, SCUO.t$gtable,
          nrow=2, ncol=3, widths=c(1,1,1.35, 1,1,1.35),
          heights=c(1,1.3))
dev.off()

#plot(simmap.janus.nuc.exons, ftype='off')
#tiplabels(getStates(simmap.janus.nuc.exons, type='tips'), frame = 'none', cex=0.5)

#loop SCUO anovas on individual exons
{
# system.time(
# ANOVs<-multi_anov(exon_stats=nuc.scuo.separate, tree = consensus.all.timetree$phy, models = tip.data.molstats$exon_models, sims = 10000)
# )
#saveRDS(ANOVs, file = './RDS/SCUO.anova.RDS')
ANOVs<- readRDS(file="./RDS/SCUO.anova.RDS")
}

# Ps<-list()
# for(i in 1:length(ANOVs)){
#   Ps[[i]]<-tmp[[i]]$Pf
# }
# names(Ps)<-names(tmp)
# Ps<-unlist(Ps)

# #Ps<-p.adjust(unlist(Ps), method="BH")
# as.data.frame(sort((Ps[Ps < 0.05]))[1:20])
as.data.frame(sort((ANOVs[ANOVs < 0.1])))

#calculate stats
t(apply(nuc.scuo.separate, 1, Stats))[,c(1:2)]

#cbind(t(apply(nuc.scuo.separate, 1, Stats))[,c(1)] ,log(nuc.scuo))
#contMap(tree=consensus.all.timetree$phy, x = setNames(nuc.scuo.separate[,1],rownames(nuc.scuo.separate)), ftype="off")
#contMap(tree=consensus.all.timetree$phy, x = exp(setNames(nuc.scuo,names(exons))), ftype="off")
#contMap(tree=consensus.all.timetree$phy, x = t(apply(nuc.scuo.separate, 1, Stats))[,c(2)], ftype="off")
}

#Section 8
##########################################
####### LHT and imputation set up#########
##########################################
{
#### analysis of life history data across the model shifts ####
{
#read in
#LHT<-read.table(file="/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/unmasked/min2x/qual20_cov2_haplo_bestonly/initial_test_filters/min50bp_min10p_aligned/ALIGNED/phased/sample1haplo/phyhetnucbf/R/bird_et_al_2020/LHT_data_11_1_21.txt", header = T, sep=",")
LHT<-read.table(file="./LHT_DATA/LHT_data_11_1_21.txt", header = T, sep=",")

#cleaned up
LHT<-data.frame(species=LHT$treeID, latitude=LHT$Range.centroid.of.latitude, 
                mass=LHT$Body.mass, mean.clutch=LHT$Mean.clutch.size, 
                freshwater=LHT$Freshwater, marine=LHT$Marine, 
                terrestrial=LHT$Terrestrial, min_altitude=LHT$Minimum.altitude, 
                max_altitude=LHT$Maximum.altitude, migrant=LHT$Migrant, 
                diet=LHT$Diet, nocturnal=LHT$Nocturnal, forest=LHT$Forest, 
                gen_length=LHT$GenLength.modeled., survival=LHT$Survival, 
                breeding=LHT$Age.at.first.breeding, longevity=LHT$Max.longevity, 
                chickPC1=LHT$chickPC1, hatchlingPC1=LHT$hatchlingPC1,
                BMR.species.OTU = LHT$BMR.species.OTU,
                BMR.species.Metabolic.rate.Ryan=LHT$BMR.species.Metabolic.rate.Ryan,
                BMR.species.Body.Mass = LHT$BMR.species.Body.Mass,
                BMR.species.bmr.sd= LHT$BMR.species.bmr.sd ,
                BMR.species.n.sd= LHT$BMR.species.n.sd,
                BMR.species.se.sd = LHT$BMR.species.se.sd,
                BMR.species.cov.sd= LHT$BMR.species.cov.sd,
                BMR.species.sd.bmr.est= LHT$BMR.species.sd.bmr.est,
                BMR.species.sd.mass.est = LHT$BMR.species.sd.mass.est,
                BMR.species.se.bmr.est.1 = LHT$BMR.species.se.bmr.est.1,
                BMR.species.se.mass.est.1 = LHT$BMR.species.se.mass.est.1,
                BMR.species.se.bmr.est.10 = LHT$BMR.species.se.bmr.est.10,
                BMR.species.se.mass.est.10= LHT$BMR.species.se.mass.est.10,
                BMR.species.se.bmr.est.20 = LHT$BMR.species.se.bmr.est.20,
                BMR.species.se.mass.est.20= LHT$BMR.species.se.mass.est.20,
                BMR.species.key = LHT$BMR.species.key) #a few of these entries were experimental and not used

#resort LHT to match tree order
rownames(LHT)<-LHT$species
LHT<-LHT[consensus.all.timetree$phy$tip.label,]
#LHT$species
#LHT$bmr.species
LHT$min_altitude<-LHT$min_altitude+1 #adding 1 to all to get rid of zeros
LHT$max_altitude<-LHT$max_altitude+1 #adding 1 to all to get rid of zeros

megaLHT<-merge(x=LHT, y = tip.data.molstats, by.x="species", by.y="label")
#resort LHT to match tree order
rownames(megaLHT)<-megaLHT$species
megaLHT<-megaLHT[consensus.all.timetree$phy$tip.label,]

#manual correction for terrestrial
megaLHT$terrestrial[which(is.na(megaLHT$terrestrial))]<-"y"
#manual corection for diet
megaLHT$diet[which(is.na(megaLHT$diet))]<-c("FruiNect", "Invertebrate", "Invertebrate")
#manual correction for nocturnal
megaLHT$nocturnal[which(is.na(megaLHT$nocturnal))]<-0

#manual type change (for filtering later)
megaLHT$nocturnal<-as.integer(megaLHT$nocturnal)
#str(megaLHT)
}

#this following dataset is formatted for input into Rphylopars, includes 
#setting up the dataset, only for quantitative traits
{
#the data are log transformed
nums <- unlist(lapply(megaLHT, is.double))

phylopars.data<-megaLHT[,nums]

#take the absolute value for latitude
phylopars.data$latitude<-abs(phylopars.data$latitude)
#log transform the right varibales
phylopars.data<-cbind(log(phylopars.data[,c(1:9)]), phylopars.data[,c(10:11)], log(phylopars.data[,c(12:48)]))
#scale each column separately
phylopars.data <- as.data.frame(lapply(phylopars.data, function(x) c(scale(x)))) ### DATA ARE SCALED AFTER
#add column for species
phylopars.data$species<-megaLHT$species
#reorder the columns
phylopars.data<-cbind(species=phylopars.data$species, phylopars.data[,c(1:47)])

#clean the NaNs
phylopars.data[is.nan(phylopars.data)] <- NA

#remove columns that are all NA
phylopars.data<-phylopars.data[,which(unlist(lapply(phylopars.data, function(x) !all(is.na(x)))))]
}

# phylopars with univariate/multivariate imputation
{

#investigate using imputation with error and phenocov_list
#bm.fit.nocorr <- phylopars(trait_data = phylopars.data[,c("species", "mass", "mean.clutch", "gen_length", "survival", "breeding","longevity", "chickPC1", "latitude")], tree = consensus.all.timetree$phy, model = "BM", phylo_correlated = F)
#bm.fit.corr <- phylopars(trait_data = phylopars.data[,c("species","mass", "mean.clutch", "gen_length", "survival", "breeding","longevity", "chickPC1", "latitude")], tree = consensus.all.timetree$phy, model = "BM", phylo_correlated = T)
#ou.fit.corr <- phylopars(trait_data = phylopars.data[,c("species","mass", "mean.clutch", "gen_length", "survival", "breeding","longevity", "chickPC1", "latitude")], tree = consensus.all.timetree$phy, model = "mvOU", phylo_correlated = T)
#ou.fit.nocorr <- phylopars(trait_data = phylopars.data[,c("species","mass", "mean.clutch", "gen_length", "survival", "breeding","longevity", "chickPC1", "latitude")], tree = consensus.all.timetree$phy, model = "OU", phylo_correlated = F)

bm.fit.nocorr<-readRDS(file="./RDS/bm.fit.nocorr.RDS")
bm.fit.corr<-readRDS(file="./RDS/bm.fit.corr.RDS")
ou.fit.corr<-readRDS(file="./RDS/ou.fit.corr.RDS")
ou.fit.nocorr<-readRDS(file="./RDS/ou.fit.nocorr.RDS")

AIC(bm.fit.nocorr)
AIC(bm.fit.corr)
AIC(ou.fit.corr)

}

# pcatest<-phyl.pca(bm.fit.corr$anc_recon[1:198,], tree=consensus.all.timetree$phy, method = 'lambda')
# plot(pcatest$S[,c(2,3)], col= megaLHT$aggregate_models, pch=16)
# rownames(bm.fit.corr$anc_recon[1:198,]) == rownames(megaLHT)
}


#fix this next section depending on the results from l1ou bootstrapping
#TESTING -- DO NOT RUN
# {
# ou.fit.corr.test <- phylopars(trait_data = phylopars.data[,c("species","mass", "mean.clutch", "gen_length", "survival", "breeding","longevity", "chickPC1", "latitude")], tree = consensus.all.timetree$phy, model = "OU", phylo_correlated = T)
# 
# mvMORPH.input<-phylopars.data[,c("species","mass", "mean.clutch", "gen_length", "survival", "breeding","longevity", "chickPC1", "latitude")]
# rownames(mvMORPH.input)<- mvMORPH.input$species
# mvMORPH.input$species <- NULL
# bm.fit.mvbm<-mvBM(data = mvMORPH.input, tree = consensus.all.timetree$phy, model="BM1")
# saveRDS(bm.fit.mvbm, file="./RDS/bm.fit.mvbm.RDS")
# bm.fit.mvou<-mvOU(data = mvMORPH.input, tree = consensus.all.timetree$phy, model="OU1", method="inverse")
# saveRDS(bm.fit.mvou, file="./RDS/bm.fit.mvou.RDS")
# bm.fit.mvou<-mvOU(data = mvMORPH.input, tree = consensus.all.timetree$phy, model="OU1", method="inverse")
# }

#Section 9
############################################################
######## setting up more analyses for mvMORPH/OUwie#########
############################################################
{
# need a simmap style trees with the mapped regimes
simmap.janus.nuc.alldata<-janus_simmap(reference = time_scale, target = consensus.all)
#now add node labels. the order of dat on the tree should reflect the 
#model labels from janus
simmap.janus.nuc.alldata$node.label<-getStates(simmap.janus.nuc.alldata, type = "nodes")

simmap.janus.nuc.exons<-janus_simmap(reference = time_scale, target = consensus.exons)
#now add node labels. the order of dat on the tree should reflect the 
#model labels from janus
simmap.janus.nuc.exons$node.label<-getStates(simmap.janus.nuc.exons, type = "nodes")

simmap.janus.nuc.introns<-janus_simmap(reference = time_scale, target = consensus.introns)
#now add node labels. the order of dat on the tree should reflect the 
#model labels from janus
simmap.janus.nuc.introns$node.label<-getStates(simmap.janus.nuc.introns, type = "nodes")


simmap.janus.nuc.utrs<-janus_simmap(reference = time_scale, target = consensus.utrs)
#now add node labels. the order of dat on the tree should reflect the 
#model labels from janus
simmap.janus.nuc.utrs$node.label<-getStates(simmap.janus.nuc.utrs, type = "nodes")


simmap.janus.mtdna.all<-janus_simmap(reference=time_scale, target = consensus.mtdnas.all)
#now add node labels. the order of dat on the tree should reflect the 
#model labels from janus
simmap.janus.mtdna.all$node.label<-getStates(simmap.janus.mtdna.all, type = "nodes")


simmap.janus.mtdna.proteins<-janus_simmap(reference=time_scale, target = consensus.mtdnas.proteins)
#now add node labels. the order of dat on the tree should reflect the 
#model labels from janus
simmap.janus.mtdna.proteins$node.label<-getStates(simmap.janus.mtdna.proteins, type = "nodes")


simmap.janus.mtdna.rrnas<-janus_simmap(reference=time_scale, target = consensus.mtdnas.rRNAs)
#now add node labels. the order of dat on the tree should reflect the 
#model labels from janus
simmap.janus.mtdna.rrnas$node.label<-getStates(simmap.janus.mtdna.rrnas, type = "nodes")
#plot(simmap.janus.mtdna.rrnas)


###generate aggretate shift input in simmap format
simmap.janus.nuc.aggregate<-consensus.all.timetree$phy
#plot(simmap.janus.aggregate, cex=0.001)

#set everything to state zero
simmap.janus.nuc.aggregate <- paintSubTree(tree = simmap.janus.nuc.aggregate, node = 199, state = as.character(0))
#plot(simmap.janus.aggregate, ftype="off")

#define the maps for aggretage nuclear data signal
remaps<-con_prop_logit[con_prop_logit$uncex.merged==1,]$node.all

for(i in 1:length(remaps)){
  simmap.janus.nuc.aggregate <- paintSubTree(tree = simmap.janus.nuc.aggregate, node = remaps[[i]], state = as.character(i), anc.state = as.character(i), stem=T)
}

plot(simmap.janus.nuc.aggregate, ftype="off", colors=setNames(c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928'), as.character(seq(0,11))))
nodelabels(node=remaps)


#11 shifts, 12 models total

check.identify(phy = simmap.janus.nuc.aggregate, data=data.frame(names=rownames(LHT),LHT$mass), simmap.tree = T)
#The regime optima are identifiable. 
#[1] 1


# ###generate aggretate shift input in simmap format (simplifed for OUwie)
{# simmap.janus.nuc.aggregate.simplified<-consensus.all.timetree$phy
# #plot(simmap.janus.aggregate, cex=0.001)
# 
# #set everything to state zero
# simmap.janus.nuc.aggregate.simplified <- paintSubTree(tree = simmap.janus.nuc.aggregate.simplified, node = 199, state = as.character(as.english(0)))
# #plot(simmap.janus.aggregate, ftype="off")
# 
# #define the maps for aggretage nuclear data signal
# remaps<-con_prop_logit[con_prop_logit$uncex.merged==1,]$node.all
# remaps<-remaps[-c(1, 9)]
# 
# for(i in 1:length(remaps)){
#   simmap.janus.nuc.aggregate.simplified <- paintSubTree(tree = simmap.janus.nuc.aggregate.simplified, node = remaps[[i]], state = as.character(as.english(i)), anc.state = as.character(as.english(i)), stem=T)
# }
}
{
simmap.janus.nuc.aggregate.simplified<-paintSubTree(simmap.janus.nuc.exons, node = 199, state = 'root_gallo')
#simmap.janus.nuc.aggregate.simplified<-paintSubTree(simmap.janus.nuc.aggregate.simplified, node = 389, state = 'notopaleognathae', anc.state = 'notopaleognathae', stem=T)
simmap.janus.nuc.aggregate.simplified<-paintSubTree(simmap.janus.nuc.aggregate.simplified, node = 393, state = 'tinamiformes', anc.state = 'tinamiformes', stem=T)
simmap.janus.nuc.aggregate.simplified<-paintSubTree(simmap.janus.nuc.aggregate.simplified, node = 390, state = 'tinamiformes_sister', anc.state = 'tinamiformes_sister', stem=T)
simmap.janus.nuc.aggregate.simplified<-paintSubTree(simmap.janus.nuc.aggregate.simplified, node = 203, state = 'otidae_sister', anc.state = 'otidae_sister', stem=T)
simmap.janus.nuc.aggregate.simplified<-paintSubTree(simmap.janus.nuc.aggregate.simplified, node = 363, state = 'columbea', anc.state = 'columbea', stem=T)
simmap.janus.nuc.aggregate.simplified<-paintSubTree(simmap.janus.nuc.aggregate.simplified, node = 344, state = 'otidae', anc.state = 'otidae', stem=T)
simmap.janus.nuc.aggregate.simplified<-paintSubTree(simmap.janus.nuc.aggregate.simplified, node = 297, state = 'aeqornithes', anc.state = 'aeqornithes', stem=T)
simmap.janus.nuc.aggregate.simplified<-paintSubTree(simmap.janus.nuc.aggregate.simplified, node = 269, state = 'coraciimorphae', anc.state = 'coraciimorphae', stem=T)
simmap.janus.nuc.aggregate.simplified<-paintSubTree(simmap.janus.nuc.aggregate.simplified, node = 256, state = 'psittaciformes', anc.state = 'psittaciformes', stem=T)
simmap.janus.nuc.aggregate.simplified<-paintSubTree(simmap.janus.nuc.aggregate.simplified, node = 242, state = 'passeri', anc.state = 'passeri', stem=T)

#nodelabels()
plot(simmap.janus.nuc.aggregate.simplified, ftype="off", colors=setNames(sample(rainbow(10)), unique(getStates(simmap.janus.nuc.aggregate.simplified, type = 'tips'))))
#tiplabels(text = getStates(simmap.janus.nuc.aggregate.simplified, type = 'tips'), frame='none')

}


{
  simmap.janus.nuc.aggregate.simplified.exons<-paintSubTree(simmap.janus.nuc.exons, node = 199, state = 'root_gallo')
  #simmap.janus.nuc.aggregate.simplified<-paintSubTree(simmap.janus.nuc.aggregate.simplified, node = 389, state = 'notopaleognathae', anc.state = 'notopaleognathae', stem=T)
  simmap.janus.nuc.aggregate.simplified.exons<-paintSubTree(simmap.janus.nuc.aggregate.simplified.exons, node = 393, state = 'tinamiformes', anc.state = 'tinamiformes', stem=T)
  simmap.janus.nuc.aggregate.simplified.exons<-paintSubTree(simmap.janus.nuc.aggregate.simplified.exons, node = 390, state = 'tinamiformes_sister', anc.state = 'tinamiformes_sister', stem=T)
  simmap.janus.nuc.aggregate.simplified.exons<-paintSubTree(simmap.janus.nuc.aggregate.simplified.exons, node = 202, state = 'passerea', anc.state = 'passerea', stem=T)
  simmap.janus.nuc.aggregate.simplified.exons<-paintSubTree(simmap.janus.nuc.aggregate.simplified.exons, node = 297, state = 'aeqornithes', anc.state = 'aeqornithes', stem=T)
  simmap.janus.nuc.aggregate.simplified.exons<-paintSubTree(simmap.janus.nuc.aggregate.simplified.exons, node = 269, state = 'coraciimorphae', anc.state = 'coraciimorphae', stem=T)
  simmap.janus.nuc.aggregate.simplified.exons<-paintSubTree(simmap.janus.nuc.aggregate.simplified.exons, node = 256, state = 'psittaciformes', anc.state = 'psittaciformes', stem=T)
  simmap.janus.nuc.aggregate.simplified.exons<-paintSubTree(simmap.janus.nuc.aggregate.simplified.exons, node = 242, state = 'passeri', anc.state = 'passeri', stem=T)
  
  #nodelabels()
  plot(simmap.janus.nuc.aggregate.simplified.exons, ftype="off", colors=setNames(sample(rainbow(10)), unique(getStates(simmap.janus.nuc.aggregate.simplified.exons, type = 'tips'))))
  #tiplabels(text = getStates(simmap.janus.nuc.aggregate.simplified, type = 'tips'), frame='none')
  
}



{
  simmap.janus.nuc.aggregate.simplified2<-paintSubTree(simmap.janus.nuc.exons, node = 199, state = 'root_struthio')
  #simmap.janus.nuc.aggregate.simplified<-paintSubTree(simmap.janus.nuc.aggregate.simplified, node = 389, state = 'notopaleognathae', anc.state = 'notopaleognathae', stem=T)
  simmap.janus.nuc.aggregate.simplified2<-paintSubTree(simmap.janus.nuc.aggregate.simplified2, node = 200, state = 'neognathae_gallo', anc.state = 'neognathae_gallo', stem=T)
  simmap.janus.nuc.aggregate.simplified2<-paintSubTree(simmap.janus.nuc.aggregate.simplified2, node = 393, state = 'tinamiformes', anc.state = 'tinamiformes', stem=T)
  simmap.janus.nuc.aggregate.simplified2<-paintSubTree(simmap.janus.nuc.aggregate.simplified2, node = 390, state = 'tinamiformes_sister', anc.state = 'tinamiformes_sister', stem=T)
  simmap.janus.nuc.aggregate.simplified2<-paintSubTree(simmap.janus.nuc.aggregate.simplified2, node = 203, state = 'otidae_sister', anc.state = 'otidae_sister', stem=T)
  simmap.janus.nuc.aggregate.simplified2<-paintSubTree(simmap.janus.nuc.aggregate.simplified2, node = 363, state = 'columbea', anc.state = 'columbea', stem=T)
  simmap.janus.nuc.aggregate.simplified2<-paintSubTree(simmap.janus.nuc.aggregate.simplified2, node = 344, state = 'otidae', anc.state = 'otidae', stem=T)
  simmap.janus.nuc.aggregate.simplified2<-paintSubTree(simmap.janus.nuc.aggregate.simplified2, node = 297, state = 'aeqornithes', anc.state = 'aeqornithes', stem=T)
  simmap.janus.nuc.aggregate.simplified2<-paintSubTree(simmap.janus.nuc.aggregate.simplified2, node = 269, state = 'coraciimorphae', anc.state = 'coraciimorphae', stem=T)
  simmap.janus.nuc.aggregate.simplified2<-paintSubTree(simmap.janus.nuc.aggregate.simplified2, node = 256, state = 'psittaciformes', anc.state = 'psittaciformes', stem=T)
  simmap.janus.nuc.aggregate.simplified2<-paintSubTree(simmap.janus.nuc.aggregate.simplified2, node = 242, state = 'passeri', anc.state = 'passeri', stem=T)
  
  #nodelabels()
  plot(simmap.janus.nuc.aggregate.simplified2, ftype="off", colors=setNames(sample(rainbow(11)), unique(getStates(simmap.janus.nuc.aggregate.simplified2, type = 'tips'))))
  #tiplabels(text = getStates(simmap.janus.nuc.aggregate.simplified, type = 'tips'), frame='none')
  
}


#generate aggregate map for nuclear + mtDNA

###generate aggretate shift input in simmap format
simmap.janus.all.aggregate<-consensus.all.timetree$phy
#plot(simmap.janus.aggregate, cex=0.001)

#set everything to state zero
simmap.janus.all.aggregate <- paintSubTree(tree = simmap.janus.all.aggregate, node = 199, state = as.character(0))
#plot(simmap.janus.aggregate, ftype="off")

#define the maps for aggregate data signal in nuc + mtDNA
remaps<-con_prop_logit[con_prop_logit$uncex.merged.mtdnas==1,]$node.all

for(i in 1:length(remaps)){
  simmap.janus.all.aggregate <- paintSubTree(tree = simmap.janus.all.aggregate, node = remaps[[i]], state = as.character(i), anc.state = as.character(i), stem=T)
}

#12 shifts, 13 models  
plot(simmap.janus.all.aggregate, ftype="off", colors=setNames(c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928', '#8dd3c7'), as.character(seq(0,12))), lwd=3)
nodelabels(node=remaps)


#generate null model with one map
simmap.janus.null<-consensus.all.timetree$phy
#plot(simmap.janus.aggregate, cex=0.001)

#set everything to state zero
simmap.janus.null <- paintSubTree(tree = simmap.janus.null, node = 199, state = as.character(0))
#plot(simmap.janus.aggregate, ftype="off")

####
check.identify(phy = simmap.janus.all.aggregate, data=data.frame(names=rownames(LHT),LHT$mass), simmap.tree = T)
#The regime optima are identifiable. 
#[1] 1


# 
# require(phyloch)
# 
# par(mfrow=c(1,5))
# #plot_jeans_time_rect_simmap(simmap.janus.nuc.alldata, title = "all", width=1.5)
# plot_jeans_time_rect_simmap(simmap.janus.nuc.exons, title = "exons", width=2.0, out.line=F, box=F)
# plot_jeans_time_rect_simmap(simmap.janus.nuc.introns, title = "introns", width=2.0, out.line=F, box=F)
# plot_jeans_time_rect_simmap(simmap.janus.nuc.utrs, title = "utrs", width=2.0, out.line=F, box=F)
# plot_jeans_time_rect_simmap(simmap.janus.mtdna.all, title = "mtdna.all", width=2.0, out.line=F, box=F)
# #plot_jeans_time_rect_simmap(simmap.janus.mtdna.proteins, title = "mtdna.prot", width=2.0, out.line=F, box=F)
# #plot_jeans_time_rect_simmap(simmap.janus.mtdna.rrnas, title = "mtdna.rrnas", width=1.5)
# #plot_jeans_time_rect_simmap(simmap.janus.nuc.aggregate, seed=1, title = "merged", width=1.5)
# plot_jeans_time_rect_simmap(simmap.janus.all.aggregate, seed=15, title = "merged", width=2.5, out.line=F, box=T)
# 
# dev.off()


#I am duplicating the protein/rrna dataset because the shift configuration is identical and the rRNA dataset is missing a few taxa
map_list<-c(simmap.janus.nuc.alldata, simmap.janus.nuc.exons, simmap.janus.nuc.introns, simmap.janus.nuc.utrs, simmap.janus.mtdna.all, simmap.janus.mtdna.proteins, simmap.janus.mtdna.proteins, simmap.janus.nuc.aggregate, simmap.janus.all.aggregate, simmap.janus.null)
names(map_list)<- c("simmap.janus.nuc.alldata", "simmap.janus.nuc.exons", "simmap.janus.nuc.introns", "simmap.janus.nuc.utrs", "simmap.janus.mtdna.all", "simmap.janus.mtdna.proteins", "simmap.janus.mtdna.rrnas", "simmap.janus.nuc.aggregate", "simmap.janus.all.aggregate", "simmap.janus.null")


#checking how many possible nodes are candidates
sum(as.numeric(count_candidates(simmap.janus.nuc.alldata) > 3)) -1 #-1 for root
#sum(as.numeric(count_candidates(simmap.janus.alldata) > 3)[199:395]-as.numeric(count_candidates(simmap.janus.alldata) > 4)[199:395])
#plot(as.phylo(simmap.janus.nuc.alldata), cex=0.2, no.margin=T)
#nodelabels(text=c(as.numeric(count_candidates(simmap.janus.alldata) > 3)[199:395]), frame = "none", adj = c(1,1))
#nodelabels(text=as.numeric(count_candidates(simmap.janus.alldata) > 3)[199:395]-as.numeric(count_candidates(simmap.janus.alldata) > 4)[199:395], frame = "none", adj = c(1,1))
# 
# cbind(which(((count_candidates(simmap.janus.nuc.alldata) > 3))[199:395])+198, )
# 
# require(BioGeoBEARS)
# tmp<-as.data.frame(prt(t=as.phylo(simmap.janus.nuc.alldata))[c(199:395),])
# tmp<-(tmp[,c(1, 10)][which(((count_candidates(simmap.janus.nuc.alldata) > 3))[199:395]),])
# 
# min(tmp$time_bp)
# 
# 86.90747-18.30677 #68.6007 for min > 4
# 86.90747-9.468076 # 77.43939 for min > 3
# 
# nodelabels(node=tmp$node)
}

#Section 10
##############################
####### running OUwie#########
##############################
{
#testing OUwie_looper (commented out)
{
# #use the datasets that have univariate BM imputation instead of phylopars.data
# imputed<-data.frame(species= phylopars.data$species, as.data.frame(bm.fit.nocorr$anc_recon)[1:198,])
# 
# OUwie.bms.mass <- OUwie_looper(maps=map_list, modlist=rep("BMS",10), df=imputed, name="mass")
# OUwie.bms.clutch <- OUwie_looper(maps=map_list, modlist=rep("BMS",10), df=imputed, name="mean.clutch")
# OUwie.bms.gen_length <- OUwie_looper(maps=map_list, modlist=rep("BMS",10), df=imputed, name="gen_length")
# OUwie.bms.chickPC1 <- OUwie_looper(maps=map_list, modlist=rep("BMS",10), df=imputed, name="chickPC1")
# OUwie.bms.nuc_div <- OUwie_looper(maps=map_list, modlist=rep("BMS",10), df=imputed, name="nuc_div")
# OUwie.bms.nuc.encprime <- OUwie_looper(maps=map_list, modlist=rep("BMS",10), df=imputed, name="nuc.encprime")
# OUwie.bms.nuc.scuo <- OUwie_looper(maps=map_list, modlist=rep("BMS",10), df=imputed, name="nuc.scuo")
# OUwie.bms.nuc.exonGC <- OUwie_looper(maps=map_list, modlist=rep("BMS",10), df=imputed, name="nuc.exonGC")
# OUwie.bms.nuc.intronGC <- OUwie_looper(maps=map_list, modlist=rep("BMS",10), df=imputed, name="nuc.intronGC")
# OUwie.bms.nuc.utrGC <- OUwie_looper(maps=map_list, modlist=rep("BMS",10), df=imputed, name="nuc.utrGC")
# 
# #OUWIE DREDGE?
# #investigate this
# #OUwie.bms.mass.dredge<-OUwie.dredge(phy=map_list$simmap.janus.null, data = imputed[,c(1,3)], criterion="mBIC", shift.max=1, algorithm="three.point")
# 
# 
# aicw(subListExtract(L=OUwie.bms.mass, name="AIC", simplify=T)[c(2,3,4,5,6,7,10)])
# aicw(subListExtract(L=OUwie.bms.clutch, name="AIC", simplify=T)[c(2,3,4,5,6,7,10)])
# aicw(subListExtract(L=OUwie.bms.gen_length, name="AIC", simplify=T)[c(2,3,4,5,6,7,10)])
# aicw(subListExtract(L=OUwie.bms.chickPC1, name="AIC", simplify=T)[c(2,3,4,5,6,7,10)])
# aicw(subListExtract(L=OUwie.bms.nuc_div, name="AIC", simplify=T)[c(2,3,4,5,6,7,10)])
# aicw(subListExtract(L=OUwie.bms.nuc.encprime, name="AIC", simplify=T)[c(2,3,4,5,6,7,10)])
# aicw(subListExtract(L=OUwie.bms.nuc.scuo, name="AIC", simplify=T)[c(2,3,4,5,6,7,10)])
# aicw(subListExtract(L=OUwie.bms.nuc.exonGC, name="AIC", simplify=T)[c(2,3,4,5,6,7,10)])
# aicw(subListExtract(L=OUwie.bms.nuc.intronGC, name="AIC", simplify=T)[c(2,3,4,5,6,7,10)])
# aicw(subListExtract(L=OUwie.bms.nuc.utrGC, name="AIC", simplify=T)[c(2,3,4,5,6,7,10)])
# 
# 
# OUwie.oum.mass <- OUwie_looper(maps=map_list, modlist=c(rep("OUM", 9), "OU1"), df=imputed, name="mass")
# OUwie.oum.clutch <- OUwie_looper(maps=map_list, modlist=c(rep("OUM", 9), "OU1"), df=imputed, name="mean.clutch")
# OUwie.oum.gen_length <- OUwie_looper(maps=map_list, modlist=c(rep("OUM", 9), "OU1"), df=imputed, name="gen_length")
# OUwie.oum.chickPC1 <- OUwie_looper(maps=map_list, modlist=c(rep("OUM", 9), "OU1"), df=imputed, name="chickPC1")
# OUwie.oum.nuc_div <- OUwie_looper(maps=map_list, modlist=c(rep("OUM", 9), "OU1"), df=imputed, name="nuc_div")
# OUwie.oum.nuc.encprime <- OUwie_looper(maps=map_list, modlist=c(rep("OUM", 9), "OU1"), df=imputed, name="nuc.encprime")
# OUwie.oum.nuc.scuo <- OUwie_looper(maps=map_list, modlist=c(rep("OUM", 9), "OU1"), df=imputed, name="nuc.scuo")
# OUwie.oum.nuc.exonGC <- OUwie_looper(maps=map_list, modlist=c(rep("OUM", 9), "OU1"), df=imputed, name="nuc.exonGC")
# OUwie.oum.nuc.intronGC <- OUwie_looper(maps=map_list, modlist=c(rep("OUM", 9), "OU1"), df=imputed, name="nuc.intronGC")
# OUwie.oum.nuc.utrGC <- OUwie_looper(maps=map_list, modlist=c(rep("OUM", 9), "OU1"), df=imputed, name="nuc.utrGC")
}

#fitting OUM to measure theta across molecular regimes (simplified tree)
{

#cbind(rownames(megaLHT), megaLHT$mass)
#rownames(megaLHT) == simmap.janus.nuc.aggregate.simplified$tip.label
bm.fit.corr.reorder<-as.data.frame(megaLHT)
bm.fit.corr.reorder<- bm.fit.corr.reorder[simmap.janus.nuc.aggregate.simplified$tip.label,]

#set up the data frame for input to OUwie
set.seed(1) #for reproducible jitter
OUwie_data <-
  data.frame(
    Genus_species = rownames(bm.fit.corr.reorder),
    Reg = getStates(simmap.janus.nuc.aggregate.simplified, type =
                      "tips"),
    mass = log(bm.fit.corr.reorder[1:198, ][, c(3)]),
    chickPC1 = jitter(bm.fit.corr.reorder[1:198, ][, c(18)]),
    mserr = rep(0.1, 198),
    mserr_2 = rep(0.05, 198)
  )

#first, reorder the input data so that it matches the tip order in the SIMMAP tree
completed_dataset<-as.data.frame(bm.fit.corr$anc_recon[1:198,])
completed_dataset.reorder<- completed_dataset[simmap.janus.nuc.aggregate.simplified$tip.label,]
OUwie_data.extra <-
  data.frame(
    Genus_species = rownames(completed_dataset.reorder),
    Reg = getStates(simmap.janus.nuc.aggregate.simplified, type = "tips"),
    clutch = (completed_dataset.reorder[1:198, ][, c(2)]),
    gen_length = (completed_dataset.reorder[1:198, ][, c(3)]),
    survival = (completed_dataset.reorder[1:198, ][, c(4)]),
    breeding = (completed_dataset.reorder[1:198, ][, c(5)]),
    longevity = (completed_dataset.reorder[1:198, ][, c(6)]),
    latitude = (completed_dataset.reorder[1:198, ][, c(8)])
  )

tip.data.molstats.reorder<-tip.data.molstats
tip.data.molstats.reorder<-tip.data.molstats.reorder[simmap.janus.nuc.exons$tip.label,]
OUwie_data.codons <-
  data.frame(
    Genus_species = simmap.janus.nuc.exons$tip.label,
    Reg = getStates(simmap.janus.nuc.exons, type = "tips"),
    log.enc = log(tip.data.molstats.reorder$nuc.enc),
    log.encprime = log(tip.data.molstats.reorder$nuc.encprime),
    log.scuo = log(tip.data.molstats.reorder$nuc.scuo)
  )

#fit the shifting theta model to body mass
# mass.OUM <- OUwie(
#   phy = simmap.janus.nuc.aggregate.simplified,
#   data = OUwie_data[, c(1, 2, 3, 5)],
#   model = 'OUM',
#   simmap.tree = T,
#   diagn = T,
#   check.identify = T,
#   algorithm = "invert",
#   scaleHeight = F,
#   mserr = "known",
#   get.root.theta = F
# )
# saveRDS(mass.OUM, file="./RDS/mass.OUM.RDS")
mass.OUM<- readRDS(file="./RDS/mass.OUM.RDS")

#try with parametric bootstrapping
# OUwie_boots.mass <-
#   OUwie.boot(
#     phy = simmap.janus.nuc.aggregate.simplified,
#     data = OUwie_data[, c(1, 2, 3, 5)],
#     model = "OUM",
#     nboot = 100,
#     alpha = mass.OUM$solution[1, ],
#     sigma.sq = mass.OUM$solution[2, ],
#     theta = mass.OUM$theta[, 1],
#     theta0 = mass.OUM$theta[1, 1],
#     algorithm = "invert",
#     simmap.tree = T,
#     quiet = F,
#     mserr = "none"
#   )
# saveRDS(OUwie_boots.mass, file="./RDS/OUwie_boots.mass.RDS")
OUwie_boots.mass<-readRDS(file="./RDS/OUwie_boots.mass.RDS")

OUwie_boots.mass <- OUwie_boots.mass[,c(21:30)] #subset the thetas
OUwie_boots.mass <- OUwie_boots.mass %>% as.data.frame(.) %>%
  tidyr::pivot_longer(cols = everything()) %>% 
  mutate_at("name", stringr::str_replace, "theta_", "") %>%
  rename(Reg = name) %>% 
  mutate(name = rep(sort(rep(c(1:10), 10)), 10)) %>%
  mutate(type = rep("theta", length(.[,1]))) %>%
  select(Reg, name, value, type)

tmp<-readRDS(file="./RDS/OUwie_boots.mass.RDS")
HDInterval::hdi(log(2)/tmp[,1]) #get the HPD 1/2 life
median(log(2)/tmp[,1])

#fit the shifting theta model to chickPC1
# chickPC1.OUM <- OUwie(
#   phy = simmap.janus.nuc.aggregate.simplified,
#   data = OUwie_data[, c(1, 2, 4)],
#   model = 'OUM',
#   simmap.tree = T,
#   diagn = T,
#   check.identify = T,
#   algorithm = "invert",
#   scaleHeight = F,
#   mserr = "none"
# )
# saveRDS(chickPC1.OUM, file="./RDS/chickPC1.OUM.RDS")
chickPC1.OUM<- readRDS(file="./RDS/chickPC1.OUM.RDS")

#try with parametric bootstrapping
# OUwie_boots.chickPC1 <-
#   OUwie.boot(
#     phy = simmap.janus.nuc.aggregate.simplified,
#     data = OUwie_data[, c(1, 2, 4, 5)],
#     model = "OUM",
#     nboot = 100,
#     alpha = chickPC1.OUM$solution[1, ],
#     sigma.sq = chickPC1.OUM$solution[2, ],
#     theta = chickPC1.OUM$theta[, 1],
#     theta0 = chickPC1.OUM$theta[1, 1],
#     algorithm = "invert",
#     simmap.tree = T,
#     quiet = F,
#     mserr = "none"
#   )
# saveRDS(OUwie_boots.chickPC1, file="./RDS/OUwie_boots.chickPC1.RDS")
OUwie_boots.chickPC1<- readRDS(file="./RDS/OUwie_boots.chickPC1.RDS")

OUwie_boots.chickPC1 <- OUwie_boots.chickPC1[,c(21:30)] #subset the thetas
OUwie_boots.chickPC1 <- OUwie_boots.chickPC1 %>% as.data.frame(.) %>%
  tidyr::pivot_longer(cols = everything()) %>% 
  mutate_at("name", stringr::str_replace, "theta_", "") %>%
  rename(Reg = name) %>% 
  mutate(name = rep(sort(rep(c(1:10), 10)), 10)) %>%
  mutate(type = rep("theta", length(.[,1]))) %>%
  select(Reg, name, value, type)

tmp<-readRDS(file="./RDS/OUwie_boots.chickPC1.RDS")
HDInterval::hdi(log(2)/tmp[,1]) #get the HPD 1/2 life
median(log(2)/tmp[,1])


###
#fit the shifting theta model to clutch size
# clutch.OUM <- OUwie(
#   phy = simmap.janus.nuc.aggregate.simplified,
#   data = OUwie_data.extra[, c(1, 2, 3)],
#   model = 'OUM',
#   simmap.tree = T,
#   diagn = T,
#   check.identify = T,
#   algorithm = "invert",
#   scaleHeight = F,
#   mserr = "none"
# )
# saveRDS(clutch.OUM, file="./RDS/clutch.OUM.RDS")
clutch.OUM<- readRDS(file="./RDS/chickPC1.OUM.RDS")

#try with parametric bootstrapping
# OUwie_boots.clutch <-
#   OUwie.boot(
#     phy = simmap.janus.nuc.aggregate.simplified,
#     data = OUwie_data.extra[, c(1, 2, 3)],
#     model = "OUM",
#     nboot = 100,
#     alpha = clutch.OUM$solution[1, ],
#     sigma.sq = clutch.OUM$solution[2, ],
#     theta = clutch.OUM$theta[, 1],
#     theta0 = clutch.OUM$theta[1, 1],
#     algorithm = "invert",
#     simmap.tree = T,
#     quiet = F,
#     mserr = "none"
#   )
# saveRDS(OUwie_boots.clutch, file="./RDS/OUwie_boots.clutch.RDS")
OUwie_boots.clutch <- readRDS(file="./RDS/OUwie_boots.clutch.RDS")

OUwie_boots.clutch <- OUwie_boots.clutch[,c(21:30)] #subset the thetas
OUwie_boots.clutch <- OUwie_boots.clutch %>% as.data.frame(.) %>%
  tidyr::pivot_longer(cols = everything()) %>%
  mutate_at("name", stringr::str_replace, "theta_", "") %>%
  rename(Reg = name) %>%
  mutate(name = rep(sort(rep(c(1:10), 10)), 10)) %>%
  mutate(type = rep("theta", length(.[,1]))) %>%
  select(Reg, name, value, type)

tmp<-readRDS(file="./RDS/OUwie_boots.clutch.RDS")
HDInterval::hdi(log(2)/tmp[,1]) #get the HPD 1/2 life
median(log(2)/tmp[,1])


###
#fit the shifting theta model to gen_length size
# gen_length.OUM <- OUwie(
#   phy = simmap.janus.nuc.aggregate.simplified,
#   data = OUwie_data.extra[, c(1, 2, 4)],
#   model = 'OUM',
#   simmap.tree = T,
#   diagn = T,
#   check.identify = T,
#   algorithm = "invert",
#   scaleHeight = F,
#   mserr = "none"#,
#   #opts = list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="500", "ftol_abs"=0.001)
# )
# saveRDS(gen_length.OUM, file="./RDS/gen_length.OUM.RDS")
gen_length.OUM <- readRDS(file="./RDS/gen_length.OUM.RDS")
### MAY BE AT SADDLE POINT ### -- but looks OK

#try with parametric bootstrapping
# OUwie_boots.gen_length <-
#   OUwie.boot(
#     phy = simmap.janus.nuc.aggregate.simplified,
#     data = OUwie_data.extra[, c(1, 2, 4)],
#     model = "OUM",
#     nboot = 100,
#     alpha = gen_length.OUM$solution[1, ],
#     sigma.sq = gen_length.OUM$solution[2, ],
#     theta = gen_length.OUM$theta[, 1],
#     theta0 = gen_length.OUM$theta[1, 1],
#     algorithm = "invert",
#     simmap.tree = T,
#     quiet = F,
#     mserr = "none"
#   )
# saveRDS(OUwie_boots.gen_length, file="./RDS/OUwie_boots.gen_length.RDS")
OUwie_boots.gen_length <- readRDS(file="./RDS/OUwie_boots.gen_length.RDS")

OUwie_boots.gen_length <- OUwie_boots.gen_length[,c(21:30)] #subset the thetas
OUwie_boots.gen_length <- OUwie_boots.gen_length %>% as.data.frame(.) %>%
  tidyr::pivot_longer(cols = everything()) %>%
  mutate_at("name", stringr::str_replace, "theta_", "") %>%
  rename(Reg = name) %>%
  mutate(name = rep(sort(rep(c(1:10), 10)), 10)) %>%
  mutate(type = rep("theta", length(.[,1]))) %>%
  select(Reg, name, value, type)

tmp<-readRDS(file="./RDS/OUwie_boots.gen_length.RDS")
HDInterval::hdi(log(2)/tmp[,1]) #get the HPD 1/2 life
median(log(2)/tmp[,1])


###
#fit the shifting theta model to survival
# survival.OUM <- OUwie(
#   phy = simmap.janus.nuc.aggregate.simplified,
#   data = OUwie_data.extra[, c(1, 2, 5)],
#   model = 'OUM',
#   simmap.tree = T,
#   diagn = T,
#   check.identify = T,
#   algorithm = "invert",
#   scaleHeight = F,
#   mserr = "none"
# )
# saveRDS(survival.OUM, file="./RDS/survival.OUM.RDS")
survival.OUM<- readRDS(file="./RDS/survival.OUM.RDS")

#try with parametric bootstrapping
# OUwie_boots.survival <-
#   OUwie.boot(
#     phy = simmap.janus.nuc.aggregate.simplified,
#     data = OUwie_data.extra[, c(1, 2, 5)],
#     model = "OUM",
#     nboot = 100,
#     alpha = survival.OUM$solution[1, ],
#     sigma.sq = survival.OUM$solution[2, ],
#     theta = survival.OUM$theta[, 1],
#     theta0 = survival.OUM$theta[1, 1],
#     algorithm = "invert",
#     simmap.tree = T,
#     quiet = F,
#     mserr = "none"
#   )
# saveRDS(OUwie_boots.survival, file="./RDS/OUwie_boots.survival.RDS")
OUwie_boots.survival <- readRDS(file="./RDS/OUwie_boots.survival.RDS")

OUwie_boots.survival <- OUwie_boots.survival[,c(21:30)] #subset the thetas
OUwie_boots.survival <- OUwie_boots.survival %>% as.data.frame(.) %>%
  tidyr::pivot_longer(cols = everything()) %>%
  mutate_at("name", stringr::str_replace, "theta_", "") %>%
  rename(Reg = name) %>%
  mutate(name = rep(sort(rep(c(1:10), 10)), 10)) %>%
  mutate(type = rep("theta", length(.[,1]))) %>%
  select(Reg, name, value, type)

tmp<-readRDS(file="./RDS/OUwie_boots.survival.RDS")
HDInterval::hdi(log(2)/tmp[,1]) #get the HPD 1/2 life
median(log(2)/tmp[,1])


###
#fit the shifting theta model to breeding
# breeding.OUM <- OUwie(
#   phy = simmap.janus.nuc.aggregate.simplified,
#   data = OUwie_data.extra[, c(1, 2, 6)],
#   model = 'OUM',
#   simmap.tree = T,
#   diagn = T,
#   check.identify = T,
#   algorithm = "invert",
#   scaleHeight = F,
#   mserr = "none"
# )
# saveRDS(breeding.OUM, file="./RDS/breeding.OUM.RDS")
breeding.OUM<- readRDS(file="./RDS/breeding.OUM.RDS")

#try with parametric bootstrapping
# OUwie_boots.breeding <-
#   OUwie.boot(
#     phy = simmap.janus.nuc.aggregate.simplified,
#     data = OUwie_data.extra[, c(1, 2, 6)],
#     model = "OUM",
#     nboot = 100,
#     alpha = breeding.OUM$solution[1, ],
#     sigma.sq = breeding.OUM$solution[2, ],
#     theta = breeding.OUM$theta[, 1],
#     theta0 = breeding.OUM$theta[1, 1],
#     algorithm = "invert",
#     simmap.tree = T,
#     quiet = F,
#     mserr = "none"
#   )
# saveRDS(OUwie_boots.breeding, file="./RDS/OUwie_boots.breeding.RDS")
OUwie_boots.breeding <- readRDS(file="./RDS/OUwie_boots.breeding.RDS")

OUwie_boots.breeding <- OUwie_boots.breeding[,c(21:30)] #subset the thetas
OUwie_boots.breeding <- OUwie_boots.breeding %>% as.data.frame(.) %>%
  tidyr::pivot_longer(cols = everything()) %>%
  mutate_at("name", stringr::str_replace, "theta_", "") %>%
  rename(Reg = name) %>%
  mutate(name = rep(sort(rep(c(1:10), 10)), 10)) %>%
  mutate(type = rep("theta", length(.[,1]))) %>%
  select(Reg, name, value, type)

tmp<-readRDS(file="./RDS/OUwie_boots.breeding.RDS")
HDInterval::hdi(log(2)/tmp[,1]) #get the HPD 1/2 life
median(log(2)/tmp[,1])


###
#fit the shifting theta model to longevity
# longevity.OUM <- OUwie(
#   phy = simmap.janus.nuc.aggregate.simplified,
#   data = OUwie_data.extra[, c(1, 2, 7)],
#   model = 'OUM',
#   simmap.tree = T,
#   diagn = T,
#   check.identify = T,
#   algorithm = "invert",
#   scaleHeight = F,
#   mserr = "none"
# )
# saveRDS(longevity.OUM, file="./RDS/longevity.OUM.RDS")
longevity.OUM<- readRDS(file="./RDS/longevity.OUM.RDS")

#try with parametric bootstrapping
# OUwie_boots.longevity <-
#   OUwie.boot(
#     phy = simmap.janus.nuc.aggregate.simplified,
#     data = OUwie_data.extra[, c(1, 2, 7)],
#     model = "OUM",
#     nboot = 100,
#     alpha = longevity.OUM$solution[1, ],
#     sigma.sq = longevity.OUM$solution[2, ],
#     theta = longevity.OUM$theta[, 1],
#     theta0 = longevity.OUM$theta[1, 1],
#     algorithm = "invert",
#     simmap.tree = T,
#     quiet = F,
#     mserr = "none"
#   )
# saveRDS(OUwie_boots.longevity, file="./RDS/OUwie_boots.longevity.RDS")
OUwie_boots.longevity <- readRDS(file="./RDS/OUwie_boots.longevity.RDS")

OUwie_boots.longevity <- OUwie_boots.longevity[,c(21:30)] #subset the thetas
OUwie_boots.longevity <- OUwie_boots.longevity %>% as.data.frame(.) %>%
  tidyr::pivot_longer(cols = everything()) %>%
  mutate_at("name", stringr::str_replace, "theta_", "") %>%
  rename(Reg = name) %>%
  mutate(name = rep(sort(rep(c(1:10), 10)), 10)) %>%
  mutate(type = rep("theta", length(.[,1]))) %>%
  select(Reg, name, value, type)

tmp<-readRDS(file="./RDS/OUwie_boots.longevity.RDS")
HDInterval::hdi(log(2)/tmp[,1]) #get the HPD 1/2 life
median(log(2)/tmp[,1])

###
#fit the shifting theta model to latitude
# latitude.OUM <- OUwie(
#   phy = simmap.janus.nuc.aggregate.simplified,
#   data = OUwie_data.extra[, c(1, 2, 8)],
#   model = 'OUM',
#   simmap.tree = T,
#   diagn = T,
#   check.identify = T,
#   algorithm = "invert",
#   scaleHeight = F,
#   mserr = "none"
# )
# saveRDS(latitude.OUM, file="./RDS/latitude.OUM.RDS")
latitude.OUM<- readRDS(file="./RDS/latitude.OUM.RDS")

#try with parametric bootstrapping
# OUwie_boots.latitude <-
#   OUwie.boot(
#     phy = simmap.janus.nuc.aggregate.simplified,
#     data = OUwie_data.extra[, c(1, 2, 8)],
#     model = "OUM",
#     nboot = 100,
#     alpha = latitude.OUM$solution[1, ],
#     sigma.sq = latitude.OUM$solution[2, ],
#     theta = latitude.OUM$theta[, 1],
#     theta0 = latitude.OUM$theta[1, 1],
#     algorithm = "invert",
#     simmap.tree = T,
#     quiet = F,
#     mserr = "none"
#   )
# saveRDS(OUwie_boots.latitude, file="./RDS/OUwie_boots.latitude.RDS")
OUwie_boots.latitude <- readRDS(file="./RDS/OUwie_boots.latitude.RDS")

OUwie_boots.latitude <- OUwie_boots.latitude[,c(21:30)] #subset the thetas
OUwie_boots.latitude <- OUwie_boots.latitude %>% as.data.frame(.) %>%
  tidyr::pivot_longer(cols = everything()) %>%
  mutate_at("name", stringr::str_replace, "theta_", "") %>%
  rename(Reg = name) %>%
  mutate(name = rep(sort(rep(c(1:10), 10)), 10)) %>%
  mutate(type = rep("theta", length(.[,1]))) %>%
  select(Reg, name, value, type)

tmp<-readRDS(file="./RDS/OUwie_boots.latitude.RDS")
HDInterval::hdi(log(2)/tmp[,1]) #get the HPD 1/2 life
median(log(2)/tmp[,1])

}
#fitting OUM to measure theta across molecular regimes (simplified tree)
{
  
  #cbind(rownames(megaLHT), megaLHT$mass)
  #rownames(megaLHT) == simmap.janus.nuc.aggregate.simplified$tip.label
  bm.fit.corr.reorder2<-as.data.frame(megaLHT)
  bm.fit.corr.reorder2<- bm.fit.corr.reorder[simmap.janus.nuc.aggregate.simplified2$tip.label,]
  
  #set up the data frame for input to OUwie
  OUwie_data2 <-
    data.frame(
      Genus_species = rownames(bm.fit.corr.reorder2),
      Reg = getStates(simmap.janus.nuc.aggregate.simplified2, type =
                        "tips"),
      mass = log(bm.fit.corr.reorder2[1:198, ][, c(3)]), #1
      chickPC1 = jitter(bm.fit.corr.reorder2[1:198, ][, c(18)]), #7
      mserr = rep(0.1, 198),
      mserr_2 = rep(0.05, 198)
    )
  

  #first, reorder the input data so that it matches the tip order in the SIMMAP tree
  completed_dataset2<-as.data.frame(bm.fit.corr$anc_recon[1:198,])
  completed_dataset.reorder2<- completed_dataset2[simmap.janus.nuc.aggregate.simplified2$tip.label,]
  OUwie_data.extra2 <-
    data.frame(
      Genus_species = rownames(completed_dataset.reorder2),
      Reg = getStates(simmap.janus.nuc.aggregate.simplified2, type = "tips"),
      clutch = (completed_dataset.reorder2[1:198, ][, c(2)]),
      gen_length = (completed_dataset.reorder2[1:198, ][, c(3)]),
      survival = (completed_dataset.reorder2[1:198, ][, c(4)]),
      breeding = (completed_dataset.reorder2[1:198, ][, c(5)]),
      longevity = (completed_dataset.reorder2[1:198, ][, c(6)]),
      latitude = (completed_dataset.reorder2[1:198, ][, c(8)])
    )
  
#commented here because it should already exist from the prior runs
  # tip.data.molstats.reorder<-tip.data.molstats
  # tip.data.molstats.reorder<-tip.data.molstats.reorder[simmap.janus.nuc.exons$tip.label,]
  # OUwie_data.codons <-
  #   data.frame(
  #     Genus_species = simmap.janus.nuc.exons$tip.label,
  #     Reg = getStates(simmap.janus.nuc.exons, type = "tips"),
  #     log.enc = log(tip.data.molstats.reorder$nuc.enc),
  #     log.encprime = log(tip.data.molstats.reorder$nuc.encprime),
  #     log.scuo = log(tip.data.molstats.reorder$nuc.scuo)
  #   )
  
  
  #fit the shifting theta model to body mass
  # mass2.OUM <- OUwie(
  #   phy = simmap.janus.nuc.aggregate.simplified2,
  #   data = OUwie_data2[, c(1, 2, 3, 5)],
  #   model = 'OUM',
  #   simmap.tree = T,
  #   diagn = T,
  #   check.identify = T,
  #   algorithm = "invert",
  #   scaleHeight = F,
  #   mserr = "known",
  #   get.root.theta = F
  # )
  # saveRDS(mass2.OUM, file="./RDS/mass2.OUM.RDS")
  mass2.OUM<- readRDS(file="./RDS/mass2.OUM.RDS")
  
  #try with parametric bootstrapping
  # OUwie_boots2.mass <-
  #   OUwie.boot(
  #     phy = simmap.janus.nuc.aggregate.simplified2,
  #     data = OUwie_data2[, c(1, 2, 3, 5)],
  #     model = "OUM",
  #     nboot = 100,
  #     alpha = mass2.OUM$solution[1, ],
  #     sigma.sq = mass2.OUM$solution[2, ],
  #     theta = mass2.OUM$theta[, 1],
  #     theta0 = mass2.OUM$theta[1, 1],
  #     algorithm = "invert",
  #     simmap.tree = T,
  #     quiet = F,
  #     mserr = "none"
  #   )
  # saveRDS(OUwie_boots2.mass, file="./RDS/OUwie_boots2.mass.RDS")
  OUwie_boots2.mass <-readRDS(file="./RDS/OUwie_boots.mass2.RDS")

  OUwie_boots2.mass <- OUwie_boots2.mass[,c(23:33)] #subset the thetas (for 2)
  OUwie_boots2.mass <- OUwie_boots2.mass %>% as.data.frame(.) %>%
    tidyr::pivot_longer(cols = everything()) %>% 
    mutate_at("name", stringr::str_replace, "theta_", "") %>%
    rename(Reg = name) %>% 
    mutate(name = rep(sort(rep(c(1:11), 10)), 10)) %>% 
    mutate(type = rep("theta", length(.[,1]))) %>%
    select(Reg, name, value, type)
  
  tmp<-readRDS(file="./RDS/OUwie_boots2.mass.RDS")
  HDInterval::hdi(log(2)/tmp[,1]) #get the HPD 1/2 life
  median(log(2)/tmp[,1])
  
  
  #fit the shifting theta model to chickPC1
  # chickPC1_2.OUM <- OUwie(
  #   phy = simmap.janus.nuc.aggregate.simplified2,
  #   data = OUwie_data2[, c(1, 2, 4)],
  #   model = 'OUM',
  #   simmap.tree = T,
  #   diagn = T,
  #   check.identify = T,
  #   algorithm = "invert",
  #   scaleHeight = F,
  #   mserr = "none"
  # )
  # saveRDS(chickPC1_2.OUM, file="./RDS/chickPC1_2.OUM.RDS")
  chickPC1_2.OUM<- readRDS(file="./RDS/chickPC1_2.OUM.RDS")
  
  #try with parametric bootstrapping
  # OUwie_boots2.chickPC1 <-
  #   OUwie.boot(
  #     phy = simmap.janus.nuc.aggregate.simplified2,
  #     data = OUwie_data2[, c(1, 2, 4, 5)],
  #     model = "OUM",
  #     nboot = 100,
  #     alpha = chickPC1_2.OUM$solution[1, ],
  #     sigma.sq = chickPC1_2.OUM$solution[2, ],
  #     theta = chickPC1_2.OUM$theta[, 1],
  #     theta0 = chickPC1_2.OUM$theta[1, 1],
  #     algorithm = "invert",
  #     simmap.tree = T,
  #     quiet = F,
  #     mserr = "none"
  #   )
  # saveRDS(OUwie_boots2.chickPC1, file="./RDS/OUwie_boots2.chickPC1.RDS")
  OUwie_boots2.chickPC1<- readRDS(file="./RDS/OUwie_boots2.chickPC1.RDS")
  
  OUwie_boots2.chickPC1 <- OUwie_boots2.chickPC1[,c(23:33)] #subset the thetas
  OUwie_boots2.chickPC1 <- OUwie_boots2.chickPC1 %>% as.data.frame(.) %>%
    tidyr::pivot_longer(cols = everything()) %>% 
    mutate_at("name", stringr::str_replace, "theta_", "") %>%
    rename(Reg = name) %>% 
    mutate(name = rep(sort(rep(c(1:11), 10)), 10)) %>%
    mutate(type = rep("theta", length(.[,1]))) %>%
    select(Reg, name, value, type)
  
  tmp<-readRDS(file="./RDS/OUwie_boots2.chickPC1.RDS")
  HDInterval::hdi(log(2)/tmp[,1]) #get the HPD 1/2 life
  median(log(2)/tmp[,1])
  
  ###
  #fit the shifting theta model to clutch size
  # clutch2.OUM <- OUwie(
  #   phy = simmap.janus.nuc.aggregate.simplified2,
  #   data = OUwie_data.extra2[, c(1, 2, 3)],
  #   model = 'OUM',
  #   simmap.tree = T,
  #   diagn = T,
  #   check.identify = T,
  #   algorithm = "invert",
  #   scaleHeight = F,
  #   mserr = "none"
  # )
  # saveRDS(clutch2.OUM, file="./RDS/clutch2.OUM.RDS")
  clutch2.OUM<- readRDS(file="./RDS/clutch2.OUM.RDS")
  
  #try with parametric bootstrapping
  # OUwie_boots2.clutch <-
  #   OUwie.boot(
  #     phy = simmap.janus.nuc.aggregate.simplified2,
  #     data = OUwie_data.extra2[, c(1, 2, 3)],
  #     model = "OUM",
  #     nboot = 100,
  #     alpha = clutch2.OUM$solution[1, ],
  #     sigma.sq = clutch2.OUM$solution[2, ],
  #     theta = clutch2.OUM$theta[, 1],
  #     theta0 = clutch2.OUM$theta[1, 1],
  #     algorithm = "invert",
  #     simmap.tree = T,
  #     quiet = F,
  #     mserr = "none"
  #   )
  # saveRDS(OUwie_boots2.clutch, file="./RDS/OUwie_boots2.clutch.RDS")
  OUwie_boots2.clutch <- readRDS(file="./RDS/OUwie_boots2.clutch.RDS")
  
  OUwie_boots2.clutch <- readRDS(file="./RDS/OUwie_boots2.clutch.RDS")
  OUwie_boots2.clutch <- OUwie_boots2.clutch[,c(23:33)] #subset the thetas
  OUwie_boots2.clutch <- OUwie_boots2.clutch %>% as.data.frame(.) %>%
    tidyr::pivot_longer(cols = everything()) %>%
    mutate_at("name", stringr::str_replace, "theta_", "") %>%
    rename(Reg = name) %>%
    mutate(name = rep(sort(rep(c(1:11), 10)), 10)) %>%
    mutate(type = rep("theta", length(.[,1]))) %>%
    select(Reg, name, value, type)
  
  tmp<-readRDS(file="./RDS/OUwie_boots2.clutch.RDS")
  HDInterval::hdi(log(2)/tmp[,1]) #get the HPD 1/2 life
  median(log(2)/tmp[,1])
  
  ###
  #fit the shifting theta model to gen_length
  # gen_length2.OUM <- OUwie(
  #   phy = simmap.janus.nuc.aggregate.simplified2,
  #   data = OUwie_data.extra2[, c(1, 2, 4)],
  #   model = 'OUM',
  #   simmap.tree = T,
  #   diagn = T,
  #   check.identify = T,
  #   algorithm = "invert",
  #   scaleHeight = F,
  #   mserr = "none"
  # )
  # saveRDS(gen_length2.OUM, file="./RDS/gen_length2.OUM.RDS")
  gen_length2.OUM<- readRDS(file="./RDS/gen_length2.OUM.RDS")
  #### MAY BE AT SADDLEPOINT ### TRYING ANYWAY THOUGH
  
  #try with parametric bootstrapping
  # OUwie_boots2.gen_length <-
  #   OUwie.boot(
  #     phy = simmap.janus.nuc.aggregate.simplified2,
  #     data = OUwie_data.extra2[, c(1, 2, 4)],
  #     model = "OUM",
  #     nboot = 100,
  #     alpha = gen_length2.OUM$solution[1, ],
  #     sigma.sq = gen_length2.OUM$solution[2, ],
  #     theta = gen_length2.OUM$theta[, 1],
  #     theta0 = gen_length2.OUM$theta[1, 1],
  #     algorithm = "invert",
  #     simmap.tree = T,
  #     quiet = F,
  #     mserr = "none"
  #   )
  # saveRDS(OUwie_boots2.gen_length, file="./RDS/OUwie_boots2.gen_length.RDS")
  OUwie_boots2.gen_length <- readRDS(file="./RDS/OUwie_boots2.gen_length.RDS")
 
  OUwie_boots2.gen_length <- OUwie_boots2.gen_length[,c(23:33)] #subset the thetas
  OUwie_boots2.gen_length <- OUwie_boots2.gen_length %>% as.data.frame(.) %>%
    tidyr::pivot_longer(cols = everything()) %>%
    mutate_at("name", stringr::str_replace, "theta_", "") %>%
    rename(Reg = name) %>%
    mutate(name = rep(sort(rep(c(1:11), 10)), 10)) %>%
    mutate(type = rep("theta", length(.[,1]))) %>%
    select(Reg, name, value, type)
  
  tmp<-readRDS(file="./RDS/OUwie_boots2.gen_length.RDS")
  HDInterval::hdi(log(2)/tmp[,1]) #get the HPD 1/2 life
  median(log(2)/tmp[,1])
  
  ###
  #fit the shifting theta model to survival
  # survival2.OUM <- OUwie(
  #   phy = simmap.janus.nuc.aggregate.simplified2,
  #   data = OUwie_data.extra[, c(1, 2, 5)],
  #   model = 'OUM',
  #   simmap.tree = T,
  #   diagn = T,
  #   check.identify = T,
  #   algorithm = "invert",
  #   scaleHeight = F,
  #   mserr = "none"
  # )
  # saveRDS(survival2.OUM, file="./RDS/survival2.OUM.RDS")
  survival2.OUM<- readRDS(file="./RDS/survival2.OUM.RDS")
  
  #try with parametric bootstrapping
  # OUwie_boots2.survival <-
  #   OUwie.boot(
  #     phy = simmap.janus.nuc.aggregate.simplified2,
  #     data = OUwie_data.extra2[, c(1, 2, 5)],
  #     model = "OUM",
  #     nboot = 100,
  #     alpha = survival2.OUM$solution[1, ],
  #     sigma.sq = survival2.OUM$solution[2, ],
  #     theta = survival2.OUM$theta[, 1],
  #     theta0 = survival2.OUM$theta[1, 1],
  #     algorithm = "invert",
  #     simmap.tree = T,
  #     quiet = F,
  #     mserr = "none"
  #   )
  # saveRDS(OUwie_boots2.survival, file="./RDS/OUwie_boots2.survival.RDS")
  OUwie_boots2.survival <- readRDS(file="./RDS/OUwie_boots2.survival.RDS")
  
  OUwie_boots2.survival <- OUwie_boots2.survival[,c(23:33)] #subset the thetas
  OUwie_boots2.survival <- OUwie_boots2.survival %>% as.data.frame(.) %>%
    tidyr::pivot_longer(cols = everything()) %>%
    mutate_at("name", stringr::str_replace, "theta_", "") %>%
    rename(Reg = name) %>%
    mutate(name = rep(sort(rep(c(1:11), 10)), 10)) %>%
    mutate(type = rep("theta", length(.[,1]))) %>%
    select(Reg, name, value, type)
  
  tmp<-readRDS(file="./RDS/OUwie_boots2.survival.RDS")
  HDInterval::hdi(log(2)/tmp[,1]) #get the HPD 1/2 life
  median(log(2)/tmp[,1])
  
  ###
  #fit the shifting theta model to breeding
  # breeding2.OUM <- OUwie(
  #   phy = simmap.janus.nuc.aggregate.simplified2,
  #   data = OUwie_data.extra[, c(1, 2, 6)],
  #   model = 'OUM',
  #   simmap.tree = T,
  #   diagn = T,
  #   check.identify = T,
  #   algorithm = "invert",
  #   scaleHeight = F,
  #   mserr = "none"
  # )
  # saveRDS(breeding2.OUM, file="./RDS/breeding2.OUM.RDS")
  breeding2.OUM <- readRDS(file="./RDS/breeding2.OUM.RDS")
  
  #try with parametric bootstrapping
  # OUwie_boots2.breeding <-
  #   OUwie.boot(
  #     phy = simmap.janus.nuc.aggregate.simplified2,
  #     data = OUwie_data.extra2[, c(1, 2, 6)],
  #     model = "OUM",
  #     nboot = 100,
  #     alpha = breeding2.OUM$solution[1, ],
  #     sigma.sq = breeding2.OUM$solution[2, ],
  #     theta = breeding2.OUM$theta[, 1],
  #     theta0 = breeding2.OUM$theta[1, 1],
  #     algorithm = "invert",
  #     simmap.tree = T,
  #     quiet = F,
  #     mserr = "none"
  #   )
  # saveRDS(OUwie_boots2.breeding, file="./RDS/OUwie_boots2.breeding.RDS")
  OUwie_boots2.breeding <- readRDS(file="./RDS/OUwie_boots2.breeding.RDS")
  
  OUwie_boots2.breeding <- OUwie_boots2.breeding[,c(23:33)] #subset the thetas
  OUwie_boots2.breeding <- OUwie_boots2.breeding %>% as.data.frame(.) %>%
    tidyr::pivot_longer(cols = everything()) %>%
    mutate_at("name", stringr::str_replace, "theta_", "") %>%
    rename(Reg = name) %>%
    mutate(name = rep(sort(rep(c(1:11), 10)), 10)) %>%
    mutate(type = rep("theta", length(.[,1]))) %>%
    select(Reg, name, value, type)
  
  tmp<-readRDS(file="./RDS/OUwie_boots2.breeding.RDS")
  HDInterval::hdi(log(2)/tmp[,1]) #get the HPD 1/2 life
  median(log(2)/tmp[,1])
  
  
  ###
  #fit the shifting theta model to longevity
  # longevity2.OUM <- OUwie(
  #   phy = simmap.janus.nuc.aggregate.simplified2,
  #   data = OUwie_data.extra2[, c(1, 2, 7)],
  #   model = 'OUM',
  #   simmap.tree = T,
  #   diagn = T,
  #   check.identify = T,
  #   algorithm = "invert",
  #   scaleHeight = F,
  #   mserr = "none"
  # )
  # saveRDS(longevity2.OUM, file="./RDS/longevity2.OUM.RDS")
  longevity2.OUM<- readRDS(file="./RDS/longevity2.OUM.RDS")
  
  #try with parametric bootstrapping
  # OUwie_boots2.longevity <-
  #   OUwie.boot(
  #     phy = simmap.janus.nuc.aggregate.simplified2,
  #     data = OUwie_data.extra2[, c(1, 2, 7)],
  #     model = "OUM",
  #     nboot = 100,
  #     alpha = longevity2.OUM$solution[1, ],
  #     sigma.sq = longevity2.OUM$solution[2, ],
  #     theta = longevity2.OUM$theta[, 1],
  #     theta0 = longevity2.OUM$theta[1, 1],
  #     algorithm = "invert",
  #     simmap.tree = T,
  #     quiet = F,
  #     mserr = "none"
  #   )
  # saveRDS(OUwie_boots2.longevity, file="./RDS/OUwie_boots2.longevity.RDS")
  OUwie_boots2.longevity <- readRDS(file="./RDS/OUwie_boots2.longevity.RDS")
  
  OUwie_boots2.longevity <- OUwie_boots2.longevity[,c(23:33)] #subset the thetas
  OUwie_boots2.longevity <- OUwie_boots2.longevity %>% as.data.frame(.) %>%
    tidyr::pivot_longer(cols = everything()) %>%
    mutate_at("name", stringr::str_replace, "theta_", "") %>%
    rename(Reg = name) %>%
    mutate(name = rep(sort(rep(c(1:11), 10)), 10)) %>%
    mutate(type = rep("theta", length(.[,1]))) %>%
    select(Reg, name, value, type)
  
  tmp<-readRDS(file="./RDS/OUwie_boots2.longevity.RDS")
  HDInterval::hdi(log(2)/tmp[,1]) #get the HPD 1/2 life
  median(log(2)/tmp[,1])
  
  
  ###
  #fit the shifting theta model to latitude
  # latitude2.OUM <- OUwie(
  #   phy = simmap.janus.nuc.aggregate.simplified2,
  #   data = OUwie_data.extra2[, c(1, 2, 8)],
  #   model = 'OUM',
  #   simmap.tree = T,
  #   diagn = T,
  #   check.identify = T,
  #   algorithm = "invert",
  #   scaleHeight = F,
  #   mserr = "none"
  # )
  # saveRDS(latitude2.OUM, file="./RDS/latitude2.OUM.RDS")
  latitude2.OUM <- readRDS(file="./RDS/latitude2.OUM.RDS")
  
  #try with parametric bootstrapping
  # OUwie_boots2.latitude <-
  #   OUwie.boot(
  #     phy = simmap.janus.nuc.aggregate.simplified2,
  #     data = OUwie_data.extra2[, c(1, 2, 8)],
  #     model = "OUM",
  #     nboot = 100,
  #     alpha = latitude2.OUM$solution[1, ],
  #     sigma.sq = latitude2.OUM$solution[2, ],
  #     theta = latitude2.OUM$theta[, 1],
  #     theta0 = latitude2.OUM$theta[1, 1],
  #     algorithm = "invert",
  #     simmap.tree = T,
  #     quiet = F,
  #     mserr = "none"
  #   )
  # saveRDS(OUwie_boots2.latitude, file="./RDS/OUwie_boots2.latitude.RDS")
  OUwie_boots2.latitude <- readRDS(file="./RDS/OUwie_boots2.latitude.RDS")
  
  OUwie_boots2.latitude <- OUwie_boots2.latitude[,c(23:33)] #subset the thetas
  OUwie_boots2.latitude <- OUwie_boots2.latitude %>% as.data.frame(.) %>%
    tidyr::pivot_longer(cols = everything()) %>%
    mutate_at("name", stringr::str_replace, "theta_", "") %>%
    rename(Reg = name) %>%
    mutate(name = rep(sort(rep(c(1:11), 10)), 10)) %>%
    mutate(type = rep("theta", length(.[,1]))) %>%
    select(Reg, name, value, type)
  
  tmp<-readRDS(file="./RDS/OUwie_boots2.latitude.RDS")
  HDInterval::hdi(log(2)/tmp[,1]) #get the HPD 1/2 life
  median(log(2)/tmp[,1])
  
  
}

# #codon usage NOT YET RUN -- RETURN TO THIS AFTER FORMATTING OTHER CODE
# {
#   #fit the shifting theta model to enc
#   # exon.enc.OUM <- OUwie(
#   #   phy = simmap.janus.nuc.exons,
#   #   data = OUwie_data.codons[, c(1, 2, 3)],
#   #   model = 'OUM',
#   #   simmap.tree = T,
#   #   diagn = T,
#   #   check.identify = T,
#   #   algorithm = "invert",
#   #   scaleHeight = F,
#   #   mserr = "none"
#   # )
#   # saveRDS(exon.enc.OUM, file="./RDS/exon.enc.OUM.RDS")
#   exon.enc.OUM <- readRDS(file="./RDS/exon.enc.OUM.RDS")
#   
#   #try with parametric bootstrapping
#   # OUwie_boots.enc <-
#   #   OUwie.boot(
#   #     phy = simmap.janus.nuc.exons,
#   #     data = OUwie_data.codons[, c(1, 2, 3)],
#   #     model = "OUM",
#   #     nboot = 100,
#   #     alpha = exon.enc.OUM$solution[1, ],
#   #     sigma.sq = exon.enc.OUM$solution[2, ],
#   #     theta = exon.enc.OUM$theta[, 1],
#   #     theta0 = exon.enc.OUM$theta[1, 1],
#   #     algorithm = "invert",
#   #     simmap.tree = T,
#   #     quiet = F,
#   #     mserr = "none"
#   #   )
#   #saveRDS(OUwie_boots.enc, file="./RDS/OUwie_boots.enc.RDS")
#   OUwie_boots.enc<- readRDS(file="./RDS/OUwie_boots.enc.RDS")
#   
#   # exon.encprime.OUM <- OUwie(
#   #   phy = simmap.janus.nuc.exons,
#   #   data = OUwie_data.codons[, c(1, 2, 4)],
#   #   model = 'OUM',
#   #   simmap.tree = T,
#   #   diagn = T,
#   #   check.identify = T,
#   #   algorithm = "invert",
#   #   scaleHeight = F,
#   #   mserr = "none"
#   # )
#   # saveRDS(exon.encprime.OUM, file="./RDS/exon.encprime.OUM.RDS")
#   exon.encprime.OUM <- readRDS(file="./RDS/exon.encprime.OUM.RDS")
#   
#   #try with parametric bootstrapping
#   OUwie_boots.encprime <-
#     OUwie.boot(
#       phy = simmap.janus.nuc.exons,
#       data = OUwie_data.codons[, c(1, 2, 4)],
#       model = "OUM",
#       nboot = 100,
#       alpha = exon.encprime.OUM$solution[1, ],
#       sigma.sq = exon.encprime.OUM$solution[2, ],
#       theta = exon.encprime.OUM$theta[, 1],
#       theta0 = exon.encprime.OUM$theta[1, 1],
#       algorithm = "invert",
#       simmap.tree = T,
#       quiet = F,
#       mserr = "none"
#     )
#   saveRDS(OUwie_boots.encprime, file="./RDS/OUwie_boots.encprime.RDS")
#   OUwie_boots.encprime<- readRDS(file="./RDS/OUwie_boots.encprime.RDS")
#   
#   # exon.scuo.OUM <- OUwie(
#   #   phy = simmap.janus.nuc.exons,
#   #   data = OUwie_data.codons[, c(1, 2, 5)],
#   #   model = 'OUM',
#   #   simmap.tree = T,
#   #   diagn = T,
#   #   check.identify = T,
#   #   algorithm = "invert",
#   #   scaleHeight = F,
#   #   mserr = "none"
#   # )
#   # saveRDS(exon.scuo.OUM, file="./RDS/exon.scuo.OUM.RDS")
#   exon.scuo.OUM <- readRDS(file="./RDS/exon.scuo.OUM.RDS")
#   
#   #try with parametric bootstrapping
#   # OUwie_boots.scuo <-
#   #   OUwie.boot(
#   #     phy = simmap.janus.nuc.exons,
#   #     data = OUwie_data.codons[, c(1, 2, 5)],
#   #     model = "OUM",
#   #     nboot = 100,
#   #     alpha = exon.scuo.OUM$solution[1, ],
#   #     sigma.sq = exon.scuo.OUM$solution[2, ],
#   #     theta = exon.scuo.OUM$theta[, 1],
#   #     theta0 = exon.scuo.OUM$theta[1, 1],
#   #     algorithm = "invert",
#   #     simmap.tree = T,
#   #     quiet = F,
#   #     mserr = "none"
#   #   )
#   # saveRDS(OUwie_boots.scuo, file="./RDS/OUwie_boots.scuo.RDS")
#   OUwie_boots.scuo <- readRDS(file="./RDS/OUwie_boots.scuo.RDS")
#   
# }

#plotting setup
{
order <- c(
  "root_gallo",
  "tinamiformes",
  "tinamiformes_sister",
  "columbea",
  "otidae",
  "otidae_sister",
  "aeqornithes",
  "coraciimorphae",
  "psittaciformes",
  "passeri"
)

  
# order <- c(
#     "root_struthio",
#     "tinamiformes",
#     "tinamiformes_sister",
#     "neognathae_gallo",
#     "columbea",
#     "otidae",
#     "otidae_sister",
#     "aeqornithes",
#     "coraciimorphae",
#     "psittaciformes",
#     "passeri"
# )
  
# sims.chickPC1 <-
#   sim_plotter(
#     janustree = simmap.janus.nuc.aggregate.simplified,
#     ouwiefit = chickPC1.OUM,
#     nsim = 10000,
#     clrs = setNames(rainbow(n = 10), unique(
#       getStates(simmap.janus.nuc.aggregate.simplified, type = 'tips')
#     )),
#     width = 0.05,
#     plot = F
#   ) %>%
#   tidyr::pivot_longer(!Reg) %>%
#   mutate(Reg = factor(Reg, levels = order)) %>%
#   arrange(Reg) %>%
#   mutate(type=rep("expected", length(.[,1])))
# saveRDS(sims.chickPC1, "./RDS/sims.chickPC1.RDS")
sims.chickPC1 <- readRDS("./RDS/sims.chickPC1.RDS")

# sims.mass <-
#   sim_plotter(
#     janustree = simmap.janus.nuc.aggregate.simplified,
#     ouwiefit = mass.OUM,
#     nsim = 10000,
#     clrs = setNames(rainbow(n = 10), unique( #change back to n = 10 for regular tree
#       getStates(simmap.janus.nuc.aggregate.simplified, type = 'tips')
#     )),
#     width = 0.05,
#     plot = F
#   ) %>%
#   tidyr::pivot_longer(!Reg) %>%
#   mutate(Reg = factor(Reg, levels = order)) %>%
#   arrange(Reg) %>%
#   mutate(type=rep("expected", length(.[,1])))
# saveRDS(sims.mass, "./RDS/sims.mass.RDS")
sims.mass <- readRDS("./RDS/sims.mass.RDS")

# sims.clutch <-
#   sim_plotter(
#     janustree = simmap.janus.nuc.aggregate.simplified,
#     ouwiefit = clutch.OUM,
#     nsim = 10000,
#     clrs = setNames(rainbow(n = 10), unique(
#       getStates(simmap.janus.nuc.aggregate.simplified, type = 'tips')
#     )),
#     width = 0.05,
#     plot = F
#   ) %>%
#   tidyr::pivot_longer(!Reg) %>%
#   mutate(Reg = factor(Reg, levels = order)) %>%
#   arrange(Reg) %>%
#   mutate(type=rep("expected", length(.[,1])))
# saveRDS(sims.clutch, "./RDS/sims.clutch.RDS")
sims.clutch<- readRDS("./RDS/sims.clutch.RDS")

# sims.gen_length <-
#   sim_plotter(
#     janustree = simmap.janus.nuc.aggregate.simplified,
#     ouwiefit = gen_length.OUM,
#     nsim = 10000,
#     clrs = setNames(rainbow(n = 10), unique(
#       getStates(simmap.janus.nuc.aggregate.simplified, type = 'tips')
#     )),
#     width = 0.05,
#     plot = F
#   ) %>%
#   tidyr::pivot_longer(!Reg) %>%
#   mutate(Reg = factor(Reg, levels = order)) %>%
#   arrange(Reg) %>%
#   mutate(type=rep("expected", length(.[,1])))
# saveRDS(sims.gen_length, "./RDS/sims.gen_length.RDS")
sims.gen_length <- readRDS("./RDS/sims.gen_length.RDS")


# sims.survival <-
#   sim_plotter(
#     janustree = simmap.janus.nuc.aggregate.simplified,
#     ouwiefit = survival.OUM,
#     nsim = 10000,
#     clrs = setNames(rainbow(n = 10), unique(
#       getStates(simmap.janus.nuc.aggregate.simplified, type = 'tips')
#     )),
#     width = 0.05,
#     plot = F
#   ) %>%
#   tidyr::pivot_longer(!Reg) %>%
#   mutate(Reg = factor(Reg, levels = order)) %>%
#   arrange(Reg) %>%
#   mutate(type=rep("expected", length(.[,1])))
# saveRDS(sims.survival, "./RDS/sims.survival.RDS")
sims.survival<- readRDS("./RDS/sims.survival.RDS")

# sims.breeding <-
#   sim_plotter(
#     janustree = simmap.janus.nuc.aggregate.simplified,
#     ouwiefit = breeding.OUM,
#     nsim = 10000,
#     clrs = setNames(rainbow(n = 10), unique(
#       getStates(simmap.janus.nuc.aggregate.simplified, type = 'tips')
#     )),
#     width = 0.05,
#     plot = F
#   ) %>%
#   tidyr::pivot_longer(!Reg) %>%
#   mutate(Reg = factor(Reg, levels = order)) %>%
#   arrange(Reg) %>%
#   mutate(type=rep("expected", length(.[,1])))
# saveRDS(sims.breeding, "./RDS/sims.breeding.RDS")
sims.breeding<- readRDS('./RDS/sims.breeding.RDS')


# sims.longevity <-
#   sim_plotter(
#     janustree = simmap.janus.nuc.aggregate.simplified,
#     ouwiefit = longevity.OUM,
#     nsim = 10000,
#     clrs = setNames(rainbow(n = 10), unique(
#       getStates(simmap.janus.nuc.aggregate.simplified, type = 'tips')
#     )),
#     width = 0.05,
#     plot = F
#   ) %>%
#   tidyr::pivot_longer(!Reg) %>%
#   mutate(Reg = factor(Reg, levels = order)) %>%
#   arrange(Reg) %>%
#   mutate(type=rep("expected", length(.[,1])))
# saveRDS(sims.longevity, "./RDS/sims.longevity.RDS")
sims.longevity<- readRDS("./RDS/sims.longevity.RDS")

# sims.latitude <-
#   sim_plotter(
#     janustree = simmap.janus.nuc.aggregate.simplified,
#     ouwiefit = latitude.OUM,
#     nsim = 10000,
#     clrs = setNames(rainbow(n = 10), unique(
#       getStates(simmap.janus.nuc.aggregate.simplified, type = 'tips')
#     )),
#     width = 0.05,
#     plot = F
#   ) %>%
#   tidyr::pivot_longer(!Reg) %>%
#   mutate(Reg = factor(Reg, levels = order)) %>%
#   arrange(Reg) %>%
#   mutate(type=rep("expected", length(.[,1])))
# saveRDS(sims.latitude, "./RDS/sims.latitude.RDS")
sims.latitude <- readRDS("./RDS/sims.latitude.RDS")


# sims.encprime <-
#   sim_plotter(
#     janustree = simmap.janus.nuc.exons,
#     ouwiefit = exon.encprime.OUM,
#     nsim = 1000,
#     clrs = setNames(rainbow(n = 8), unique(
#       getStates(simmap.janus.nuc.aggregate.simplified, type = 'tips')
#     )),
#     width = 0.05,
#     plot = T
#   ) %>%
#   tidyr::pivot_longer(!Reg) %>%
#   mutate(Reg = factor(Reg, levels = order)) %>%
#   arrange(Reg) %>%
#   mutate(type=rep("theta", length(.[,1])))


}

#more plotting setup
{
Plot_data.mass <- rbind(
  OUwie_data[, c(2, 3)] %>%
    tidyr::pivot_longer(!Reg) %>%
    mutate(Reg = factor(Reg, levels = order)) %>%
    arrange(Reg) %>%
    mutate(type = rep("data", length(.[, 1]))),
  sims.mass, OUwie_boots.mass)

#Plot_data.mass$value <- Plot_data.mass$value * attributes(scale(log(megaLHT$mass)))[[3]] #unscale
#Plot_data.mass$value <- Plot_data.mass$value + attributes(scale(log(megaLHT$mass)))[[2]] #uncenter
#masslabs <- ((c(-2, 0, 2, 4)) * attributes(scale(log(megaLHT$mass)))[[3]]) + attributes(scale(log(megaLHT$mass)))[[2]]
#max(Plot_data.mass[Plot_data.mass$type=='data',]$value)
  
#colors for regimes
regimes<-c("black", "#ddcc77", "#aa4499","#6699cc", "#117733", "#999933", "#332288", "#88ccee", "#888888", "#882255")
#regimes<-c("black", "#ddcc77", "#aa4499",'#661100',"#6699cc", "#117733", "#999933", "#332288", "#88ccee", "#888888", "#882255")
#second row is alt test with 11 groups

Plot_data.chickPC1 <- rbind(
  OUwie_data[, c(2, 4)] %>%
    tidyr::pivot_longer(!Reg) %>%
    mutate(Reg = factor(Reg, levels = order)) %>%
    arrange(Reg) %>%
    mutate(type = rep("data", length(.[, 1]))),
  sims.chickPC1, OUwie_boots.chickPC1)

Plot_data.clutch <- rbind(
  OUwie_data.extra[, c(2, 3)] %>%
    tidyr::pivot_longer(!Reg) %>%
    mutate(Reg = factor(Reg, levels = order)) %>%
    arrange(Reg) %>%
    mutate(type = rep("data", length(.[, 1]))),
  sims.clutch, OUwie_boots.clutch)

Plot_data.gen_length <- rbind(
  OUwie_data.extra[, c(2, 4)] %>%
    tidyr::pivot_longer(!Reg) %>%
    mutate(Reg = factor(Reg, levels = order)) %>%
    arrange(Reg) %>%
    mutate(type = rep("data", length(.[, 1]))),
  sims.gen_length, OUwie_boots.gen_length)

Plot_data.survival <- rbind(
  OUwie_data.extra[, c(2, 5)] %>%
    tidyr::pivot_longer(!Reg) %>%
    mutate(Reg = factor(Reg, levels = order)) %>%
    arrange(Reg) %>%
    mutate(type = rep("data", length(.[, 1]))),
  sims.survival, OUwie_boots.survival)


Plot_data.breeding <- rbind(
  OUwie_data.extra[, c(2, 6)] %>%
    tidyr::pivot_longer(!Reg) %>%
    mutate(Reg = factor(Reg, levels = order)) %>%
    arrange(Reg) %>%
    mutate(type = rep("data", length(.[, 1]))),
  sims.breeding, OUwie_boots.breeding)


Plot_data.longevity <- rbind(
  OUwie_data.extra[, c(2, 7)] %>%
    tidyr::pivot_longer(!Reg) %>%
    mutate(Reg = factor(Reg, levels = order)) %>%
    arrange(Reg) %>%
    mutate(type = rep("data", length(.[, 1]))),
  sims.longevity, OUwie_boots.longevity)


Plot_data.latitude <- rbind(
  OUwie_data.extra[, c(2, 8)] %>%
    tidyr::pivot_longer(!Reg) %>%
    mutate(Reg = factor(Reg, levels = order)) %>%
    arrange(Reg) %>%
    mutate(type = rep("data", length(.[, 1]))),
  sims.latitude, OUwie_boots.latitude)


}

#generate ggplot objects 
{
#require(ggpubr)
massplot <-
  {ggplot(Plot_data.mass,
           aes(x = value,
               y = Reg,
               fill = Reg,
               color= Reg,
               stroke=0.01, height=..ndensity..)) +
    #height = ..density..)) +
    #plotting theta
    geom_density_ridges(color="grey",
      quantile_lines = F,
      quantile_fun = estimate_mode,
      data = dplyr::filter(Plot_data.mass, type == "expected"),
      #stat = "density",
      scale = 0.85,
      size = 0.01,
      lty = 1,
      rel_min_height = 0.001,
      fill = c("grey"),
      #viridis(2)[1],#"#fc8d62",
      alpha = 0.25,
      #0.4,
      bandwidth = 0.35
    ) + #scale_fill_manual(values = regimes, labels = NULL) +
    #plotting expected
    geom_density_ridges(
      quantile_lines = TRUE,
      quantile_fun = estimate_mode,
      data = dplyr::filter(Plot_data.mass, type == "theta"),
      #stat = "density",
      scale = 0.9,
      size = 0.75,
      rel_min_height = 0.001,
      #fill = c("red"),#viridis(2)[1],#"#fc8d62",
      alpha = 0.6,
      #0.4,
      bandwidth = 0.7
    ) + scale_fill_manual(values = regimes, labels = NULL) + scale_color_manual(values=regimes) + 
    #plotting data
    # geom_density_ridges(
    #   quantile_lines = F,
    #   quantile_fun = estimate_mode,
    #   data = dplyr::filter(Plot_data.mass, type == "data"),
    #   #stat = "density",
    #   scale = 0.9,
    #   size = 0.5,
    #   rel_min_height = 0.01,
    #   fill = "white",
    #   #viridis(2)[2],#"#8da0cb",
  #   alpha = 0.75,
  #   bandwidth = 0.65
  # ) +
  geom_density_ridges(size=0, color="black", fill="white",
    data = dplyr::filter(Plot_data.mass, type == "data"),
    jittered_points = TRUE,
    position = position_raincloud(width=0.75, height=0.05, ygap = -0.025, seed=5), scale=0, point_alpha=0.5, point_size=1.5, point_shape=21
    #position = position_points_jitter(width = 0.1, height = 0, yoffset = -0.025),
    #point_shape = '|', point_size = 1, point_alpha = 0.75, alpha = 0, scale=0, bandwidth = 0.5,
  ) + 
  theme_ridges() +
    theme(
      #axis.text.x=element_blank(), #remove x axis labels
      #axis.ticks.x=element_blank(), #remove x axis ticks
      axis.text.y = element_blank(),
      #remove y axis labels
      #axis.ticks.y=element_blank()  #remove y axis ticks
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(size=8.5),
      #axis.text.x = element_text(vjust = 10),
      legend.position = "none"
    ) + grids(axis="x", linetype = "dotted") + grids(axis="y", linetype = "solid", size=1.25) +
    scale_x_continuous(breaks = c(2.5, 5, 7.5, 10, 12.5), limits = c(0.5, 13)) +
    scale_y_discrete(expand = c(0.015, 0))
  }
  
chickpc1plot <-
  {ggplot(Plot_data.chickPC1,
           aes(x = value,
               y = Reg,
               fill = Reg, 
               color= Reg,
               stroke=0.01, height=..ndensity..)) +
    #height = ..density..)) +
    #plotting theta
    geom_density_ridges(color="grey",
      quantile_lines = F,
      quantile_fun = estimate_mode,
      data = dplyr::filter(Plot_data.chickPC1, type == "expected"),
      #stat = "density",
      scale = 0.9,
      size = 0.01,
      lty = 1,
      rel_min_height = 0.001,
      fill = c("grey"),
      #viridis(2)[1],#"#fc8d62",
      alpha = 0.25,
      #0.4,
      bandwidth = 0.35
    ) + #scale_fill_manual(values = regimes, labels = NULL) +
    #plotting expected
    geom_density_ridges(
      quantile_lines = TRUE,
      quantile_fun = estimate_mode,
      data = dplyr::filter(Plot_data.chickPC1, type == "theta"),
      #stat = "density",
      scale = 0.9,
      size = 0.75,
      rel_min_height = 0.001,
      #fill = c("red"),#viridis(2)[1],#"#fc8d62",
      alpha = 0.6,
      #0.4,
      bandwidth = 0.7
    ) + scale_fill_manual(values = regimes, labels = NULL) + scale_color_manual(values=regimes) + 
    #plotting data
    # geom_density_ridges(
    #   quantile_lines = F,
    #   quantile_fun = estimate_mode,
    #   data = dplyr::filter(Plot_data.chickPC1, type == "data"),
    #   #stat = "density",
    #   scale = 0.9,
    #   size = 0.5,
    #   rel_min_height = 0.01,
    #   fill = "white",
    #   #viridis(2)[2],#"#8da0cb",
  #   alpha = 0.75,
  #   bandwidth = 0.65
  # ) +
  geom_density_ridges(size=0, color="black", fill="white",
    data = dplyr::filter(Plot_data.chickPC1, type == "data"),
    jittered_points = TRUE,
    position = position_raincloud(width=0.75, height=0.05, ygap = -0.025, seed=1), scale=0, point_alpha=0.5, point_size=1.5, point_shape=21,
    #position = position_points_jitter(width = 0.1, height = 0, yoffset = -0.025),
    #point_shape = '|', point_size = 1, point_alpha = 0.75, alpha = 0, scale=0, bandwidth = 0.5
  ) + #scale_fill_manual(values = regimes, labels = NULL) #+ scale_color_manual(values="black")
  theme_ridges() +
    theme(
      #axis.text.x=element_blank(), #remove x axis labels
      #axis.ticks.x=element_blank(), #remove x axis ticks
      #axis.text.y=element_blank(),  #remove y axis labels
      #axis.ticks.y=element_blank()  #remove y axis ticks
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(size=8.5),
      legend.position = "none"
    ) + grids(axis="x", linetype = "dotted") + grids(axis="y", linetype = "solid", size=1.25) +
    scale_x_continuous(breaks = c(-4, -2, 0, 2, 4), limits = c(-5, 5)) +
    scale_y_discrete(expand = c(0.015, 0))
    }

clutchplot <-
  {ggplot(Plot_data.clutch,
         aes(x = value,
             y = Reg,
             fill = Reg, 
             color= Reg,
             stroke=0.01, height=..ndensity..)) +
  #height = ..density..)) +
  #plotting theta
  geom_density_ridges(color="grey",
                      quantile_lines = F,
                      quantile_fun = estimate_mode,
                      data = dplyr::filter(Plot_data.clutch, type == "expected"),
                      #stat = "density",
                      scale = 0.9,
                      size = 0.01,
                      lty = 1,
                      rel_min_height = 0.001,
                      fill = c("grey"),
                      #viridis(2)[1],#"#fc8d62",
                      alpha = 0.25,
                      #0.4,
                      bandwidth = 0.35
  ) + #scale_fill_manual(values = regimes, labels = NULL) +
  #plotting expected
  geom_density_ridges(
    quantile_lines = TRUE,
    quantile_fun = estimate_mode,
    data = dplyr::filter(Plot_data.clutch, type == "theta"),
    #stat = "density",
    scale = 0.9,
    size = 0.75,
    rel_min_height = 0.001,
    #fill = c("red"),#viridis(2)[1],#"#fc8d62",
    alpha = 0.6,
    #0.4,
    bandwidth = 0.7
  ) + scale_fill_manual(values = regimes, labels = NULL) + scale_color_manual(values=regimes) + 
  #plotting data
  # geom_density_ridges(
  #   quantile_lines = F,
  #   quantile_fun = estimate_mode,
  #   data = dplyr::filter(Plot_data.chickPC1, type == "data"),
  #   #stat = "density",
  #   scale = 0.9,
  #   size = 0.5,
  #   rel_min_height = 0.01,
  #   fill = "white",
  #   #viridis(2)[2],#"#8da0cb",
#   alpha = 0.75,
#   bandwidth = 0.65
# ) +
geom_density_ridges(size=0, color="black", fill="white",
                    data = dplyr::filter(Plot_data.clutch, type == "data"),
                    jittered_points = TRUE,
                    position = position_raincloud(width=0.75, height=0.05, ygap = -0.025, seed=1), scale=0, point_alpha=0.5, point_size=1.5, point_shape=21,
                    #position = position_points_jitter(width = 0.1, height = 0, yoffset = -0.025),
                    #point_shape = '|', point_size = 1, point_alpha = 0.75, alpha = 0, scale=0, bandwidth = 0.5
) + #scale_fill_manual(values = regimes, labels = NULL) #+ scale_color_manual(values="black")
  theme_ridges() +
      theme(
        #axis.text.x=element_blank(), #remove x axis labels
        #axis.ticks.x=element_blank(), #remove x axis ticks
        axis.text.y = element_blank(),
        #remove y axis labels
        #axis.ticks.y=element_blank()  #remove y axis ticks
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=8.5),
        #axis.text.x = element_text(vjust = 10),
        legend.position = "none"
      ) + grids(axis="x", linetype = "dotted") + grids(axis="y", linetype = "solid", size=1.25) +
  scale_x_continuous(breaks = c(-4, -2, 0, 2, 4), limits = c(-5, 5)) +
  scale_y_discrete(expand = c(0.015, 0))
  }

gen_length_plot <-
  {ggplot(Plot_data.gen_length,
         aes(x = value,
             y = Reg,
             fill = Reg, 
             color= Reg,
             stroke=0.01, height=..ndensity..)) +
  #height = ..density..)) +
  #plotting theta
  geom_density_ridges(color="grey",
                      quantile_lines = F,
                      quantile_fun = estimate_mode,
                      data = dplyr::filter(Plot_data.gen_length, type == "expected"),
                      #stat = "density",
                      scale = 0.9,
                      size = 0.01,
                      lty = 1,
                      rel_min_height = 0.001,
                      fill = c("grey"),
                      #viridis(2)[1],#"#fc8d62",
                      alpha = 0.25,
                      #0.4,
                      bandwidth = 0.35
  ) + #scale_fill_manual(values = regimes, labels = NULL) +
  #plotting expected
  geom_density_ridges(
    quantile_lines = TRUE,
    quantile_fun = estimate_mode,
    data = dplyr::filter(Plot_data.gen_length, type == "theta"),
    #stat = "density",
    scale = 0.9,
    size = 0.75,
    rel_min_height = 0.001,
    #fill = c("red"),#viridis(2)[1],#"#fc8d62",
    alpha = 0.6,
    #0.4,
    bandwidth = 0.7
  ) + scale_fill_manual(values = regimes, labels = NULL) + scale_color_manual(values=regimes) + 
  #plotting data
  # geom_density_ridges(
  #   quantile_lines = F,
  #   quantile_fun = estimate_mode,
  #   data = dplyr::filter(Plot_data.chickPC1, type == "data"),
  #   #stat = "density",
  #   scale = 0.9,
  #   size = 0.5,
  #   rel_min_height = 0.01,
  #   fill = "white",
  #   #viridis(2)[2],#"#8da0cb",
#   alpha = 0.75,
#   bandwidth = 0.65
# ) +
geom_density_ridges(size=0, color="black", fill="white",
                    data = dplyr::filter(Plot_data.gen_length, type == "data"),
                    jittered_points = TRUE,
                    position = position_raincloud(width=0.75, height=0.05, ygap = -0.025, seed=1), scale=0, point_alpha=0.5, point_size=1.5, point_shape=21,
                    #position = position_points_jitter(width = 0.1, height = 0, yoffset = -0.025),
                    #point_shape = '|', point_size = 1, point_alpha = 0.75, alpha = 0, scale=0, bandwidth = 0.5
) + #scale_fill_manual(values = regimes, labels = NULL) #+ scale_color_manual(values="black")
  theme_ridges() +
      theme(
        #axis.text.x=element_blank(), #remove x axis labels
        #axis.ticks.x=element_blank(), #remove x axis ticks
        axis.text.y = element_blank(),
        #remove y axis labels
        #axis.ticks.y=element_blank()  #remove y axis ticks
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=8.5),
        #axis.text.x = element_text(vjust = 10),
        legend.position = "none"
      ) + grids(axis="x", linetype = "dotted") + grids(axis="y", linetype = "solid", size=1.25) +
  scale_x_continuous(breaks = c(-4, -2, 0, 2, 4), limits = c(-5, 5)) +
  scale_y_discrete(expand = c(0.015, 0))
  }

survival_plot <-
  {ggplot(Plot_data.survival,
          aes(x = value,
              y = Reg,
              fill = Reg, 
              color= Reg,
              stroke=0.01, height=..ndensity..)) +
      #height = ..density..)) +
      #plotting theta
      geom_density_ridges(color="grey",
                          quantile_lines = F,
                          quantile_fun = estimate_mode,
                          data = dplyr::filter(Plot_data.survival, type == "expected"),
                          #stat = "density",
                          scale = 0.9,
                          size = 0.01,
                          lty = 1,
                          rel_min_height = 0.001,
                          fill = c("grey"),
                          #viridis(2)[1],#"#fc8d62",
                          alpha = 0.25,
                          #0.4,
                          bandwidth = 0.35
      ) + #scale_fill_manual(values = regimes, labels = NULL) +
      #plotting expected
      geom_density_ridges(
        quantile_lines = TRUE,
        quantile_fun = estimate_mode,
        data = dplyr::filter(Plot_data.survival, type == "theta"),
        #stat = "density",
        scale = 0.9,
        size = 0.75,
        rel_min_height = 0.001,
        #fill = c("red"),#viridis(2)[1],#"#fc8d62",
        alpha = 0.6,
        #0.4,
        bandwidth = 0.7
      ) + scale_fill_manual(values = regimes, labels = NULL) + scale_color_manual(values=regimes) + 
      #plotting data
      # geom_density_ridges(
      #   quantile_lines = F,
      #   quantile_fun = estimate_mode,
      #   data = dplyr::filter(Plot_data.chickPC1, type == "data"),
      #   #stat = "density",
      #   scale = 0.9,
      #   size = 0.5,
      #   rel_min_height = 0.01,
      #   fill = "white",
      #   #viridis(2)[2],#"#8da0cb",
    #   alpha = 0.75,
    #   bandwidth = 0.65
    # ) +
    geom_density_ridges(size=0, color="black", fill="white",
                        data = dplyr::filter(Plot_data.survival, type == "data"),
                        jittered_points = TRUE,
                        position = position_raincloud(width=0.75, height=0.05, ygap = -0.025, seed=1), scale=0, point_alpha=0.5, point_size=1.5, point_shape=21,
                        #position = position_points_jitter(width = 0.1, height = 0, yoffset = -0.025),
                        #point_shape = '|', point_size = 1, point_alpha = 0.75, alpha = 0, scale=0, bandwidth = 0.5
    ) + #scale_fill_manual(values = regimes, labels = NULL) #+ scale_color_manual(values="black")
      theme_ridges() +
      theme(
        #axis.text.x=element_blank(), #remove x axis labels
        #axis.ticks.x=element_blank(), #remove x axis ticks
        axis.text.y = element_blank(),
        #remove y axis labels
        #axis.ticks.y=element_blank()  #remove y axis ticks
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=8.5),
        #axis.text.x = element_text(vjust = 10),
        legend.position = "none"
      ) + grids(axis="x", linetype = "dotted") + grids(axis="y", linetype = "solid", size=1.25) +
      scale_x_continuous(breaks = c(-4, -2, 0, 2, 4), limits = c(-5, 5)) +
      scale_y_discrete(expand = c(0.015, 0))
  }

breeding_plot <-
  {ggplot(Plot_data.breeding,
          aes(x = value,
              y = Reg,
              fill = Reg, 
              color= Reg,
              stroke=0.01, height=..ndensity..)) +
      #height = ..density..)) +
      #plotting theta
      geom_density_ridges(color="grey",
                          quantile_lines = F,
                          quantile_fun = estimate_mode,
                          data = dplyr::filter(Plot_data.breeding, type == "expected"),
                          #stat = "density",
                          scale = 0.9,
                          size = 0.01,
                          lty = 1,
                          rel_min_height = 0.001,
                          fill = c("grey"),
                          #viridis(2)[1],#"#fc8d62",
                          alpha = 0.25,
                          #0.4,
                          bandwidth = 0.35
      ) + #scale_fill_manual(values = regimes, labels = NULL) +
      #plotting expected
      geom_density_ridges(
        quantile_lines = TRUE,
        quantile_fun = estimate_mode,
        data = dplyr::filter(Plot_data.breeding, type == "theta"),
        #stat = "density",
        scale = 0.9,
        size = 0.75,
        rel_min_height = 0.001,
        #fill = c("red"),#viridis(2)[1],#"#fc8d62",
        alpha = 0.6,
        #0.4,
        bandwidth = 0.7
      ) + scale_fill_manual(values = regimes, labels = NULL) + scale_color_manual(values=regimes) + 
      #plotting data
      # geom_density_ridges(
      #   quantile_lines = F,
      #   quantile_fun = estimate_mode,
      #   data = dplyr::filter(Plot_data.chickPC1, type == "data"),
      #   #stat = "density",
      #   scale = 0.9,
      #   size = 0.5,
      #   rel_min_height = 0.01,
      #   fill = "white",
      #   #viridis(2)[2],#"#8da0cb",
    #   alpha = 0.75,
    #   bandwidth = 0.65
    # ) +
    geom_density_ridges(size=0, color="black", fill="white",
                        data = dplyr::filter(Plot_data.breeding, type == "data"),
                        jittered_points = TRUE,
                        position = position_raincloud(width=0.75, height=0.05, ygap = -0.025, seed=1), scale=0, point_alpha=0.5, point_size=1.5, point_shape=21,
                        #position = position_points_jitter(width = 0.1, height = 0, yoffset = -0.025),
                        #point_shape = '|', point_size = 1, point_alpha = 0.75, alpha = 0, scale=0, bandwidth = 0.5
    ) + #scale_fill_manual(values = regimes, labels = NULL) #+ scale_color_manual(values="black")
      theme_ridges() +
      theme(
        #axis.text.x=element_blank(), #remove x axis labels
        #axis.ticks.x=element_blank(), #remove x axis ticks
        axis.text.y = element_blank(),
        #remove y axis labels
        #axis.ticks.y=element_blank()  #remove y axis ticks
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=8.5),
        #axis.text.x = element_text(vjust = 10),
        legend.position = "none"
      ) + grids(axis="x", linetype = "dotted") + grids(axis="y", linetype = "solid", size=1.25) +
      scale_x_continuous(breaks = c(-4, -2, 0, 2, 4), limits = c(-5, 5)) +
      scale_y_discrete(expand = c(0.015, 0))
  }

longevity_plot <-
  {ggplot(Plot_data.longevity,
          aes(x = value,
              y = Reg,
              fill = Reg, 
              color= Reg,
              stroke=0.01, height=..ndensity..)) +
      #height = ..density..)) +
      #plotting theta
      geom_density_ridges(color="grey",
                          quantile_lines = F,
                          quantile_fun = estimate_mode,
                          data = dplyr::filter(Plot_data.longevity, type == "expected"),
                          #stat = "density",
                          scale = 0.9,
                          size = 0.01,
                          lty = 1,
                          rel_min_height = 0.001,
                          fill = c("grey"),
                          #viridis(2)[1],#"#fc8d62",
                          alpha = 0.25,
                          #0.4,
                          bandwidth = 0.35
      ) + #scale_fill_manual(values = regimes, labels = NULL) +
      #plotting expected
      geom_density_ridges(
        quantile_lines = TRUE,
        quantile_fun = estimate_mode,
        data = dplyr::filter(Plot_data.longevity, type == "theta"),
        #stat = "density",
        scale = 0.9,
        size = 0.75,
        rel_min_height = 0.001,
        #fill = c("red"),#viridis(2)[1],#"#fc8d62",
        alpha = 0.6,
        #0.4,
        bandwidth = 0.7
      ) + scale_fill_manual(values = regimes, labels = NULL) + scale_color_manual(values=regimes) + 
      #plotting data
      # geom_density_ridges(
      #   quantile_lines = F,
      #   quantile_fun = estimate_mode,
      #   data = dplyr::filter(Plot_data.chickPC1, type == "data"),
      #   #stat = "density",
      #   scale = 0.9,
      #   size = 0.5,
      #   rel_min_height = 0.01,
      #   fill = "white",
      #   #viridis(2)[2],#"#8da0cb",
    #   alpha = 0.75,
    #   bandwidth = 0.65
    # ) +
    geom_density_ridges(size=0, color="black", fill="white",
                        data = dplyr::filter(Plot_data.longevity, type == "data"),
                        jittered_points = TRUE,
                        position = position_raincloud(width=0.75, height=0.05, ygap = -0.025, seed=1), scale=0, point_alpha=0.5, point_size=1.5, point_shape=21,
                        #position = position_points_jitter(width = 0.1, height = 0, yoffset = -0.025),
                        #point_shape = '|', point_size = 1, point_alpha = 0.75, alpha = 0, scale=0, bandwidth = 0.5
    ) + #scale_fill_manual(values = regimes, labels = NULL) #+ scale_color_manual(values="black")
      theme_ridges() +
      theme(
        #axis.text.x=element_blank(), #remove x axis labels
        #axis.ticks.x=element_blank(), #remove x axis ticks
        axis.text.y = element_blank(),
        #remove y axis labels
        #axis.ticks.y=element_blank()  #remove y axis ticks
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=8.5),
        #axis.text.x = element_text(vjust = 10),
        legend.position = "none"
      ) + grids(axis="x", linetype = "dotted") + grids(axis="y", linetype = "solid", size=1.25) +
      scale_x_continuous(breaks = c(-4, -2, 0, 2, 4), limits = c(-5, 5)) +
      scale_y_discrete(expand = c(0.015, 0))
  }

latitude_plot <-
  {ggplot(Plot_data.latitude,
          aes(x = value,
              y = Reg,
              fill = Reg, 
              color= Reg,
              stroke=0.01, height=..ndensity..)) +
      #height = ..density..)) +
      #plotting theta
      geom_density_ridges(color="grey",
                          quantile_lines = F,
                          quantile_fun = estimate_mode,
                          data = dplyr::filter(Plot_data.latitude, type == "expected"),
                          #stat = "density",
                          scale = 0.9,
                          size = 0.01,
                          lty = 1,
                          rel_min_height = 0.001,
                          fill = c("grey"),
                          #viridis(2)[1],#"#fc8d62",
                          alpha = 0.25,
                          #0.4,
                          bandwidth = 0.35
      ) + #scale_fill_manual(values = regimes, labels = NULL) +
      #plotting expected
      geom_density_ridges(
        quantile_lines = TRUE,
        quantile_fun = estimate_mode,
        data = dplyr::filter(Plot_data.latitude, type == "theta"),
        #stat = "density",
        scale = 0.9,
        size = 0.75,
        rel_min_height = 0.001,
        #fill = c("red"),#viridis(2)[1],#"#fc8d62",
        alpha = 0.6,
        #0.4,
        bandwidth = 0.7
      ) + scale_fill_manual(values = regimes, labels = NULL) + scale_color_manual(values=regimes) + 
      #plotting data
      # geom_density_ridges(
      #   quantile_lines = F,
      #   quantile_fun = estimate_mode,
      #   data = dplyr::filter(Plot_data.chickPC1, type == "data"),
      #   #stat = "density",
      #   scale = 0.9,
      #   size = 0.5,
      #   rel_min_height = 0.01,
      #   fill = "white",
      #   #viridis(2)[2],#"#8da0cb",
    #   alpha = 0.75,
    #   bandwidth = 0.65
    # ) +
    geom_density_ridges(size=0, color="black", fill="white",
                        data = dplyr::filter(Plot_data.latitude, type == "data"),
                        jittered_points = TRUE,
                        position = position_raincloud(width=0.75, height=0.05, ygap = -0.025, seed=1), scale=0, point_alpha=0.5, point_size=1.5, point_shape=21,
                        #position = position_points_jitter(width = 0.1, height = 0, yoffset = -0.025),
                        #point_shape = '|', point_size = 1, point_alpha = 0.75, alpha = 0, scale=0, bandwidth = 0.5
    ) + #scale_fill_manual(values = regimes, labels = NULL) #+ scale_color_manual(values="black")
      theme_ridges() +
      theme(
        #axis.text.x=element_blank(), #remove x axis labels
        #axis.ticks.x=element_blank(), #remove x axis ticks
        axis.text.y = element_blank(),
        #remove y axis labels
        #axis.ticks.y=element_blank()  #remove y axis ticks
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=8.5),
        #axis.text.x = element_text(vjust = 10),
        legend.position = "none"
      ) + grids(axis="x", linetype = "dotted") + grids(axis="y", linetype = "solid", size=1.25) +
      scale_x_continuous(breaks = c(-4, -2, 0, 2, 4), limits = c(-5, 5)) +
      scale_y_discrete(expand = c(0.015, 0))
  }

}

#plotting notes/experimental/ commented out
{
# ggplot(Plot_data.mass, aes(x=value, y=Reg, fill=Reg)) +
#   geom_density_ridges(quantile_lines = TRUE, quantiles=2, data=dplyr::filter(Plot_data.mass, type=="theta"), scale = 0.9, size = 0.25, rel_min_height = 0.01, fill="red", alpha=0.4) +
#   geom_density_ridges(quantile_lines = TRUE, quantiles=2, data=dplyr::filter(Plot_data.mass, type=="data"), scale = 0.9, size = 0.25, rel_min_height = 0.01, fill="grey", alpha=0.9) +
#   theme_ridges()
# 
# 
# chickpc1plot<-ggplot(Plot_data.chickPC1, aes(x=value, y=Reg, fill=Reg)) +
#   geom_density_ridges(data=dplyr::filter(Plot_data.chickPC1, type=="theta"), scale = 0.9, size = 0.25, rel_min_height = 0.01, fill="red", alpha=0.4) +
#   geom_density_ridges(data=dplyr::filter(Plot_data.chickPC1, type=="data"), scale = 0.9, size = 0.25, rel_min_height = 0.01, fill="grey", alpha=0.9) +
#   theme_ridges()
  
  # 
  # # library
  # library(ggridges)
  # library(ggplot2)
  # library(viridis)
  # library(hrbrthemes)
  # 
  # # Plot
  # ggplot(sims.chickPC1, aes(x = `value`, y = `Reg`, fill = ..x..)) +
  #   geom_density_ridges_gradient(scale = 4, show.legend = F) +
  #   scale_fill_viridis(name = "Temp. [F]", option = "C") +
  #   labs(title = 'Body mass Theta') +
  #   theme_ipsum() +
  #   theme(
  #     legend.position="none",
  #     panel.spacing = unit(0.1, "lines"),
  #     strip.text.x = element_text(size = 8)
  #   )

}
}

#plotting
pdf(file="LH_shift_plots.pdf", height=11.5, width=7)
ggarrange(chickpc1plot, massplot, labels=c("ChickPC1", "Mass"),
          widths=c(1.7,1.0))
dev.off()

#plotting supplemental shifts
pdf(file="LH_shift_plots_supp.pdf", width=11*2, height=8.5*1.75)
ggarrange(chickpc1plot, massplot, clutchplot, gen_length_plot, survival_plot, breeding_plot,longevity_plot, latitude_plot,
          labels=c("ChickPC1", "Mass", "Clutch Size", "Generation Length", "Survival", "Age First Breeding", "Longevity", "Latitude"),
          widths=c(1.6,1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0), nrow=1)
dev.off()


#plotting setup for simplified 2 (copied from above, and modified)
#plotting setup
{

  order2 <- c(
      "root_struthio",
      "tinamiformes",
      "tinamiformes_sister",
      "neognathae_gallo",
      "columbea",
      "otidae",
      "otidae_sister",
      "aeqornithes",
      "coraciimorphae",
      "psittaciformes",
      "passeri"
  )
  
  # sims2.chickPC1 <-
  #   sim_plotter(
  #     janustree = simmap.janus.nuc.aggregate.simplified2,
  #     ouwiefit = chickPC1_2.OUM,
  #     nsim = 10000,
  #     clrs = setNames(rainbow(n = 11), unique(
  #       getStates(simmap.janus.nuc.aggregate.simplified2, type = 'tips')
  #     )),
  #     width = 0.05,
  #     plot = F
  #   ) %>%
  #   tidyr::pivot_longer(!Reg) %>%
  #   mutate(Reg = factor(Reg, levels = order2)) %>%
  #   arrange(Reg) %>%
  #   mutate(type=rep("expected", length(.[,1])))
  # saveRDS(sims2.chickPC1, "./RDS/sims2.chickPC1.RDS")
  sims2.chickPC1 <- readRDS("./RDS/sims2.chickPC1.RDS")
  
  # sims2.mass <-
  #   sim_plotter(
  #     janustree = simmap.janus.nuc.aggregate.simplified2,
  #     ouwiefit = mass2.OUM,
  #     nsim = 10000,
  #     clrs = setNames(rainbow(n = 11), unique( #change back to n = 10 for regular tree
  #       getStates(simmap.janus.nuc.aggregate.simplified2, type = 'tips')
  #     )),
  #     width = 0.05,
  #     plot = F
  #   ) %>%
  #   tidyr::pivot_longer(!Reg) %>%
  #   mutate(Reg = factor(Reg, levels = order2)) %>%
  #   arrange(Reg) %>%
  #   mutate(type=rep("expected", length(.[,1])))
  # saveRDS(sims2.mass, "./RDS/sims2.mass.RDS")
  sims2.mass <- readRDS("./RDS/sims2.mass.RDS")
  
  # sims2.clutch <-
  #   sim_plotter(
  #     janustree = simmap.janus.nuc.aggregate.simplified2,
  #     ouwiefit = clutch2.OUM,
  #     nsim = 10000,
  #     clrs = setNames(rainbow(n = 11), unique(
  #       getStates(simmap.janus.nuc.aggregate.simplified2, type = 'tips')
  #     )),
  #     width = 0.05,
  #     plot = F
  #   ) %>%
  #   tidyr::pivot_longer(!Reg) %>%
  #   mutate(Reg = factor(Reg, levels = order2)) %>%
  #   arrange(Reg) %>%
  #   mutate(type=rep("expected", length(.[,1])))
  # saveRDS(sims2.clutch, "./RDS/sims2.clutch.RDS")
  sims2.clutch<- readRDS("./RDS/sims2.clutch.RDS")
  
  # sims2.gen_length <-
  #   sim_plotter(
  #     janustree = simmap.janus.nuc.aggregate.simplified2,
  #     ouwiefit = gen_length2.OUM,
  #     nsim = 10000,
  #     clrs = setNames(rainbow(n = 11), unique(
  #       getStates(simmap.janus.nuc.aggregate.simplified2, type = 'tips')
  #     )),
  #     width = 0.05,
  #     plot = F
  #   ) %>%
  #   tidyr::pivot_longer(!Reg) %>%
  #   mutate(Reg = factor(Reg, levels = order2)) %>%
  #   arrange(Reg) %>%
  #   mutate(type=rep("expected", length(.[,1])))
  # saveRDS(sims2.gen_length, "./RDS/sims2.gen_length.RDS")
  sims2.gen_length <- readRDS("./RDS/sims2.gen_length.RDS")
  
  
  # sims2.survival <-
  #   sim_plotter(
  #     janustree = simmap.janus.nuc.aggregate.simplified2,
  #     ouwiefit = survival2.OUM,
  #     nsim = 10000,
  #     clrs = setNames(rainbow(n = 11), unique(
  #       getStates(simmap.janus.nuc.aggregate.simplified2, type = 'tips')
  #     )),
  #     width = 0.05,
  #     plot = F
  #   ) %>%
  #   tidyr::pivot_longer(!Reg) %>%
  #   mutate(Reg = factor(Reg, levels = order2)) %>%
  #   arrange(Reg) %>%
  #   mutate(type=rep("expected", length(.[,1])))
  # saveRDS(sims2.survival, "./RDS/sims2.survival.RDS")
  sims2.survival<- readRDS("./RDS/sims2.survival.RDS")
  
  # sims2.breeding <-
  #   sim_plotter(
  #     janustree = simmap.janus.nuc.aggregate.simplified2,
  #     ouwiefit = breeding2.OUM,
  #     nsim = 10000,
  #     clrs = setNames(rainbow(n = 11), unique(
  #       getStates(simmap.janus.nuc.aggregate.simplified2, type = 'tips')
  #     )),
  #     width = 0.05,
  #     plot = F
  #   ) %>%
  #   tidyr::pivot_longer(!Reg) %>%
  #   mutate(Reg = factor(Reg, levels = order2)) %>%
  #   arrange(Reg) %>%
  #   mutate(type=rep("expected", length(.[,1])))
  # saveRDS(sims2.breeding, "./RDS/sims2.breeding.RDS")
  sims2.breeding<- readRDS('./RDS/sims2.breeding.RDS')
  
  
  # sims2.longevity <-
  #   sim_plotter(
  #     janustree = simmap.janus.nuc.aggregate.simplified2,
  #     ouwiefit = longevity2.OUM,
  #     nsim = 10000,
  #     clrs = setNames(rainbow(n = 11), unique(
  #       getStates(simmap.janus.nuc.aggregate.simplified2, type = 'tips')
  #     )),
  #     width = 0.05,
  #     plot = F
  #   ) %>%
  #   tidyr::pivot_longer(!Reg) %>%
  #   mutate(Reg = factor(Reg, levels = order2)) %>%
  #   arrange(Reg) %>%
  #   mutate(type=rep("expected", length(.[,1])))
  # saveRDS(sims2.longevity, "./RDS/sims2.longevity.RDS")
  sims2.longevity<- readRDS("./RDS/sims2.longevity.RDS")
  
  # sims2.latitude <-
  #   sim_plotter(
  #     janustree = simmap.janus.nuc.aggregate.simplified2,
  #     ouwiefit = latitude2.OUM,
  #     nsim = 10000,
  #     clrs = setNames(rainbow(n = 11), unique(
  #       getStates(simmap.janus.nuc.aggregate.simplified2, type = 'tips')
  #     )),
  #     width = 0.05,
  #     plot = F
  #   ) %>%
  #   tidyr::pivot_longer(!Reg) %>%
  #   mutate(Reg = factor(Reg, levels = order2)) %>%
  #   arrange(Reg) %>%
  #   mutate(type=rep("expected", length(.[,1])))
  # saveRDS(sims2.latitude, "./RDS/sims2.latitude.RDS")
  sims2.latitude <- readRDS("./RDS/sims2.latitude.RDS")
  
  # sims.encprime <-
  #   sim_plotter(
  #     janustree = simmap.janus.nuc.exons,
  #     ouwiefit = exon.encprime.OUM,
  #     nsim = 1000,
  #     clrs = setNames(rainbow(n = 8), unique(
  #       getStates(simmap.janus.nuc.aggregate.simplified, type = 'tips')
  #     )),
  #     width = 0.05,
  #     plot = T
  #   ) %>%
  #   tidyr::pivot_longer(!Reg) %>%
  #   mutate(Reg = factor(Reg, levels = order)) %>%
  #   arrange(Reg) %>%
  #   mutate(type=rep("theta", length(.[,1])))
  
  
}

#more plotting setup
{
  Plot_data2.mass <- rbind(
    OUwie_data2[, c(2, 3)] %>%
      tidyr::pivot_longer(!Reg) %>%
      mutate(Reg = factor(Reg, levels = order2)) %>%
      arrange(Reg) %>%
      mutate(type = rep("data", length(.[, 1]))),
    sims2.mass, OUwie_boots2.mass)
  
  #Plot_data.mass$value <- Plot_data.mass$value * attributes(scale(log(megaLHT$mass)))[[3]] #unscale
  #Plot_data.mass$value <- Plot_data.mass$value + attributes(scale(log(megaLHT$mass)))[[2]] #uncenter
  #masslabs <- ((c(-2, 0, 2, 4)) * attributes(scale(log(megaLHT$mass)))[[3]]) + attributes(scale(log(megaLHT$mass)))[[2]]
  #max(Plot_data.mass[Plot_data.mass$type=='data',]$value)
  
  #colors for regimes
  #regimes<-c("black", "#ddcc77", "#aa4499","#6699cc", "#117733", "#999933", "#332288", "#88ccee", "#888888", "#882255")
  regimes2<-c("black", "#ddcc77", "#aa4499",'#661100',"#6699cc", "#117733", "#999933", "#332288", "#88ccee", "#888888", "#882255")
  #second row is alt test with 11 groups
  
  Plot_data2.chickPC1 <- rbind(
    OUwie_data2[, c(2, 4)] %>%
      tidyr::pivot_longer(!Reg) %>%
      mutate(Reg = factor(Reg, levels = order2)) %>%
      arrange(Reg) %>%
      mutate(type = rep("data", length(.[, 1]))),
    sims2.chickPC1, OUwie_boots2.chickPC1)
  
  Plot_data2.clutch <- rbind(
    OUwie_data.extra2[, c(2, 3)] %>%
      tidyr::pivot_longer(!Reg) %>%
      mutate(Reg = factor(Reg, levels = order2)) %>%
      arrange(Reg) %>%
      mutate(type = rep("data", length(.[, 1]))),
    sims2.clutch, OUwie_boots2.clutch)
  
  Plot_data2.gen_length <- rbind(
    OUwie_data.extra2[, c(2, 4)] %>%
      tidyr::pivot_longer(!Reg) %>%
      mutate(Reg = factor(Reg, levels = order2)) %>%
      arrange(Reg) %>%
      mutate(type = rep("data", length(.[, 1]))),
    sims2.gen_length, OUwie_boots2.gen_length)
  
  Plot_data2.survival <- rbind(
    OUwie_data.extra2[, c(2, 5)] %>%
      tidyr::pivot_longer(!Reg) %>%
      mutate(Reg = factor(Reg, levels = order2)) %>%
      arrange(Reg) %>%
      mutate(type = rep("data", length(.[, 1]))),
    sims2.survival, OUwie_boots2.survival)
  
  Plot_data2.breeding <- rbind(
    OUwie_data.extra2[, c(2, 6)] %>%
      tidyr::pivot_longer(!Reg) %>%
      mutate(Reg = factor(Reg, levels = order2)) %>%
      arrange(Reg) %>%
      mutate(type = rep("data", length(.[, 1]))),
    sims2.breeding, OUwie_boots2.breeding)
  
  Plot_data2.longevity <- rbind(
    OUwie_data.extra2[, c(2, 7)] %>%
      tidyr::pivot_longer(!Reg) %>%
      mutate(Reg = factor(Reg, levels = order2)) %>%
      arrange(Reg) %>%
      mutate(type = rep("data", length(.[, 1]))),
    sims2.longevity, OUwie_boots2.longevity)
  
  Plot_data2.latitude <- rbind(
    OUwie_data.extra2[, c(2, 8)] %>%
      tidyr::pivot_longer(!Reg) %>%
      mutate(Reg = factor(Reg, levels = order2)) %>%
      arrange(Reg) %>%
      mutate(type = rep("data", length(.[, 1]))),
    sims2.latitude, OUwie_boots2.latitude)
  
}

#generate ggplot objects
{
  massplot2 <-
    {ggplot(Plot_data2.mass,
            aes(x = value,
                y = Reg,
                fill = Reg,
                color= Reg,
                stroke=0.01, height=..ndensity..)) +
        #height = ..density..)) +
        #plotting theta
        geom_density_ridges(color="grey",
                            quantile_lines = F,
                            quantile_fun = estimate_mode,
                            data = dplyr::filter(Plot_data2.mass, type == "expected"),
                            #stat = "density",
                            scale = 0.85,
                            size = 0.01,
                            lty = 1,
                            rel_min_height = 0.001,
                            fill = c("grey"),
                            #viridis(2)[1],#"#fc8d62",
                            alpha = 0.25,
                            #0.4,
                            bandwidth = 0.35
        ) + #scale_fill_manual(values = regimes, labels = NULL) +
        #plotting expected
        geom_density_ridges(
          quantile_lines = TRUE,
          quantile_fun = estimate_mode,
          data = dplyr::filter(Plot_data2.mass, type == "theta"),
          #stat = "density",
          scale = 0.9,
          size = 0.75,
          rel_min_height = 0.001,
          #fill = c("red"),#viridis(2)[1],#"#fc8d62",
          alpha = 0.6,
          #0.4,
          bandwidth = 0.7
        ) + scale_fill_manual(values = regimes2, labels = NULL) + scale_color_manual(values=regimes2) + 
        #plotting data
        # geom_density_ridges(
        #   quantile_lines = F,
        #   quantile_fun = estimate_mode,
        #   data = dplyr::filter(Plot_data.mass, type == "data"),
        #   #stat = "density",
        #   scale = 0.9,
        #   size = 0.5,
        #   rel_min_height = 0.01,
        #   fill = "white",
        #   #viridis(2)[2],#"#8da0cb",
      #   alpha = 0.75,
      #   bandwidth = 0.65
      # ) +
      geom_density_ridges(size=0, color="black", fill="white",
                          data = dplyr::filter(Plot_data2.mass, type == "data"),
                          jittered_points = TRUE,
                          position = position_raincloud(width=0.75, height=0.05, ygap = -0.025, seed=5), scale=0, point_alpha=0.5, point_size=1.5, point_shape=21
                          #position = position_points_jitter(width = 0.1, height = 0, yoffset = -0.025),
                          #point_shape = '|', point_size = 1, point_alpha = 0.75, alpha = 0, scale=0, bandwidth = 0.5,
      ) + 
        theme_ridges() +
        theme(
          #axis.text.x=element_blank(), #remove x axis labels
          #axis.ticks.x=element_blank(), #remove x axis ticks
          axis.text.y = element_blank(),
          #remove y axis labels
          #axis.ticks.y=element_blank()  #remove y axis ticks
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size=8.5),
          #axis.text.x = element_text(vjust = 10),
          legend.position = "none"
        ) + grids(axis="x", linetype = "dotted") + grids(axis="y", linetype = "solid", size=1.25) +
        scale_x_continuous(breaks = c(2.5, 5, 7.5, 10, 12.5), limits = c(0.5, 13)) +
        scale_y_discrete(expand = c(0.015, 0))
    }
  
  chickpc1plot2 <-
    {ggplot(Plot_data2.chickPC1,
            aes(x = value,
                y = Reg,
                fill = Reg, 
                color= Reg,
                stroke=0.01, height=..ndensity..)) +
        #height = ..density..)) +
        #plotting theta
        geom_density_ridges(color="grey",
                            quantile_lines = F,
                            quantile_fun = estimate_mode,
                            data = dplyr::filter(Plot_data2.chickPC1, type == "expected"),
                            #stat = "density",
                            scale = 0.9,
                            size = 0.01,
                            lty = 1,
                            rel_min_height = 0.001,
                            fill = c("grey"),
                            #viridis(2)[1],#"#fc8d62",
                            alpha = 0.25,
                            #0.4,
                            bandwidth = 0.35
        ) + #scale_fill_manual(values = regimes, labels = NULL) +
        #plotting expected
        geom_density_ridges(
          quantile_lines = TRUE,
          quantile_fun = estimate_mode,
          data = dplyr::filter(Plot_data2.chickPC1, type == "theta"),
          #stat = "density",
          scale = 0.9,
          size = 0.75,
          rel_min_height = 0.001,
          #fill = c("red"),#viridis(2)[1],#"#fc8d62",
          alpha = 0.6,
          #0.4,
          bandwidth = 0.7
        ) + scale_fill_manual(values = regimes2, labels = NULL) + scale_color_manual(values=regimes2) + 
        #plotting data
        # geom_density_ridges(
        #   quantile_lines = F,
        #   quantile_fun = estimate_mode,
        #   data = dplyr::filter(Plot_data.chickPC1, type == "data"),
        #   #stat = "density",
        #   scale = 0.9,
        #   size = 0.5,
        #   rel_min_height = 0.01,
        #   fill = "white",
        #   #viridis(2)[2],#"#8da0cb",
      #   alpha = 0.75,
      #   bandwidth = 0.65
      # ) +
      geom_density_ridges(size=0, color="black", fill="white",
                          data = dplyr::filter(Plot_data2.chickPC1, type == "data"),
                          jittered_points = TRUE,
                          position = position_raincloud(width=0.75, height=0.05, ygap = -0.025, seed=1), scale=0, point_alpha=0.5, point_size=1.5, point_shape=21,
                          #position = position_points_jitter(width = 0.1, height = 0, yoffset = -0.025),
                          #point_shape = '|', point_size = 1, point_alpha = 0.75, alpha = 0, scale=0, bandwidth = 0.5
      ) + #scale_fill_manual(values = regimes, labels = NULL) #+ scale_color_manual(values="black")
        theme_ridges() +
        theme(
          #axis.text.x=element_blank(), #remove x axis labels
          #axis.ticks.x=element_blank(), #remove x axis ticks
          #axis.text.y=element_blank(),  #remove y axis labels
          #axis.ticks.y=element_blank()  #remove y axis ticks
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size=8.5),
          legend.position = "none"
        ) + grids(axis="x", linetype = "dotted") + grids(axis="y", linetype = "solid", size=1.25) +
        scale_x_continuous(breaks = c(-4, -2, 0, 2, 4), limits = c(-5, 5)) +
        scale_y_discrete(expand = c(0.015, 0))
    }
  
  clutchplot2 <-
    {ggplot(Plot_data2.clutch,
            aes(x = value,
                y = Reg,
                fill = Reg, 
                color= Reg,
                stroke=0.01, height=..ndensity..)) +
        #height = ..density..)) +
        #plotting theta
        geom_density_ridges(color="grey",
                            quantile_lines = F,
                            quantile_fun = estimate_mode,
                            data = dplyr::filter(Plot_data2.clutch, type == "expected"),
                            #stat = "density",
                            scale = 0.9,
                            size = 0.01,
                            lty = 1,
                            rel_min_height = 0.001,
                            fill = c("grey"),
                            #viridis(2)[1],#"#fc8d62",
                            alpha = 0.25,
                            #0.4,
                            bandwidth = 0.35
        ) + #scale_fill_manual(values = regimes, labels = NULL) +
        #plotting expected
        geom_density_ridges(
          quantile_lines = TRUE,
          quantile_fun = estimate_mode,
          data = dplyr::filter(Plot_data2.clutch, type == "theta"),
          #stat = "density",
          scale = 0.9,
          size = 0.75,
          rel_min_height = 0.001,
          #fill = c("red"),#viridis(2)[1],#"#fc8d62",
          alpha = 0.6,
          #0.4,
          bandwidth = 0.7
        ) + scale_fill_manual(values = regimes2, labels = NULL) + scale_color_manual(values=regimes2) + 
        #plotting data
        # geom_density_ridges(
        #   quantile_lines = F,
        #   quantile_fun = estimate_mode,
        #   data = dplyr::filter(Plot_data.chickPC1, type == "data"),
        #   #stat = "density",
        #   scale = 0.9,
        #   size = 0.5,
        #   rel_min_height = 0.01,
        #   fill = "white",
        #   #viridis(2)[2],#"#8da0cb",
      #   alpha = 0.75,
      #   bandwidth = 0.65
      # ) +
      geom_density_ridges(size=0, color="black", fill="white",
                          data = dplyr::filter(Plot_data2.clutch, type == "data"),
                          jittered_points = TRUE,
                          position = position_raincloud(width=0.75, height=0.05, ygap = -0.025, seed=1), scale=0, point_alpha=0.5, point_size=1.5, point_shape=21,
                          #position = position_points_jitter(width = 0.1, height = 0, yoffset = -0.025),
                          #point_shape = '|', point_size = 1, point_alpha = 0.75, alpha = 0, scale=0, bandwidth = 0.5
      ) + #scale_fill_manual(values = regimes, labels = NULL) #+ scale_color_manual(values="black")
        theme_ridges() +
        theme(
          #axis.text.x=element_blank(), #remove x axis labels
          #axis.ticks.x=element_blank(), #remove x axis ticks
          axis.text.y = element_blank(),
          #remove y axis labels
          #axis.ticks.y=element_blank()  #remove y axis ticks
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size=8.5),
          #axis.text.x = element_text(vjust = 10),
          legend.position = "none"
        ) + grids(axis="x", linetype = "dotted") + grids(axis="y", linetype = "solid", size=1.25) +
        scale_x_continuous(breaks = c(-4, -2, 0, 2, 4), limits = c(-5, 5)) +
        scale_y_discrete(expand = c(0.015, 0))
    }
  
  gen_length_plot2 <-
    {ggplot(Plot_data2.gen_length,
            aes(x = value,
                y = Reg,
                fill = Reg, 
                color= Reg,
                stroke=0.01, height=..ndensity..)) +
        #height = ..density..)) +
        #plotting theta
        geom_density_ridges(color="grey",
                            quantile_lines = F,
                            quantile_fun = estimate_mode,
                            data = dplyr::filter(Plot_data2.gen_length, type == "expected"),
                            #stat = "density",
                            scale = 0.9,
                            size = 0.01,
                            lty = 1,
                            rel_min_height = 0.001,
                            fill = c("grey"),
                            #viridis(2)[1],#"#fc8d62",
                            alpha = 0.25,
                            #0.4,
                            bandwidth = 0.35
        ) + #scale_fill_manual(values = regimes, labels = NULL) +
        #plotting expected
        geom_density_ridges(
          quantile_lines = TRUE,
          quantile_fun = estimate_mode,
          data = dplyr::filter(Plot_data2.gen_length, type == "theta"),
          #stat = "density",
          scale = 0.9,
          size = 0.75,
          rel_min_height = 0.001,
          #fill = c("red"),#viridis(2)[1],#"#fc8d62",
          alpha = 0.6,
          #0.4,
          bandwidth = 0.7
        ) + scale_fill_manual(values = regimes2, labels = NULL) + scale_color_manual(values=regimes2) + 
        #plotting data
        # geom_density_ridges(
        #   quantile_lines = F,
        #   quantile_fun = estimate_mode,
        #   data = dplyr::filter(Plot_data.chickPC1, type == "data"),
        #   #stat = "density",
        #   scale = 0.9,
        #   size = 0.5,
        #   rel_min_height = 0.01,
        #   fill = "white",
        #   #viridis(2)[2],#"#8da0cb",
      #   alpha = 0.75,
      #   bandwidth = 0.65
      # ) +
      geom_density_ridges(size=0, color="black", fill="white",
                          data = dplyr::filter(Plot_data2.gen_length, type == "data"),
                          jittered_points = TRUE,
                          position = position_raincloud(width=0.75, height=0.05, ygap = -0.025, seed=1), scale=0, point_alpha=0.5, point_size=1.5, point_shape=21,
                          #position = position_points_jitter(width = 0.1, height = 0, yoffset = -0.025),
                          #point_shape = '|', point_size = 1, point_alpha = 0.75, alpha = 0, scale=0, bandwidth = 0.5
      ) + #scale_fill_manual(values = regimes, labels = NULL) #+ scale_color_manual(values="black")
        theme_ridges() +
        theme(
          #axis.text.x=element_blank(), #remove x axis labels
          #axis.ticks.x=element_blank(), #remove x axis ticks
          axis.text.y = element_blank(),
          #remove y axis labels
          #axis.ticks.y=element_blank()  #remove y axis ticks
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size=8.5),
          #axis.text.x = element_text(vjust = 10),
          legend.position = "none"
        ) + grids(axis="x", linetype = "dotted") + grids(axis="y", linetype = "solid", size=1.25) +
        scale_x_continuous(breaks = c(-4, -2, 0, 2, 4), limits = c(-5, 5)) +
        scale_y_discrete(expand = c(0.015, 0))
    }
  
  survival_plot2 <-
    {ggplot(Plot_data2.survival,
            aes(x = value,
                y = Reg,
                fill = Reg, 
                color= Reg,
                stroke=0.01, height=..ndensity..)) +
        #height = ..density..)) +
        #plotting theta
        geom_density_ridges(color="grey",
                            quantile_lines = F,
                            quantile_fun = estimate_mode,
                            data = dplyr::filter(Plot_data2.survival, type == "expected"),
                            #stat = "density",
                            scale = 0.9,
                            size = 0.01,
                            lty = 1,
                            rel_min_height = 0.001,
                            fill = c("grey"),
                            #viridis(2)[1],#"#fc8d62",
                            alpha = 0.25,
                            #0.4,
                            bandwidth = 0.35
        ) + #scale_fill_manual(values = regimes, labels = NULL) +
        #plotting expected
        geom_density_ridges(
          quantile_lines = TRUE,
          quantile_fun = estimate_mode,
          data = dplyr::filter(Plot_data2.survival, type == "theta"),
          #stat = "density",
          scale = 0.9,
          size = 0.75,
          rel_min_height = 0.001,
          #fill = c("red"),#viridis(2)[1],#"#fc8d62",
          alpha = 0.6,
          #0.4,
          bandwidth = 0.7
        ) + scale_fill_manual(values = regimes2, labels = NULL) + scale_color_manual(values=regimes2) + 
        #plotting data
        # geom_density_ridges(
        #   quantile_lines = F,
        #   quantile_fun = estimate_mode,
        #   data = dplyr::filter(Plot_data.chickPC1, type == "data"),
        #   #stat = "density",
        #   scale = 0.9,
        #   size = 0.5,
        #   rel_min_height = 0.01,
        #   fill = "white",
        #   #viridis(2)[2],#"#8da0cb",
      #   alpha = 0.75,
      #   bandwidth = 0.65
      # ) +
      geom_density_ridges(size=0, color="black", fill="white",
                          data = dplyr::filter(Plot_data2.survival, type == "data"),
                          jittered_points = TRUE,
                          position = position_raincloud(width=0.75, height=0.05, ygap = -0.025, seed=1), scale=0, point_alpha=0.5, point_size=1.5, point_shape=21,
                          #position = position_points_jitter(width = 0.1, height = 0, yoffset = -0.025),
                          #point_shape = '|', point_size = 1, point_alpha = 0.75, alpha = 0, scale=0, bandwidth = 0.5
      ) + #scale_fill_manual(values = regimes, labels = NULL) #+ scale_color_manual(values="black")
        theme_ridges() +
        theme(
          #axis.text.x=element_blank(), #remove x axis labels
          #axis.ticks.x=element_blank(), #remove x axis ticks
          axis.text.y = element_blank(),
          #remove y axis labels
          #axis.ticks.y=element_blank()  #remove y axis ticks
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size=8.5),
          #axis.text.x = element_text(vjust = 10),
          legend.position = "none"
        ) + grids(axis="x", linetype = "dotted") + grids(axis="y", linetype = "solid", size=1.25) +
        scale_x_continuous(breaks = c(-4, -2, 0, 2, 4), limits = c(-5, 5)) +
        scale_y_discrete(expand = c(0.015, 0))
    }
  
  breeding_plot2 <-
    {ggplot(Plot_data2.breeding,
            aes(x = value,
                y = Reg,
                fill = Reg, 
                color= Reg,
                stroke=0.01, height=..ndensity..)) +
        #height = ..density..)) +
        #plotting theta
        geom_density_ridges(color="grey",
                            quantile_lines = F,
                            quantile_fun = estimate_mode,
                            data = dplyr::filter(Plot_data2.breeding, type == "expected"),
                            #stat = "density",
                            scale = 0.9,
                            size = 0.01,
                            lty = 1,
                            rel_min_height = 0.001,
                            fill = c("grey"),
                            #viridis(2)[1],#"#fc8d62",
                            alpha = 0.25,
                            #0.4,
                            bandwidth = 0.35
        ) + #scale_fill_manual(values = regimes, labels = NULL) +
        #plotting expected
        geom_density_ridges(
          quantile_lines = TRUE,
          quantile_fun = estimate_mode,
          data = dplyr::filter(Plot_data2.breeding, type == "theta"),
          #stat = "density",
          scale = 0.9,
          size = 0.75,
          rel_min_height = 0.001,
          #fill = c("red"),#viridis(2)[1],#"#fc8d62",
          alpha = 0.6,
          #0.4,
          bandwidth = 0.7
        ) + scale_fill_manual(values = regimes2, labels = NULL) + scale_color_manual(values=regimes2) + 
        #plotting data
        # geom_density_ridges(
        #   quantile_lines = F,
        #   quantile_fun = estimate_mode,
        #   data = dplyr::filter(Plot_data.chickPC1, type == "data"),
        #   #stat = "density",
        #   scale = 0.9,
        #   size = 0.5,
        #   rel_min_height = 0.01,
        #   fill = "white",
        #   #viridis(2)[2],#"#8da0cb",
      #   alpha = 0.75,
      #   bandwidth = 0.65
      # ) +
      geom_density_ridges(size=0, color="black", fill="white",
                          data = dplyr::filter(Plot_data2.breeding, type == "data"),
                          jittered_points = TRUE,
                          position = position_raincloud(width=0.75, height=0.05, ygap = -0.025, seed=1), scale=0, point_alpha=0.5, point_size=1.5, point_shape=21,
                          #position = position_points_jitter(width = 0.1, height = 0, yoffset = -0.025),
                          #point_shape = '|', point_size = 1, point_alpha = 0.75, alpha = 0, scale=0, bandwidth = 0.5
      ) + #scale_fill_manual(values = regimes, labels = NULL) #+ scale_color_manual(values="black")
        theme_ridges() +
        theme(
          #axis.text.x=element_blank(), #remove x axis labels
          #axis.ticks.x=element_blank(), #remove x axis ticks
          axis.text.y = element_blank(),
          #remove y axis labels
          #axis.ticks.y=element_blank()  #remove y axis ticks
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size=8.5),
          #axis.text.x = element_text(vjust = 10),
          legend.position = "none"
        ) + grids(axis="x", linetype = "dotted") + grids(axis="y", linetype = "solid", size=1.25) +
        scale_x_continuous(breaks = c(-4, -2, 0, 2, 4), limits = c(-5, 5)) +
        scale_y_discrete(expand = c(0.015, 0))
    }
  
  longevity_plot2 <-
    {ggplot(Plot_data2.longevity,
            aes(x = value,
                y = Reg,
                fill = Reg, 
                color= Reg,
                stroke=0.01, height=..ndensity..)) +
        #height = ..density..)) +
        #plotting theta
        geom_density_ridges(color="grey",
                            quantile_lines = F,
                            quantile_fun = estimate_mode,
                            data = dplyr::filter(Plot_data2.longevity, type == "expected"),
                            #stat = "density",
                            scale = 0.9,
                            size = 0.01,
                            lty = 1,
                            rel_min_height = 0.001,
                            fill = c("grey"),
                            #viridis(2)[1],#"#fc8d62",
                            alpha = 0.25,
                            #0.4,
                            bandwidth = 0.35
        ) + #scale_fill_manual(values = regimes, labels = NULL) +
        #plotting expected
        geom_density_ridges(
          quantile_lines = TRUE,
          quantile_fun = estimate_mode,
          data = dplyr::filter(Plot_data2.longevity, type == "theta"),
          #stat = "density",
          scale = 0.9,
          size = 0.75,
          rel_min_height = 0.001,
          #fill = c("red"),#viridis(2)[1],#"#fc8d62",
          alpha = 0.6,
          #0.4,
          bandwidth = 0.7
        ) + scale_fill_manual(values = regimes2, labels = NULL) + scale_color_manual(values=regimes2) + 
        #plotting data
        # geom_density_ridges(
        #   quantile_lines = F,
        #   quantile_fun = estimate_mode,
        #   data = dplyr::filter(Plot_data.chickPC1, type == "data"),
        #   #stat = "density",
        #   scale = 0.9,
        #   size = 0.5,
        #   rel_min_height = 0.01,
        #   fill = "white",
        #   #viridis(2)[2],#"#8da0cb",
      #   alpha = 0.75,
      #   bandwidth = 0.65
      # ) +
      geom_density_ridges(size=0, color="black", fill="white",
                          data = dplyr::filter(Plot_data2.longevity, type == "data"),
                          jittered_points = TRUE,
                          position = position_raincloud(width=0.75, height=0.05, ygap = -0.025, seed=1), scale=0, point_alpha=0.5, point_size=1.5, point_shape=21,
                          #position = position_points_jitter(width = 0.1, height = 0, yoffset = -0.025),
                          #point_shape = '|', point_size = 1, point_alpha = 0.75, alpha = 0, scale=0, bandwidth = 0.5
      ) + #scale_fill_manual(values = regimes, labels = NULL) #+ scale_color_manual(values="black")
        theme_ridges() +
        theme(
          #axis.text.x=element_blank(), #remove x axis labels
          #axis.ticks.x=element_blank(), #remove x axis ticks
          axis.text.y = element_blank(),
          #remove y axis labels
          #axis.ticks.y=element_blank()  #remove y axis ticks
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size=8.5),
          #axis.text.x = element_text(vjust = 10),
          legend.position = "none"
        ) + grids(axis="x", linetype = "dotted") + grids(axis="y", linetype = "solid", size=1.25) +
        scale_x_continuous(breaks = c(-4, -2, 0, 2, 4), limits = c(-5, 5)) +
        scale_y_discrete(expand = c(0.015, 0))
    }
  
  latitude_plot2 <-
    {ggplot(Plot_data2.latitude,
            aes(x = value,
                y = Reg,
                fill = Reg, 
                color= Reg,
                stroke=0.01, height=..ndensity..)) +
        #height = ..density..)) +
        #plotting theta
        geom_density_ridges(color="grey",
                            quantile_lines = F,
                            quantile_fun = estimate_mode,
                            data = dplyr::filter(Plot_data2.latitude, type == "expected"),
                            #stat = "density",
                            scale = 0.9,
                            size = 0.01,
                            lty = 1,
                            rel_min_height = 0.001,
                            fill = c("grey"),
                            #viridis(2)[1],#"#fc8d62",
                            alpha = 0.25,
                            #0.4,
                            bandwidth = 0.35
        ) + #scale_fill_manual(values = regimes, labels = NULL) +
        #plotting expected
        geom_density_ridges(
          quantile_lines = TRUE,
          quantile_fun = estimate_mode,
          data = dplyr::filter(Plot_data2.latitude, type == "theta"),
          #stat = "density",
          scale = 0.9,
          size = 0.75,
          rel_min_height = 0.001,
          #fill = c("red"),#viridis(2)[1],#"#fc8d62",
          alpha = 0.6,
          #0.4,
          bandwidth = 0.7
        ) + scale_fill_manual(values = regimes2, labels = NULL) + scale_color_manual(values=regimes2) + 
        #plotting data
        # geom_density_ridges(
        #   quantile_lines = F,
        #   quantile_fun = estimate_mode,
        #   data = dplyr::filter(Plot_data.chickPC1, type == "data"),
        #   #stat = "density",
        #   scale = 0.9,
        #   size = 0.5,
        #   rel_min_height = 0.01,
        #   fill = "white",
        #   #viridis(2)[2],#"#8da0cb",
      #   alpha = 0.75,
      #   bandwidth = 0.65
      # ) +
      geom_density_ridges(size=0, color="black", fill="white",
                          data = dplyr::filter(Plot_data2.latitude, type == "data"),
                          jittered_points = TRUE,
                          position = position_raincloud(width=0.75, height=0.05, ygap = -0.025, seed=1), scale=0, point_alpha=0.5, point_size=1.5, point_shape=21,
                          #position = position_points_jitter(width = 0.1, height = 0, yoffset = -0.025),
                          #point_shape = '|', point_size = 1, point_alpha = 0.75, alpha = 0, scale=0, bandwidth = 0.5
      ) + #scale_fill_manual(values = regimes, labels = NULL) #+ scale_color_manual(values="black")
        theme_ridges() +
        theme(
          #axis.text.x=element_blank(), #remove x axis labels
          #axis.ticks.x=element_blank(), #remove x axis ticks
          axis.text.y = element_blank(),
          #remove y axis labels
          #axis.ticks.y=element_blank()  #remove y axis ticks
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size=8.5),
          #axis.text.x = element_text(vjust = 10),
          legend.position = "none"
        ) + grids(axis="x", linetype = "dotted") + grids(axis="y", linetype = "solid", size=1.25) +
        scale_x_continuous(breaks = c(-4, -2, 0, 2, 4), limits = c(-5, 5)) +
        scale_y_discrete(expand = c(0.015, 0))
    }
  
}

#plotting
#require(ggpubr)
pdf(file="LH_shift_plots2.pdf", height=11.5, width=7)
ggarrange(chickpc1plot2, massplot2, labels=c("ChickPC1", "Mass"),
          widths=c(1.7,1.0))
dev.off()

#plotting supplemental shifts
pdf(file="LH_shift_plots_supp2.pdf", width=11*2, height=8.5*1.75)
ggarrange(chickpc1plot2, massplot2, clutchplot2, gen_length_plot2, survival_plot2, breeding_plot2,longevity_plot2, latitude_plot2,
          labels=c("ChickPC1", "Mass", "Clutch Size", "Generation Length", "Survival", "Age First Breeding", "Longevity", "Latitude"),
          widths=c(1.6,1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0), nrow=1)
dev.off()


#Section 11
#########################################################
#testing base compositional plots/janus parameter plots##
#processing equilibrium base composition from janus######
#set up datasets for raw data############################
#########################################################
{
nuc.exonACGT <- nuc.exonACGT [consensus.exons.data$label,]
nuc.exonACGT <- cbind(consensus.exons.data$exon_models, nuc.exonACGT)
colnames(nuc.exonACGT) <- c( "group", "A", "C", "G", "T")
#nuc.exonACGT[,1]<-seq(1:198)
nuc.exonACGT<-aggregate(nuc.exonACGT, by= list(nuc.exonACGT[,1]), FUN=mean)
nuc.exonACGT$Group.1<-NULL


nuc.intronACGT <- nuc.intronACGT [consensus.introns.data$label,]
nuc.intronACGT <- cbind(consensus.introns.data$intron_models, nuc.intronACGT)
colnames(nuc.intronACGT) <- c( "group", "A", "C", "G", "T")
#nuc.exonACGT[,1]<-seq(1:198)
nuc.intronACGT<-aggregate(nuc.intronACGT, by= list(nuc.intronACGT[,1]), FUN=mean)
nuc.intronACGT$Group.1<-NULL

nuc.utrACGT <- nuc.utrACGT [consensus.utrs.data$label,]
nuc.utrACGT <- cbind(consensus.utrs.data$utr_models, nuc.utrACGT)
colnames(nuc.utrACGT) <- c( "group", "A", "C", "G", "T")
#nuc.exonACGT[,1]<-seq(1:198)
nuc.utrACGT<-aggregate(nuc.utrACGT, by= list(nuc.utrACGT[,1]), FUN=mean)
nuc.utrACGT$Group.1<-NULL

mtdna.ACGT <- mtdna.ACGT [consensus.mtdnas.all.data$label,]
mtdna.ACGT <- cbind(consensus.mtdnas.all.data$mtdna_all_models, mtdna.ACGT)
colnames(mtdna.ACGT) <- c( "group", "A", "C", "G", "T")
#nuc.exonACGT[,1]<-seq(1:198)
mtdna.ACGT<-aggregate(mtdna.ACGT, by= list(mtdna.ACGT[,1]), FUN=mean)
mtdna.ACGT$Group.1<-NULL


mtdna.ACGT.eq<-baseFreq_plotSet(consensus.mtdnas.all.data)
nuc.exonACGT.eq<-baseFreq_plotSet(consensus.exons.data)
nuc.utrACGT.eq<-baseFreq_plotSet(consensus.utrs.data)
nuc.intronACGT.eq<-baseFreq_plotSet(consensus.introns.data)

}

#barplot(t(nuc.exonACGT.eq[,c(2:5)]), beside = F, col = rainbow(4))
#consensus.exons.data$exon_models
#tmp<-data.frame(group=as.character(c(1,2,3,4,5,6,7,8)), tmp)

#exon plots
{
nuc.exonACGT.data.plot<-ggradar(
  nuc.exonACGT,
  grid.min = 0.05,
  grid.mid = 0.25,
  grid.max = 0.45,
  #centre.y=0.2,
  group.line.width = 0.5, 
  group.point.size = 1,
  background.circle.colour = "white",
  gridline.mid.colour = "grey",
  legend.position = "left",
  values.radar = c("05%", "25%", "45%"), plot.legend=T
)

nuc.exonACGT.eq.plot<- ggradar(
  #rbind(tmp, nuc.exonACGT)[c(1,9)+0,],
  nuc.exonACGT.eq,
  grid.min = 0.05,
  grid.mid = 0.25,
  grid.max = 0.45,
  group.line.width = 0.5, 
  group.point.size = 1,
  background.circle.colour = "white",
  gridline.mid.colour = "grey",
  legend.position = "none",
  values.radar = c("05%", "25%", "45%")
)
}

#intron plots
{
  nuc.intronACGT.data.plot<-ggradar(
    nuc.intronACGT,
    grid.min = 0.05,
    grid.mid = 0.25,
    grid.max = 0.45,
    #centre.y=0.2,
    group.line.width = 0.5, 
    group.point.size = 1,
    background.circle.colour = "white",
    gridline.mid.colour = "grey",
    legend.position = "left",
    values.radar = c("05%", "25%", "45%"), plot.legend=T
  )
  
  
  
  nuc.intronACGT.eq.plot<- ggradar(
    #rbind(tmp, nuc.exonACGT)[c(1,9)+0,],
    nuc.intronACGT.eq,
    grid.min = 0.15,
    grid.mid = 0.25,
    grid.max = 0.35,
    #centre.y = 0.15,
    group.line.width = 0.5, 
    group.point.size = 1,
    background.circle.colour = "white",
    gridline.mid.colour = "grey",
    legend.position = "left",
    values.radar = c("15%", "25%", "35%")
  )
  
}

#utr plots
{
  nuc.utrACGT.data.plot<-ggradar(
    nuc.utrACGT,
    grid.min = 0.05,
    grid.mid = 0.25,
    grid.max = 0.45,
    #centre.y=0.2,
    group.line.width = 0.5, 
    group.point.size = 1,
    background.circle.colour = "white",
    gridline.mid.colour = "grey",
    legend.position = "left",
    values.radar = c("05%", "25%", "45%"), plot.legend=T
  )
  
  
  
  nuc.utrACGT.eq.plot<- ggradar(
    #rbind(tmp, nuc.exonACGT)[c(1,9)+0,],
    nuc.utrACGT.eq,
    grid.min = 0.05,
    grid.mid = 0.25,
    grid.max = 0.45,
    group.line.width = 0.5, 
    group.point.size = 1,
    background.circle.colour = "white",
    gridline.mid.colour = "grey",
    legend.position = "left",
    values.radar = c("05%", "25%", "45%")
  )
  
}

#mtdna plots
{
  mtdna.ACGT.data.plot<-ggradar(
    mtdna.ACGT,
    grid.min = 0.05,
    grid.mid = 0.25,
    grid.max = 0.45,
    #centre.y=0.2,
    group.line.width = 0.5, 
    group.point.size = 1,
    background.circle.colour = "white",
    gridline.mid.colour = "grey",
    legend.position = "left",
    values.radar = c("05%", "25%", "45%"), plot.legend=T
  )
  
  
  
  mtdna.ACGT.eq.plot<- ggradar(
    #rbind(tmp, nuc.exonACGT)[c(1,9)+0,],
    mtdna.ACGT.eq,
    grid.min = 0.05,
    grid.mid = 0.25,
    grid.max = 0.45,
    group.line.width = 0.5, 
    group.point.size = 1,
    background.circle.colour = "white",
    gridline.mid.colour = "grey",
    legend.position = "left",
    values.radar = c("05%", "25%", "45%")
  )
  
}

# this code uses the ggradar function but not using any more
# ggarrange(nuc.exonACGT.data.plot, 
#           nuc.exonACGT.eq.plot,
#           
#           nuc.intronACGT.data.plot,
#           nuc.intronACGT.eq.plot,
#           
#           nuc.utrACGT.data.plot,
#           nuc.utrACGT.eq.plot,
#           
#           mtdna.ACGT.data.plot,
#           mtdna.ACGT.eq.plot, 
#           ncol=2, nrow=4,
#           common.legend=F,
#           legend="left"
#           )
# 
# 
# ggarrange(
#           nuc.exonACGT.eq.plot,
#           nuc.intronACGT.eq.plot,
#           nuc.utrACGT.eq.plot,
#           mtdna.ACGT.eq.plot, 
#           ncol=4, nrow=1,
#           common.legend=T,
#           legend="bottom"
# )


#reference plots for checking colors
# pdf(file='exons_ref2.pdf',width=5, height=15)
# plot_janus_time_rect(reference = time_scale, target = consensus.exons)
# tiplabels(getStates(simmap.janus.nuc.exons, type = 'tips'), frame="none", cex=0.5)
# dev.off()
# 
# pdf(file='introns_ref2.pdf',width=5, height=15)
# plot_janus_time_rect(reference = time_scale, target = consensus.introns)
# tiplabels(getStates(simmap.janus.nuc.introns, type = 'tips'), frame="none", cex=0.5)
# dev.off()
# 
# pdf(file='utrs_ref2.pdf',width=5, height=15)
# plot_janus_time_rect(reference = time_scale, target = consensus.utrs)
# tiplabels(getStates(simmap.janus.nuc.utrs, type = 'tips'), frame="none", cex=0.5)
# dev.off()
# 
# pdf(file='mtdna_ref2.pdf',width=5, height=15)
# plot_janus_time_rect(reference = time_scale, target = consensus.mtdnas.all)
# tiplabels(getStates(simmap.janus.mtdna.all, type = 'tips'), frame="none", cex=0.5)
# dev.off()


#testing radarchart function instead

#exons
{
nuc.exonACGT.eq.rad<-nuc.exonACGT.eq[, 2:5]
rownames(nuc.exonACGT.eq.rad)<-nuc.exonACGT.eq[,1]
set.seed(seed=5)
radExonCols<-sample(brewer.pal(n=8, name="Paired"))
radarchart(
  rbind(rep(0.4, 4) , rep(0.2, 4), nuc.exonACGT.eq.rad),
  axistype = 1,
  plty = 1,
  plwd = 2,
  seg = 2,
  #custom the grid
  cglcol = "black",
  cglty = 3,
  axislabcol = "black",
  caxislabels = seq(0.2, 0.4, 0.1),
  cglwd = 0.8,
  #custom labels
  vlcex = 2,
  #pdensity=0.5,
  #palcex = 0.1,
  pcol = radExonCols
)

# Add a legend
legend(x=0.7, y=1.2, legend = rownames(nuc.exonACGT.eq.rad), bty = "n", pch=20 , text.col = "grey", cex=1.2, pt.cex=3, col=radExonCols)
}

#introns
{
  nuc.intronACGT.eq.rad<-nuc.intronACGT.eq[, 2:5]
  rownames(nuc.intronACGT.eq.rad)<-nuc.intronACGT.eq[,1]
  set.seed(seed=5)
  radIntronCols<-sample(brewer.pal(n=6, name="Paired"))
  radarchart(
    rbind(rep(0.35, 4) , rep(0.2, 4), nuc.intronACGT.eq.rad),
    axistype = 1,
    plty = 1,
    plwd = 2,
    seg = 2,
    #custom the grid
    cglcol = "black",
    cglty = 3,
    axislabcol = "black",
    caxislabels = seq(0.2, 0.35, length.out=3),
    cglwd = 0.8,
    #custom labels
    vlcex = 2,
    #pdensity=0.5,
    #palcex = 0.1,
    pcol = radIntronCols
  )
  
  # Add a legend
  legend(x=0.7, y=1, legend = rownames(nuc.intronACGT.eq.rad), bty = "n", pch=20 , text.col = "grey", cex=1.2, pt.cex=3, col=radIntronCols)
}

#utrs
{
  nuc.utrACGT.eq.rad<-nuc.utrACGT.eq[, 2:5]
  rownames(nuc.utrACGT.eq.rad)<-nuc.utrACGT.eq[,1]
  set.seed(seed=5)
  radUtrCols<-sample(brewer.pal(n=4, name="Paired"))
  radarchart(
    rbind(rep(0.3, 4) , rep(0.2, 4), nuc.utrACGT.eq.rad),
    axistype = 1,
    plty = 1,
    plwd = 2,
    seg = 2,
    #custom the grid
    cglcol = "black",
    cglty = 3,
    axislabcol = "black",
    caxislabels = seq(0.2, 0.3, length.out=3),
    cglwd = 0.8,
    #custom labels
    vlcex = 2,
    #pdensity=0.5,
    #palcex = 0.1,
    pcol = radUtrCols
  )
  
  # Add a legend
  legend(x=0.7, y=1, legend = rownames(nuc.utrACGT.eq.rad), bty = "n", pch=20 , text.col = "grey", cex=1.2, pt.cex=3, col=radUtrCols)
}

#mtdnas
{
  mtdna.ACGT.eq.rad<-mtdna.ACGT.eq[, 2:5]
  rownames(mtdna.ACGT.eq.rad)<-mtdna.ACGT.eq[,1]
  set.seed(seed=5)
  radMtdnaCols<-sample(brewer.pal(n=3, name="Paired"))
  radarchart(
    rbind(rep(0.4, 4) , rep(0.0, 4), mtdna.ACGT.eq.rad),
    axistype = 1,
    plty = 1,
    plwd = 2,
    seg = 4,
    #custom the grid
    cglcol = "black",
    cglty = 3,
    axislabcol = "black",
    caxislabels = seq(0.0, 0.4, length.out=5),
    cglwd = 0.8,
    #custom labels
    vlcex = 2,
    #pdensity=0.5,
    #palcex = 0.1,
    pcol = radMtdnaCols
  )
  
  # Add a legend
  legend(x=0.7, y=1, legend = rownames(mtdna.ACGT.eq.rad), bty = "n", pch=20 , text.col = "grey", cex=1.2, pt.cex=3, col=radMtdnaCols)
}

#plotting the simmap tree with radar plots
#for some reasons these generate an error when plotting interactively, but plot fine to PDF

#exon radar plot
{
  pdf(file='exons_radarplot.pdf',width=8.5, height=11)
  plot_janus_time_rect_radar(reference = time_scale, target = consensus.exons)
  #par(new=T, mar=c(0, 0, 20, 20))
  subplot(fun={
    nuc.exonACGT.eq.rad<-nuc.exonACGT.eq[, 2:5]
    rownames(nuc.exonACGT.eq.rad)<-nuc.exonACGT.eq[,1]
    set.seed(seed=5)
    radExonCols<-sample(brewer.pal(n=8, name="Paired"))
    radarchart(
      rbind(rep(0.4, 4) , rep(0.2, 4), nuc.exonACGT.eq.rad),
      axistype = 1,
      plty = 1,
      plwd = 1.5,
      seg = 2,
      #custom the grid
      cglcol = "black",
      cglty = 3,
      axislabcol = "black",
      caxislabels = seq(0.2, 0.4, 0.1),
      cglwd = 0.5,
      #custom labels
      vlcex = 1.5,
      calcex = 0.5,
      palcex = 0.5,
      #pdensity=0.5,
      #palcex = 0.1,
      pcol = radExonCols
    )
    
    # Add a legend
    legend(x=-1.2, y=2.9, legend = rownames(nuc.exonACGT.eq.rad), bty = "n", pch=20 , text.col = "grey", cex=2.0, pt.cex=5, col=radExonCols, xpd=T)
  }, x = -70, y= 80, size=c(4.75,4.75))
  #tiplabels(unlist(lapply(strsplit(consensus.exons@phylo$tip.label, split="_"), `[`, 1)), frame="none", bg="none", cex=0.5, adj=0)
  #tiplabels(getStates(simmap.janus.nuc.exons, type = 'tips'), frame="none", cex=0.5)
  text("Exon signal", x=-100, y=3, cex=3)
  mtext(text="Supplementary Figure 7a")
  dev.off()
}

#intron radar plot
{
  pdf(file='introns_radarplot.pdf',width=8.5, height=11)
  plot_janus_time_rect_radar(reference = time_scale, target = consensus.introns)
  #par(new=T, mar=c(0, 0, 20, 20))
  subplot(fun={
    nuc.intronACGT.eq.rad<-nuc.intronACGT.eq[, 2:5]
    rownames(nuc.intronACGT.eq.rad)<-nuc.intronACGT.eq[,1]
    set.seed(seed=5)
    radIntronCols<-sample(brewer.pal(n=6, name="Paired"))
    radarchart(
      rbind(rep(0.35, 4) , rep(0.2, 4), nuc.intronACGT.eq.rad),
      axistype = 1,
      plty = 1,
      plwd = 1.5,
      seg = 2,
      #custom the grid
      cglcol = "black",
      cglty = 3,
      axislabcol = "black",
      caxislabels = seq(0.2, 0.35, length.out=3),
      cglwd = 0.5,
      #custom labels
      vlcex = 1.5,
      calcex = 0.5,
      palcex = 0.5,
      #pdensity=0.5,
      #palcex = 0.1,
      pcol = radIntronCols
    )
    
    # Add a legend
    legend(x=-1.2, y=2.9, legend = rownames(nuc.intronACGT.eq.rad), bty = "n", pch=20 , text.col = "grey", cex=2.0, pt.cex=5, col=radIntronCols, xpd=T)
  }, x = -70, y= 80, size=c(4.75,4.75))
  #tiplabels(unlist(lapply(strsplit(consensus.exons@phylo$tip.label, split="_"), `[`, 1)), frame="none", bg="none", cex=0.5, adj=0)
  #tiplabels(getStates(simmap.janus.nuc.exons, type = 'tips'), frame="none", cex=0.5)
  text("Intron signal", x=-100, y=3, cex=3)
  mtext(text="Supplementary Figure 7b")
  dev.off()
}

#utr radar plot
{
  pdf(file='utrs_radarplot.pdf',width=8.5, height=11)
  plot_janus_time_rect_radar(reference = time_scale, target = consensus.utrs)
  #par(new=T, mar=c(0, 0, 20, 20))
  subplot(fun={
    nuc.utrACGT.eq.rad<-nuc.utrACGT.eq[, 2:5]
    rownames(nuc.utrACGT.eq.rad)<-nuc.utrACGT.eq[,1]
    set.seed(seed=5)
    radUtrCols<-sample(brewer.pal(n=4, name="Paired"))
    radarchart(
      rbind(rep(0.3, 4) , rep(0.2, 4), nuc.utrACGT.eq.rad),
      axistype = 1,
      plty = 1,
      plwd = 1.5,
      seg = 2,
      #custom the grid
      cglcol = "black",
      cglty = 3,
      axislabcol = "black",
      caxislabels = seq(0.2, 0.3, length.out=3),
      cglwd = 0.5,
      #custom labels
      vlcex = 1.5,
      calcex = 0.5,
      palcex = 0.5,
      #pdensity=0.5,
      #palcex = 0.1,
      pcol = radUtrCols
    )
    
    # Add a legend
    legend(x=-1.2, y=2.9, legend = rownames(nuc.utrACGT.eq.rad), bty = "n", pch=20 , text.col = "grey", cex=2.0, pt.cex=5, col=radUtrCols, xpd=T)
  }, x = -70, y= 80, size=c(4.75,4.75))
  #tiplabels(unlist(lapply(strsplit(consensus.exons@phylo$tip.label, split="_"), `[`, 1)), frame="none", bg="none", cex=0.5, adj=0)
  #tiplabels(getStates(simmap.janus.nuc.exons, type = 'tips'), frame="none", cex=0.5)
  text("UTR signal", x=-100, y=3, cex=3)
  mtext(text="Supplementary Figure 7c")
  dev.off()
}

#mtdna radar plot
{
  pdf(file='mtdnas_radarplot.pdf',width=8.5, height=11)
  plot_janus_time_rect_radar(reference = time_scale, target = consensus.mtdnas.all)
  #par(new=T, mar=c(0, 0, 20, 20))
  subplot(fun={
    mtdna.ACGT.eq.rad<-mtdna.ACGT.eq[, 2:5]
    rownames(mtdna.ACGT.eq.rad)<-mtdna.ACGT.eq[,1]
    set.seed(seed=5)
    radMtdnaCols<-sample(brewer.pal(n=3, name="Paired"))
    radarchart(
      rbind(rep(0.4, 4) , rep(0.0, 4), mtdna.ACGT.eq.rad),
      axistype = 1,
      plty = 1,
      plwd = 1.5,
      seg = 4,
      #custom the grid
      cglcol = "black",
      cglty = 3,
      axislabcol = "black",
      caxislabels = seq(0.0, 0.4, length.out=5),
      cglwd = 0.5,
      #custom labels
      vlcex = 1.5,
      calcex = 0.5,
      palcex = 0.5,
      #pdensity=0.5,
      #palcex = 0.1,
      pcol = radMtdnaCols
    )
    
    # Add a legend
    legend(x=-1.2, y=2.9, legend = rownames(mtdna.ACGT.eq.rad), bty = "n", pch=20 , text.col = "grey", cex=2.0, pt.cex=5, col=radMtdnaCols, xpd=T)
  }, x = -70, y= 80, size=c(4.75,4.75))
  #tiplabels(unlist(lapply(strsplit(consensus.exons@phylo$tip.label, split="_"), `[`, 1)), frame="none", bg="none", cex=0.5, adj=0)
  #tiplabels(getStates(simmap.janus.nuc.exons, type = 'tips'), frame="none", cex=0.5)
  text("mtDNA signal", x=-100, y=3, cex=3)
  mtext(text="Supplementary Figure 7d")
  dev.off()
}


#base comp shifts
pdf(file="bf_shifts.pdf", width=12,height=3.5)
layout(matrix(c(1,2,3,4), nrow=1), widths=c(1.15,1,0.7,0.60))
baseFreq_shiftplot(nuc.exonACGT[,c("group", "A", "T", "G", "C")], nuc.exonACGT.eq[,c("group", "A", "T", "G", "C")], title="Exon")
baseFreq_shiftplot(nuc.intronACGT[,c("group", "A", "T", "G", "C")], nuc.intronACGT.eq[,c("group", "A", "T", "G", "C")], title="Intron")
baseFreq_shiftplot(nuc.utrACGT[,c("group", "A", "T", "G", "C")], nuc.utrACGT.eq[,c("group", "A", "T", "G", "C")], title="UTR")
baseFreq_shiftplot(mtdna.ACGT[,c("group", "A", "T", "G", "C")], mtdna.ACGT.eq[,c("group", "A", "T", "G", "C")], title="mtDNA")
dev.off()


#Section 12
#########################################################
#setting up datasets for analysis of metabolic scaling ##
#########################################################
{
#metabolic scaling test
{#simmap.janus.nuc.alldata

# plot(log(megaLHT$BMR.species.Metabolic.rate.Ryan/megaLHT$BMR.species.Body.Mass) ~ log(megaLHT$BMR.species.Body.Mass))
# plot(log(megaLHT$BMR.species.Metabolic.rate.Ryan) ~ log(megaLHT$BMR.species.Body.Mass))

scaling_data <- cbind(mass=megaLHT$BMR.species.Body.Mass, bmr=megaLHT$BMR.species.Metabolic.rate.Ryan)
scaling_data <- log(scaling_data)
rownames(scaling_data)<- megaLHT$species
#scale
#scaling_data.scaled <- as.data.frame(lapply(as.data.frame(scaling_data), function(x) c(scale(x))))
#rownames(scaling_data.scaled)<- megaLHT$species
#scaling_data$species<-rownames(scaling_data)
#scaling_data<-scaling_data[,c(3,1,2)]

#errors
#scaling_data_errors<-cbind(mass=megaLHT$BMR.species.se.mass.est.10, bmr=megaLHT$BMR.species.se.bmr.est.10)
#scaling_data_errors<-log(scaling_data_errors)^2 #squared standard errors in log space

#taking the 'mean' sd for birds / sqrt(10) for mass
#taking the value from the allometry paper for bmr
scaling_data_errors<-cbind(mass=rep(0.055, length(rownames(scaling_data))), bmr=rep(0.1, length(rownames(scaling_data))))
scaling_data_errors<-(scaling_data_errors)^2 #squared standard errors in log space
#scaling_data_errors <- as.data.frame(lapply(as.data.frame(scaling_data_errors), function(x) c(scale(x))))
rownames(scaling_data_errors)<- megaLHT$species
}
#which regime map to use for imputation? test them mvBM
{
#janus.null.mvBM<-mvBM(tree=simmap.janus.null, data = scaling_data, error = scaling_data_errors, model="BM1", diagnostic=FALSE, echo=FALSE, method='rpf', scale.height = F)
#saveRDS(janus.null.mvBM, "./RDS/janus.null.mvBM.RDS")
janus.null.mvBM<-readRDS("./RDS/janus.null.mvBM.RDS")

#all.aggregte.mvBM<-mvBM(tree=simmap.janus.all.aggregate, data = scaling_data, error = scaling_data_errors, model="BMM", diagnostic=FALSE, echo=FALSE, method='rpf', scale.height = F)
#saveRDS(all.aggregte.mvBM, "./RDS/all.aggregte.mvBM.RDS")
all.aggregte.mvBM<-readRDS("./RDS/all.aggregte.mvBM.RDS")

#mtdna.all.mvBM<-mvBM(tree=simmap.janus.mtdna.all, data = scaling_data, error = scaling_data_errors, model="BMM", diagnostic=FALSE, echo=FALSE, method='rpf', scale.height = F)
#saveRDS(mtdna.all.mvBM, "./RDS/mtdna.all.mvBM.RDS")
mtdna.all.mvBM<-readRDS("./RDS/mtdna.all.mvBM.RDS")

#mtdna.proteins.mvBM<-mvBM(tree=simmap.janus.mtdna.proteins, data = scaling_data, error = scaling_data_errors, model="BMM", diagnostic=FALSE, echo=FALSE, method='rpf', scale.height = F)
#saveRDS(mtdna.proteins.mvBM, "./RDS/mtdna.proteins.mvBM.RDS")
mtdna.proteins.mvBM<- readRDS("./RDS/mtdna.proteins.mvBM.RDS")

#match proteins bc the rrna signal is identical
#mtdna.rrnas.mvBM<-mvBM(tree=simmap.janus.mtdna.proteins, data = scaling_data, error = scaling_data_errors, model="BMM", diagnostic=FALSE, echo=FALSE, method='rpf', scale.height = F)
#saveRDS(mtdna.rrnas.mvBM, "./RDS/mtdna.rrnas.mvBM.RDS")
mtdna.rrnas.mvBM<-readRDS("./RDS/mtdna.rrnas.mvBM.RDS")

#nuc.aggregate.mvBM<-mvBM(tree=simmap.janus.nuc.aggregate, data = scaling_data, error = scaling_data_errors, model="BMM", diagnostic=FALSE, echo=FALSE, method='rpf', scale.height = F)
#saveRDS(nuc.aggregate.mvBM, "./RDS/nuc.aggregate.mvBM.RDS")
nuc.aggregate.mvBM<-readRDS("./RDS/nuc.aggregate.mvBM.RDS")

#nuc.alldata.mvBM<-mvBM(tree=simmap.janus.nuc.alldata, data = scaling_data, error = scaling_data_errors, model="BMM", diagnostic=FALSE, echo=FALSE, method='rpf', scale.height = F)
#saveRDS(nuc.alldata.mvBM, "./RDS/nuc.alldata.mvBM.RDS")
nuc.alldata.mvBM <- readRDS("./RDS/nuc.alldata.mvBM.RDS")

#nuc.exons.mvBM<-mvBM(tree=simmap.janus.nuc.exons, data = scaling_data, error = scaling_data_errors, model="BMM", diagnostic=FALSE, echo=FALSE, method='rpf', scale.height = F)
#saveRDS(nuc.exons.mvBM, "./RDS/nuc.exons.mvBM.RDS")
nuc.exons.mvBM<- readRDS("./RDS/nuc.exons.mvBM.RDS")

#nuc.introns.mvBM<-mvBM(tree=simmap.janus.nuc.introns, data = scaling_data, error = scaling_data_errors, model="BMM", diagnostic=FALSE, echo=FALSE, method='rpf', scale.height = F)
#saveRDS(nuc.introns.mvBM, "./RDS/nuc.introns.mvBM.RDS")
nuc.introns.mvBM<-readRDS("./RDS/nuc.introns.mvBM.RDS")

#nuc.utrs.mvBM<-mvBM(tree=simmap.janus.nuc.utrs, data = scaling_data, error = scaling_data_errors, model="BMM", diagnostic=FALSE, echo=FALSE, method='rpf', scale.height = F)
#saveRDS(nuc.utrs.mvBM, "./RDS/nuc.utrs.mvBM.RDS")
nuc.utrs.mvBM<-readRDS("./RDS/nuc.utrs.mvBM.RDS")

results<-list(janus.null.mvBM, all.aggregte.mvBM, mtdna.all.mvBM, mtdna.proteins.mvBM, 
              mtdna.rrnas.mvBM, nuc.aggregate.mvBM, nuc.alldata.mvBM, nuc.exons.mvBM,
              nuc.introns.mvBM, nuc.utrs.mvBM)

aicw(results)

#warnings ok
janus.null.mvBM.imputed <- (estim(tree=simmap.janus.null, data = scaling_data, object = janus.null.mvBM, error = scaling_data_errors))

all.aggregte.mvBM.imputed <-(estim(tree=simmap.janus.all.aggregate, data = scaling_data, object = all.aggregte.mvBM, error = scaling_data_errors))

nuc.aggregate.mvBM.imputed <-(estim(tree=simmap.janus.nuc.aggregate, data = scaling_data, object = nuc.aggregate.mvBM, error = scaling_data_errors))

summary(lm((nuc.aggregate.mvBM.imputed$estimates[,2][janus.null.mvBM.imputed$NA_index-198] ~ all.aggregte.mvBM.imputed$estimates[,2][janus.null.mvBM.imputed$NA_index-198])))

janus.null.mvBM.imputed$NA_index
}


#now trying to fit a variety of alternative models
# # Fitting the models
# # BM1 - (Equal rate matrix)
# model_1<-mvBM(tree=simmap.janus.nuc.alldata, data = scaling_data, error = scaling_data_errors, model="BM1", diagnostic=FALSE, echo=FALSE, method='inverse', scale.height = F)
# 
# # BMM - (Proportional rate matrices)
# model_2<-mvBM(tree=simmap.janus.nuc.alldata, data = scaling_data, error = scaling_data_errors, model="BMM", param=list(constraint="proportional"), diagnostic=FALSE, echo=FALSE, method='inverse', scale.height = F)
# 
# # BMM - (Shared eigenvectors between rate matrices)
# model_3<-mvBM(tree=simmap.janus.nuc.alldata, data = scaling_data, error = scaling_data_errors, model="BMM", param=list(constraint="shared"), diagnostic=FALSE, echo=FALSE, method='inverse', scale.height = F)
# 
# # BMM - (Similar correlations between rate matrices)
# model_4<-mvBM(tree=simmap.janus.nuc.alldata, data = scaling_data, error = scaling_data_errors, model="BMM", param=list(constraint="correlation"), diagnostic=FALSE, echo=FALSE, method='inverse', scale.height = F)
# 
# # BMM - (Similar variances between rate matrices)
# model_5<-mvBM(tree=simmap.janus.nuc.alldata, data = scaling_data, error = scaling_data_errors, model="BMM", param=list(constraint="variance"), diagnostic=FALSE, echo=FALSE, method='inverse', scale.height = F)
# 
# # BMM - (Independent rate matrices)
# model_6<-mvBM(tree=simmap.janus.nuc.alldata, data = scaling_data, error = scaling_data_errors, model="BMM", diagnostic=FALSE, echo=FALSE, method='inverse', scale.height = F)
# 
# #try increasing
# control = list(maxit = 20000)
# 
# LRT(model_6,model_5)
# LRT(model_6,model_4)
# LRT(model_6,model_3)
# LRT(model_6,model_2)
# LRT(model_6,model_1)
# 
# 
# results <- list(model_1,model_2,model_3, model_4, model_5,model_6)
# aicw(results)
# 
# cov2cor(model_6$sigma[,,"zero"])
# cov2cor(model_6$sigma[,,"one"])
# cov2cor(model_6$sigma[,,"two"])
# cov2cor(model_6$sigma[,,"three"])
# cov2cor(model_6$sigma[,,"four"])
# cov2cor(model_6$sigma[,,"five"])
# cov2cor(model_6$sigma[,,"six"])
# 
# mvbm.imputed<-(estim(tree=simmap.janus.nuc.alldata, data = scaling_data, object = model_6, error = scaling_data_errors))
# 
# 
# 
# bmm.fit.scaled<-mvMORPH::mvBM(tree=simmap.janus.nuc.alldata, data = scaling_data.scaled, error = scaling_data_errors, method='rpf')
# #oum.fit<-mvOU(tree=simmap.janus.nuc.alldata, data = scaling_data.scaled, error = scaling_data_errors, method='rpf')
# bmm.fit<-mvMORPH::mvBM(tree=simmap.janus.nuc.alldata, data = scaling_data, error = scaling_data_errors, method='rpf')
# 
# 
# plot(scaling_data.scaled,asp=1, cex=0.8, main="bmm")
# ellipse(c(0,0), bmm.fit.scaled$sigma[,,"zero"], sqrt(qchisq(.95,2)), col="#a6cee3")
# #plot(scaling_data.scaled,asp=1, cex=0.8, main="bmm", ylim=c(-5,5))
# ellipse(c(0,0), bmm.fit.scaled$sigma[,,"one"], sqrt(qchisq(.95,2)), col="#1f78b4")
# #plot(scaling_data.scaled,asp=1, cex=0.8, main="bmm", ylim=c(-5,5))
# ellipse(c(0,0), bmm.fit.scaled$sigma[,,"two"], sqrt(qchisq(.95,2)), col="#b2df8a")
# #plot(scaling_data.scaled,asp=1, cex=0.8, main="bmm", ylim=c(-5,5))
# ellipse(c(0,0), bmm.fit.scaled$sigma[,,"three"], sqrt(qchisq(.95,2)), col="#33a02c")
# #plot(scaling_data.scaled,asp=1, cex=0.8, main="bmm", ylim=c(-5,5))
# ellipse(c(0,0), bmm.fit.scaled$sigma[,,"four"], sqrt(qchisq(.95,2)), col="#fb9a99")
# #plot(scaling_data.scaled,asp=1, cex=0.8, main="bmm", ylim=c(-5,5))
# ellipse(c(0,0), bmm.fit.scaled$sigma[,,"five"], sqrt(qchisq(.95,2)), col="#e31a1c")
# #plot(scaling_data.scaled,asp=1, cex=0.8, main="bmm", ylim=c(-5,5))
# ellipse(c(0,0), bmm.fit.scaled$sigma[,,"six"], sqrt(qchisq(.95,2)), col="#fdbf6f")
# 
# plot(scaling_data.scaled,asp=1, cex=0.8, main="bmm", ylim=c(-5,5))
# ellipse(c(0,0), bmm.fit$sigma[,,"seven"], sqrt(qchisq(.95,2)), col="#ff7f00")
# plot(scaling_data.scaled,asp=1, cex=0.8, main="bmm", ylim=c(-5,5))
# ellipse(c(0,0), bmm.fit$sigma[,,"eight"], sqrt(qchisq(.95,2)), col="#cab2d6")
# 
# plot(simmap.janus.nuc.alldata)
# 
# 


# 
# #trying gls.ancova https://smaerslab.com/gls-ancova/
# 
# mvbm.imputed
# require(evomap)
# Y<-"bmr"
# X<-"mass"
# data<-mvbm.imputed$estimates[,c(which(colnames(mvbm.imputed$estimates)==Y),which(colnames(mvbm.imputed$estimates)==X)),drop=F]
# tree<-as.phylo(simmap.janus.nuc.alldata)
# tree$node.label<-NULL
# tree<-treedata(tree,data,sort=T,warnings=T)$phy
# data<-as.data.frame(treedata(tree,data,sort=T,warnings=T)$data)   #match the data to the tree
# colnames(data)<-c("Dependent","Independent")     #Note that naming the variables 'Dependent' and 'Independent' serves only to standardize the procedure below across different data sets.
# rownames(data)<-rownames(mvbm.imputed$estimates)
# data<-data[tree$tip.label,]
# 
# statezero<-names(getStates(simmap.janus.nuc.alldata, type = 'tips')[getStates(simmap.janus.nuc.alldata, type = 'tips') == "zero"])
# stateone<-names(getStates(simmap.janus.nuc.alldata, type = 'tips')[getStates(simmap.janus.nuc.alldata, type = 'tips') == "one"])
# statetwo<-names(getStates(simmap.janus.nuc.alldata, type = 'tips')[getStates(simmap.janus.nuc.alldata, type = 'tips') == "two"])
# statethree<-names(getStates(simmap.janus.nuc.alldata, type = 'tips')[getStates(simmap.janus.nuc.alldata, type = 'tips') == "three"])
# statefour<-names(getStates(simmap.janus.nuc.alldata, type = 'tips')[getStates(simmap.janus.nuc.alldata, type = 'tips') == "four"])
# statefive<-names(getStates(simmap.janus.nuc.alldata, type = 'tips')[getStates(simmap.janus.nuc.alldata, type = 'tips') == "five"])
# statesix<-names(getStates(simmap.janus.nuc.alldata, type = 'tips')[getStates(simmap.janus.nuc.alldata, type = 'tips') == "six"])
# 
# is_tip <- tree$edge[,2] <= length(tree$tip.label)
# ordered_tips <- tree$edge[is_tip, 2]
# ordered_tips<-setNames(ordered_tips, tree$tip.label)
# statezero<-ordered_tips[statezero]
# stateone<-ordered_tips[stateone]
# statetwo<-ordered_tips[statetwo]
# statethree<-ordered_tips[statethree]
# statefour<-ordered_tips[statefour]
# statefive<-ordered_tips[statefive]
# statesix<-ordered_tips[statesix]
# 
# 
# layout(t(matrix(c(1,2,
#                   1,3,
#                   1,4,
#                   1,5,
#                   1,6,
#                   1,7), nrow=2)))
# 
# #.pardefault <- par()
# #par(mfrow=c(1,2))
# #par(mar=c(5.1, 4.1, 4.1, 2.1))
# plot(simmap.janus.nuc.alldata, fsize=0.0001)
# #tiplabels(getStates(simmap.janus.nuc.alldata, type = 'tips'), cex=0.5)
# 
# par(mar=c(5.1, 4.1, 4.1, 2.1))
# 
# plot(data$Dependent~data$Independent,col="white",pch=19,xlab="log body mass", ylab="log BMR",cex.lab=1)
# #pGLS.plotGrade("Dependent","Independent", data,tree,model="BM", group=statezero, linecol=colorspace::adjust_transparency(col='#F5C710', alpha = 0.85), pchcol='black', pchbg=colorspace::adjust_transparency(col='#F5C710', alpha = 0.75),linelwd=5, pchlwd=0.5,cex=2.5,pch=21)
# pGLS.plotGrade("Dependent","Independent",data,tree,model="BM",group=stateone,col="#61D04F",lwd=3,cex=0.25,pch=19)
# pGLS.plotGrade("Dependent","Independent",data,tree,model="BM",group=statetwo,col="#CD0BBC",lwd=3,cex=0.25,pch=19)
# pGLS.plotGrade("Dependent","Independent",data,tree,model="BM",group=statethree,col="#28E2E5",lwd=3,cex=0.25,pch=19)
# pGLS.plotGrade("Dependent","Independent",data,tree,model="BM",group=statefour,col="#DF536B",lwd=3,cex=0.25,pch=19)
# pGLS.plotGrade("Dependent","Independent",data,tree,model="BM",group=statefive,col="black",lwd=3,cex=0.25,pch=19)
# pGLS.plotGrade("Dependent","Independent",data,tree,model="BM",group=statesix,col="#2297E6",lwd=3,cex=0.25,pch=19)
# 
# # par(mar=c(5.1, 4.1, 4.1, 2.1))
# # plot(data$Dependent~data$Independent,col="white",pch=19,xlab="log body mass", ylab="log BMR",cex.lab=1)
# # pGLS.plotGrade("Dependent","Independent",data,tree,model="lambda",group=statezero,col="#F5C710",lwd=3,cex=0.25,pch=19)
# # pGLS.plotGrade("Dependent","Independent",data,tree,model="lambda",group=stateone,col="#61D04F",lwd=3,cex=0.25,pch=19)
# # pGLS.plotGrade("Dependent","Independent",data,tree,model="lambda",group=statetwo,col="#CD0BBC",lwd=3,cex=0.25,pch=19)
# # pGLS.plotGrade("Dependent","Independent",data,tree,model="lambda",group=statethree,col="#28E2E5",lwd=3,cex=0.25,pch=19)
# # pGLS.plotGrade("Dependent","Independent",data,tree,model="lambda",group=statefour,col="#DF536B",lwd=3,cex=0.25,pch=19)
# # pGLS.plotGrade("Dependent","Independent",data,tree,model="lambda",group=statefive,col="black",lwd=3,cex=0.25,pch=19)
# # pGLS.plotGrade("Dependent","Independent",data,tree,model="lambda",group=statesix,col="#2297E6",lwd=3,cex=0.25,pch=19)
# # 
# # dev.off()
# # 
# 
# 
# layout(t(matrix(c(1,2,
#                   1,3,
#                   1,4,
#                   1,5,
#                   1,6,
#                   1,7), nrow=2)))
# 
# par(mar=c(5.1, 4.1, 4.1, 2.1))
# plot(simmap.janus.nuc.alldata, fsize=0.0001)
# 
# 
# 
# 
# #trying gls.ancova https://smaerslab.com/gls-ancova/
# 
# #zero to four
# par(mar=c(4, 4, 0, 3))
# plot(data$Dependent~data$Independent,pch=21,xlab="log body mass", ylab="log BMR",cex.lab=1, col='black', bg=colorspace::adjust_transparency(col='grey', alpha = 0.85), lwd=0.5, cex=0.75)#, asp=1)
# #plot(data$Dependent~data$Independent,pch=21,xlab="log body mass", ylab="log BMR",cex.lab=1, col='white', bg=colorspace::adjust_transparency(col='white', alpha = 0.85), lwd=0.5, cex=2.5, asp=1)
# pGLS.plotGrade("Dependent","Independent", data, tree,model="BM", group=statezero, linecol=colorspace::adjust_transparency(col='#F5C710', alpha = 0.85), pchcol='black', pchbg=colorspace::adjust_transparency(col='#F5C710', alpha = 0.75),linelwd=5, pchlwd=0.5,cex=2.5,pch=21)
# #ellipse(c(mean(data[names(statezero),][,2]),mean(data[names(statezero),][,1])), bmm.fit$sigma[,,"zero"], sqrt(qchisq(.95,2)), col="#F5C710")
# pGLS.plotGrade("Dependent","Independent", data,tree,model="BM", group=statefour, linecol=colorspace::adjust_transparency(col='#DF536B', alpha = 0.85), pchcol='black', pchbg=colorspace::adjust_transparency(col='#DF536B', alpha = 0.75),linelwd=5, pchlwd=0.5,cex=2.5,pch=21)
# #ellipse(c(mean(data[names(statefour),][,2]),mean(data[names(statefour),][,1])), bmm.fit$sigma[,,"four"], sqrt(qchisq(.95,2)), col="#DF536B")
# 
# 
# #zero to six
# plot(data$Dependent~data$Independent,pch=21,xlab="log body mass", ylab="log BMR",cex.lab=1, col='black', bg=colorspace::adjust_transparency(col='grey', alpha = 0.85), lwd=0.5, cex=0.75)#, asp=1)
# pGLS.plotGrade("Dependent","Independent", data,tree,model="BM", group=statezero, linecol=colorspace::adjust_transparency(col='#F5C710', alpha = 0.85), pchcol='black', pchbg=colorspace::adjust_transparency(col='#F5C710', alpha = 0.75),linelwd=5, pchlwd=0.5,cex=2.5,pch=21)
# pGLS.plotGrade("Dependent","Independent", data,tree,model="BM", group=statesix, linecol=colorspace::adjust_transparency(col='#2297E6', alpha = 0.85), pchcol='black', pchbg=colorspace::adjust_transparency(col='#2297E6', alpha = 0.75),linelwd=5, pchlwd=0.5,cex=2.5,pch=21)
# 
# #six to two
# plot(data$Dependent~data$Independent,pch=21,xlab="log body mass", ylab="log BMR",cex.lab=1, col='black', bg=colorspace::adjust_transparency(col='grey', alpha = 0.85), lwd=0.5, cex=0.75)#, asp=1)
# pGLS.plotGrade("Dependent","Independent", data,tree,model="BM", group=statesix, linecol=colorspace::adjust_transparency(col='#2297E6', alpha = 0.85), pchcol='black', pchbg=colorspace::adjust_transparency(col='#2297E6', alpha = 0.75),linelwd=5, pchlwd=0.5,cex=2.5,pch=21)
# pGLS.plotGrade("Dependent","Independent", data,tree,model="BM", group=statetwo, linecol=colorspace::adjust_transparency(col='#CD0BBC', alpha = 0.85), pchcol='black', pchbg=colorspace::adjust_transparency(col='#CD0BBC', alpha = 0.75),linelwd=5, pchlwd=0.5,cex=2.5,pch=21)
# 
# #six to one
# plot(data$Dependent~data$Independent,pch=21,xlab="log body mass", ylab="log BMR",cex.lab=1, col='black', bg=colorspace::adjust_transparency(col='grey', alpha = 0.85), lwd=0.5, cex=0.75)#, asp=1)
# pGLS.plotGrade("Dependent","Independent", data,tree,model="BM", group=statesix, linecol=colorspace::adjust_transparency(col='#2297E6', alpha = 0.85), pchcol='black', pchbg=colorspace::adjust_transparency(col='#2297E6', alpha = 0.75),linelwd=5, pchlwd=0.5,cex=2.5,pch=21)
# pGLS.plotGrade("Dependent","Independent", data,tree,model="BM", group=stateone, linecol=colorspace::adjust_transparency(col='#61D04F', alpha = 0.85), pchcol='black', pchbg=colorspace::adjust_transparency(col='#61D04F', alpha = 0.75),linelwd=5, pchlwd=0.5,cex=2.5,pch=21)
# 
# #six to five
# plot(data$Dependent~data$Independent,pch=21,xlab="log body mass", ylab="log BMR",cex.lab=1, col='black', bg=colorspace::adjust_transparency(col='grey', alpha = 0.85), lwd=0.5, cex=0.75)#, asp=1)
# pGLS.plotGrade("Dependent","Independent", data,tree,model="BM", group=statesix, linecol=colorspace::adjust_transparency(col='#2297E6', alpha = 0.85), pchcol='black', pchbg=colorspace::adjust_transparency(col='#2297E6', alpha = 0.75),linelwd=5, pchlwd=0.5,cex=2.5,pch=21)
# pGLS.plotGrade("Dependent","Independent", data,tree,model="BM", group=statefive, linecol=colorspace::adjust_transparency(col='black', alpha = 0.85), pchcol='black', pchbg=colorspace::adjust_transparency(col='black', alpha = 0.75),linelwd=5, pchlwd=0.5,cex=2.5,pch=21)
# 
# #six to three
# plot(data$Dependent~data$Independent,pch=21,xlab="log body mass", ylab="log BMR",cex.lab=1, col='black', bg=colorspace::adjust_transparency(col='grey', alpha = 0.85), lwd=0.5, cex=0.75)#, asp=1)
# pGLS.plotGrade("Dependent","Independent", data,tree,model="BM", group=statesix, linecol=colorspace::adjust_transparency(col='#2297E6', alpha = 0.85), pchcol='black', pchbg=colorspace::adjust_transparency(col='#2297E6', alpha = 0.75),linelwd=5, pchlwd=0.5,cex=2.5,pch=21)
# pGLS.plotGrade("Dependent","Independent", data,tree,model="BM", group=statethree, linecol=colorspace::adjust_transparency(col='#28E2E5', alpha = 0.85), pchcol='black', pchbg=colorspace::adjust_transparency(col='#28E2E5', alpha = 0.75),linelwd=5, pchlwd=0.5,cex=2.5,pch=21)
# 
# dev.off()
# 
# 
# 
# plot(scaling_data,asp=1, cex=0.8, main="bmm")
# ellipse(c(mean(scaling_data[,1]), mean(scaling_data[,2], na.rm=T)), bmm.fit$sigma[,,"zero"], sqrt(qchisq(.95,2)), col="#a6cee3")
# 
# 
# 
# 
# grpS<-as.factor(getStates(simmap.janus.nuc.alldata, type = 'tips'))
# grpI<-as.factor(getStates(simmap.janus.nuc.alldata, type = 'tips'))
# 
# #running pANCOVA
# Model<-model.matrix(as.formula(Dependent~Independent),data)
# 
# #(1) Differences in slopes, holding intercept constant: 
# Model_S<-model.matrix(as.formula(Dependent~grpS:Independent),data) 
# 
# #(2) Differences in intercept, holding slopes constant: 
# Model_I<-model.matrix(as.formula(Dependent~grpI + Independent),data) 
# 
# #(3) Differences in slopes and differences in intercept: 
# Model_SI<-model.matrix(as.formula(Dependent~grpI + grpS:Independent),data)
# 
# #(1) Differences in slopes, holding intercept constant:
# gls.ancova(Dependent~Independent,vcv(tree),Model,Model_S)
# #vcv(rescale(x=tree, model="lambda", 0.99)
# #(2) Differences in intercept, holding slopes constant:
# gls.ancova(Dependent~Independent,vcv(tree),Model,Model_I)
# 
# #(3) Differences in slopes and differences in intercept:
# gls.ancova(Dependent~Independent,vcv.phylo(tree),Model,Model_SI)
# 
# gls.ancova(Dependent~Independent,vcv(tree),Model_I,Model_SI)
# 
# ?vcv
# rescale(x=tree, "OU", 0.5)
# 

}

#Section 13
###########################################################
#setting up data export for random forest classification ##
###########################################################
{
#data export for random forest test
{
  
  #saveRDS(bm.fit.corr$anc_recon[1:198,][,c(1:8)], file="./RDS/RF_test_data.RDS")
  #add in the "aggregate shift signal" to megaLHT
  megaLHT$aggregate_models<-getStates(simmap.janus.all.aggregate, type="tips")
  #saveRDS(megaLHT, file="./RDS/megaLHT.RDS")
  
  #generate distances
  
{
  distances<-as.dist(cophenetic.phylo(as.phylo(simmap.janus.all.aggregate)), diag=T)
  #saveRDS(distances, file="./RDS/dist.RDS")
  #adephylo:::distTips(as.phylo(simmap.janus.all.aggregate))
  #distances<-cophenetic.phylo(as.phylo(simmap.janus.all.aggregate))
  distances.gene<-as.dist(cophenetic.phylo(as.phylo(consensus.all))[labels(distances),labels(distances)], diag=T)
  
  #export tree
  #saveRDS(as.phylo(simmap.janus.all.aggregate), file="./RDS/RF_timetree.RDS")
  
  pcoa.dist<-ape::pcoa(distances)
  #rownames(pcoa.dist$values)<-labels(distances)
  pcoa.dist.gene<-ape::pcoa(distances.gene)
  #plot(scale(pcoa.dist$vectors[,c(1,2)]))

}
#save objects
#saveRDS(pcoa.dist, file="./RDS/pcoa.dist.RDS")
#saveRDS(pcoa.dist.gene, file="./RDS/pcoa.dist.gene.RDS")
}


#plot(log(megaLHT$mass) ~ log(megaLHT$nuc.scuo))
}

#Section 14
###################
#setting up l1ou ##
###################

{
require(l1ou)
l1ou_test<-adjust_data(as.phylo(simmap.janus.null), bm.fit.corr$anc_recon[1:198,][,c(1:8)])

#simulate null data for l1ou
library(ratematrix)
set.seed(5)
l1ou_test_null<-adjust_data(as.phylo(simmap.janus.null), simRatematrix(tree = consensus.all.timetree$phy, vcv = bm.fit.corr$pars$phylocov))

##NULL BOOTSTRAP VCV SHOULD BE GLOBAL OU, not only BM???###
### REVISIT THIS ####


## use parallel computing to accelerate the computation 
edges <- edge_indices_N(tree=l1ou_test$tree, min=4) #add 14

#unconstrained models
{
# eModel.par.unconstrained.aic <- estimate_shift_configuration(l1ou_test$tree, 
#                                                              l1ou_test$Y, 
#                                                              nCores=50, 
#                                                              quiet=F, 
#                                                              candid.edges = c(edges), 
#                                                              edge.length.threshold = 7.066136e-05, 
#                                                              max.nShifts = 60, 
#                                                              criterion="AIC")
# saveRDS(eModel.par.unconstrained.aic, file="./RDS/eModel.par.unconstrained.aic.RDS")
eModel.par.unconstrained.aic <- readRDS(file="./RDS/eModel.par.unconstrained.aic.RDS")

# eModel.par.unconstrained.bic <- estimate_shift_configuration(l1ou_test$tree, 
#                                                              l1ou_test$Y, 
#                                                              nCores=50, 
#                                                              quiet=F, 
#                                                              candid.edges = c(edges), 
#                                                              edge.length.threshold = 7.066136e-05, 
#                                                              max.nShifts = 60, 
#                                                              criterion="BIC")
# saveRDS(eModel.par.unconstrained.bic, file="./RDS/eModel.par.unconstrained.bic.RDS")
eModel.par.unconstrained.bic <- readRDS(file="./RDS/eModel.par.unconstrained.bic.RDS")


# eModel.par.unconstrained.pbic <- estimate_shift_configuration(l1ou_test$tree, 
#                                                               l1ou_test$Y, 
#                                                               nCores=50, 
#                                                               quiet=F, 
#                                                               candid.edges = c(edges), 
#                                                               edge.length.threshold = 7.066136e-05, 
#                                                               max.nShifts = 60, 
#                                                               criterion="pBIC")
# saveRDS(eModel.par.unconstrained.pbic, file="./RDS/eModel.par.unconstrained.pbic.RDS")
eModel.par.unconstrained.pbic <- readRDS("./RDS/eModel.par.unconstrained.pbic.RDS")


# eModel.par.unconstrained.aicc <- estimate_shift_configuration(l1ou_test$tree, 
#                                                               l1ou_test$Y, 
#                                                               nCores=50, 
#                                                               quiet=F, 
#                                                               candid.edges = c(edges), 
#                                                               edge.length.threshold = 7.066136e-05, 
#                                                               max.nShifts = 60, 
#                                                               criterion="AICc")
# saveRDS(eModel.par.unconstrained.aicc, file="./RDS/eModel.par.unconstrained.aicc.RDS")
eModel.par.unconstrained.aicc <- readRDS(file="./RDS/eModel.par.unconstrained.aicc.RDS")

# eModel.par.unconstrained.pBICess <- estimate_shift_configuration(l1ou_test$tree, 
#                                                                  l1ou_test$Y, 
#                                                                  nCores=50, 
#                                                                  quiet=F, 
#                                                                  candid.edges = c(edges), 
#                                                                  edge.length.threshold = 7.066136e-05, 
#                                                                  max.nShifts = 60, 
#                                                                  criterion="pBICess")
# saveRDS(eModel.par.unconstrained.pBICess, file="./RDS/eModel.par.unconstrained.pBICess.RDS")
eModel.par.unconstrained.pBICess <- readRDS(file="./RDS/eModel.par.unconstrained.pBICess.RDS")

# eModel.par.unconstrained.mBIC <- estimate_shift_configuration(l1ou_test$tree, 
#                                                               l1ou_test$Y, 
#                                                               nCores=50, 
#                                                               quiet=F, 
#                                                               candid.edges = c(edges), 
#                                                               edge.length.threshold = 7.066136e-05, 
#                                                               max.nShifts = 60, 
#                                                               criterion="mBIC")
# saveRDS(eModel.par.unconstrained.mBIC, file="./RDS/eModel.par.unconstrained.mBIC.RDS")
eModel.par.unconstrained.mBIC <- readRDS(file="./RDS/eModel.par.unconstrained.mBIC.RDS")
}

#these two lines generate a result that is very similar to janus, but it is quite constrained
edges <- edge_indices_nodes(tree=l1ou_test$tree, nodes=con_prop_logit[con_prop_logit$uncex.merged.mtdnas==1,]$node.all)

#constrained models
{
# eModel.par.constrained.aic <- estimate_shift_configuration(l1ou_test$tree, 
#                                                            l1ou_test$Y, 
#                                                            nCores=50, 
#                                                            quiet=F, 
#                                                            candid.edges = c(350, 349, edges), 
#                                                            edge.length.threshold = 7.066136e-05, 
#                                                            max.nShifts = 12, 
#                                                            criterion="AIC")
# saveRDS(eModel.par.constrained.aic, file="./RDS/eModel.par.constrained.aic.RDS")
eModel.par.constrained.aic <- readRDS(file="./RDS/eModel.par.constrained.aic.RDS")


# eModel.par.constrained.bic <- estimate_shift_configuration(l1ou_test$tree,
#                                                            l1ou_test$Y, 
#                                                            nCores=50, 
#                                                            quiet=F, 
#                                                            candid.edges = c(350, 349, edges), 
#                                                            edge.length.threshold = 7.066136e-05, 
#                                                            max.nShifts = 12, 
#                                                            criterion="BIC")
# saveRDS(eModel.par.constrained.bic, file="./RDS/eModel.par.constrained.bic.RDS")
eModel.par.constrained.bic <- readRDS(file="./RDS/eModel.par.constrained.bic.RDS")
# 
# eModel.par.constrained.pbic <- estimate_shift_configuration(l1ou_test$tree, 
#                                                             l1ou_test$Y, 
#                                                             nCores=50, 
#                                                             quiet=F, 
#                                                             candid.edges = c(350, 349, edges), 
#                                                             edge.length.threshold = 7.066136e-05, 
#                                                             max.nShifts = 12, 
#                                                             criterion="pBIC")
# saveRDS(eModel.par.constrained.pbic, file="./RDS/eModel.par.constrained.pbic.RDS")
eModel.par.constrained.pbic <- readRDS(file="./RDS/eModel.par.constrained.pbic.RDS")

# eModel.par.constrained.aicc <- estimate_shift_configuration(l1ou_test$tree, 
#                                                             l1ou_test$Y, 
#                                                             nCores=50, 
#                                                             quiet=F, 
#                                                             candid.edges = c(350, 349, edges), 
#                                                             edge.length.threshold = 7.066136e-05, 
#                                                             max.nShifts = 12, 
#                                                             criterion="AICc")
# saveRDS(eModel.par.constrained.aicc, file="./RDS/eModel.par.constrained.aicc.RDS")
eModel.par.constrained.aicc <- readRDS(file="./RDS/eModel.par.constrained.aicc.RDS")


# eModel.par.constrained.pBICess <- estimate_shift_configuration(l1ou_test$tree, 
#                                                                l1ou_test$Y, 
#                                                                nCores=50, 
#                                                                quiet=F, 
#                                                                candid.edges = c(350, 349, edges), 
#                                                                edge.length.threshold = 7.066136e-05, 
#                                                                max.nShifts = 12, 
#                                                                criterion="pBICess")
# saveRDS(eModel.par.constrained.pBICess, file="./RDS/eModel.par.constrained.pBICess.RDS")
eModel.par.constrained.pBICess <- readRDS(file="./RDS/eModel.par.constrained.pBICess.RDS")

# 
# eModel.par.constrained.mBIC <- estimate_shift_configuration(l1ou_test$tree, 
#                                                             l1ou_test$Y, 
#                                                             nCores=50, 
#                                                             quiet=F, 
#                                                             candid.edges = c(350, 349, edges), 
#                                                             edge.length.threshold = 7.066136e-05, 
#                                                             max.nShifts = 12, 
#                                                             criterion="mBIC")
# saveRDS(eModel.par.constrained.mBIC, file="./RDS/eModel.par.constrained.mBIC.RDS")
eModel.par.constrained.mBIC <- readRDS(file="./RDS/eModel.par.constrained.mBIC.RDS")
}

#testing null fits

## use parallel computing to accelerate the computation 
edges <- edge_indices_N(tree=l1ou_test$tree, min=4) #add 14
{
# eModel.par.unconstrained.aic.null <- estimate_shift_configuration(l1ou_test_null$tree,
#                                                              l1ou_test_null$Y,
#                                                              nCores=50,
#                                                              quiet=F,
#                                                              candid.edges = c(edges),
#                                                              edge.length.threshold = 7.066136e-05,
#                                                              max.nShifts = 60,
#                                                              criterion="AIC")
# saveRDS(eModel.par.unconstrained.aic.null, file="./RDS/eModel.par.unconstrained.aic.null.RDS")
eModel.par.unconstrained.aic.null <- readRDS(file="./RDS/eModel.par.unconstrained.aic.null.RDS")

# eModel.par.unconstrained.bic.null <- estimate_shift_configuration(l1ou_test_null$tree,
#                                                              l1ou_test_null$Y,
#                                                              nCores=50,
#                                                              quiet=F,
#                                                              candid.edges = c(edges),
#                                                              edge.length.threshold = 7.066136e-05,
#                                                              max.nShifts = 60,
#                                                              criterion="BIC")
# saveRDS(eModel.par.unconstrained.bic.null, file="./RDS/eModel.par.unconstrained.bic.null.RDS")
eModel.par.unconstrained.bic.null <- readRDS(file="./RDS/eModel.par.unconstrained.bic.null.RDS")


# eModel.par.unconstrained.pbic.null <- estimate_shift_configuration(l1ou_test_null$tree,
#                                                               l1ou_test_null$Y,
#                                                               nCores=50,
#                                                               quiet=F,
#                                                               candid.edges = c(edges),
#                                                               edge.length.threshold = 7.066136e-05,
#                                                               max.nShifts = 60,
#                                                               criterion="pBIC")
# saveRDS(eModel.par.unconstrained.pbic.null, file="./RDS/eModel.par.unconstrained.pbic.null.RDS")
eModel.par.unconstrained.pbic.null <- readRDS("./RDS/eModel.par.unconstrained.pbic.null.RDS")


# eModel.par.unconstrained.aicc.null <- estimate_shift_configuration(l1ou_test_null$tree,
#                                                               l1ou_test_null$Y,
#                                                               nCores=50,
#                                                               quiet=F,
#                                                               candid.edges = c(edges),
#                                                               edge.length.threshold = 7.066136e-05,
#                                                               max.nShifts = 60,
#                                                               criterion="AICc")
# saveRDS(eModel.par.unconstrained.aicc.null, file="./RDS/eModel.par.unconstrained.aicc.null.RDS")
eModel.par.unconstrained.aicc.null <- readRDS(file="./RDS/eModel.par.unconstrained.aicc.null.RDS")

# eModel.par.unconstrained.pBICess.null <- estimate_shift_configuration(l1ou_test_null$tree,
#                                                                  l1ou_test_null$Y,
#                                                                  nCores=50,
#                                                                  quiet=F,
#                                                                  candid.edges = c(edges),
#                                                                  edge.length.threshold = 7.066136e-05,
#                                                                  max.nShifts = 60,
#                                                                  criterion="pBICess")
# saveRDS(eModel.par.unconstrained.pBICess.null, file="./RDS/eModel.par.unconstrained.pBICess.null.RDS")
eModel.par.unconstrained.pBICess.null <- readRDS(file="./RDS/eModel.par.unconstrained.pBICess.null.RDS")

# eModel.par.unconstrained.mBIC.null <- estimate_shift_configuration(l1ou_test_null$tree,
#                                                               l1ou_test_null$Y,
#                                                               nCores=50,
#                                                               quiet=F,
#                                                               candid.edges = c(edges),
#                                                               edge.length.threshold = 7.066136e-05,
#                                                               max.nShifts = 60,
#                                                               criterion="mBIC")
# saveRDS(eModel.par.unconstrained.mBIC.null, file="./RDS/eModel.par.unconstrained.mBIC.null.RDS")
eModel.par.unconstrained.mBIC.null <- readRDS(file="./RDS/eModel.par.unconstrained.mBIC.null.RDS")


#these two lines generate a result that is very similar to janus, but it is quite constrained
edges <- edge_indices_nodes(tree=l1ou_test$tree, nodes=con_prop_logit[con_prop_logit$uncex.merged.mtdnas==1,]$node.all)

# eModel.par.constrained.aic.null <- estimate_shift_configuration(l1ou_test_null$tree,
#                                                            l1ou_test_null$Y,
#                                                            nCores=50,
#                                                            quiet=F,
#                                                            candid.edges = c(350, 349, edges),
#                                                            edge.length.threshold = 7.066136e-05,
#                                                            max.nShifts = 12,
#                                                            criterion="AIC")
# saveRDS(eModel.par.constrained.aic.null, file="./RDS/eModel.par.constrained.aic.null.RDS")
eModel.par.constrained.aic.null <- readRDS(file="./RDS/eModel.par.constrained.aic.null.RDS")


# eModel.par.constrained.bic.null <- estimate_shift_configuration(l1ou_test_null$tree,
#                                                            l1ou_test_null$Y,
#                                                            nCores=50,
#                                                            quiet=F,
#                                                            candid.edges = c(350, 349, edges),
#                                                            edge.length.threshold = 7.066136e-05,
#                                                            max.nShifts = 12,
#                                                            criterion="BIC")
# saveRDS(eModel.par.constrained.bic.null, file="./RDS/eModel.par.constrained.bic.null.RDS")
eModel.par.constrained.bic.null <- readRDS(file="./RDS/eModel.par.constrained.bic.null.RDS")

# eModel.par.constrained.pbic.null <- estimate_shift_configuration(l1ou_test_null$tree,
#                                                             l1ou_test_null$Y,
#                                                             nCores=50,
#                                                             quiet=F,
#                                                             candid.edges = c(350, 349, edges),
#                                                             edge.length.threshold = 7.066136e-05,
#                                                             max.nShifts = 12,
#                                                             criterion="pBIC")
# saveRDS(eModel.par.constrained.pbic.null, file="./RDS/eModel.par.constrained.pbic.null.RDS")
eModel.par.constrained.pbic.null <- readRDS(file="./RDS/eModel.par.constrained.pbic.null.RDS")

# eModel.par.constrained.aicc.null <- estimate_shift_configuration(l1ou_test_null$tree,
#                                                             l1ou_test_null$Y,
#                                                             nCores=50,
#                                                             quiet=F,
#                                                             candid.edges = c(350, 349, edges),
#                                                             edge.length.threshold = 7.066136e-05,
#                                                             max.nShifts = 12,
#                                                             criterion="AICc")
# saveRDS(eModel.par.constrained.aicc.null, file="./RDS/eModel.par.constrained.aicc.null.RDS")
eModel.par.constrained.aicc.null <- readRDS(file="./RDS/eModel.par.constrained.aicc.null.RDS")

# eModel.par.constrained.pBICess.null <- estimate_shift_configuration(l1ou_test_null$tree,
#                                                                l1ou_test_null$Y,
#                                                                nCores=50,
#                                                                quiet=F,
#                                                                candid.edges = c(350, 349, edges),
#                                                                edge.length.threshold = 7.066136e-05,
#                                                                max.nShifts = 12,
#                                                                criterion="pBICess")
# saveRDS(eModel.par.constrained.pBICess.null, file="./RDS/eModel.par.constrained.pBICess.null.RDS")
eModel.par.constrained.pBICess.null <- readRDS(file="./RDS/eModel.par.constrained.pBICess.null.RDS")

# eModel.par.constrained.mBIC.null <- estimate_shift_configuration(l1ou_test_null$tree,
#                                                             l1ou_test_null$Y,
#                                                             nCores=50,
#                                                             quiet=F,
#                                                             candid.edges = c(350, 349, edges),
#                                                             edge.length.threshold = 7.066136e-05,
#                                                             max.nShifts = 12,
#                                                             criterion="mBIC")
# saveRDS(eModel.par.constrained.mBIC.null, file="./RDS/eModel.par.constrained.mBIC.null.RDS")
eModel.par.constrained.mBIC.null <- readRDS(file="./RDS/eModel.par.constrained.mBIC.null.RDS")
}

#running bootstrap models
{
# fit.boot.bic.constrained<-l1ou_bootstrap_support(eModel.par.constrained.bic, nItrs = 100, multicore = T, nCores = 50, quietly = F)
# saveRDS(fit.boot.bic.constrained, file="./RDS/fit.boot.bic.constrained.RDS")
fit.boot.bic.constrained<- readRDS(file="./RDS/fit.boot.bic.constrained.RDS")

# fit.boot.aicc.constrained<-l1ou_bootstrap_support(eModel.par.constrained.aicc, nItrs = 100, multicore = T, nCores = 50, quietly = F)
# saveRDS(fit.boot.aicc.constrained, file="./RDS/fit.boot.aicc.constrained.RDS")
fit.boot.aicc.constrained <- readRDS(file="./RDS/fit.boot.aicc.constrained.RDS")
#this only completed 94 replicates -- investigate further? 

# fit.boot.aic.constrained<-l1ou_bootstrap_support(eModel.par.constrained.aic, nItrs = 100, multicore = T, nCores = 50, quietly = F)
# saveRDS(fit.boot.aic.constrained, file="./RDS/fit.boot.aic.constrained.RDS")
fit.boot.aic.constrained <-readRDS("./RDS/fit.boot.aic.constrained.RDS")

#fit.boot.pbic.constrained<-l1ou_bootstrap_support(eModel.par.constrained.pbic, nItrs = 100, multicore = T, nCores = 30, quietly = F)
#saveRDS(fit.boot.pbic.constrained, file="./RDS/fit.boot.pbic.constrained.RDS")
fit.boot.pbic.constrained<- readRDS("./RDS/fit.boot.pbic.constrained.RDS")

#fit.boot.pbicess.constrained<-l1ou_bootstrap_support(eModel.par.constrained.pBICess, nItrs = 100, multicore = T, nCores = 30, quietly = F)
#saveRDS(fit.boot.pbicess.constrained, file="./RDS/fit.boot.pbicess.constrained.RDS")
fit.boot.pbicess.constrained<- readRDS("./RDS/fit.boot.pbicess.constrained.RDS")

#fit.boot.mbic.constrained<-l1ou_bootstrap_support(eModel.par.constrained.mBIC, nItrs = 100, multicore = T, nCores = 30, quietly = F)
#saveRDS(fit.boot.mbic.constrained, file="./RDS/fit.boot.mbic.constrained.RDS")
fit.boot.mbic.constrained<- readRDS("./RDS/fit.boot.mbic.constrained.RDS")
}
## running bootstrap models for unconstrained analysis
{
#fit.boot.aicc.unconstrained<-l1ou_bootstrap_support(eModel.par.unconstrained.aicc, nItrs = 100, multicore = T, nCores = 40, quietly = F)
#saveRDS(fit.boot.aicc.unconstrained, file="./RDS/fit.boot.aicc.unconstrained.RDS")
fit.boot.aicc.unconstrained <- readRDS(file="./RDS/fit.boot.aicc.unconstrained.RDS")
#only completed 44 replicates -- need to run more

#fit.boot.pbic.unconstrained<-l1ou_bootstrap_support(eModel.par.unconstrained.pbic, nItrs = 100, multicore = T, nCores = 40, quietly = F)
#saveRDS(fit.boot.pbic.unconstrained, file="./RDS/fit.boot.pbic.unconstrained.RDS")
fit.boot.pbic.unconstrained<- readRDS("./RDS/fit.boot.pbic.unconstrained.RDS")
#only completed 54 replicates -- need to try this again 

#running bootstrap models with AICc but boots with pBIC (as suggested by authors)
eModel.par.unconstrained.aicc.mod<-eModel.par.unconstrained.aicc
eModel.par.unconstrained.aicc.mod$l1ou.options$criterion <- "pBIC"
eModel.par.constrained.aicc.mod<-eModel.par.constrained.aicc
eModel.par.constrained.aicc.mod$l1ou.options$criterion <- "pBIC"

#fit.boot.aicc.unconstrained.mod<-l1ou_bootstrap_support(eModel.par.unconstrained.aicc.mod, nItrs = 100, multicore = T, nCores = 30, quietly = F)
#saveRDS(fit.boot.aicc.unconstrained.mod, file="./RDS/fit.boot.aicc.unconstrained.mod.RDS")
fit.boot.aicc.unconstrained.mod <- readRDS(file="./RDS/fit.boot.aicc.unconstrained.mod.RDS")

#fit.boot.aicc.constrained.mod<-l1ou_bootstrap_support(eModel.par.constrained.aicc.mod, nItrs = 100, multicore = T, nCores = 30, quietly = F)
#saveRDS(fit.boot.aicc.constrained.mod, file="./RDS/fit.boot.aicc.constrained.mod")
fit.boot.aicc.constrained.mod<- readRDS("./RDS/fit.boot.aicc.constrained.mod.RDS")
}

####Testing null bootstrap#####
{
#fit.boot.aicc.constrained.null<-l1ou_bootstrap_support(eModel.par.constrained.aicc.null, nItrs = 100, multicore = T, nCores = 50, quietly = T)
#saveRDS(fit.boot.aicc.constrained.null, file="./RDS/fit.boot.aicc.constrained.null.RDS")
fit.boot.aicc.constrained.null<-readRDS(file="./RDS/fit.boot.aicc.constrained.null.RDS")


#this only completed 94 replicates -- investigate further?

# system.time(
#   boots.aicc<-l1ou_nullboot(seed = 5, 
#                             input_tree = consensus.all.timetree$phy, 
#                             boots = 500, 
#                             input_vcv =  bm.fit.corr$pars$phylocov, 
#                             ncores = 50, 
#                             IC = "AICc", 
#                             considered_edges = c(350, 349, edges), 
#                             threshold = 7.066136e-05, 
#                             max.nshifts = 12,
#                             timelim=100)
# )
# saveRDS(boots.aicc,file="./RDS/boots.aicc.RDS")
boots.aicc<- readRDS(file="./RDS/boots.aicc.RDS")

# system.time(
#   boots.pBIC<-l1ou_nullboot(seed = 5, 
#                             input_tree = consensus.all.timetree$phy, 
#                             boots = 500, 
#                             input_vcv =  bm.fit.corr$pars$phylocov, 
#                             ncores = 50, 
#                             IC = "pBIC", 
#                             considered_edges = c(350, 349, edges), 
#                             threshold = 7.066136e-05, 
#                             max.nshifts = 12,
#                             timelim=100)
# )
# saveRDS(boots.pBIC, file="./RDS/boots.pBIC.RDS")
boots.pBIC<-readRDS(file="./RDS/boots.pBIC.RDS")

}

}


#testing LHTs for EB bias
{
##testing mass
BM.mass<-fitContinuous(phy = simmap.janus.null, dat=l1ou_test$Y[,1], model='BM')
EB.mass<-fitContinuous(phy = simmap.janus.null, dat=l1ou_test$Y[,1], model='EB', bounds=list(a=c(-10,0)))
OU.mass<-fitContinuous(phy = simmap.janus.null, dat=l1ou_test$Y[,1], model='OU', bounds=list(alpha=c(0,100)))
aicw(c("BM"=AIC(BM.mass),"EB"=AIC(EB.mass), "OU"=AIC(OU.mass)))

#confirming OU
OU.mass.OUwie<-OUwie(
  phy = simmap.janus.nuc.aggregate.simplified2,
  data = data.frame(
    "Genus_species" = names(l1ou_test$Y[, 1]),
    "Reg" = rep(0, length(names(l1ou_test$Y[, 1]))),
    "X" = l1ou_test$Y[, 1]
  )
  ,
  model = "OU1",
  simmap.tree = F
)

##testing mean.clutch
BM.mc<-fitContinuous(phy = simmap.janus.null, dat=l1ou_test$Y[,2], model='BM')
EB.mc<-fitContinuous(phy = simmap.janus.null, dat=l1ou_test$Y[,2], model='EB', bounds=list(a=c(-10,0)))
OU.mc<-fitContinuous(phy = simmap.janus.null, dat=l1ou_test$Y[,2], model='OU', bounds=list(alpha=c(0,100)))
aicw(c("BM"=AIC(BM.mc),"EB"=AIC(EB.mc), "OU"=AIC(OU.mc)))

#confirming OU
OU.mc.OUwie<-OUwie(
  phy = simmap.janus.nuc.aggregate.simplified2,
  data = data.frame(
    "Genus_species" = names(l1ou_test$Y[, 2]),
    "Reg" = rep(0, length(names(l1ou_test$Y[, 2]))),
    "X" = l1ou_test$Y[, 2],
    #"mserr" = rep(0.1, length(names(l1ou_test$Y[, 2])))
  )
  ,
  model = "OU1",
  simmap.tree = F
)

##testing gen_length
BM.gen<-fitContinuous(phy = simmap.janus.null, dat=l1ou_test$Y[,3], model='BM')
EB.gen<-fitContinuous(phy = simmap.janus.null, dat=l1ou_test$Y[,3], model='EB', bounds=list(a=c(-10,0)))
OU.gen<-fitContinuous(phy = simmap.janus.null, dat=l1ou_test$Y[,3], model='OU', bounds=list(alpha=c(0,100)))
aicw(c("BM"=AIC(BM.gen),"EB"=AIC(EB.gen), "OU"=AIC(OU.gen)))

#confirming OU
OU.gen.OUwie<-OUwie(
  phy = simmap.janus.nuc.aggregate.simplified2,
  data = data.frame(
    "Genus_species" = names(l1ou_test$Y[, 3]),
    "Reg" = rep(0, length(names(l1ou_test$Y[, 3]))),
    "X" = l1ou_test$Y[, 3]
  )
  ,
  model = "OU1",
  simmap.tree = F
)

##testing survival
BM.surv<-fitContinuous(phy = simmap.janus.null, dat=l1ou_test$Y[,4], model='BM')
EB.surv<-fitContinuous(phy = simmap.janus.null, dat=l1ou_test$Y[,4], model='EB', bounds=list(a=c(-10,0)))
OU.surv<-fitContinuous(phy = simmap.janus.null, dat=l1ou_test$Y[,4], model='OU', bounds=list(alpha=c(0,100)))
aicw(c("BM"=AIC(BM.surv),"EB"=AIC(EB.surv), "OU"=AIC(OU.surv)))

#confirming OU
OU.surv.OUwie<-OUwie(
  phy = simmap.janus.nuc.aggregate.simplified2,
  data = data.frame(
    "Genus_species" = names(l1ou_test$Y[, 4]),
    "Reg" = rep(0, length(names(l1ou_test$Y[, 4]))),
    "X" = l1ou_test$Y[, 4]
  )
  ,
  model = "OU1",
  simmap.tree = F
)


##testing breeding
BM.breed<-fitContinuous(phy = simmap.janus.null, dat=l1ou_test$Y[,5], model='BM')
EB.breed<-fitContinuous(phy = simmap.janus.null, dat=l1ou_test$Y[,5], model='EB', bounds=list(a=c(-10,0)))
OU.breed<-fitContinuous(phy = simmap.janus.null, dat=l1ou_test$Y[,5], model='OU', bounds=list(alpha=c(0,100)))
aicw(c("BM"=AIC(BM.breed),"EB"=AIC(EB.breed), "OU"=AIC(OU.breed)))

#confirming OU
OU.breed.OUwie<-OUwie(
  phy = simmap.janus.nuc.aggregate.simplified2,
  data = data.frame(
    "Genus_species" = names(l1ou_test$Y[, 5]),
    "Reg" = rep(0, length(names(l1ou_test$Y[, 5]))),
    "X" = l1ou_test$Y[, 5]
  )
  ,
  model = "OU1",
  simmap.tree = F
)


##testing longevity
BM.lg<-fitContinuous(phy = simmap.janus.null, dat=l1ou_test$Y[,6], model='BM')
EB.lg<-fitContinuous(phy = simmap.janus.null, dat=l1ou_test$Y[,6], model='EB', bounds=list(a=c(-10,0)))
OU.lg<-fitContinuous(phy = simmap.janus.null, dat=l1ou_test$Y[,6], model='OU', bounds=list(alpha=c(0,100)))
aicw(c("BM"=AIC(BM.lg),"EB"=AIC(EB.lg), "OU"=AIC(OU.lg)))

#confirming OU
OU.lg.OUwie<-OUwie(
  phy = simmap.janus.nuc.aggregate.simplified2,
  data = data.frame(
    "Genus_species" = names(l1ou_test$Y[, 6]),
    "Reg" = rep(0, length(names(l1ou_test$Y[, 6]))),
    "X" = l1ou_test$Y[, 6]
  )
  ,
  model = "OU1",
  simmap.tree = F
)


##testing chick PC1
BM.pc1<-fitContinuous(phy = simmap.janus.null, dat=l1ou_test$Y[,7], model='BM')
EB.pc1<-fitContinuous(phy = simmap.janus.null, dat=l1ou_test$Y[,7], model='EB', bounds=list(a=c(-10,0)))
OU.pc1<-fitContinuous(phy = simmap.janus.null, dat=l1ou_test$Y[,7], model='OU', bounds=list(alpha=c(0,100)))
aicw(c("BM"=AIC(BM.pc1),"EB"=AIC(EB.pc1), "OU"=AIC(OU.pc1)))

#confirming OU
OU.pc1.OUwie<-OUwie(
  phy = simmap.janus.nuc.aggregate.simplified2,
  data = data.frame(
    "Genus_species" = names(l1ou_test$Y[, 7]),
    "Reg" = rep(0, length(names(l1ou_test$Y[, 7]))),
    "X" = l1ou_test$Y[, 7]
  )
  ,
  model = "OU1",
  simmap.tree = F
)


##testing latitude
BM.lat<-fitContinuous(phy = simmap.janus.null, dat=l1ou_test$Y[,8], model='BM')
EB.lat<-fitContinuous(phy = simmap.janus.null, dat=l1ou_test$Y[,8], model='EB', bounds=list(a=c(-10,0)))
OU.lat<-fitContinuous(phy = simmap.janus.null, dat=l1ou_test$Y[,8], model='OU', bounds=list(alpha=c(0,100)))
aicw(c("BM"=AIC(BM.lat),"EB"=AIC(EB.lat), "OU"=AIC(OU.lat)))

#confirming OU
OU.lat.OUwie<-OUwie(
  phy = simmap.janus.nuc.aggregate.simplified2,
  data = data.frame(
    "Genus_species" = names(l1ou_test$Y[, 8]),
    "Reg" = rep(0, length(names(l1ou_test$Y[, 8]))),
    "X" = l1ou_test$Y[, 8]
  )
  ,
  model = "OU1",
  simmap.tree = F
)


}
#no bias toward EB

#processing
#experimental bootstraps (NOT YET RUN)
require("R.utils")
# {
# system.time(
#   boots.aicc.mvmorph<-l1ou_nullboot.mvmorph(seed = 5, 
#                                             input_tree = consensus.all.timetree$phy, 
#                                             boots = 500, 
#                                             model_fit = bm.fit.mvbm, 
#                                             ncores = 50, 
#                                             IC = "AICc", 
#                                             considered_edges = c(350, 349, edges), 
#                                             threshold = 7.066136e-05, 
#                                             max.nshifts = 12,
#                                             timelim=100)
# )
# saveRDS(boots.aicc,file="./RDS/boots.aicc.mvmorph.RDS")
# 
# 
# system.time(
#   boots.pBIC.mvmorph<-l1ou_nullboot.mvmorph(seed = 5, 
#                                             input_tree = consensus.all.timetree$phy, 
#                                             boots = 500, 
#                                             model_fit = bm.fit.mvbm, 
#                                             ncores = 50, 
#                                             IC = "pBIC", 
#                                             considered_edges = c(350, 349, edges), 
#                                             threshold = 7.066136e-05, 
#                                             max.nshifts = 12,
#                                             timelim=100)
# )
# saveRDS(boots.pBIC, file="./RDS/boots.pBIC.mvmorph.RDS")
# 
# }



#Section 15
###################
#processing l1ou ##
###################
{
#process null bootstrapping
{#l1ou_nullboot_process(boots_output = boots.pBIC)[-c(6,7)]
edges <- edge_indices_nodes(tree=l1ou_test$tree, nodes=con_prop_logit[con_prop_logit$uncex.merged.mtdnas==1,]$node.all)


edge_names<-c("sister_tinamiformes",
"tinamiformes",
"notopaleognathae",
"aeqornithes",
"coraciimorphae",
"passeri",
"psittaciformes",
"sister_otidae",
"otidae",
"passerea",
"columbea",
"neognathae")

par(mfrow=c(1,2), mar=c(7, 4.1, 4.1, 2.1))
barplot(
  rbind(
    l1ou_nullboot_process(boots_output = boots.aicc)[-c(6,7,15)],
    setNames(fit.boot.aicc.constrained$detection.rate[edges],edges)),
  beside=T, ylim=c(0,1.0), main= "AICc", cex.names = 0.8, names.arg=edge_names, las=2)


barplot(
  rbind(
    l1ou_nullboot_process(boots_output = boots.pBIC)[-c(6,7,15)],
    setNames(fit.boot.pbic.constrained$detection.rate[edges],edge_names)), 
  beside=T, ylim=c(0,1.0), main = "pBIC", cex.names=0.8, names.arg=edge_names, las=2)
}

#plot(l1ou_test$tree, cex=0.0001, no.margin=T)
#edgelabels(edge=edges)



#plotting support for shifts on the aggregate shift tree

#generate relative summaries 

{
aicc.summary<-rbind(
  (l1ou_nullboot_process(boots_output = boots.aicc)[-c(6,7)]),
  c(setNames(fit.boot.aicc.constrained$detection.rate[edges],edges), n=length(fit.boot.aicc.constrained$all.shifts)))

pbic.summary<-rbind(
  (l1ou_nullboot_process(boots_output = boots.pBIC)[-c(6,7)]),
  c(setNames(fit.boot.pbic.constrained$detection.rate[edges],edges), n=length(fit.boot.pbic.constrained$all.shifts)))


#proportional support
l1ou.support.aicc<- data.frame(null=aicc.summary[1,][1:12]/(aicc.summary[1,][1:12]+aicc.summary[2,][1:12]),
                          alt = aicc.summary[2,][1:12]/(aicc.summary[1,][1:12]+aicc.summary[2,][1:12]))


#proportional support
l1ou.support.pbic<- data.frame(null=pbic.summary[1,][1:12]/(pbic.summary[1,][1:12]+pbic.summary[2,][1:12]),
                               alt = pbic.summary[2,][1:12]/(pbic.summary[1,][1:12]+pbic.summary[2,][1:12]))

}
#summary(l1ou.support.pbic$alt)

# #notes on which edge reflect which clade#
# edge 15 - notopaleognathae
# edge 14 - tinamiformes
# edge 13 - sister to tinamiformes
# edge 393 - neognathate
# edge 390 - columbea
# edge 388 - otidae 
# edge 389 - passerea
# edge 387 - remainder of neoaves
# edge 199 - aequornithes 
# edge 265 - coraciimorphae
# edge 374 - psittaciformes
# edge 370 - passeri

require(stargazer)
fisherloop(aicc.summary[c(2,1),])
fisherloop(pbic.summary[c(2,1),])

#set up matricies for heatmaps
{
#exon, intron, utr, mtdna
r1<- matrix(c(0,0,1,0), 2,2)
r2<- matrix(c(1,0,0,1), 2,2)
r3<- matrix(c(1,1,0,0), 2,2)
r4<- matrix(c(0,0,0,1), 2,2)
r5<- matrix(c(0,1,0,0), 2,2)
r6<- matrix(c(1,0,0,0), 2,2)
r7<- matrix(c(1,0,0,0), 2,2)
r8<- matrix(c(0,1,1,0), 2,2)
r9<- matrix(c(1,0,0,0), 2,2)
r10<- matrix(c(1,1,1,0), 2,2)
r11<- matrix(c(1,0,0,0), 2,2)
r12<- matrix(c(1,0,0,0), 2,2)

r.comb<- rbind(unlist(as.list(r1)),
               unlist(as.list(r2)),
               unlist(as.list(r3)),
               unlist(as.list(r4)),
               unlist(as.list(r5)),
               unlist(as.list(r6)),
               unlist(as.list(r7)),
               unlist(as.list(r8)),
               unlist(as.list(r9)),
               unlist(as.list(r10)),
               unlist(as.list(r11)),
               unlist(as.list(r12)))
}
dev.off()
#plot aggregate tree with proportaional support vals relative to the null false positive rate
}


pdf(file="Figure1.pdf", width=2.5, height=8)
{

#set up the plotting area and plot the tree
{
#window<-par()
#saveRDS(window, file="./RDS/Figure1_coordinate_par.RDS")
window<-readRDS("./RDS/Figure1_coordinate_par.RDS")
require(phyloch)

par(lwd=1)
plot_janus_time_rect_simmap(simmap.janus.all.aggregate, seed=11, title = "merged", width=1.5, out.line=F, box=F)
}

#set up the candidate edges and marks
{
jan.candidates<-as.numeric(names(node_indices_N(tree=simmap.janus.all.aggregate, min = 4))[
  node_indices_N(tree=consensus.all.timetree$phy, min=4) >=4])[-1]

par(lwd=0.4)
#nodelabels(node=jan.candidates, pch = 19, bg='black', cex=0.5)
nodelabels(node=jan.candidates, pch = 22, bg='white', cex=0.7)
}

#set up the hits and heatmaps
{
jan.hits<-con_prop_logit$node.all[as.logical(con_prop_logit$uncex.merged.mtdnas)]

#set up the color matrix for janus hits

#datacols<-matrix(tmp, ncol=4)
colfn<-function(data){
  tmp<<-tail(data,length(data)-1)
  #print(length(tmp))
  return(tail(data,length(data)-1)[1])
}

s <- 2
par(xpd = TRUE)

par(lwd=0.5)
#add heatmap for regime 1
x<-getcoords(simmap.janus.all.aggregate, node=jan.hits[10])[1]
y<-getcoords(simmap.janus.all.aggregate, node=jan.hits[10])[2]
tmp<-c("tmp","white","white","red","white")
for (h in s * (c(-.5, .5)-0.5))
  for (v in s * (c(-.5, .5)-0.5))
    rect(x + h, y + v, x + h + s, y + v + s, col=colfn(tmp))

#add heatmap for regime 2
x<-getcoords(simmap.janus.all.aggregate, node=jan.hits[12])[1]
y<-getcoords(simmap.janus.all.aggregate, node=jan.hits[12])[2]
tmp<-c("tmp","red","white","white","red")
for (h in s * (c(-.5, .5)-0.5))
  for (v in s * (c(-.5, .5)-0.5))
    rect(x + h, y + v, x + h + s, y + v + s, col=colfn(tmp))

#add heatmap for regime 3
x<-getcoords(simmap.janus.all.aggregate, node=jan.hits[11])[1]
y<-getcoords(simmap.janus.all.aggregate, node=jan.hits[11])[2]
tmp<-c("tmp","red","red","white","white")
for (h in s * (c(-.5, .5)-0.5))
  for (v in s * (c(-.5, .5)-0.5))
    rect(x + h, y + v, x + h + s, y + v + s, col=colfn(tmp))

#add heatmap for regime 4
x<-getcoords(simmap.janus.all.aggregate, node=jan.hits[1])[1]
y<-getcoords(simmap.janus.all.aggregate, node=jan.hits[1])[2]
tmp<-c("tmp","white","white","white", "red")
for (h in s * (c(-.5, .5)-0.5))
  for (v in s * (c(-.5, .5)-0.5))
    rect(x + h, y + v, x + h + s, y + v + s, col=colfn(tmp))

#add heatmap for regime 5
x<-getcoords(simmap.janus.all.aggregate, node=jan.hits[9])[1]
y<-getcoords(simmap.janus.all.aggregate, node=jan.hits[9])[2]
tmp<-c("tmp","white","red","white","white")
for (h in s * (c(-.5, .5)-0.5))
  for (v in s * (c(-.5, .5)-0.5))
    rect(x + h, y + v, x + h + s, y + v + s, col=colfn(tmp))

#add heatmap for regime 6
x<-getcoords(simmap.janus.all.aggregate, node=jan.hits[2])[1]
y<-getcoords(simmap.janus.all.aggregate, node=jan.hits[2])[2]
tmp<-c("tmp","red","white","white","white")
for (h in s * (c(-.5, .5)-0.5))
  for (v in s * (c(-.5, .5)-0.5))
    rect(x + h, y + v, x + h + s, y + v + s, col=colfn(tmp))

#add heatmap for regime 7
x<-getcoords(simmap.janus.all.aggregate, node=jan.hits[8])[1]
y<-getcoords(simmap.janus.all.aggregate, node=jan.hits[8])[2]
tmp<-c("tmp","white","red","white","white")
for (h in s * (c(-.5, .5)-0.5))
  for (v in s * (c(-.5, .5)-0.5))
    rect(x + h, y + v, x + h + s, y + v + s, col=colfn(tmp))

#add heatmap for regime 8
x<-getcoords(simmap.janus.all.aggregate, node=jan.hits[3])[1]
y<-getcoords(simmap.janus.all.aggregate, node=jan.hits[3])[2]
tmp<-c("tmp","white","red","red","white")
for (h in s * (c(-.5, .5)-0.5))
  for (v in s * (c(-.5, .5)-0.5))
    rect(x + h, y + v, x + h + s, y + v + s, col=colfn(tmp))

#add heatmap for regime 9
x<-getcoords(simmap.janus.all.aggregate, node=jan.hits[7])[1]
y<-getcoords(simmap.janus.all.aggregate, node=jan.hits[7])[2]
tmp<-c("tmp","red","white","white","white")
for (h in s * (c(-.5, .5)-0.5))
  for (v in s * (c(-.5, .5)-0.5))
    rect(x + h, y + v, x + h + s, y + v + s, col=colfn(tmp))

#add heatmap for regime 10
x<-getcoords(simmap.janus.all.aggregate, node=jan.hits[6])[1]
y<-getcoords(simmap.janus.all.aggregate, node=jan.hits[6])[2]
tmp<-c("tmp","red","red","red","white")
for (h in s * (c(-.5, .5)-0.5))
  for (v in s * (c(-.5, .5)-0.5))
    rect(x + h, y + v, x + h + s, y + v + s, col=colfn(tmp))

#add heatmap for regime 11
x<-getcoords(simmap.janus.all.aggregate, node=jan.hits[5])[1]
y<-getcoords(simmap.janus.all.aggregate, node=jan.hits[5])[2]
tmp<-c("tmp","red","white","white","white")
for (h in s * (c(-.5, .5)-0.5))
  for (v in s * (c(-.5, .5)-0.5))
    rect(x + h, y + v, x + h + s, y + v + s, col=colfn(tmp))

#add heatmap for regime 12
x<-getcoords(simmap.janus.all.aggregate, node=jan.hits[4])[1]
y<-getcoords(simmap.janus.all.aggregate, node=jan.hits[4])[2]
tmp<-c("tmp","red","white","white","white")
for (h in s * (c(-.5, .5)-0.5))
  for (v in s * (c(-.5, .5)-0.5))
    rect(x + h, y + v, x + h + s, y + v + s, col=colfn(tmp))
}

#add legend 
{
x<-8
y<-16
s=4
tmp<-c("tmp","red","white","white","white")
for (h in s * (c(-.5, .5)-0.5))
  for (v in s * (c(-.5, .5)-0.5))
    rect(x + h, y + v, x + h + s, y + v + s, col=colfn(tmp))

text(x= 12, y=22, "mtDNA", cex=0.45, font=2)
text(x= 4, y=22, "Introns", cex=0.45, font=2)
text(x= 12, y=9.5, "UTRs", cex=0.45, font=2)
text(x= 4, y=9.5, "Exons", cex=0.45, font=2)
}

#add node markers
{
mapplots:::add.pie(z=c(35,65), x=7.8, y=35, radius = 5, 
                   labels = c("FP", "P"), col = c("grey", make.transparent("white", 0.0)),
                    cex=0.5, label.dist = 0.35)

edgelabels(edge = c(380, 387, 379, 189, 121, 68, 98, 4, 287, 3, 326, 1), pie=l1ou.support.aicc, piecol=c("grey", "white"), cex=1.25, adj=-3)
edgelabels(edge = c(380, 387, 379, 189, 121, 68, 98, 4, 287, 3, 326, 1), pie=l1ou.support.pbic, piecol=c("grey", "white"), cex=1.25, adj=-10.5)


##posterior probabilities for allometries
pp<-l1ou.support.aicc
pp$null<-NULL
pp$alt<-NULL
pp$pp<-c(0.17, 0.12, 0, 0.10, 0.99, 0, 0.14, 0.37, 0.15, 0, 0.45, 0)
pp$pp_cols<-c(0.17, 0.12, 0.1, 0.10, 0.99, 0.1, 0.14, 0.37, 0.15, 0.1, 0.45, 0.1)
pp$alt<-1-pp$pp
#pp$col<-hcl.colors(100, rev=T)[pp$pp_cols*100]#"black"
#pp$col<-heat.colors(100, rev=T)[pp$pp_cols*100]#"black"
pp$col<-"black"

pp$col[c(3, 6, 10, 12)]<-"#00000000"

#edgelabels(edge = c(380, 387, 379, 189, 121, 68, 98, 4, 287, 3, 326, 1), pie=pp, piecol=c("red", "grey"), cex=0.9, adj=-3)
pp.edges<-c(380, 387, 379, 189, 121, 68, 98, 4, 287, 3, 326, 1)
for(i in 1:length(pp.edges)){
  par(new=T)
  #par(lwd=2)
  par(lwd=(pp$pp[i]+0.1)*2)
  edgelabels(edge = pp.edges[i], pch=1, cex=2, adj=-3, col=pp$col[i])
}


}
# 
# plot(l1ou_test$tree, cex=0.001,no.margin=T)
# edgelabels(edge = c(13,14,15,199,265, 370, 374, 387, 388, 389, 390, 393))

#edgelabels(edge = c(387, 380, 379, 1, 326, 287, 3, 4, 189, 121, 98, 68), text=as.character(c(2, 3, 1, 4, 5, 7, 6, 8, 9, 10, 11, 12)), cex=0.5, srt=270, frame="none", font=2, adj=c(0.6,0.2), col=make.transparent('black', 0.9))
#lapply(locator(), function(x) x + c(-1,1))

}
dev.off()

#code for plotting references
# 
# pdf(file="exons_ref.pdf", width=5, height=15)
# plot(simmap.janus.nuc.exons, ftype="off")
# tiplabels(getStates(simmap.janus.nuc.exons, type="tips"), frame='none', cex=0.5)
# dev.off()
# pdf(file="introns_ref.pdf", width=5, height=15)
# plot(simmap.janus.nuc.introns, ftype="off")
# tiplabels(getStates(simmap.janus.nuc.introns, type="tips"), frame='none', cex=0.5)
# dev.off()
# pdf(file="utrs_ref.pdf", width=5, height=15)
# plot(simmap.janus.nuc.utrs, ftype="off")
# tiplabels(getStates(simmap.janus.nuc.utrs, type="tips"), frame='none', cex=0.5)
# dev.off()
# pdf(file="mtdna_ref.pdf", width=5, height=15)
# plot(simmap.janus.mtdna.all, ftype="off")
# tiplabels(getStates(simmap.janus.mtdna.all, type="tips"), frame='none', cex=0.5)
# dev.off()

#generates warnings but works


#Section 16
#########################################################
###more testing with bmr/mass allometry or curvature ####
#########################################################
{
#create bmr dataset for plotting the rjMCMC results
nomiss<-nuc.aggregate.mvBM.imputed$estimates
#saveRDS(nomiss, file="./RDS/nomiss.RDS")
nomiss<-readRDS(file="./RDS/nomiss.RDS")


# #what is the variance in bmrs within each regime
bmr_vars<-data.frame(state=(getStates(simmap.janus.all.aggregate, type="tips")), nuc.aggregate.mvBM.imputed$estimates)
bmr_vars$mass2<-bmr_vars$mass^2
# bmr_vars$mass <- NULL
# bmr_vars$bmr <- (bmr_vars$bmr)
# bmr_vars <- bmr_vars %>% group_by(state) %>% summarise_at(vars(bmr), list(name = var))
# bmr_vars$name <- exp(bmr_vars$name)
# mean(bmr_vars$name, na.rm=T)
# sd(bmr_vars$name, na.rm=T)

dev.off()


#lin<-lm(bmr~mass, data=bmr_vars)
require(nlme)
lin<-gls(bmr ~ mass, correlation = corMartins(1, phy = simmap.janus.all.aggregate),
         data = bmr_vars, method = "ML")
#quad<-lm(bmr~mass+mass2, data=bmr_vars)
quad<-gls(bmr ~ mass + mass2, correlation = corMartins(1, phy = simmap.janus.all.aggregate),
          data = bmr_vars, method = "ML")

aic.w(c(AIC(quad), AIC(lin)))

####
lin.pgls <- phylolm::phylolm(bmr ~ mass, phy = simmap.janus.all.aggregate, data = bmr_vars, model = "OUfixedRoot", upper.bound=10)
quad.pgls <- phylolm::phylolm(bmr ~ mass+mass2, phy = simmap.janus.all.aggregate, data = bmr_vars, model = "OUfixedRoot", upper.bound=10)


AIC(lin)
AIC(quad) ## indicates curvature model is a better fit

plot(bmr~mass, data=bmr_vars)
abline(lin, col='red')
curve(quad$coefficient[3]*x^2 + quad$coefficient[2]*x + quad$coefficient[1], add=T, col="blue")

AIC(lin)
AIC(quad)
BIC(lin)
BIC(quad)

# #checking the rr2 values
# # This also works for models fit with nlme::gls()
# z.f.gls <- nlme::gls(y_pgls ~ x_trait, data = d, correlation = ape::corPagel(1, phy), method = "ML")
# z.x.gls <- nlme::gls(y_pgls ~ 1, data = d, correlation = ape::corPagel(1, phy), method = "ML")
# R2(mod = z.f.gls, mod.r = z.v.lm)
# ##    R2_lik  R2_resid   R2_pred 
# ## 0.3826794 0.4854591 0.4599150
# R2(mod = z.f.gls)



}


#reference plot for l1ou
require(ape)
plot(
l1ou_test$tree, cex=0.1, no.margin=T)
edgelabels(edge= c(13,14,199,265,374,387,388, 390))
# 
# root
# 13
# 14
# 199
# 265
# 374
# 387
# 388
# 390


rootstate<-l1ou_test$tree$tip.label[
  !l1ou_test$tree$tip.label %in% 
    c(tips_from_edge(tree= l1ou_test$tree, edge=13),
      tips_from_edge(tree= l1ou_test$tree, edge=14),
      tips_from_edge(tree= l1ou_test$tree, edge=199),
      tips_from_edge(tree= l1ou_test$tree, edge=265),
      tips_from_edge(tree= l1ou_test$tree, edge=374),
      tips_from_edge(tree= l1ou_test$tree, edge=387),
      tips_from_edge(tree= l1ou_test$tree, edge=388),
      tips_from_edge(tree= l1ou_test$tree, edge=390))]

nested<-tips_from_edge(tree= l1ou_test$tree, edge=387)[
  !tips_from_edge(tree= l1ou_test$tree, edge=387) %in% 
    c(tips_from_edge(tree= l1ou_test$tree, edge=199),
      tips_from_edge(tree= l1ou_test$tree, edge=265),
      tips_from_edge(tree= l1ou_test$tree, edge=374))]



#experimenting with plotting bmr/mass for each section (not used)
#root to 13
{plot(nomiss, pch=21, bg=make.transparent('grey',0.5), lwd=0.5, bty='n', xlim=c(0, 17), ylim=c(-4.5, 12))
x.lim<-range(nomiss[,1])
y.lim<-range(nomiss[,2])

#root
points(nomiss[rootstate,], pch=21, bg=make.transparent('black',0.5), lwd=0.5, bty='n', xlim=x.lim, ylim=y.lim)
curve(expr=-3.397 + x*	0.690, add=T, col="black")

#root to 13
plot(nomiss[rootstate,], pch=21, bg=make.transparent('grey',0.5), lwd=0.5, bty='n', xlim=x.lim, ylim=y.lim)
curve(expr=-3.397 + x*	0.690, add=T, col="grey")
points(nomiss[tips_from_edge(tree= l1ou_test$tree, edge=13),], pch=21, bg=make.transparent('red',0.5), lwd=0.5)
curve(expr=-4.281 + x*	0.763, add=T, col=rcartocolor:::carto_pal(n = 12, "Safe")[1])

#root to 14
plot(nomiss[rootstate,], pch=21, bg=make.transparent('grey',0.5), lwd=0.5, bty='n', xlim=x.lim, ylim=y.lim)
curve(expr=-3.397 + x*	0.690, add=T, col="grey")
points(nomiss[tips_from_edge(tree= l1ou_test$tree, edge=14),], pch=21, bg=make.transparent('red',0.5), lwd=0.5)
curve(expr=-3.877 + x*	0.756, add=T, col=rcartocolor:::carto_pal(n = 12, "Safe")[2])

#root to 390
plot(nomiss[rootstate,], pch=21, bg=make.transparent('grey',0.5), lwd=0.5, bty='n', xlim=x.lim, ylim=y.lim)
curve(expr=-3.397 + x*	0.690, add=T, col="grey")
points(nomiss[tips_from_edge(tree= l1ou_test$tree, edge=390),], pch=21, bg=make.transparent('red',0.5), lwd=0.5)
curve(expr=-4.261 + x*	0.837, add=T, col=rcartocolor:::carto_pal(n = 12, "Safe")[3])

#root to 388
plot(nomiss[rootstate,], pch=21, bg=make.transparent('grey',0.5), lwd=0.5, bty='n', xlim=x.lim, ylim=y.lim)
curve(expr=-3.397 + x*	0.690, add=T, col="grey")
points(nomiss[tips_from_edge(tree= l1ou_test$tree, edge=388),], pch=21, bg=make.transparent('red',0.5), lwd=0.5)
curve(expr=-3.462 + x*	0.681, add=T, col=rcartocolor:::carto_pal(n = 12, "Safe")[4])

#root to 387
plot(nomiss[rootstate,], pch=21, bg=make.transparent('grey',0.5), lwd=0.5, bty='n', xlim=x.lim, ylim=y.lim)
curve(expr=-3.397 + x*	0.690, add=T, col="grey")
points(nomiss[tips_from_edge(tree= l1ou_test$tree, edge=387),], pch=21, bg=make.transparent('red',0.5), lwd=0.5)
curve(expr=-3.260 + x*	0.682, add=T, col=rcartocolor:::carto_pal(n = 12, "Safe")[5])

#387 to 199
plot(nomiss[nested,], pch=21, bg=make.transparent('grey',0.5), lwd=0.5, bty='n', xlim=x.lim, ylim=y.lim)
curve(expr=-3.260 + x*	0.682, add=T, col="grey")
points(nomiss[tips_from_edge(tree= l1ou_test$tree, edge=199),], pch=21, bg=make.transparent('red',0.5), lwd=0.5)
curve(expr=-3.510 + x*	0.722, add=T, col=rcartocolor:::carto_pal(n = 12, "Safe")[6])

#387 to 265
plot(nomiss[nested,], pch=21, bg=make.transparent('grey',0.5), lwd=0.5, bty='n', xlim=x.lim, ylim=y.lim)
curve(expr=-3.260 + x*	0.682, add=T, col="grey")
points(nomiss[tips_from_edge(tree= l1ou_test$tree, edge=265),], pch=21, bg=make.transparent('red',0.5), lwd=0.5)
curve(expr=-3.866 + x*	0.748, add=T, col=rcartocolor:::carto_pal(n = 12, "Safe")[7])

#387 to 374
plot(nomiss[nested,], pch=21, bg=make.transparent('grey',0.5), lwd=0.5, bty='n', xlim=x.lim, ylim=y.lim)
curve(expr=-3.397 + x*	0.690, add=T, col="grey")
points(nomiss[tips_from_edge(tree= l1ou_test$tree, edge=374),], pch=21, bg=make.transparent('red',0.5), lwd=0.5)
curve(expr=-3.841 + x*	0.748, add=T, col=rcartocolor:::carto_pal(n = 12, "Safe")[8])

}


#plotting as segments
{plot(nomiss, pch=21, bg=make.transparent('grey',0.001), lwd=0.5, bty='n')

#root
points(nomiss[rootstate,], pch=21, bg=make.transparent('red',0.5), lwd=0.5, bty='n', xlim=x.lim, ylim=y.lim)
curve(expr=-3.397 + x*	0.690, add=T, col="red", from=min(nomiss[rootstate,][,1]), to=max(nomiss[rootstate,][,1]))

#root to 13
points(nomiss[tips_from_edge(tree= l1ou_test$tree, edge=13),], pch=21, bg=make.transparent('red',0.5), lwd=0.5)
curve(expr=-4.281 + x*	0.763, add=T, col="red", from=min(nomiss[tips_from_edge(tree= l1ou_test$tree, edge=13),][,1]), to=max(nomiss[tips_from_edge(tree= l1ou_test$tree, edge=13),][,1]))

#root to 14
points(nomiss[tips_from_edge(tree= l1ou_test$tree, edge=14),], pch=21, bg=make.transparent('red',0.5), lwd=0.5)
curve(expr=-3.877 + x*	0.756, add=T, col="red", from=min(nomiss[tips_from_edge(tree= l1ou_test$tree, edge=14),][,1]), to=max(nomiss[tips_from_edge(tree= l1ou_test$tree, edge=14),][,1]))

#root to 390
points(nomiss[tips_from_edge(tree= l1ou_test$tree, edge=390),], pch=21, bg=make.transparent('red',0.5), lwd=0.5)
curve(expr=-4.261 + x*	0.837, add=T, col="red", from=min(nomiss[tips_from_edge(tree= l1ou_test$tree, edge=390),][,1]), to=max(nomiss[tips_from_edge(tree= l1ou_test$tree, edge=390),][,1]))

#root to 388
points(nomiss[tips_from_edge(tree= l1ou_test$tree, edge=388),], pch=21, bg=make.transparent('red',0.5), lwd=0.5)
curve(expr=-3.462 + x*	0.681, add=T, col="red", from=min(nomiss[tips_from_edge(tree= l1ou_test$tree, edge=388),][,1]), to=max(nomiss[tips_from_edge(tree= l1ou_test$tree, edge=388),][,1]))

#root to 387
points(nomiss[tips_from_edge(tree= l1ou_test$tree, edge=387),], pch=21, bg=make.transparent('red',0.5), lwd=0.5)
curve(expr=-3.260 + x*	0.682, add=T, col="red", from=min(nomiss[tips_from_edge(tree= l1ou_test$tree, edge=387),][,1]), to=max(nomiss[tips_from_edge(tree= l1ou_test$tree, edge=387),][,1]))

#387 to 199
points(nomiss[tips_from_edge(tree= l1ou_test$tree, edge=199),], pch=21, bg=make.transparent('red',0.5), lwd=0.5)
curve(expr=-3.510 + x*	0.722, add=T, col="red", from=min(nomiss[tips_from_edge(tree= l1ou_test$tree, edge=199),][,1]), to=max(nomiss[tips_from_edge(tree= l1ou_test$tree, edge=199),][,1]))

#387 to 265
points(nomiss[tips_from_edge(tree= l1ou_test$tree, edge=265),], pch=21, bg=make.transparent('red',0.5), lwd=0.5)
curve(expr=-3.866 + x*	0.748, add=T, col="red", from=min(nomiss[tips_from_edge(tree= l1ou_test$tree, edge=265),][,1]), to=max(nomiss[tips_from_edge(tree= l1ou_test$tree, edge=265),][,1]))

#387 to 374
points(nomiss[tips_from_edge(tree= l1ou_test$tree, edge=374),], pch=21, bg=make.transparent('red',0.5), lwd=0.5)
curve(expr=-3.841 + x*	0.748, add=T, col="red", from=min(nomiss[tips_from_edge(tree= l1ou_test$tree, edge=374),][,1]), to=max(nomiss[tips_from_edge(tree= l1ou_test$tree, edge=374),][,1]))

}


#### experimenting with l1ou plots #### 
# {
# nEdges <- Nedge(l1ou_test$tree) # total number of edges
# ew <- rep(1,nEdges)  # to set default edge width of 1
# ew[eModel.par$shift.configuration] <- 3   # to widen edges with a shift 
# #edge.width=ew
# par(mfrow=c(1,2))
# #plot(eModel.par.unconstrained.aicc, cex=0.001, label.offset=0.02, asterisk=F, edge.shift.ann=F, plot.bar=F)
# plot(eModel.par.constrained.aicc, cex=0.001, label.offset=0.02, asterisk=F, edge.shift.ann=F, plot.bar=F)
# edges <- edge_indices_nodes(tree=l1ou_test$tree, nodes=con_prop_logit[con_prop_logit$uncex.merged.mtdnas==1,]$node.all)
# edgelabels(edge=edges)
# 
# 
# edgelabels(edge = edges[-3] , pie=fit.boot.aicc.constrained$detection.rate[edges[-3]], cex=0.4, adj = c(0.49,2))
# edgelabels(edge = edges[-3] , pie=fit.boot.aicc.constrained.mod$detection.rate[edges[-3]], cex=c(0.4), adj = c(0.51, 2))
# #edgelabels(edge = edges[-3] , pie=fit.boot.aicc.unconstrained$detection.rate[edges[-3]], cex=0.4, adj = c(0.49,-2))
# #edgelabels(edge = edges[-3] , pie=fit.boot.aicc.unconstrained.mod$detection.rate[edges[-3]], cex=c(0.4), adj = c(0.51, -2))
# 
# ##
# edgelabels(edge = edges[-3] , pie=fit.boot.aicc.constrained$detection.rate[edges[-3]], cex=0.4, adj = c(0.54,2))
# edgelabels(edge = edges[-3] , pie=fit.boot.pbic.constrained$detection.rate[edges[-3]], cex=0.4, adj = c(0.56,2))
# #edgelabels(edge = edges[-3] , pie=fit.boot.aicc.unconstrained$detection.rate[edges[-3]], cex=c(0.4), adj = c(0.54, -2))
# #edgelabels(edge = edges[-3] , pie=fit.boot.pbic.unconstrained$detection.rate[edges[-3]], cex=c(0.4), adj = c(0.56, -2))
# 
# 
# #plotting directly on the model shift tree
# plot(simmap.janus.all.aggregate, ftype="off", colors=setNames(c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928', '#8dd3c7'), as.character(seq(0,12))))
# 
# edges1 <- edge_indices_nodes(tree=l1ou_test$tree, nodes=con_prop_logit[con_prop_logit$uncex.merged.mtdnas==1,]$node.all)
# edges2 <- edge_indices_nodes(tree=simmap.janus.all.aggregate, nodes=con_prop_logit[con_prop_logit$uncex.merged.mtdnas==1,]$node.all)
# 
# edgelabels(edge = edges2 , pie=fit.boot.aicc.constrained$detection.rate[edges1], cex=0.4, adj = c(-1,2.5))
# edgelabels(edge = edges2 , pie=fit.boot.aicc.constrained.mod$detection.rate[edges1], cex=c(0.4), adj = c(0.8, 2.5))
# edgelabels(edge = edges2 , pie=fit.boot.aicc.unconstrained$detection.rate[edges1], cex=0.4, adj = c(-1,-1.9))
# edgelabels(edge = edges2 , pie=fit.boot.aicc.unconstrained.mod$detection.rate[edges1], cex=c(0.4), adj = c(0.8, -1.9))
# 
# }
# 



#alternative reference plots

par(mfrow=c(1,8))
#plot(simmap.janus.nuc.alldata, ftype="off")
plot(simmap.janus.nuc.exons, ftype="off")
plot(simmap.janus.nuc.introns, ftype="off")
plot(simmap.janus.nuc.utrs, ftype="off")
plot(simmap.janus.mtdna.all, ftype="off")
plot(simmap.janus.mtdna.proteins, ftype="off")
plot(simmap.janus.mtdna.rrnas, ftype="off")
plot(simmap.janus.nuc.aggregate, ftype="off", colors=setNames(c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928'), as.character(seq(0,11))))
plot(lwd=4,direction="upwards",simmap.janus.all.aggregate, ftype="off", colors=setNames(c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928', '#8dd3c7'), as.character(seq(0,12))))

dev.off()


# 
# #data for sean
# tree<-l1ou_test$tree
# edge_indices_nodes(tree=l1ou_test$tree, nodes=con_prop_logit[con_prop_logit$uncex.merged.mtdnas==1,]$node.all)

#uncex.merged.mtdnas (the most complex shift configuration, intersection of all analyses)
#edge IDs
#13  14  15 199 265 370 374 387 388 389 390 393

#uncex.merged (intersection of nuclear datasets)
#edge IDs
#13  14  15 199 265 370 374 387 388 389 390

#uncex.allnucdata (all nuclear data concatenated)
#edgeIDs
#14 199 265 370 374 389

#uncex.exons (all exons concatenated)
#edge IDs
#13  14 199 265 370 374 389

#uncex.introns (all introns concatenated)
#edge IDs
#13 265 387 388 390

#uncex.utrs (all utrs concatenated)
#edge IDs
#15 265 387

#uncex.mtDNAs.all (all mtdna data concatenated)
#edge IDs
#14 393

#mtDNAs.proteins, mtdna.rrnas (mtdna rrna or protein separately, one shift)
#edge IDs
#393


#plot(l1ou_test$tree, cex=0.001)
#edgelabels(edge=edges)


##testing patterns of life history integration with mvMORPH/RPANDA
#not in the paper, just experimenting

{
require(RPANDA)
require(mvMORPH)

LHT.mvgls.bmm <-
  mvgls(l1ou_test$Y ~ 1,
        tree = simmap.janus.all.aggregate,
        model = 'BMM',
        error = T)

LHT.mvgls.OU <-
  mvgls(l1ou_test$Y ~ 1,
        tree = simmap.janus.all.aggregate,
        model = 'OU',
        error = T)

grps<-factor(getStates(simmap.janus.all.aggregate, type='tips'))
dat.grp<-data.frame(l1ou_test$Y, groups=grps)
LHT.mvgls.bmm.grp <-
  mvgls(as.matrix(dat.grp[,1:8]) ~ groups, data=dat.grp,
        tree = simmap.janus.all.aggregate,
        model = 'BMM',
        error = T, method='LL')

manova.gls(LHT.mvgls.bmm.grp, type='III')

glht.test.p<-pairwise.glh(LHT.mvgls.bmm.grp, verbose=F, term='groups', adjust='none', test='Pillai')
glht.test.w<-pairwise.glh(LHT.mvgls.bmm.grp, verbose=F, term='groups', adjust='none', test='Wilks')
create_distance_matrix(glht.test.p, param="pvalue")


###
grps<-factor(getStates(simmap.janus.nuc.aggregate.simplified, type='tips'))
dat.grp<-data.frame(l1ou_test$Y, groups=grps)
LHT.mvgls.bmm.grp <-
  mvgls(as.matrix(dat.grp[,1:8]) ~ groups, data=dat.grp,
        tree = simmap.janus.nuc.aggregate.simplified,
        model = 'BM',
        error = T, method='LL')

manova.gls(LHT.mvgls.bmm.grp)

glht.test.p<-pairwise.glh(LHT.mvgls.bmm.grp, verbose=F, term='groups', adjust='none', test='Pillai')
glht.test.w<-pairwise.glh(LHT.mvgls.bmm.grp, verbose=F, term='groups', adjust='none', test='Wilks')

pval.dist.p<-as.dist(create_distance_matrix(glht.test.p, param="pvalue"))
stat.dist.p<-as.dist(create_distance_matrix(glht.test.p, param="stat"))
pval.dist.w<-as.dist(create_distance_matrix(glht.test.w, param="pvalue"))
stat.dist.w<-as.dist(create_distance_matrix(glht.test.w, param="stat"))


require(qgraph)
qgraph(stat.dist.p, vsize=10, minimum=0.1)

qgraph(stat.dist.w, vsize=10)
qgraph(pval.dist.p, vsize=10)

##testing scOU models for theta
errormat<-as.matrix(set_all_values(l1ou_test$Y, 0.1^2))

LHT.mvBM.bmm <-
  mvBM(data = l1ou_test$Y,
       tree = simmap.janus.all.aggregate,
       model = 'BMM', error=errormat)


#fit the scOU model (eg. PhyloEM)
#input data in this data frame are log transformed, scaled, and then imputed under mvBM
#this means that the variances are not all 1 after imputation...hmm
LHT.mvOU.OUM <-
  mvOU(
    data = l1ou_test$Y,
    tree = simmap.janus.nuc.aggregate.simplified2,
    model = 'OUM',
    param = list(decomp = "equaldiagonal", root = "stationary")
  )

#saveRDS(LHT.mvOU.OUM, file='./RDS/LHT.mvOU.OUM.RDS')
LHT.mvOU.OUM<-readRDS('./RDS/LHT.mvOU.OUM.RDS')

qgraph(LHT.mvOU.OUM$sigma, graph='cor')
qgraph(LHT.mvOU.OUM$sigma, graph='pcor')

LHT.mvOU.OUM.error <-
  mvOU(
    data = l1ou_test$Y,
    tree = simmap.janus.nuc.aggregate.simplified2,
    model = 'OUM',
    param = list(decomp = "equaldiagonal", root = "stationary"),
    error=errormat
  )

#saveRDS(LHT.mvOU.OUM.error, file='./RDS/LHT.mvOU.OUM.error.RDS')
LHT.mvOU.OUM.error<-readRDS('./RDS/LHT.mvOU.OUM.error.RDS')

phylopars.data.mvOU<-phylopars.data[,c("species", "mass", "mean.clutch", "gen_length", "survival", "breeding","longevity", "chickPC1", "latitude")]
rownames(phylopars.data.mvOU)<-phylopars.data.mvOU$species
phylopars.data.mvOU$species<-NULL

LHT.mvOU.OUM.missing<-
  mvOU(
    data = phylopars.data.mvOU ,
    tree = simmap.janus.nuc.aggregate.simplified2,
    model = 'OUM',
    param = list(decomp = "equaldiagonal", root = "stationary")
  )

#saveRDS(LHT.mvOU.OUM.missing, file='./RDS/LHT.mvOU.OUM.missing.RDS')
LHT.mvOU.OUM.missing<-readRDS(file='./RDS/LHT.mvOU.OUM.missing.RDS')

LHT.mvOU.OUM.missing.error<-
  mvOU(
    data = phylopars.data.mvOU,
    tree = simmap.janus.nuc.aggregate.simplified2,
    model = 'OUM',
    param = list(decomp = "equaldiagonal", root = "stationary"),
    error=errormat
  )

#saveRDS(LHT.mvOU.OUM.missing.error, file='./RDS/LHT.mvOU.OUM.missing.error.RDS')
LHT.mvOU.OUM.missing.error<-readRDS(file='./RDS/LHT.mvOU.OUM.missing.error.RDS')

thetas<-do.call(cbind,
lapply(colnames(LHT.mvOU.OUM.missing.error$theta), function(x) get_dup_thetas(tree=simmap.janus.nuc.aggregate.simplified2, fitted = LHT.mvOU.OUM.missing.error, trait=x)),
)
colnames(thetas)<-colnames(LHT.mvOU.OUM.missing.error$theta)

thetas<-thetas[,rev(c('chickPC1', 'mass', 'mean.clutch', 'longevity', 'breeding', 'latitude', 'gen_length','survival'))]
thetas.norm<-thetas/max(thetas,na.rm=TRUE)
thetas.zeromin<-apply(thetas, 2, function(x) x-min(x))
thetas.zeromin.norm<-thetas.zeromin/max(thetas.zeromin,na.rm=TRUE)

thetas.pca<-phyl.pca(tree=simmap.janus.nuc.aggregate.simplified2, Y = thetas)
thetas.pca.zeromin<-apply(thetas.pca$S, 2, function(x) x-min(x))
thetas.norm<-thetas/max(thetas,na.rm=TRUE)

require(phytools)
require(plotrix)


pdf(file='test.pdf', width=10, height=10)
plot_phylo_with_bars(tree=simmap.janus.nuc.aggregate.simplified2, thetas.zeromin, mapcolors=setNames(rainbow(11), unique(getStates(simmap.janus.nuc.aggregate.simplified2, type = 'tips'))), part=0.5, colors=rev(RColorBrewer::brewer.pal(n=8,"Spectral")))
dev.off()

plot_phylo_with_bars(tree=simmap.janus.nuc.aggregate.simplified2, thetas.pca.zeromin[,c("PC1", "PC2")], mapcolors=setNames(rainbow(11), unique(getStates(simmap.janus.nuc.aggregate.simplified2, type = 'tips'))), part=0.5, colors=rev(RColorBrewer::brewer.pal(n=8,"Spectral")))
phylo.heatmap(tree = simmap.janus.nuc.aggregate.simplified2, X = thetas.zeromin, standardize=F, legend=F, ftype='off')
}