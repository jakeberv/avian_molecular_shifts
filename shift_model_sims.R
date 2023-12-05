#script for generating input templates for iqtree to simulate shift models
#collab between Jacob Berv and ChatGPT4

library(ape)
library(phytools)

# Function to identify edge indices based on a minimum clade size
edge_indices_N <- function(tree, min_clade_size) {
  pruned <- tree
  theNodes <- length(pruned$tip.label) + 1:pruned$Nnode
  results <- numeric(length(theNodes))
  for(i in seq_along(theNodes)){
    temp <- extract.clade(pruned, theNodes[i])
    results[i] <- length(temp$tip.label)
  }
  names(results) <- theNodes
  edgetable <- pruned$edge
  rownames(edgetable) <- seq_len(nrow(edgetable))
  eligible_edges <- as.numeric(rownames(edgetable)[as.numeric(edgetable[,2]) %in% names(results[results >= min_clade_size])])
  cat("Eligible edge indices for annotation based on minimum clade size:\n", eligible_edges, "\n\n")
  return(eligible_edges)
}

# Function to convert edge indices to node indices
node_indices_edge <- function(tree, edge_indices) {
  edgetable <- tree$edge
  node_indices <- edgetable[edge_indices, 2]
  cat("Node indices corresponding to the eligible edges:\n", node_indices, "\n\n")
  return(node_indices)
}

# Function to randomly select nodes for annotation
select_nodes_for_annotation <- function(candidate_nodes, reps) {
  if (length(candidate_nodes) < reps) {
    stop("The number of reps is greater than the number of candidate nodes.")
  }
  selected_nodes <- sample(candidate_nodes, reps, replace = FALSE)
  cat("Randomly selected nodes for annotation:\n", selected_nodes, "\n\n")
  return(selected_nodes)
}

# Function to annotate selected nodes with a given model
annotate_nodes <- function(tree, selected_nodes, model) {
  tree$node.label <- rep("", tree$Nnode)  # Initialize all node labels to empty
  for (node in selected_nodes) {
    # Ensure node index is corrected for the internal node label indexing
    corrected_index <- node - length(tree$tip.label)
    tree$node.label[corrected_index] <- model
    cat(sprintf("Node %d (internal index %d) annotated with model %s\n", node, corrected_index, model))
  }
  return(tree)
}

# Main function to annotate branches and write to Newick format
annotate_branches <- function(input_tree, model, reps, min_clade_size) {
  if (is.character(input_tree)) {
    tree <- read.tree(text = input_tree)
  } else if (inherits(input_tree, "phylo")) {
    tree <- input_tree
  } else {
    stop("Invalid input: input_tree must be a Newick string or a phylo object.")
  }
  
  # Identify candidate edges and nodes
  edge_indices <- edge_indices_N(tree, min_clade_size)
  candidate_nodes <- node_indices_edge(tree, edge_indices)
  print(candidate_nodes)
  # Select nodes and annotate them
  selected_nodes <- select_nodes_for_annotation(candidate_nodes, reps)
  annotated_tree <- annotate_nodes(tree, selected_nodes, model)
  
  # Use write.tree to output the annotated tree in Newick format
  newick_with_labels <- write.tree(annotated_tree, file = "", digits = 10)
  cat("Annotated Newick string with write.tree:\n", newick_with_labels, "\n\n")
  
  # Plot the tree with the selected nodes highlighted
  plot_selected_nodes(annotated_tree, selected_nodes)
  return(newick_with_labels)
}

# Function to plot the tree with selected nodes highlighted
plot_selected_nodes <- function(tree, selected_nodes) {
  plot.phylo(tree, show.node.label = TRUE, cex = 1, no.margin=T, show.tip.label=F)
  nodelabels(frame = "none", node = selected_nodes, bg = "red", cex = 0.8)
}


write_vector_to_files <- function(vector, prefix) {
  for (i in seq_along(vector)) {
    file_name <- sprintf("%s_%d.tre", prefix, i)
    writeLines(as.character(vector[i]), file_name)
  }
}

# Simulate a tree using pbtree from the ape package
set.seed(123) # Set a seed for reproducibility
tree <- pbtree(n = 100)
cat("Simulated tree in Newick format:\n", write.tree(tree), "\n\n")

# Example usage
min_clade_size <- 4
reps <- 1
model <- "[&model=HKY{2.0}+GC]"

# Run the annotation process
annotated_newick <- annotate_branches(tree, model, reps, min_clade_size)
cat("Final annotated Newick string:\n", annotated_newick, "\n")




#read in trees for preparation
exontree<-read.tree(file="/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/unmasked/min2x/qual20_cov2_haplo_bestonly/initial_test_filters/min50bp_min10p_aligned/ALIGNED/phased/sample1haplo/phyhetnucbf/backbone_constraint/janus/8952e31d/gamma/redo_m3/exons_MFP_MERGE_MRL3_constraint_G_UE_UL_M3/exons_MFP_MERGE_MRL3_constraint.rooted.treefile")

introntree<-read.tree(file="/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/unmasked/min2x/qual20_cov2_haplo_bestonly/initial_test_filters/min50bp_min10p_aligned/ALIGNED/phased/sample1haplo/phyhetnucbf/backbone_constraint/janus/8952e31d/gamma/redo_m3/introns_MFP_MERGE_MRL3_constraint_G_UE_UL_M3/introns_MFP_MERGE_MRL3_constraint.rooted.treefile")

utrtree<-read.tree(file="/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/unmasked/min2x/qual20_cov2_haplo_bestonly/initial_test_filters/min50bp_min10p_aligned/ALIGNED/phased/sample1haplo/phyhetnucbf/backbone_constraint/janus/8952e31d/gamma/redo_m3/utrs_MFP_MERGE_MRL3_constraint_G_UE_UL_M3/utrs_MFP_MERGE_MRL3_constraint.rooted.treefile")

mtdnatree<-read.tree(file="/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/mtdnas/janus/concat_all/mtDNA_all_MRL3_constraint.janus.tre")

set.seed(123)
# Run the annotation process
exon_base <- rep(annotate_branches(exontree, model, reps, min_clade_size), times=10)
intron_base <- rep(annotate_branches(introntree, model, reps, min_clade_size), times=10)
utr_base <- rep(annotate_branches(utrtree, model, reps, min_clade_size), times=10)
mtdna_base <-rep(annotate_branches(mtdnatree, model, reps, min_clade_size), times=10)

# Assuming your vector is named 'exon_base'
write_vector_to_files(exon_base, "exon")
write_vector_to_files(intron_base, "intron")
write_vector_to_files(utr_base, "utr")
write_vector_to_files(mtdna_base, "mtdna")

## 2 derived shifts
set.seed(123)
# Run the annotation process
exon_base_2 <- rep(annotate_branches(exontree, model, reps=2, min_clade_size), times=10)
intron_base_2 <- rep(annotate_branches(introntree, model, reps=2, min_clade_size), times=10)
utr_base_2 <- rep(annotate_branches(utrtree, model, reps=2, min_clade_size), times=10)
mtdna_base_2 <-rep(annotate_branches(mtdnatree, model, reps=2, min_clade_size), times=10)

# Assuming your vector is named 'exon_base'
write_vector_to_files(exon_base_2, "exon")
write_vector_to_files(intron_base_2, "intron")
write_vector_to_files(utr_base_2, "utr")
write_vector_to_files(mtdna_base_2, "mtdna")


#three derived shifts
set.seed(123)
# Run the annotation process
exon_base_3 <- rep(annotate_branches(exontree, model, reps=3, min_clade_size), times=10)
intron_base_3 <- rep(annotate_branches(introntree, model, reps=3, min_clade_size), times=10)
utr_base_3 <- rep(annotate_branches(utrtree, model, reps=3, min_clade_size), times=10)
mtdna_base_3 <-rep(annotate_branches(mtdnatree, model, reps=3, min_clade_size), times=10)

# Assuming your vector is named 'exon_base'
write_vector_to_files(exon_base_3, "exon")
write_vector_to_files(intron_base_3, "intron")
write_vector_to_files(utr_base_3, "utr")
write_vector_to_files(mtdna_base_3, "mtdna")

