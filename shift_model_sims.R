# Script for generating input templates for iqtree to simulate shift models
# Collaboration between Jacob Berv and ChatGPT4 with modifications for nested or independent shifts

library(ape)
library(phytools)
require(Biostrings)
require(phylotate)
library(DirichletReg)

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

# Function to randomly select nodes for annotation with option for nested or independent shifts
select_nodes_for_annotation <- function(tree, candidate_nodes, reps, nested = FALSE, buffer = 1) {
  if (length(candidate_nodes) < reps) {
    stop("The number of reps is greater than the number of candidate nodes.")
  }
  
  selected_nodes <- c()  # To keep track of nodes selected for shifts
  excluded_nodes <- c()  # For independent shifts
  first_selected_node <- NULL  # To store the first selected node in nested shifts
  nested_pool <- NULL    # To keep track of eligible nodes for nested shifts after the first selection
  
  cat("Starting shifts selection...\n")
  for (i in 1:reps) {
    if (nested) {
      if (i == 1) {
        # For the first node in nested mode, any node is eligible
        first_selected_node <- sample(candidate_nodes, 1, replace = FALSE)
        selected_nodes <- c(selected_nodes, first_selected_node)
        nested_pool <- union(get_descendants(tree, first_selected_node), get_ancestors(tree, first_selected_node))
        cat("First node selected, nested pool defined: ", nested_pool, "\n")
      } else {
        # Apply buffer and lineage constraints for subsequent selections in nested mode
        valid_candidates <- intersect(candidate_nodes, nested_pool)
        buffer_candidates <- vector()
        for (node in valid_candidates) {
          path <- nodepath(tree, from = first_selected_node, to = node)
          if (length(path) - 1 >= buffer) {
            buffer_candidates <- c(buffer_candidates, node)
          }
        }
        
        if (length(buffer_candidates) == 0) {
          stop("No more candidate nodes left to select for shifts considering the buffer.")
        }
        
        selected <- sample(buffer_candidates, 1, replace = FALSE)
        selected_nodes <- c(selected_nodes, selected)
        cat("Selected node for shift (nested): ", selected, "\n")
      }
    } else {
      # Independent shifts (original behavior)
      candidate_nodes <- setdiff(candidate_nodes, excluded_nodes)
      
      if (length(candidate_nodes) == 0) {
        stop("No more candidate nodes left to select for shifts.")
      }
      
      selected <- sample(candidate_nodes, 1, replace = FALSE)
      selected_nodes <- c(selected_nodes, selected)
      
      # Exclude the selected node and its descendants
      descendants_and_ancestors <- union(get_descendants(tree, selected), get_ancestors(tree, selected))
      excluded_nodes <- union(excluded_nodes, descendants_and_ancestors)
      cat("Selected node for shift (independent): ", selected, "\n")
    }
  }
  
  cat("Final randomly selected nodes for annotation: ", selected_nodes, "\n\n")
  return(selected_nodes)
}

# Helper function to find all ancestor nodes up to the root for a given node
get_ancestors <- function(tree, node) {
  ancestors <- c()
  current_node <- node
  while(current_node %in% tree$edge[,2]) {
    parent <- tree$edge[tree$edge[,2] == current_node, 1]
    if (length(parent) == 1) {
      ancestors <- c(ancestors, parent)
      current_node <- parent
    } else {
      break  # If we reach a node that has no parent, i.e., the root node.
    }
  }
  return(ancestors)
}

# Helper function to find all descendant nodes for a given node
get_descendants <- function(tree, node) {
  descendants <- c()
  nodes_to_visit <- node
  while(length(nodes_to_visit) > 0) {
    current_node <- nodes_to_visit[1]
    nodes_to_visit <- nodes_to_visit[-1]
    children <- tree$edge[tree$edge[,1] == current_node, 2]
    descendants <- c(descendants, children)
    nodes_to_visit <- c(nodes_to_visit, children)
  }
  return(descendants)
}

# Main function to annotate branches and write to Newick format with nested or independent option
annotate_branches <- function(input_tree, models, reps, min_clade_size, nested = FALSE, annotate_tips = FALSE, buffer = 1) {
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
  
  repeat {
    # Select nodes for annotation
    selected_nodes <- select_nodes_for_annotation(tree, candidate_nodes, reps, nested, buffer)
    annotated_tree <- annotate_nodes(tree, selected_nodes, models, nested)
    newick_with_labels <- write.tree(annotated_tree, file = "", digits = 10)
    
    if (annotate_tips) {
      annotation_results <- annotate_tips_based_on_parents(newick_with_labels)
      annotated_newick <- annotation_results$annotated_newick
      tip_states <- annotation_results$tip_states
      
      # Count the number of tips in each state
      state_counts <- table(unlist(tip_states))
      
      # Check if all states have at least the minimum number of tips
      if (all(state_counts >= min_clade_size)) {
        break
      } else {
        cat("Re-running selection process to meet minimum clade size criterion\n")
      }
    } else {
      break
    }
  }
  
  # Plot the tree with the selected nodes highlighted
  plot_selected_nodes(annotated_tree, selected_nodes)
  
  return(newick_with_labels)
}

# Function to annotate selected nodes with a given model, including their descendants
annotate_nodes <- function(tree, selected_nodes, models, nested) {
  if (is.null(tree$node.label)) {
    tree$node.label <- rep("", tree$Nnode)
  }
  
  for (i in seq_along(selected_nodes)) {
    node <- selected_nodes[i]
    model <- models[i]
    
    if (nested) {
      if (i > 1) {
        if (is_descendant(tree, node, selected_nodes[1:(i-1)])) {
          tree <- apply_model_to_subtree(tree, node, model, overwrite = TRUE)
        } else if (is_ancestor(tree, node, selected_nodes[1:(i-1)])) {
          tree <- apply_model_to_subtree(tree, node, model, overwrite = FALSE)
        }
      } else {
        tree <- apply_model_to_node(tree, node, model)
      }
    } else {
      tree <- apply_model_to_node(tree, node, model)
    }
  }
  return(tree)
}

# Helper function to check if a node is a descendant of any of the given nodes
is_descendant <- function(tree, node, potential_ancestors) {
  descendants_of_ancestors <- unlist(lapply(potential_ancestors, function(anc) {
    get_descendants(tree, anc)
  }))
  return(node %in% descendants_of_ancestors)
}

# Helper function to check if a node is an ancestor of any of the given nodes
is_ancestor <- function(tree, node, potential_descendants) {
  ancestors_of_descendants <- unlist(lapply(potential_descendants, function(desc) {
    get_ancestors(tree, desc)
  }))
  return(node %in% ancestors_of_descendants)
}

# Function to apply model to node and its descendants with option to overwrite
apply_model_to_subtree <- function(tree, node, model, overwrite = TRUE) {
  corrected_index <- node - length(tree$tip.label)
  if (corrected_index > 0 && corrected_index <= length(tree$node.label)) {
    if (overwrite || tree$node.label[corrected_index] == "") {
      tree$node.label[corrected_index] <- model
    }
    descendants <- get_descendants(tree, node)
    for (descendant in descendants) {
      corrected_descendant_index <- descendant - length(tree$tip.label)
      if (corrected_descendant_index > 0 && corrected_descendant_index <= length(tree$node.label)) {
        if (overwrite || tree$node.label[corrected_descendant_index] == "") {
          tree$node.label[corrected_descendant_index] <- model
        }
      }
    }
  }
  return(tree)
}

# Helper function to apply model to a node and its descendants
apply_model_to_node <- function(tree, node, model) {
  corrected_index <- node - length(tree$tip.label)
  if (corrected_index > 0 && corrected_index <= length(tree$node.label)) {
    tree$node.label[corrected_index] <- model
    descendants <- get_descendants(tree, node)
    for (descendant in descendants) {
      corrected_descendant_index <- descendant - length(tree$tip.label)
      if (corrected_descendant_index > 0 && corrected_descendant_index <= length(tree$node.label)) {
        tree$node.label[corrected_descendant_index] <- model
      }
    }
  }
  return(tree)
}

# Function to plot the tree with selected nodes highlighted
plot_selected_nodes <- function(tree, selected_nodes) {
  plot.phylo(tree, show.node.label = TRUE, cex = 0.5, no.margin=T, show.tip.label=F)
  nodelabels(frame = "none", node = selected_nodes, col = "red", cex = 1)
}

write_vector_to_files <- function(vector, prefix) {
  for (i in seq_along(vector)) {
    file_name <- sprintf("%s_%d.tre", prefix, i)
    writeLines(as.character(vector[i]), file_name)
  }
}

write_clean_vector_to_files <- function(vector, prefix) {
  for (i in seq_along(vector)) {
    file_name <- sprintf("%s_%d.tre", prefix, i)
    writeLines(as.character(write.tree(read.tree(text=vector[[i]]))), file_name)
  }
}

#getParent is from https://github.com/DomBennett/MoreTreeTools/blob/master/R/get-methods.R
getParent <- function (tree, node=NULL, tips=NULL, edges=NULL) {
  if (!is.null (node) & length (node) == 1) {
    if (!is.numeric (node)) {
      stop ('Node must be numeric')
    }
    if (node > getSize (tree) + tree$Nnode) {
      stop ('Node not in tree')
    }
    if ((node == getSize (tree) + 1) & is.rooted (tree)) {
      # if node is root, return it
      return (node)
    }
    return (tree$edge[tree$edge[ ,2] == node, 1])
  } else if (!is.null (tips)) {
    if (is.character (tips)) {
      # if tips are labels
      edges <- match (match (tips, tree$tip.label), tree$edge[,2])
    } else {
      # ... else they're numbers
      edges <- match (tips, tree$edge[,2])
    }
  } else if (!is.null (node)) {
    edges <- which (tree$edge[ ,2] %in% node)
  } else if (!is.null (edges)) {
    if (is.character (edges) & !is.null (tree$edge.label)) {
      # assume they are labels
      edges <- match (edges, tree$edge.label)
    }
  } else {
    stop ('Must provide either edges, tips or nodes argument')
  }
  end.nodes <- tree$edge[edges, 1]
  term.node <- length (tree$tip.label) + 1
  while (TRUE){
    if (sum (end.nodes[1] == end.nodes) == length (end.nodes)){
      break
    }
    end.nodes <- sort (end.nodes, TRUE)
    start.node <- end.nodes[1]
    edge <- match (start.node, tree$edge[,2])
    end.node <- tree$edge[edge,1]
    edges <- c(edges, edge)
    end.nodes <- c(end.nodes[!end.nodes %in% start.node], end.node)
  }
  return (end.nodes[1])
}

#function to generate tip annotations
annotate_tips_based_on_parents <- function(annotated_newick) {
  cat("Parsing Annotated Newick string into a tree with annotations...\n")
  
  tree <- parse_annotated(annotated_newick, format = "newick")
  tip_states <- vector("list", length(tree$tip.label))
  cat("Modifying annotations for tips based on their parent nodes...\n")
  
  for (tip_label in tree$tip.label) {
    parent_node <- getParent(tree, tips = tip_label)
    if (!is.na(tree$node.comment[parent_node]) && nzchar(tree$node.comment[parent_node])) {
      annotation_str <- gsub("[[:space:]]*\\[&", "", tree$node.comment[parent_node])
      annotation_str <- gsub("\\]", "", annotation_str)
      pattern <- sprintf("(%s:[^,;)]+)", tip_label)
      replacement <- sprintf("\\1[%s]", annotation_str)
      annotated_newick <- gsub(pattern, replacement, annotated_newick, perl = TRUE)
      tip_states[[tip_label]] <- annotation_str
      cat(sprintf("Tip %s annotated with %s\n", tip_label, annotation_str))
    } else {
      cat(sprintf("No annotation for parent node of tip %s\n", tip_label))
    }
  }
  
  cat("Final annotated Newick string:\n")
  cat(annotated_newick, "\n\n")
  return(list(annotated_newick = annotated_newick, tip_states = tip_states))
}

# Function to generate a nucleotide frequency string
generate_nucleotide_freq_string <- function(alpha = c(2, 2, 2, 2)) {
  # Generate a sample
  sample <- rdirichlet(1, alpha)
  
  # Convert the sample to a formatted string
  freq_string <- paste0("[&model=HKY{2.0}+F{", 
                        paste(round(sample, 2), collapse = "/"), 
                        "}]")
  
  return(freq_string)
}


### 
#independent shift models


# Simulate a tree using pbtree from the ape package
set.seed(123) # Set a seed for reproducibility
tree <- pbtree(n = 100)
cat("Simulated tree in Newick format:\n", write.tree(tree), "\n\n")

# Example usage
min_clade_size <- 4
reps <- 1
#models <- c("[&model=HKY{2.0}+GC]", "[&model=GTR{2.0}+GC]")
#models <- c("[HKY{2.0}]", "[GTR]")
models <- c("[&model=HKY{2.0}+F{0.35/0.15/0.15/0.35}]")

nested_shifts <- F # Set to TRUE if you want nested shifts

# Run the annotation process with the option for nested or independent shifts
annotated_newick <- annotate_branches(tree, models, reps, min_clade_size, nested_shifts, annotate_tips=T)


#read in trees for preparation
exontree<-read.tree(file="/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/unmasked/min2x/qual20_cov2_haplo_bestonly/initial_test_filters/min50bp_min10p_aligned/ALIGNED/phased/sample1haplo/phyhetnucbf/backbone_constraint/janus/8952e31d/gamma/redo_m3/exons_MFP_MERGE_MRL3_constraint_G_UE_UL_M3/exons_MFP_MERGE_MRL3_constraint.rooted.treefile")
exontree$tip.label <- gsub("_1", "", exontree$tip.label)
exontree$node.label<-NULL

introntree<-read.tree(file="/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/unmasked/min2x/qual20_cov2_haplo_bestonly/initial_test_filters/min50bp_min10p_aligned/ALIGNED/phased/sample1haplo/phyhetnucbf/backbone_constraint/janus/8952e31d/gamma/redo_m3/introns_MFP_MERGE_MRL3_constraint_G_UE_UL_M3/introns_MFP_MERGE_MRL3_constraint.rooted.treefile")
introntree$tip.label <- gsub("_1", "", introntree$tip.label)
introntree$node.label<-NULL

utrtree<-read.tree(file="/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/unmasked/min2x/qual20_cov2_haplo_bestonly/initial_test_filters/min50bp_min10p_aligned/ALIGNED/phased/sample1haplo/phyhetnucbf/backbone_constraint/janus/8952e31d/gamma/redo_m3/utrs_MFP_MERGE_MRL3_constraint_G_UE_UL_M3/utrs_MFP_MERGE_MRL3_constraint.rooted.treefile")
utrtree$tip.label <- gsub("_1", "", utrtree$tip.label)
utrtree$node.label<-NULL

mtdnatree<-read.tree(file="/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/mtdnas/janus/concat_all/mtDNA_all_MRL3_constraint.janus.tre")
mtdnatree$tip.label <- gsub("_1", "", mtdnatree$tip.label)
mtdnatree$node.label<-NULL


set.seed(123)

# Run the annotation process
exon_base <- list()
for(i in 1:10){exon_base[[i]]<-annotate_branches(exontree, models, reps, min_clade_size, annotate_tips = T)}

# Run the annotation process
intron_base <- list()
for(i in 1:10){intron_base[[i]]<-annotate_branches(introntree, models, reps, min_clade_size, annotate_tips = T)}

# Run the annotation process
utr_base <- list()
for(i in 1:10){utr_base[[i]]<-annotate_branches(utrtree, models, reps, min_clade_size, annotate_tips = T)}

# Run the annotation process
mtdna_base <- list()
for(i in 1:10){mtdna_base[[i]]<-annotate_branches(mtdnatree, models, reps, min_clade_size, annotate_tips = T)}

# Assuming your vector is named 'exon_base'
write_vector_to_files(exon_base, "exon")
write_vector_to_files(intron_base, "intron")
write_vector_to_files(utr_base, "utr")
write_vector_to_files(mtdna_base, "mtdna")


# Assuming your vector is named 'exon_base'
write_clean_vector_to_files(exon_base, "exon")
write_clean_vector_to_files(intron_base, "intron")
write_clean_vector_to_files(utr_base, "utr")
write_clean_vector_to_files(mtdna_base, "mtdna")


##testing indicates things look good

#function to visulize the output, to check that it is working
readFastaAndPlot <- function(fasta_file) {
  # Read the FASTA file
  sequences <- readDNAStringSet(fasta_file)
  names(sequences)<- gsub(" ", "", names(sequences))
  
  # Create a named vector of sequences
  seq_vector <- sapply(sequences, as.character)
  #names(seq_vector) == exon_1$tip.label
  
  
  # Prepare data for plotting
  base_counts <- matrix(0, nrow = length(seq_vector), ncol = 4)
  colnames(base_counts) <- c("A", "T", "G", "C")
  rownames(base_counts) <- names(seq_vector)
  
  # Calculate base frequencies for each sequence
  for (i in seq_along(seq_vector)) {
    seq_bases <- unlist(strsplit(seq_vector[i], ""))
    for (base in colnames(base_counts)) {
      base_counts[i, base] <- sum(seq_bases == base) / length(seq_bases)
    }
  }
  
  # Plotting
  bp <- barplot(t(base_counts), main = "Base Frequencies in FASTA Sequences", xlab = "Sequence ID", ylab = "Frequency", col = rainbow(4), beside = FALSE, legend = colnames(base_counts), horiz=T, las=2, cex.names=0.15)
  
}
readFastaAndGenerateGCTable <- function(fasta_file, phylo) {
  # Read the FASTA file
  sequences <- readDNAStringSet(fasta_file)
  names(sequences) <- gsub(" ", "", names(sequences))
  
  # Create a named vector of sequences
  seq_vector <- sapply(sequences, as.character)
  names(seq_vector) <- names(sequences)
  
  # Reorder the sequences according to the phylo object
  reordered_seq_vector <- seq_vector[phylo$tip.label]
  
  # Check for any missing sequences
  if (any(is.na(reordered_seq_vector))) {
    stop("Some sequences in the phylogeny are not present in the FASTA file.")
  }
  
  # Prepare data for GC content calculation
  gc_content <- numeric(length(reordered_seq_vector))
  names(gc_content) <- phylo$tip.label
  
  # Calculate GC content for each sequence
  for (i in seq_along(reordered_seq_vector)) {
    seq_bases <- unlist(strsplit(reordered_seq_vector[i], ""))
    gc_count <- sum(seq_bases %in% c("G", "C"))
    gc_content[i] <- gc_count / length(seq_bases) * 100
  }
  
  # Generate a table of %GC content
  gc_content_table <- data.frame(Sequence_ID = names(gc_content), Percent_GC = gc_content)
  return(gc_content_table)
}

exon_1<-read.tree(file="/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/simulated_alignments/false_negatives/1_shift/exon_1_1.tre")
readFastaAndPlot(fasta_file = '/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/simulated_alignments/false_negatives/1_shift/EXONS3_1.fa')
readFastaAndGenerateGCTable(fasta_file = '/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/simulated_alignments/false_negatives/1_shift/EXONS3_1.fa', phylo = exon_1)


exon_2<-read.tree(file="/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/simulated_alignments/false_negatives/2_shift/exon_1_1.tre")
readFastaAndPlot(fasta_file = '/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/simulated_alignments/false_negatives/2_shift/EXONS1_1.fa')
readFastaAndGenerateGCTable(fasta_file = '/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/simulated_alignments/false_negatives/2_shift/EXONS1_1.fa', phylo = exon_2)


exon_3<-read.tree(file="/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/simulated_alignments/false_negatives/1_shift/exon_3.tre")
readFastaAndPlot(fasta_file = '/Applications/Phylogenetics/iqtree-2.2.2.6-MacOSX/bin/EXONS3_1.fa')
readFastaAndGenerateGCTable(fasta_file = '/Applications/Phylogenetics/iqtree-2.2.2.6-MacOSX/bin/EXONS3_1.fa', phylo = exon_3)


#checking names
test<-readDNAStringSet(filepath = '/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/simulated_alignments/false_negatives/1_shift/EXONS1_1.fa')
names(test)
cbind(exon_1$tip.label, exon_1$tip.label %in% names(test), names(test))
#cbind(exon_2$tip.label, exon_2$tip.label %in% names(test), names(test))



### nested shift models


# Simulate a tree using pbtree from the ape package
set.seed(123) # Set a seed for reproducibility
tree <- pbtree(n = 200)
cat("Simulated tree in Newick format:\n", write.tree(tree), "\n\n")

# Example usage
min_clade_size <- 4
reps <- 2
#models <- c("[&model=HKY{2.0}+GC]", "[&model=GTR{2.0}+GC]")
#models <- c("[HKY{2.0}]", "[GTR]")
#models <- c("[A]", "[B]")
models <- c("[&model=HKY{2.0}+F{0.35/0.15/0.15/0.35}]", "[&model=HKY{2.0}+F{0.15/0.35/0.35/0.15}]")
#models<-c(generate_nucleotide_freq_string(), generate_nucleotide_freq_string())


nested_shifts <- T # Set to TRUE if you want nested shifts

# Run the annotation process with the option for nested or independent shifts
annotated_newick <- annotate_branches(tree, models, reps, min_clade_size, nested_shifts, annotate_tips=T)


#read in trees for preparation
exontree<-read.tree(file="/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/unmasked/min2x/qual20_cov2_haplo_bestonly/initial_test_filters/min50bp_min10p_aligned/ALIGNED/phased/sample1haplo/phyhetnucbf/backbone_constraint/janus/8952e31d/gamma/redo_m3/exons_MFP_MERGE_MRL3_constraint_G_UE_UL_M3/exons_MFP_MERGE_MRL3_constraint.rooted.treefile")
exontree$tip.label <- gsub("_1", "", exontree$tip.label)
exontree$node.label<-NULL

introntree<-read.tree(file="/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/unmasked/min2x/qual20_cov2_haplo_bestonly/initial_test_filters/min50bp_min10p_aligned/ALIGNED/phased/sample1haplo/phyhetnucbf/backbone_constraint/janus/8952e31d/gamma/redo_m3/introns_MFP_MERGE_MRL3_constraint_G_UE_UL_M3/introns_MFP_MERGE_MRL3_constraint.rooted.treefile")
introntree$tip.label <- gsub("_1", "", introntree$tip.label)
introntree$node.label<-NULL

utrtree<-read.tree(file="/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/unmasked/min2x/qual20_cov2_haplo_bestonly/initial_test_filters/min50bp_min10p_aligned/ALIGNED/phased/sample1haplo/phyhetnucbf/backbone_constraint/janus/8952e31d/gamma/redo_m3/utrs_MFP_MERGE_MRL3_constraint_G_UE_UL_M3/utrs_MFP_MERGE_MRL3_constraint.rooted.treefile")
utrtree$tip.label <- gsub("_1", "", utrtree$tip.label)
utrtree$node.label<-NULL

mtdnatree<-read.tree(file="/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/mtdnas/janus/concat_all/mtDNA_all_MRL3_constraint.janus.tre")
mtdnatree$tip.label <- gsub("_1", "", mtdnatree$tip.label)
mtdnatree$node.label<-NULL


set.seed(123)

# Run the annotation process
exon_base <- list()
for(i in 1:10){exon_base[[i]]<-annotate_branches(exontree, models, reps, min_clade_size, nested_shifts, annotate_tips=T)}
#annotate_branches(exontree, models, reps, min_clade_size, nested_shifts, annotate_tips=T)

# Run the annotation process
intron_base <- list()
for(i in 1:10){intron_base[[i]]<-annotate_branches(introntree, models, reps, min_clade_size, nested_shifts, annotate_tips=T)}

# Run the annotation process
utr_base <- list()
for(i in 1:10){utr_base[[i]]<-annotate_branches(utrtree, models, reps, min_clade_size, nested_shifts, annotate_tips=T)}

# Run the annotation process
mtdna_base <- list()
for(i in 1:10){mtdna_base[[i]]<-annotate_branches(mtdnatree, models, reps, min_clade_size, nested_shifts, annotate_tips=T)}

# Assuming your vector is named 'exon_base'
write_vector_to_files(exon_base, "exon")
write_vector_to_files(intron_base, "intron")
write_vector_to_files(utr_base, "utr")
write_vector_to_files(mtdna_base, "mtdna")


# Assuming your vector is named 'exon_base'
write_clean_vector_to_files(exon_base, "exon")
write_clean_vector_to_files(intron_base, "intron")
write_clean_vector_to_files(utr_base, "utr")
write_clean_vector_to_files(mtdna_base, "mtdna")


##testing indicates things look good

exon_1<-read.tree(file="/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/simulated_alignments/false_negatives/1_shift/exon_3_1.tre")
readFastaAndPlot(fasta_file = '/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/simulated_alignments/false_negatives/1_shift/EXONS3_1.fa')
readFastaAndGenerateGCTable(fasta_file = '/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/simulated_alignments/false_negatives/1_shift/EXONS3_1.fa', phylo = exon_1)


exon_2<-read.tree(file="/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/simulated_alignments/false_negatives/2_shift/exon_1_1.tre")
readFastaAndPlot(fasta_file = '/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/simulated_alignments/false_negatives/2_shift/EXONS1_1.fa')
readFastaAndGenerateGCTable(fasta_file = '/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/simulated_alignments/false_negatives/2_shift/EXONS1_1.fa', phylo = exon_2)


#checking names
exon_1<-read.tree(file="/Users/cotinga/temp/2shift_temp/exon_1_1.tre")
test<-readDNAStringSet(filepath = "/Users/cotinga/temp/2shift_temp/EXONS1_1.fa")
names(test)
cbind(exon_1$tip.label, exon_1$tip.label %in% names(test), names(test))

cbind(exon_2$tip.label, exon_2$tip.label %in% names(test), names(test))



##consider adding a flag to the shift node selection that excludes nodes immediately parent or daughter to the initial round one sample

