# Script for generating input templates for iqtree to simulate shift models
# Collaboration between Jacob Berv and ChatGPT4 with modifications for nested or independent shifts

library(ape)
library(phytools)
require(Biostrings)
require(phylotate)
library(DirichletReg)
library(parallel)
library(pbmcapply)

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
select_nodes_for_annotation <- function(tree, candidate_nodes, reps, nested = FALSE, buffer = 1, derived=FALSE) {
  # Check for derived = TRUE with reps = 2
  if (derived && reps != 2) {
    stop("The 'derived' option can only be used when 'reps' is equal to 2.")
  }
  
  if (length(candidate_nodes) < reps) {
    stop("The number of reps is greater than the number of candidate nodes.")
  }
  
  selected_nodes <- c()  # To keep track of nodes selected for shifts
  excluded_nodes <- c()  # For independent shifts
  
  cat("Starting shifts selection...\n")
  
  if (nested) {
    retry_count <- 0
    max_retries <- 100  # Maximum number of retries
    
    while (TRUE) {
      if (length(selected_nodes) < reps) {
        if (length(selected_nodes) == 0) {
          selected <- sample(candidate_nodes, 1, replace = FALSE)
          selected_nodes <- c(selected_nodes, selected)
          cat("First node selected for shift (nested): ", selected, "\n")
        }
        
        # Modification to nested pool calculation
        if (derived) {
          nested_pool <- intersect(get_descendants(tree, selected_nodes[1]), candidate_nodes)
        } else {
          nested_pool <- intersect(union(get_descendants(tree, selected_nodes[1]), get_ancestors(tree, selected_nodes[1])), candidate_nodes)
        }
        
        cat("Nested pool for buffer estimation: ", nested_pool, "\n")
        
        buffer_candidates <- vector()
        for (node in nested_pool) {
          path <- nodepath(tree, from = selected_nodes[1], to = node)
          cat("Node: ", node, " Path Length: ", length(path) - 1, " (Buffer: ", buffer, ")\n")
          if (length(path) - 1 >= buffer) {
            buffer_candidates <- c(buffer_candidates, node)
          }
        }
        
        cat("Buffer candidates: ", buffer_candidates, "\n")
        
        if (length(buffer_candidates) == 1) {
          selected <- buffer_candidates  # Auto-select the single candidate
          selected_nodes <- c(selected_nodes, selected)
          cat("Only one buffer candidate available. Auto-selected node for shift (nested): ", selected, "\n")
        } else if (length(buffer_candidates) > 1) {
          selected <- sample(buffer_candidates, 1, replace = FALSE)
          selected_nodes <- c(selected_nodes, selected)
          cat("Selected node for shift (nested): ", selected, "\n")
        } else {
          retry_count <- retry_count + 1
          if (retry_count <= max_retries) {
            cat("Retry ", retry_count, "/", max_retries, ": Resetting selection process\n")
            selected_nodes <- c()  # Resetting for a new attempt
            next
          } else {
            stop("Maximum number of retries reached. Unable to find suitable nodes for nested shifts.")
          }
        }
      } else {
        return(selected_nodes)
      }
    }
  } else {
    # Logic for independent shifts (nested = FALSE)
    for (i in 1:reps) {
      buffer_candidates <- setdiff(candidate_nodes, excluded_nodes)
      
      if (i > 1) {
        temp_buffer_candidates <- vector()
        for (node in buffer_candidates) {
          is_buffer_compliant <- TRUE
          for (selected_node in selected_nodes) {
            path <- nodepath(tree, from = selected_node, to = node)
            cat("Checking buffer for Node: ", node, " against Selected Node: ", selected_node, " Path Length: ", length(path) - 1, " (Buffer: ", buffer, ")\n")
            if (length(path) - 1 < buffer) {
              is_buffer_compliant <- FALSE
              break
            }
          }
          if (is_buffer_compliant) {
            temp_buffer_candidates <- c(temp_buffer_candidates, node)
          }
        }
        buffer_candidates <- temp_buffer_candidates
      }
      
      cat("Buffer candidates for independent shift: ", buffer_candidates, "\n")
      if (length(buffer_candidates) == 0) {
        cat("No buffer-compliant candidates found for independent shifts.\n")
        next
      }
      
      selected <- sample(buffer_candidates, 1, replace = FALSE)
      selected_nodes <- c(selected_nodes, selected)
      
      # Update excluded_nodes with lineage of the selected node
      descendants_and_ancestors <- union(get_descendants(tree, selected), get_ancestors(tree, selected))
      excluded_nodes <- union(excluded_nodes, c(descendants_and_ancestors, selected))
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
annotate_branches <- function(input_tree, models, reps, min_clade_size, nested = FALSE, annotate_tips = FALSE, buffer = 1, max_clade_size, derived) {
  if (is.character(input_tree)) {
    tree <- read.tree(text = input_tree)
  } else if (inherits(input_tree, "phylo")) {
    tree <- input_tree
  } else {
    stop("Invalid input: input_tree must be a Newick string or a phylo object.")
  }
  
  edge_indices <- edge_indices_N(tree, min_clade_size)
  candidate_nodes <- node_indices_edge(tree, edge_indices)
  print(candidate_nodes)
  
  max_retries <- 100  # Limit the number of retries to prevent infinite loops
  retry_count <- 0
  previous_selected_nodes <- NULL
  
  repeat {
    selected_nodes <- select_nodes_for_annotation(tree, candidate_nodes, reps, nested, buffer, derived)
    # Check if the new selection is different from the previous one
    if (!is.null(previous_selected_nodes) && all(selected_nodes %in% previous_selected_nodes)) {
      cat("No new nodes selected, re-running selection process...\n")
      if (retry_count >= max_retries) {
        stop("Maximum retry limit reached. Unable to find suitable nodes.")
      }
      retry_count <- retry_count + 1
    }
    previous_selected_nodes <- selected_nodes
    
    annotated_tree <- annotate_nodes(tree, selected_nodes, models, nested)
    newick_with_labels <- write.tree(annotated_tree, file = "", digits = 10)
    
    if (annotate_tips) {
      annotation_results <- annotate_tips_based_on_parents(newick_with_labels)
      annotated_newick <- annotation_results$annotated_newick
      tip_states <- annotation_results$tip_states
      
      state_counts <- table(unlist(tip_states))
      if (all(state_counts >= min_clade_size & state_counts < max_clade_size)) {
        break
      } else {
        cat("Re-running selection process to meet clade size criteria\n")
        if (retry_count >= max_retries) {
          stop("Maximum retry limit reached. Unable to find suitable nodes.")
        }
        retry_count <- retry_count + 1
      }
    } else {
      break
    }
  }
  
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

# Function to generate a nucleotide frequency string with samples > 0
generate_nucleotide_freq_string <- function(alpha = c(1, 1, 1, 1)) {
  repeat {
    # Generate a sample
    sample <- rdirichlet(1, alpha)
    
    # Check if all elements are greater than 0.01 and their sum equals 1
    if (all(sample > 0.01) && sum(round(sample, 2)) == 1) {
      # Convert the sample to a formatted string
      freq_string <- paste0("[&model=HKY{2.0}+F{", 
                            paste(round(sample, 2), collapse = "/"), 
                            "}]")
      return(freq_string)
    }
    # If conditions are not met, repeat the sampling
  }
}

#generat model strings following dirichlet sampling
models <- function(n) {
  sapply(1:n, function(x) generate_nucleotide_freq_string())
}

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
processDirectoryForGCTable <- function(directory, phylo) {
  # List all .fa files in the directory
  fasta_files <- list.files(directory, pattern = "\\.fa$", full.names = TRUE)
  
  # Initialize a list to store data frames
  all_gc_content <- list()
  
  # Loop over the FASTA files
  for (fasta_file in fasta_files) {
    gc_content_table <- readFastaAndGenerateGCTable(fasta_file, phylo)
    col_name <- basename(fasta_file)
    colnames(gc_content_table)[2] <- col_name  # Rename Percent_GC column to file name
    all_gc_content[[col_name]] <- gc_content_table
  }
  
  # Combine all data frames by Sequence_ID
  combined_gc_content <- Reduce(function(x, y) merge(x, y, by = "Sequence_ID", all = TRUE), all_gc_content)
  
  return(combined_gc_content)
}


##parsing the output (go version)
parse_files_hmsj <- function(directory, expected, exclude_zero_lines = FALSE, group_by_reps = T, group_by_type = T) {
  # Get list of files with 'stderror' in their names
  file_names <- list.files(path = directory, pattern = "stderror", full.names = TRUE)
  
  # Initialize a data frame to store results
  results <- data.frame(file_name = character(), prefix = character(), type = character(), line_count = integer(), stringsAsFactors = FALSE)
  
  # Process each file
  for (file in file_names) {
    # Read the file content
    lines <- readLines(file)
    
    # Find the index of the line that contains 'Final models'
    final_model_index <- which(grepl("Final models", lines))
    
    # Initialize count
    count <- 0
    
    # Check if the next line contains '-------'
    if (length(final_model_index) > 0 && (final_model_index + 1) <= length(lines)) {
      if (grepl("-------", lines[final_model_index + 1])) {
        # Count the number of lines after 'Final models' and '-------'
        count <- length(lines) - (final_model_index + 1)
      }
    }
    
    # Extract the prefix and type from the file name
    file_prefix <- ifelse(group_by_reps, strsplit(basename(file), "_")[[1]][1], NA)
    type <- ifelse(group_by_type, gsub("[0-9]+", "", file_prefix), NA)
    
    # Add results to the data frame if they meet the condition
    if (!exclude_zero_lines || count > 0) {
      results <- rbind(results, data.frame(file_name = basename(file), prefix = file_prefix, type = type, line_count = count))
    }
  }
  
  # Calculate success rate, false_neg, false_pos and return results based on group_by_reps and group_by_type flags
  if (group_by_reps) {
    if (group_by_type) {
      # Group by type and summarize line counts
      type_grouped_results <- aggregate(line_count ~ type, data = results, FUN = sum)
      # Calculate success rate
      type_grouped_results$success_rate <- type_grouped_results$line_count / (expected * sapply(type_grouped_results$type, function(x) sum(results$type == x)))
      # Calculate false negative rate
      type_grouped_results$false_neg <- pmax(0, 1 - type_grouped_results$success_rate)
      # Calculate false positive rate
      type_grouped_results$false_pos <- pmax(type_grouped_results$success_rate - 1, 0)
      return(type_grouped_results[, c("type", "line_count", "success_rate", "false_neg", "false_pos")])
    } else {
      # Group by prefix and summarize line counts
      grouped_results <- aggregate(line_count ~ prefix, data = results, FUN = sum)
      # Calculate success rate
      grouped_results$success_rate <- grouped_results$line_count / (expected * sapply(grouped_results$prefix, function(x) sum(results$prefix == x)))
      # Calculate false negative rate
      grouped_results$false_neg <- pmax(0, 1 - grouped_results$success_rate)
      # Calculate false positive rate
      grouped_results$false_pos <- pmax(grouped_results$success_rate - 1, 0)
      return(grouped_results[, c("prefix", "line_count", "success_rate", "false_neg", "false_pos")])
    }
  } else {
    # Add success rate column
    results$success_rate <- results$line_count / expected
    # Add false negative rate column
    results$false_neg <- pmax(0, 1 - results$success_rate)
    # Add false positive rate column
    results$false_pos <- pmax(results$success_rate - 1, 0)
    # Return results without grouping
    return(results[, c("file_name", "line_count", "success_rate", "false_neg", "false_pos")])
  }
}

##parsing the output (c version)
parse_files_hmshj <- function(directory, expected, exclude_zero_lines = FALSE, group_by_reps = T, group_by_type = T) {
  # Get list of files with 'stderror' in their names
  file_names <- list.files(path = directory, pattern = "stderror", full.names = TRUE)
  
  # Initialize a data frame to store results
  results <- data.frame(file_name = character(), prefix = character(), type = character(), line_count = integer(), stringsAsFactors = FALSE)
  
  # Process each file
  for (file in file_names) {
    # Read the file content
    lines <- readLines(file)
    
    # Find the index of the line that contains 'Final models'
    final_model_index <- which(grepl("final:", lines))
    
    # Initialize count
    count <- 0
    
    # Check if the next line contains '-------'
    if (length(final_model_index) > 0 && (final_model_index + 1) <= length(lines)) {
      if (grepl("rm: ", lines[final_model_index + 1])) {
        # Count the number of lines after 'Final models' and '-------'
        count <- length(lines) - (final_model_index + 1)
      }
    }
    
    # Extract the prefix and type from the file name
    file_prefix <- ifelse(group_by_reps, strsplit(basename(file), "_")[[1]][1], NA)
    type <- ifelse(group_by_type, gsub("[0-9]+", "", file_prefix), NA)
    
    # Add results to the data frame if they meet the condition
    if (!exclude_zero_lines || count > 0) {
      results <- rbind(results, data.frame(file_name = basename(file), prefix = file_prefix, type = type, line_count = count))
    }
  }
  
  # Calculate success rate, false_neg, false_pos and return results based on group_by_reps and group_by_type flags
  if (group_by_reps) {
    if (group_by_type) {
      # Group by type and summarize line counts
      type_grouped_results <- aggregate(line_count ~ type, data = results, FUN = sum)
      # Calculate success rate
      type_grouped_results$success_rate <- type_grouped_results$line_count / (expected * sapply(type_grouped_results$type, function(x) sum(results$type == x)))
      # Calculate false negative rate
      type_grouped_results$false_neg <- pmax(0, 1 - type_grouped_results$success_rate)
      # Calculate false positive rate
      type_grouped_results$false_pos <- pmax(type_grouped_results$success_rate - 1, 0)
      return(type_grouped_results[, c("type", "line_count", "success_rate", "false_neg", "false_pos")])
    } else {
      # Group by prefix and summarize line counts
      grouped_results <- aggregate(line_count ~ prefix, data = results, FUN = sum)
      # Calculate success rate
      grouped_results$success_rate <- grouped_results$line_count / (expected * sapply(grouped_results$prefix, function(x) sum(results$prefix == x)))
      # Calculate false negative rate
      grouped_results$false_neg <- pmax(0, 1 - grouped_results$success_rate)
      # Calculate false positive rate
      grouped_results$false_pos <- pmax(grouped_results$success_rate - 1, 0)
      return(grouped_results[, c("prefix", "line_count", "success_rate", "false_neg", "false_pos")])
    }
  } else {
    # Add success rate column
    results$success_rate <- results$line_count / expected
    # Add false negative rate column
    results$false_neg <- pmax(0, 1 - results$success_rate)
    # Add false positive rate column
    results$false_pos <- pmax(results$success_rate - 1, 0)
    # Return results without grouping
    return(results[, c("file_name", "line_count", "success_rate", "false_neg", "false_pos")])
  }
}

#corrected parsing function includig CI calculation and test > 0
parse_files_hmshj_CI_t_test <- function(directory, expected, exclude_zero_lines = FALSE, group_by_reps = TRUE, group_by_type = TRUE, confidence_level = 0.95, p_value_correction = "none") {
  z_value <- qnorm(1 - (1 - confidence_level) / 2)
  file_names <- list.files(path = directory, pattern = "stderror", full.names = TRUE)
  results <- data.frame(file_name = character(), prefix = character(), type = character(), line_count = integer(), variance = numeric(), stringsAsFactors = FALSE)
  
  # Process each file
  for (file in file_names) {
    lines <- readLines(file)
    final_model_index <- which(grepl("final:", lines))
    count <- 0
    if (length(final_model_index) > 0 && (final_model_index + 1) <= length(lines)) {
      if (grepl("rm: ", lines[final_model_index + 1])) {
        count <- length(lines) - (final_model_index + 1)
      }
    }
    file_prefix <- ifelse(group_by_reps, strsplit(basename(file), "_")[[1]][1], NA)
    type <- ifelse(group_by_type, gsub("[0-9]+", "", file_prefix), NA)
    file_variance <- (count - expected)^2
    results <- rbind(results, data.frame(file_name = basename(file), prefix = file_prefix, type = type, line_count = count, variance = file_variance))
  }
  
  # Aggregate line counts and calculate mean variance by type
  grouped_results <- aggregate(cbind(line_count, variance) ~ type, data = results, FUN = sum)
  variance_means <- aggregate(variance ~ type, data = results, FUN = mean)
  grouped_results$variance_mean <- variance_means$variance
  grouped_results$expected_lines <- expected * sapply(grouped_results$type, function(x) sum(results$type == x))
  grouped_results$actual_sample_size <- sapply(grouped_results$type, function(x) sum(results$type == x))
  
  # Define a helper function for SE and CI calculations
  calculate_se_ci <- function(rate, variance, n) {
    if (is.na(rate) || is.na(variance) || variance == 0 || is.na(n) || n == 0 || rate == 0) {
      return(list(se = NA, ci_lower = NA, ci_upper = NA))
    }
    se <- sqrt(variance / n)
    ci_lower <- pmax(0, rate - z_value * se)
    ci_upper <- rate + z_value * se
    return(list(se = se, ci_lower = ci_lower, ci_upper = ci_upper))
  }
  
  # Initialize columns for rates, SE, CI, test statistics, and p-values
  grouped_results$success_rate <- numeric(nrow(grouped_results))
  grouped_results$false_neg <- numeric(nrow(grouped_results))
  grouped_results$false_neg_se <- numeric(nrow(grouped_results))
  grouped_results$false_neg_ci_lower <- numeric(nrow(grouped_results))
  grouped_results$false_neg_ci_upper <- numeric(nrow(grouped_results))
  grouped_results$false_neg_test_statistic <- numeric(nrow(grouped_results))
  grouped_results$false_neg_p_value <- numeric(nrow(grouped_results))
  grouped_results$false_pos <- numeric(nrow(grouped_results))
  grouped_results$false_pos_se <- numeric(nrow(grouped_results))
  grouped_results$false_pos_ci_lower <- numeric(nrow(grouped_results))
  grouped_results$false_pos_ci_upper <- numeric(nrow(grouped_results))
  grouped_results$false_pos_test_statistic <- numeric(nrow(grouped_results))
  grouped_results$false_pos_p_value <- numeric(nrow(grouped_results))
  
  # Process each group
  for (i in 1:nrow(grouped_results)) {
    success_rate <- grouped_results$line_count[i] / grouped_results$expected_lines[i]
    grouped_results$success_rate[i] <- success_rate
    grouped_results$false_neg[i] <- pmax(0, 1 - success_rate)
    grouped_results$false_pos[i] <- pmax(0, success_rate - 1)
    
    # Calculate SE and CI for false negatives
    fn_metrics <- calculate_se_ci(grouped_results$false_neg[i], grouped_results$variance_mean[i], grouped_results$actual_sample_size[i])
    grouped_results$false_neg_se[i] <- fn_metrics$se
    grouped_results$false_neg_ci_lower[i] <- fn_metrics$ci_lower
    grouped_results$false_neg_ci_upper[i] <- fn_metrics$ci_upper
    
    # Calculate SE and CI for false positives
    fp_metrics <- calculate_se_ci(grouped_results$false_pos[i], grouped_results$variance_mean[i], grouped_results$actual_sample_size[i])
    grouped_results$false_pos_se[i] <- fp_metrics$se
    grouped_results$false_pos_ci_lower[i] <- fp_metrics$ci_lower
    grouped_results$false_pos_ci_upper[i] <- fp_metrics$ci_upper
    
    # Calculate test statistics and p-values
    # (handling NA values for cases where rate is zero)
    grouped_results$false_neg_test_statistic[i] <- if (!is.na(fn_metrics$se) && fn_metrics$se > 0) (grouped_results$false_neg[i] - 0) / fn_metrics$se else NA
    grouped_results$false_neg_p_value[i] <- if (!is.na(fn_metrics$se) && fn_metrics$se > 0) pt(grouped_results$false_neg_test_statistic[i], df = grouped_results$actual_sample_size[i] - 1, lower.tail = FALSE) else NA
    grouped_results$false_pos_test_statistic[i] <- if (!is.na(fp_metrics$se) && fp_metrics$se > 0) (grouped_results$false_pos[i] - 0) / fp_metrics$se else NA
    grouped_results$false_pos_p_value[i] <- if (!is.na(fp_metrics$se) && fp_metrics$se > 0) pt(grouped_results$false_pos_test_statistic[i], df = grouped_results$actual_sample_size[i] - 1, lower.tail = FALSE) else NA
  }
  
  # Adjust p-values for multiple comparisons if required
  if (p_value_correction == "bonferroni") {
    num_tests <- nrow(grouped_results)*2
    grouped_results$false_neg_p_value <- p.adjust(grouped_results$false_neg_p_value, method = "bonferroni", n = num_tests)
    grouped_results$false_pos_p_value <- p.adjust(grouped_results$false_pos_p_value, method = "bonferroni", n = num_tests)
  } else if (p_value_correction == "fdr") {
    num_tests <- nrow(grouped_results)*2
    grouped_results$false_neg_p_value <- p.adjust(grouped_results$false_neg_p_value, method = "fdr", n = num_tests)
    grouped_results$false_pos_p_value <- p.adjust(grouped_results$false_pos_p_value, method = "fdr", n = num_tests)
  }
  
  # Select the columns to return, including test statistics
  columns_to_return <- c("type", "line_count", "success_rate", "expected_lines", "actual_sample_size", "false_neg", "false_neg_se", "false_neg_ci_lower", "false_neg_ci_upper", "false_neg_test_statistic", "false_neg_p_value", "false_pos", "false_pos_se", "false_pos_ci_lower", "false_pos_ci_upper", "false_pos_test_statistic", "false_pos_p_value")
  return(grouped_results[, columns_to_return, drop = FALSE])
}
parse_files_hmshj_CI_z_test <- function(directory, expected, exclude_zero_lines = FALSE, group_by_reps = TRUE, group_by_type = TRUE, confidence_level = 0.95, p_value_correction = "none") {
  z_value <- qnorm(1 - (1 - confidence_level) / 2)
  file_names <- list.files(path = directory, pattern = "stderror", full.names = TRUE)
  results <- data.frame(file_name = character(), prefix = character(), type = character(), line_count = integer(), variance = numeric(), stringsAsFactors = FALSE)
  
  # Process each file
  for (file in file_names) {
    lines <- readLines(file)
    final_model_index <- which(grepl("final:", lines))
    count <- 0
    if (length(final_model_index) > 0 && (final_model_index + 1) <= length(lines)) {
      if (grepl("rm: ", lines[final_model_index + 1])) {
        count <- length(lines) - (final_model_index + 1)
      }
    }
    file_prefix <- ifelse(group_by_reps, strsplit(basename(file), "_")[[1]][1], NA)
    type <- ifelse(group_by_type, gsub("[0-9]+", "", file_prefix), NA)
    file_variance <- (count - expected)^2
    results <- rbind(results, data.frame(file_name = basename(file), prefix = file_prefix, type = type, line_count = count, variance = file_variance))
  }
  
  # Aggregate line counts and calculate mean variance by type
  grouped_results <- aggregate(cbind(line_count, variance) ~ type, data = results, FUN = sum)
  variance_means <- aggregate(variance ~ type, data = results, FUN = mean)
  grouped_results$variance_mean <- variance_means$variance
  grouped_results$expected_lines <- expected * sapply(grouped_results$type, function(x) sum(results$type == x))
  grouped_results$actual_sample_size <- sapply(grouped_results$type, function(x) sum(results$type == x))
  
  # Define a helper function for SE and CI calculations
  calculate_se_ci <- function(rate, variance, n) {
    if (is.na(rate) || is.na(variance) || variance == 0 || is.na(n) || n == 0 || rate == 0) {
      return(list(se = NA, ci_lower = NA, ci_upper = NA))
    }
    se <- sqrt(variance / n)
    ci_lower <- pmax(0, rate - z_value * se)
    ci_upper <- rate + z_value * se
    return(list(se = se, ci_lower = ci_lower, ci_upper = ci_upper))
  }
  
  # Initialize columns for rates, SE, CI, test statistics, and p-values
  grouped_results$success_rate <- numeric(nrow(grouped_results))
  grouped_results$false_neg <- numeric(nrow(grouped_results))
  grouped_results$false_neg_se <- numeric(nrow(grouped_results))
  grouped_results$false_neg_ci_lower <- numeric(nrow(grouped_results))
  grouped_results$false_neg_ci_upper <- numeric(nrow(grouped_results))
  grouped_results$false_neg_test_statistic <- numeric(nrow(grouped_results))
  grouped_results$false_neg_p_value <- numeric(nrow(grouped_results))
  grouped_results$false_pos <- numeric(nrow(grouped_results))
  grouped_results$false_pos_se <- numeric(nrow(grouped_results))
  grouped_results$false_pos_ci_lower <- numeric(nrow(grouped_results))
  grouped_results$false_pos_ci_upper <- numeric(nrow(grouped_results))
  grouped_results$false_pos_test_statistic <- numeric(nrow(grouped_results))
  grouped_results$false_pos_p_value <- numeric(nrow(grouped_results))
  
  # Process each group
  for (i in 1:nrow(grouped_results)) {
    success_rate <- grouped_results$line_count[i] / grouped_results$expected_lines[i]
    grouped_results$success_rate[i] <- success_rate
    grouped_results$false_neg[i] <- pmax(0, 1 - success_rate)
    grouped_results$false_pos[i] <- pmax(0, success_rate - 1)
    
    # Calculate SE and CI for false negatives
    fn_metrics <- calculate_se_ci(grouped_results$false_neg[i], grouped_results$variance_mean[i], grouped_results$actual_sample_size[i])
    grouped_results$false_neg_se[i] <- fn_metrics$se
    grouped_results$false_neg_ci_lower[i] <- fn_metrics$ci_lower
    grouped_results$false_neg_ci_upper[i] <- fn_metrics$ci_upper
    
    # Calculate SE and CI for false positives
    fp_metrics <- calculate_se_ci(grouped_results$false_pos[i], grouped_results$variance_mean[i], grouped_results$actual_sample_size[i])
    grouped_results$false_pos_se[i] <- fp_metrics$se
    grouped_results$false_pos_ci_lower[i] <- fp_metrics$ci_lower
    grouped_results$false_pos_ci_upper[i] <- fp_metrics$ci_upper
    
    # Calculate Z-test statistics and p-values
    grouped_results$false_neg_test_statistic[i] <- if (!is.na(fn_metrics$se) && fn_metrics$se > 0) (grouped_results$false_neg[i] - 0) / fn_metrics$se else NA
    grouped_results$false_neg_p_value[i] <- if (!is.na(fn_metrics$se) && fn_metrics$se > 0) pnorm(grouped_results$false_neg_test_statistic[i], lower.tail = FALSE) else NA
    grouped_results$false_pos_test_statistic[i] <- if (!is.na(fp_metrics$se) && fp_metrics$se > 0) (grouped_results$false_pos[i] - 0) / fp_metrics$se else NA
    grouped_results$false_pos_p_value[i] <- if (!is.na(fp_metrics$se) && fp_metrics$se > 0) pnorm(grouped_results$false_pos_test_statistic[i], lower.tail = FALSE) else NA
  }
  
  # Adjust p-values for multiple comparisons if required
  if (p_value_correction == "bonferroni") {
    num_tests <- nrow(grouped_results) * 2
    grouped_results$false_neg_p_value <- p.adjust(grouped_results$false_neg_p_value, method = "bonferroni", n = num_tests)
    grouped_results$false_pos_p_value <- p.adjust(grouped_results$false_pos_p_value, method = "bonferroni", n = num_tests)
  } else if (p_value_correction == "fdr") {
    num_tests <- nrow(grouped_results) * 2
    grouped_results$false_neg_p_value <- p.adjust(grouped_results$false_neg_p_value, method = "fdr", n = num_tests)
    grouped_results$false_pos_p_value <- p.adjust(grouped_results$false_pos_p_value, method = "fdr", n = num_tests)
  }
  
  # Select the columns to return, including test statistics
  columns_to_return <- c("type", "line_count", "success_rate", "expected_lines", "actual_sample_size", "false_neg", "false_neg_se", "false_neg_ci_lower", "false_neg_ci_upper", "false_neg_test_statistic", "false_neg_p_value", "false_pos", "false_pos_se", "false_pos_ci_lower", "false_pos_ci_upper", "false_pos_test_statistic", "false_pos_p_value")
  return(grouped_results[, columns_to_return, drop = FALSE])
}

##  plotting function
plot_rate_comparison_CI <- function(object_list, y_lim = NULL, plot_type = "ci", global_title = NULL) {
  # Check if the plot_type is valid
  if (!plot_type %in% c("ci", "se")) {
    stop("Invalid plot type. Please specify 'ci' for confidence intervals or 'se' for standard error.")
  }
  
  # Extract names from the object list
  object_names <- names(object_list)
  
  # Extract Numeric Suffixes
  suffixes <- as.numeric(gsub("\\D", "", object_names))
  
  # Create an empty list to store aggregated data
  aggregated_data <- list()
  
  # Extract and Aggregate Data including test statistics and p-values
  for (i in seq_along(object_list)) {
    obj <- object_list[[i]]
    for (type in unique(obj$type)) {
      if (!type %in% names(aggregated_data)) {
        aggregated_data[[type]] <- data.frame(suffix = numeric(), false_neg = numeric(), 
                                              false_neg_se = numeric(), false_neg_ci_lower = numeric(), false_neg_ci_upper = numeric(), 
                                              false_neg_test_statistic = numeric(), false_neg_p_value = numeric(), 
                                              false_pos = numeric(), false_pos_se = numeric(), false_pos_ci_lower = numeric(), false_pos_ci_upper = numeric(),
                                              false_pos_test_statistic = numeric(), false_pos_p_value = numeric())
      }
      type_data <- obj[obj$type == type, ]
      aggregated_data[[type]] <- rbind(aggregated_data[[type]], 
                                       data.frame(suffix = suffixes[i], 
                                                  false_neg = type_data$false_neg, 
                                                  false_neg_se = type_data$false_neg_se, false_neg_ci_lower = type_data$false_neg_ci_lower, false_neg_ci_upper = type_data$false_neg_ci_upper, 
                                                  false_neg_test_statistic = type_data$false_neg_test_statistic, false_neg_p_value = type_data$false_neg_p_value, 
                                                  false_pos = type_data$false_pos, 
                                                  false_pos_se = type_data$false_pos_se, false_pos_ci_lower = type_data$false_pos_ci_lower, false_pos_ci_upper = type_data$false_pos_ci_upper,
                                                  false_pos_test_statistic = type_data$false_pos_test_statistic, false_pos_p_value = type_data$false_pos_p_value))
    }
  }
  
  # Determine global min and max rate values if y_lim is not set
  if (is.null(y_lim)) {
    global_min_rate <- Inf
    global_max_rate <- -Inf
    for (plot_data in aggregated_data) {
      rate_min <- if (plot_type == "ci") min(plot_data$false_neg_ci_lower, plot_data$false_pos_ci_lower, na.rm = TRUE) else 0
      rate_max <- if (plot_type == "ci") max(plot_data$false_neg_ci_upper, plot_data$false_pos_ci_upper, na.rm = TRUE) else 0
      global_min_rate <- min(global_min_rate, rate_min)
      global_max_rate <- max(global_max_rate, rate_max)
    }
    y_lim <- c(global_min_rate, global_max_rate)
  }
  
  # Adjust the width of the error bar lines
  horizontal_line_width <- 0.5
  
  # Adjust the text position below the x-axis
  text_offset <- 0.2 #* (y_lim[2] - y_lim[1])  # Increased offset for more space
  text_cex <- 0.45  # Adjust text size as needed
  
  # Plotting area setup
  par(mfrow = c(2, 2), oma = c(0, 0, 3, 0))  # Set up plotting area for 4 plots with space for global title
  
  # Adjust margins to make room for text annotations
  par(mar = c(7, 4, 4, 2) + 0.1)  # Adjust bottom margin as needed
  
  for (type in names(aggregated_data)) {
    plot_data <- aggregated_data[[type]]
    if (nrow(plot_data) > 0) {
      with(plot_data, {
        # Plot false negatives and false positives
        plot(suffix, false_neg, type = "o", pch = 21, bg = "blue", col = "blue", xlab = "", ylab = "Average Rate",
             main = paste("Type:", type, "tree"), ylim = y_lim, xaxt = "n", yaxt = "n")
        
        # Add the x-axis with default labels
        axis(1, at = suffix, labels = TRUE)
        mtext("Sequence length (kbp)", side = 1, line = 3.5)  # Adjust the 'line' parameter to lower or raise the label
        axis(2)
        
        points(suffix, false_pos, type = "o", pch = 22, bg = "red", col = "red")
        
        # Add error bars based on plot_type
        y_neg_lower <- if (plot_type == "ci") false_neg_ci_lower else false_neg - false_neg_se
        y_neg_upper <- if (plot_type == "ci") false_neg_ci_upper else false_neg + false_neg_se
        y_pos_lower <- if (plot_type == "ci") false_pos_ci_lower else false_pos - false_pos_se
        y_pos_upper <- if (plot_type == "ci") false_pos_ci_upper else false_pos + false_pos_se
        
        # Add error bars for false negatives
        for (i in 1:length(suffix)) {
          lines(c(suffix[i], suffix[i]), c(y_neg_lower[i], y_neg_upper[i]), col = "blue")
          lines(c(suffix[i] - horizontal_line_width, suffix[i] + horizontal_line_width), c(y_neg_upper[i], y_neg_upper[i]), col = "blue") # Horizontal top line
          lines(c(suffix[i] - horizontal_line_width, suffix[i] + horizontal_line_width), c(y_neg_lower[i], y_neg_lower[i]), col = "blue") # Horizontal bottom line
        }
        
        # Add error bars for false positives
        for (i in 1:length(suffix)) {
          lines(c(suffix[i], suffix[i]), c(y_pos_lower[i], y_pos_upper[i]), col = "red")
          lines(c(suffix[i] - horizontal_line_width, suffix[i] + horizontal_line_width), c(y_pos_upper[i], y_pos_upper[i]), col = "red") # Horizontal top line
          lines(c(suffix[i] - horizontal_line_width, suffix[i] + horizontal_line_width), c(y_pos_lower[i], y_pos_lower[i]), col = "red") # Horizontal bottom line
        }
        
        # # Add text for test statistics and p-values below the x-axis
        # for (i in 1:length(suffix)) {
        #   text_x_pos <- suffix[i]
        #   # Annotations for false positives
        #   text(text_x_pos, y = par("usr")[3] - text_offset, labels = paste0("t: ", round(false_pos_test_statistic[i], 2), ", p: ", round(false_pos_p_value[i], 2)), cex = text_cex, adj = c(0.5, 0), xpd = TRUE, col='red')
        #   # Annotations for false negatives
        #   text(text_x_pos, y = par("usr")[3] - 1.2 * text_offset, labels = paste0("t: ", round(false_neg_test_statistic[i], 2), ", p: ", round(false_neg_p_value[i], 2)), cex = text_cex, adj = c(0.5, 0), xpd = TRUE, col='blue')
        # }
        
        # Calculate midpoints between x-positions
        midpoints <- (suffix[-length(suffix)] + suffix[-1]) / 2
        
        # Add text for test statistics and p-values below the x-axis (centered with respect to the tick marks, left-justified with respect to each other)
        for (i in 1:length(suffix)) {
          text_x_pos <- suffix[i]
          
          # Annotations for false positives (centered with respect to the tick mark, left-justified with respect to each other)
          text(text_x_pos, y = par("usr")[3] - text_offset, labels = paste0("t:", round(false_pos_test_statistic[i], 2), ", p:", round(false_pos_p_value[i], 2)), cex = text_cex, adj = c(0.5, 0), xpd = TRUE, col='red')
          
          # Annotations for false negatives (centered with respect to the tick mark, left-justified with respect to each other)
          text(text_x_pos, y = par("usr")[3] - 1.2 * text_offset, labels = paste0("t:", round(false_neg_test_statistic[i], 2), ", p:", round(false_neg_p_value[i], 2)), cex = text_cex, adj = c(0.5, 0), xpd = TRUE, col='blue')
          
          # # Add vertical line to separate annotations (except for the last one)
          # if (i < length(suffix)) {
          #   axis(side = 1, at = midpoints[i], labels = NA, col = "black", lwd = 0.1, pos = -0.275, tck=0.075)
          # }
        }
        
        
        # Add a reference line for a specific rate value
        abline(h=0.1, lty=2, col='black')
        
        # Add legend
        legend("topright", inset = 0.025, legend = c("False Neg", "False Pos"), pch = c(21, 22), pt.bg = c("blue", "red"), col = c("blue", "red"), bty='n')
        
      })
    } else {
      plot(1, type = "n", xlab = "", ylab = "", main = paste("No data for type:", type), ylim = y_lim)
    }
  }
  
  # If a global title is provided, add it to the top of the plot
  if (!is.null(global_title)) {
    mtext(global_title, side = 3, outer = TRUE, line = 1, cex = 1.5, font = 2)
  }
  
  # Reset margin changes for subsequent plots
  par(mar = c(5, 4, 4, 2) + 0.1)
}

##  plotting function with inset reference tree
plot_rate_comparison_CI_inset <- function(object_list, phylogram_list, y_lim = NULL, plot_type = "ci", global_title = NULL, nested=NULL, test_type=NULL) {
  if (!plot_type %in% c("ci", "se")) {
    stop("Invalid plot type. Please specify 'ci' for confidence intervals or 'se' for standard error.")
  }
  
  object_names <- names(object_list)
  suffixes <- as.numeric(gsub("\\D", "", object_names))
  aggregated_data <- list()
  
  for (i in seq_along(object_list)) {
    obj <- object_list[[i]]
    for (type in unique(obj$type)) {
      if (!type %in% names(aggregated_data)) {
        aggregated_data[[type]] <- data.frame(suffix = numeric(), false_neg = numeric(),
                                              false_neg_se = numeric(), false_neg_ci_lower = numeric(), false_neg_ci_upper = numeric(),
                                              false_neg_test_statistic = numeric(), false_neg_p_value = numeric(),
                                              false_pos = numeric(), false_pos_se = numeric(), false_pos_ci_lower = numeric(), false_pos_ci_upper = numeric(),
                                              false_pos_test_statistic = numeric(), false_pos_p_value = numeric())
      }
      type_data <- obj[obj$type == type, ]
      aggregated_data[[type]] <- rbind(aggregated_data[[type]],
                                       data.frame(suffix = suffixes[i],
                                                  false_neg = type_data$false_neg,
                                                  false_neg_se = type_data$false_neg_se, false_neg_ci_lower = type_data$false_neg_ci_lower, false_neg_ci_upper = type_data$false_neg_ci_upper,
                                                  false_neg_test_statistic = type_data$false_neg_test_statistic, false_neg_p_value = type_data$false_neg_p_value,
                                                  false_pos = type_data$false_pos,
                                                  false_pos_se = type_data$false_pos_se, false_pos_ci_lower = type_data$false_pos_ci_lower, false_pos_ci_upper = type_data$false_pos_ci_upper,
                                                  false_pos_test_statistic = type_data$false_pos_test_statistic, false_pos_p_value = type_data$false_pos_p_value))
    }
  }
  
  if (is.null(y_lim)) {
    global_min_rate <- Inf
    global_max_rate <- -Inf
    for (plot_data in aggregated_data) {
      rate_min <- if (plot_type == "ci") min(plot_data$false_neg_ci_lower, plot_data$false_pos_ci_lower, na.rm = TRUE) else 0
      rate_max <- if (plot_type == "ci") max(plot_data$false_neg_ci_upper, plot_data$false_pos_ci_upper, na.rm = TRUE) else 0
      global_min_rate <- min(global_min_rate, rate_min)
      global_max_rate <- max(global_max_rate, rate_max)
    }
    y_lim <- c(global_min_rate, global_max_rate)
  }
  
  horizontal_line_width <- 0.5
  text_offset <- 0.2 # Define text_offset
  text_cex <- 0.45  # Adjust text size as needed
  
  par(mfrow = c(2, 2), oma = c(0, 0, 3, 0))  # Set up plotting area for 4 plots with space for global title
  par(mar = c(7, 4, 4, 2) + 0.1)  
  
  for (type_idx in seq_along(names(aggregated_data))) {
    type <- names(aggregated_data)[type_idx]
    plot_data <- aggregated_data[[type]]
    if (nrow(plot_data) > 0) {
      with(plot_data, {
        plot(suffix, false_neg, type = "o", pch = 21, bg = "blue", col = "blue", xlab = "", ylab = "Average Rate",
             main = paste("Type:", type), ylim = y_lim, xaxt = "n", yaxt = "n", bty="n")
        
        axis(1, at = suffix, labels = TRUE)
        mtext("Sequence length (kbp)", side = 1, line = 3.5)
        axis(2)
        
        points(suffix, false_pos, type = "o", pch = 22, bg = "red", col = "red")
        
        y_neg_lower <- if (plot_type == "ci") false_neg_ci_lower else false_neg - false_neg_se
        y_neg_upper <- if (plot_type == "ci") false_neg_ci_upper else false_neg + false_neg_se
        y_pos_lower <- if (plot_type == "ci") false_pos_ci_lower else false_pos - false_pos_se
        y_pos_upper <- if (plot_type == "ci") false_pos_ci_upper else false_pos + false_pos_se
        
        for (i in 1:length(suffix)) {
          lines(c(suffix[i], suffix[i]), c(y_neg_lower[i], y_neg_upper[i]), col = "blue")
          lines(c(suffix[i] - horizontal_line_width, suffix[i] + horizontal_line_width), c(y_neg_upper[i], y_neg_upper[i]), col = "blue")
          lines(c(suffix[i] - horizontal_line_width, suffix[i] + horizontal_line_width), c(y_neg_lower[i], y_neg_lower[i]), col = "blue")
          lines(c(suffix[i], suffix[i]), c(y_pos_lower[i], y_pos_upper[i]), col = "red")
          lines(c(suffix[i] - horizontal_line_width, suffix[i] + horizontal_line_width), c(y_pos_upper[i], y_pos_upper[i]), col = "red")
          lines(c(suffix[i] - horizontal_line_width, suffix[i] + horizontal_line_width), c(y_pos_lower[i], y_pos_lower[i]), col = "red")
        }
        
        # for (i in 1:length(suffix)) {
        #   text_x_pos <- suffix[i]
        #   text(text_x_pos, y = par("usr")[3] - text_offset, labels = paste0("t: ", round(false_pos_test_statistic[i], 2), ", p: ", format.pval(false_pos_p_value[i], digits = 2)), cex = text_cex, adj = c(0.5, 0), xpd = TRUE, col='red')
        #   text(text_x_pos, y = par("usr")[3] - 1.2 * text_offset, labels = paste0("t: ", round(false_neg_test_statistic[i], 2), ", p: ", format.pval(false_neg_p_value[i], digits = 2)), cex = text_cex, adj = c(0.5, 0), xpd = TRUE, col='blue')
        # }
        
        
        # Calculate midpoints between x-positions
        midpoints <- (suffix[-length(suffix)] + suffix[-1]) / 2
        
        if(test_type=='z'){
        # Add text for test statistics and p-values below the x-axis (centered with respect to the tick marks, left-justified with respect to each other)
        for (i in 1:length(suffix)) {
          text_x_pos <- suffix[i]
          
          # Annotations for false positives (centered with respect to the tick mark, left-justified with respect to each other)
          text(text_x_pos, y = par("usr")[3] - text_offset, labels = paste0("z:", round(false_pos_test_statistic[i], 2), ", p:", round(false_pos_p_value[i], 2)), cex = text_cex, adj = c(0.5, 0), xpd = TRUE, col='red')
          
          # Annotations for false negatives (centered with respect to the tick mark, left-justified with respect to each other)
          text(text_x_pos, y = par("usr")[3] - 1.2 * text_offset, labels = paste0("z:", round(false_neg_test_statistic[i], 2), ", p:", round(false_neg_p_value[i], 2)), cex = text_cex, adj = c(0.5, 0), xpd = TRUE, col='blue')
          
          # # Add vertical line to separate annotations (except for the last one)
          # if (i < length(suffix)) {
          #   axis(side = 1, at = midpoints[i], labels = NA, col = "black", lwd = 0.1, pos = -0.275, tck=0.075)
          # }
        }
        }
        if(test_type=='t'){
          # Add text for test statistics and p-values below the x-axis (centered with respect to the tick marks, left-justified with respect to each other)
          for (i in 1:length(suffix)) {
            text_x_pos <- suffix[i]
            
            # Annotations for false positives (centered with respect to the tick mark, left-justified with respect to each other)
            text(text_x_pos, y = par("usr")[3] - text_offset, labels = paste0("t:", round(false_pos_test_statistic[i], 2), ", p:", round(false_pos_p_value[i], 2)), cex = text_cex, adj = c(0.5, 0), xpd = TRUE, col='red')
            
            # Annotations for false negatives (centered with respect to the tick mark, left-justified with respect to each other)
            text(text_x_pos, y = par("usr")[3] - 1.2 * text_offset, labels = paste0("t:", round(false_neg_test_statistic[i], 2), ", p:", round(false_neg_p_value[i], 2)), cex = text_cex, adj = c(0.5, 0), xpd = TRUE, col='blue')
            
            # # Add vertical line to separate annotations (except for the last one)
            # if (i < length(suffix)) {
            #   axis(side = 1, at = midpoints[i], labels = NA, col = "black", lwd = 0.1, pos = -0.275, tck=0.075)
            # }
          }
        }
        
        abline(h=0.1, lty=2, col='black')
        
        legend("topleft", inset = 0.025, legend = c("False Neg", "False Pos"), pch = c(21, 22), pt.bg = c("blue", "red"), col = c("blue", "red"), bty='n', cex = 0.8)
        
        # First, you need to determine the current plot limits
        usr <- par("usr")
        
        # Define the coordinates for the corners of the rectangle
        # Let's say you want the box to take up the upper right quarter of the plot
        xleft <- (usr[1] + usr[2]) / 2  # midpoint of the x-axis
        ybottom <- (usr[3] + usr[4]) / 2  # midpoint of the y-axis
        xright <- usr[2]  # right limit of the x-axis
        ytop <- usr[4]  # top limit of the y-axis
        
        # Now draw the rectangle
        par(xpd=NA)
        if(nested == F){
          rect(xleft-6.5, ybottom-0.05, xright, ytop, border="black", lwd=0.5)
        } else {
          rect(xleft-13, ybottom-0.05, xright, ytop, border="black", lwd=0.5)
        }
        par(xpd=F)
        
        # Inset phylogram plot
        if (!is.null(phylogram_list[[type_idx]])) {
          par(new=T)
          plot.phylo(ladderize(phylogram_list[[type_idx]]), show.tip.label = FALSE, direction = 'upwards', x.lim=c(-125,200), y.lim=c(-0.65,0.75), edge.width = 0.25, col='lightgrey')
          #box()
          axisPhylo(side=4, backward=F, las=2, lwd.ticks = 0.4, cex.axis = 0.65, lwd=0.4)
        }
        
        legend("topright", legend="Phylogram (substitutions/site)", bty="n", cex=0.7, text.col='black')
       
      })
    } else {
      plot(1, type = "n", xlab = "", ylab = "", main = paste("No data for type:", type), ylim = y_lim)
    }
  }
  
  if (!is.null(global_title)) {
    mtext(global_title, side = 3, outer = TRUE, line = 1, cex = 1.5, font = 2)
  }
  
  par(mar = c(5, 4, 4, 2) + 0.1)
}
plot_rate_comparison_CI_inset_fn <- function(object_list, phylogram_list, y_lim = NULL, plot_type = "ci", global_title = NULL, nested=NULL, include_false_negatives = TRUE, test_type=NULL) {
  if (!plot_type %in% c("ci", "se")) {
    stop("Invalid plot type. Please specify 'ci' for confidence intervals or 'se' for standard error.")
  }
  
  object_names <- names(object_list)
  suffixes <- as.numeric(gsub("\\D", "", object_names))
  aggregated_data <- list()
  
  for (i in seq_along(object_list)) {
    obj <- object_list[[i]]
    for (type in unique(obj$type)) {
      if (!type %in% names(aggregated_data)) {
        aggregated_data[[type]] <- data.frame(suffix = numeric(), false_neg = numeric(),
                                              false_neg_se = numeric(), false_neg_ci_lower = numeric(), false_neg_ci_upper = numeric(),
                                              false_neg_test_statistic = numeric(), false_neg_p_value = numeric(),
                                              false_pos = numeric(), false_pos_se = numeric(), false_pos_ci_lower = numeric(), false_pos_ci_upper = numeric(),
                                              false_pos_test_statistic = numeric(), false_pos_p_value = numeric())
      }
      type_data <- obj[obj$type == type, ]
      aggregated_data[[type]] <- rbind(aggregated_data[[type]],
                                       data.frame(suffix = suffixes[i],
                                                  false_neg = type_data$false_neg,
                                                  false_neg_se = type_data$false_neg_se, false_neg_ci_lower = type_data$false_neg_ci_lower, false_neg_ci_upper = type_data$false_neg_ci_upper,
                                                  false_neg_test_statistic = type_data$false_neg_test_statistic, false_neg_p_value = type_data$false_neg_p_value,
                                                  false_pos = type_data$false_pos,
                                                  false_pos_se = type_data$false_pos_se, false_pos_ci_lower = type_data$false_pos_ci_lower, false_pos_ci_upper = type_data$false_pos_ci_upper,
                                                  false_pos_test_statistic = type_data$false_pos_test_statistic, false_pos_p_value = type_data$false_pos_p_value))
    }
  }
  
  if (is.null(y_lim)) {
    global_min_rate <- Inf
    global_max_rate <- -Inf
    for (plot_data in aggregated_data) {
      rate_min <- if (plot_type == "ci") min(plot_data$false_neg_ci_lower, plot_data$false_pos_ci_lower, na.rm = TRUE) else 0
      rate_max <- if (plot_type == "ci") max(plot_data$false_neg_ci_upper, plot_data$false_pos_ci_upper, na.rm = TRUE) else 0
      global_min_rate <- min(global_min_rate, rate_min)
      global_max_rate <- max(global_max_rate, rate_max)
    }
    y_lim <- c(global_min_rate, global_max_rate)
  }
  
  horizontal_line_width <- 0.5
  text_offset <- 0.2 # Define text_offset
  text_cex <- 0.45  # Adjust text size as needed
  
  par(mfrow = c(2, 2), oma = c(0, 0, 3, 0))  # Set up plotting area for 4 plots with space for global title
  par(mar = c(7, 4, 4, 2) + 0.1)  
  
  for (type_idx in seq_along(names(aggregated_data))) {
    type <- names(aggregated_data)[type_idx]
    plot_data <- aggregated_data[[type]]
    if (nrow(plot_data) > 0) {
      with(plot_data, {
        plot(suffix, false_pos, type = "o", pch = 22, bg = "red", col = "red", xlab = "", ylab = "Average Rate",
             main = paste("Type:", type), ylim = y_lim, xaxt = "n", yaxt = "n", bty="n")
        
        axis(1, at = suffix, labels = TRUE)
        mtext("Sequence length (kbp)", side = 1, line = 3.5)
        axis(2)
        
        y_pos_lower <- if (plot_type == "ci") false_pos_ci_lower else false_pos - false_pos_se
        y_pos_upper <- if (plot_type == "ci") false_pos_ci_upper else false_pos + false_pos_se
        
        for (i in 1:length(suffix)) {
          lines(c(suffix[i], suffix[i]), c(y_pos_lower[i], y_pos_upper[i]), col = "red")
          lines(c(suffix[i] - horizontal_line_width, suffix[i] + horizontal_line_width), c(y_pos_upper[i], y_pos_upper[i]), col = "red")
          lines(c(suffix[i] - horizontal_line_width, suffix[i] + horizontal_line_width), c(y_pos_lower[i], y_pos_lower[i]), col = "red")
        }
        
        if (include_false_negatives) {
          points(suffix, false_neg, type = "o", pch = 21, bg = "blue", col = "blue")
          
          y_neg_lower <- if (plot_type == "ci") false_neg_ci_lower else false_neg - false_neg_se
          y_neg_upper <- if (plot_type == "ci") false_neg_ci_upper else false_neg + false_neg_se
          
          for (i in 1:length(suffix)) {
            lines(c(suffix[i], suffix[i]), c(y_neg_lower[i], y_neg_upper[i]), col = "blue")
            lines(c(suffix[i] - horizontal_line_width, suffix[i] + horizontal_line_width), c(y_neg_upper[i], y_neg_upper[i]), col = "blue")
            lines(c(suffix[i] - horizontal_line_width, suffix[i] + horizontal_line_width), c(y_neg_lower[i], y_neg_lower[i]), col = "blue")
          }
        }
        
        if(test_type=='z'){
          # Add text for test statistics and p-values below the x-axis (centered with respect to the tick marks, left-justified with respect to each other)
          for (i in 1:length(suffix)) {
            text_x_pos <- suffix[i]
            
            # Annotations for false positives (centered with respect to the tick mark, left-justified with respect to each other)
            text(text_x_pos, y = par("usr")[3] - text_offset, labels = paste0("z:", round(false_pos_test_statistic[i], 2), ", p:", round(false_pos_p_value[i], 2)), cex = text_cex, adj = c(0.5, 0), xpd = TRUE, col='red')
            
            # Annotations for false negatives (centered with respect to the tick mark, left-justified with respect to each other)
            #text(text_x_pos, y = par("usr")[3] - 1.2 * text_offset, labels = paste0("z:", round(false_neg_test_statistic[i], 2), ", p:", round(false_neg_p_value[i], 2)), cex = text_cex, adj = c(0.5, 0), xpd = TRUE, col='blue')
            
            # # Add vertical line to separate annotations (except for the last one)
            # if (i < length(suffix)) {
            #   axis(side = 1, at = midpoints[i], labels = NA, col = "black", lwd = 0.1, pos = -0.275, tck=0.075)
            # }
          }
        }
        if(test_type=='t'){
          # Add text for test statistics and p-values below the x-axis (centered with respect to the tick marks, left-justified with respect to each other)
          for (i in 1:length(suffix)) {
            text_x_pos <- suffix[i]
            
            # Annotations for false positives (centered with respect to the tick mark, left-justified with respect to each other)
            text(text_x_pos, y = par("usr")[3] - text_offset, labels = paste0("t:", round(false_pos_test_statistic[i], 2), ", p:", round(false_pos_p_value[i], 2)), cex = text_cex, adj = c(0.5, 0), xpd = TRUE, col='red')
            
            # Annotations for false negatives (centered with respect to the tick mark, left-justified with respect to each other)
            #text(text_x_pos, y = par("usr")[3] - 1.2 * text_offset, labels = paste0("t:", round(false_neg_test_statistic[i], 2), ", p:", round(false_neg_p_value[i], 2)), cex = text_cex, adj = c(0.5, 0), xpd = TRUE, col='blue')
            
            # # Add vertical line to separate annotations (except for the last one)
            # if (i < length(suffix)) {
            #   axis(side = 1, at = midpoints[i], labels = NA, col = "black", lwd = 0.1, pos = -0.275, tck=0.075)
            # }
          }
        }
        
        
        abline(h=0.1, lty=2, col='black')
        
        if (include_false_negatives){
        legend("topleft", inset = 0.025, legend = c("False Neg", "False Pos"), pch = c(21, 22), pt.bg = c("blue", "red"), col = c("blue", "red"), bty='n', cex = 0.8)
        } else {
          legend("topleft", inset = 0.025, legend = c("False Pos"), pch = c(22), pt.bg = c("red"), col = c("red"), bty='n', cex = 0.8)
        }
        
          
        # Define the coordinates for the corners of the rectangle
        usr <- par("usr")
        xleft <- (usr[1] + usr[2]) / 2
        ybottom <- (usr[3] + usr[4]) / 2
        xright <- usr[2]
        ytop <- usr[4]
        
        # Draw the rectangle
        par(xpd=NA)
        if(nested == F){
          rect(xleft-6.5, ybottom-0.05, xright, ytop, border="black", lwd=0.5)
        } else {
          rect(xleft-13, ybottom-0.05, xright, ytop, border="black", lwd=0.5)
        }
        par(xpd=F)
        
        # Inset phylogram plot
        if (!is.null(phylogram_list[[type_idx]])) {
          par(new=T)
          plot.phylo(ladderize(phylogram_list[[type_idx]]), show.tip.label = FALSE, direction = 'upwards', x.lim=c(-125,200), y.lim=c(-0.65,0.75), edge.width = 0.25, col='lightgrey')
          axisPhylo(side=4, backward=F, las=2, lwd.ticks = 0.4, cex.axis = 0.65, lwd=0.4)
        }
        
        legend("topright", legend="Phylogram (substitutions/site)", bty="n", cex=0.7, text.col='black')
      })
    } else {
      plot(1, type = "n", xlab = "", ylab = "", main = paste("No data for type:", type), ylim = y_lim)
    }
  }
  
  if (!is.null(global_title)) {
    mtext(global_title, side = 3, outer = TRUE, line = 1, cex = 1.5, font = 2)
  }
  
  par(mar = c(5, 4, 4, 2) + 0.1)
}


###
#independent shift models

# Simulate a tree using pbtree from the ape package
set.seed(123) # Set a seed for reproducibility
tree <- pbtree(n = 200)
cat("Simulated tree in Newick format:\n", write.tree(tree), "\n\n")

# Example usage
min_clade_size <- 4 #whaver you want
reps <- 1 #or 2, 3, 4, etc


# Run the annotation process with the option for nested or independent shifts
#annotated_newick <- annotate_branches(tree, models(reps), reps, min_clade_size, nested_shifts, annotate_tips=T, buffer=2)
annotated_newick <- annotate_branches(tree, models(reps), reps, min_clade_size=4, nested = F, annotate_tips=T, buffer=2, max_clade_size = 198, derived=F)


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

# #Run the annotation process
# exon_base <- list()
# for(i in 1:10){exon_base[[i]]<-annotate_branches(exontree, models, reps=1, min_clade_size=4, annotate_tips = T, buffer=2, max_clade_size = 198, derived=F)}
# 
# intron_base <- list()
# for(i in 1:10){intron_base[[i]]<-annotate_branches(introntree, models, reps=1, min_clade_size=4, annotate_tips = T, buffer=2, max_clade_size = 198, derived=F)}
# 
# # Run the annotation process
# utr_base <- list()
# for(i in 1:10){utr_base[[i]]<-annotate_branches(utrtree, models, reps=1, min_clade_size=4, annotate_tips = T, buffer=2, max_clade_size = 198, derived=F)}
# 
# # Run the annotation process
# mtdna_base <- list()
# for(i in 1:10){mtdna_base[[i]]<-annotate_branches(mtdnatree, models, reps=1, min_clade_size=4, annotate_tips = T, buffer=2, max_clade_size = 198, derived=F)}
# 
# 
# 
# #user specified model params (10 reps)
# reps <- 2
# set.seed(123)
# #now running in parallel
# exon_base <- pbmclapply(1:10, function(i) {
#   annotate_branches(exontree, models, reps=2, min_clade_size=4, annotate_tips = T, buffer=2, max_clade_size = 198, derived=F)
# }, mc.cores = detectCores()-1) # Number of cores to use
# 
# intron_base <- pbmclapply(1:10, function(i) {
#   annotate_branches(introntree, models, reps=2, min_clade_size=4, annotate_tips = T, buffer=2, max_clade_size = 198, derived=F)
# }, mc.cores = detectCores()-1) # Number of cores to use
# 
# utr_base <- pbmclapply(1:10, function(i) {
#   annotate_branches(utrtree, models, reps=2, min_clade_size=4, annotate_tips = T, buffer=2, max_clade_size = 198, derived=F)
# }, mc.cores = detectCores()-1) # Number of cores to use
# 
# mtdna_base <- pbmclapply(1:10, function(i) {
#   annotate_branches(mtdnatree, models, reps=2, min_clade_size=4, annotate_tips = T, buffer=2, max_clade_size = 198, derived=F)
# }, mc.cores = detectCores()-1) # Number of cores to use


### runing the annotation process in parallel

reps <- 4
set.seed(123)
#now running in parallel
exon_base <- pbmclapply(1:100, function(i) {
  annotate_branches(exontree, models(reps), reps, min_clade_size=4, annotate_tips = T, buffer=2, max_clade_size = 150, derived=F)
}, mc.cores = detectCores()-1) # Number of cores to use

intron_base <- pbmclapply(1:100, function(i) {
  annotate_branches(introntree, models(reps), reps, min_clade_size=4, annotate_tips = T, buffer=2, max_clade_size = 150, derived=F)
}, mc.cores = detectCores()-1) # Number of cores to use

utr_base <- pbmclapply(1:100, function(i) {
  annotate_branches(utrtree, models(reps), reps, min_clade_size=4, annotate_tips = T, buffer=2, max_clade_size = 150, derived=F)
}, mc.cores = detectCores()-1) # Number of cores to use

mtdna_base <- pbmclapply(1:100, function(i) {
  annotate_branches(mtdnatree, models(reps), reps, min_clade_size=4, annotate_tips = T, buffer=2, max_clade_size = 150, derived=F)
}, mc.cores = detectCores()-1) # Number of cores to use


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


### nested shift models
#models <- c("[&model=HKY{2.0}+F{0.45/0.05/0.05/0.45}]", "[&model=HKY{2.0}+F{0.45/0.05/0.05/0.45}]")

# Simulate a tree using pbtree from the ape package
set.seed(123) # Set a seed for reproducibility
tree <- pbtree(n = 20)
cat("Simulated tree in Newick format:\n", write.tree(tree), "\n\n")

# Example usage
min_clade_size <- 4
reps <- 2
#models <- c("[&model=HKY{2.0}+F{0.35/0.15/0.15/0.35}]", "[&model=HKY{2.0}+F{0.15/0.35/0.35/0.15}]")
#models <- c("[&model=HKY{2.0}+F{0.20/0.30/0.30/0.20}]", "[&model=HKY{2.0}+F{0.35/0.15/0.15/0.35}]")
#models <- c("[A]", "[B]")

nested_shifts <- T # Set to TRUE if you want nested shifts

# Run the annotation process with the option for nested or independent shifts
annotated_newick <- annotate_branches(tree, models(reps), reps, min_clade_size, nested_shifts, annotate_tips=T, buffer=2, derived = T)
#annotated_newick <- annotate_branches(tree, models, reps, min_clade_size, nested_shifts, annotate_tips=T, buffer=2, derived = T, max_clade_size = 50)


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


# Run the annotation process
#models <- c("[&model=HKY{2.0}+F{0.35/0.15/0.15/0.35}]", "[&model=HKY{2.0}+F{0.45/0.05/0.05/0.45}]")
#models <- c("[&model=HKY{2.0}+F{0.35/0.15/0.15/0.35}]", "[&model=HKY{2.0}+F{0.49/0.01/0.01/0.49}]")

# 
# # Run the annotation process
# set.seed(123)
# exon_base <- list()
# for(i in 1:100){exon_base[[i]]<-annotate_branches(exontree, models, reps=2, min_clade_size=8, nested=T, annotate_tips=T, buffer=2, derived=T, max_clade_size=150)} #length(exon_base$tip.label)/2
# 
# # Run the annotation process
# intron_base <- list()
# for(i in 1:100){intron_base[[i]]<-annotate_branches(introntree, models, reps=2, min_clade_size=8, nested=T, annotate_tips=T, buffer=2, derived=T, max_clade_size=150)}
# 
# # Run the annotation process
# utr_base <- list()
# for(i in 1:100){utr_base[[i]]<-annotate_branches(utrtree, models, reps=2, min_clade_size=8, nested=T, annotate_tips=T, buffer=2, derived=T, max_clade_size=150)}
# 
# # Run the annotation process
# mtdna_base <- list()
# for(i in 1:100){mtdna_base[[i]]<-annotate_branches(mtdnatree, models, reps=2, min_clade_size=8, nested=T, annotate_tips=T, buffer=2, derived=T, max_clade_size=150)}


##  dirichlet 100 nested -- with dirichlet root

reps <- 2
set.seed(123)
#now running in parallel
exon_base <- pbmclapply(1:100, function(i) {
  annotate_branches(exontree, models(reps), reps=2, min_clade_size=4, nested=T, annotate_tips=T, buffer=2, derived=F, max_clade_size=150)
}, mc.cores = detectCores()-1) # Number of cores to use

intron_base <- pbmclapply(1:100, function(i) {
  annotate_branches(introntree, models(reps), reps=2, min_clade_size=4, nested=T, annotate_tips=T, buffer=2, derived=F, max_clade_size=150)
}, mc.cores = detectCores()-1) # Number of cores to use

utr_base <- pbmclapply(1:100, function(i) {
  annotate_branches(utrtree, models(reps), reps=2, min_clade_size=4, nested=T, annotate_tips=T, buffer=2, derived=F, max_clade_size=150)
}, mc.cores = detectCores()-1) # Number of cores to use

mtdna_base <- pbmclapply(1:100, function(i) {
  annotate_branches(mtdnatree, models(reps), reps=2, min_clade_size=4, nested=T, annotate_tips=T, buffer=2, derived=F, max_clade_size=150)
}, mc.cores = detectCores()-1) # Number of cores to use

writeLines(models(400), 'models.txt')

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


##### FINAL RUNS BELOW THIS #####


#t test version
{
  ### dirichlet all root models, no shifts
  shift_zero_2k <- parse_files_hmshj_CI_t_test(expected=1, group_by_reps = T, group_by_type = T, directory = '/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/simulated_alignments/false_negatives/second_attempt/output/fp_random_2k/fp_random_2k', confidence_level = 0.95, p_value_correction = 'fdr')
  shift_zero_10k <- parse_files_hmshj_CI_t_test(expected=1, group_by_reps = T, group_by_type = T, directory = '/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/simulated_alignments/false_negatives/second_attempt/output/fp_random_10k/fp_attempt2_hmshj', confidence_level = 0.95, p_value_correction = 'fdr')
  shift_zero_50k <- parse_files_hmshj_CI_t_test(expected=1, group_by_reps = T, group_by_type = T, directory = '/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/simulated_alignments/false_negatives/second_attempt/output/fp_random_50k', confidence_level = 0.95, p_value_correction = 'fdr')
  
  ### dirichlet all shifts, independent
  shift_one_2k <- parse_files_hmshj_CI_t_test(expected=2, group_by_reps = T, group_by_type = T, directory = '/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/simulated_alignments/false_negatives/third_attempt/fn/independent/output/1_shift_2k', confidence_level = 0.95, p_value_correction = 'fdr')
  shift_one_10k <- parse_files_hmshj_CI_t_test(expected=2, group_by_reps = T, group_by_type = T, directory = '/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/simulated_alignments/false_negatives/third_attempt/fn/independent/output/1_shift_10k', confidence_level = 0.95, p_value_correction = 'fdr')
  shift_one_25k <- parse_files_hmshj_CI_t_test(expected=2, group_by_reps = T, group_by_type = T, directory = '/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/simulated_alignments/false_negatives/third_attempt/fn/independent/output/1_shift_25k', confidence_level = 0.95, p_value_correction = 'fdr')
  shift_one_50k <- parse_files_hmshj_CI_t_test(expected=2, group_by_reps = T, group_by_type = T, directory = '/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/simulated_alignments/false_negatives/third_attempt/fn/independent/output/1_shift_50k', confidence_level = 0.95, p_value_correction = 'fdr')
  
  shift_two_2k <- parse_files_hmshj_CI_t_test(expected=3, group_by_reps = T, group_by_type = T, directory = '/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/simulated_alignments/false_negatives/third_attempt/fn/independent/output/2_shift_2k', confidence_level = 0.95, p_value_correction = 'fdr')
  shift_two_10k <- parse_files_hmshj_CI_t_test(expected=3, group_by_reps = T, group_by_type = T, directory = '/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/simulated_alignments/false_negatives/third_attempt/fn/independent/output/2_shift_10k', confidence_level = 0.95, p_value_correction = 'fdr')
  shift_two_25k <- parse_files_hmshj_CI_t_test(expected=3, group_by_reps = T, group_by_type = T, directory = '/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/simulated_alignments/false_negatives/third_attempt/fn/independent/output/2_shift_25k', confidence_level = 0.95, p_value_correction = 'fdr')
  shift_two_50k <- parse_files_hmshj_CI_t_test(expected=3, group_by_reps = T, group_by_type = T, directory = '/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/simulated_alignments/false_negatives/third_attempt/fn/independent/output/2_shift_50k', confidence_level = 0.95, p_value_correction = 'fdr')

  shift_three_2k <- parse_files_hmshj_CI_t_test(expected=4, group_by_reps = T, group_by_type = T, directory = '/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/simulated_alignments/false_negatives/third_attempt/fn/independent/output/3_shift_2k', confidence_level = 0.95, p_value_correction = 'fdr')
  shift_three_10k <- parse_files_hmshj_CI_t_test(expected=4, group_by_reps = T, group_by_type = T, directory = '/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/simulated_alignments/false_negatives/third_attempt/fn/independent/output/3_shift_10k', confidence_level = 0.95, p_value_correction = 'fdr')
  shift_three_25k <- parse_files_hmshj_CI_t_test(expected=4, group_by_reps = T, group_by_type = T, directory = '/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/simulated_alignments/false_negatives/third_attempt/fn/independent/output/3_shift_25k', confidence_level = 0.95, p_value_correction = 'fdr')
  shift_three_50k <- parse_files_hmshj_CI_t_test(expected=4, group_by_reps = T, group_by_type = T, directory = '/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/simulated_alignments/false_negatives/third_attempt/fn/independent/output/3_shift_50k', confidence_level = 0.95, p_value_correction = 'fdr')

  shift_four_2k <- parse_files_hmshj_CI_t_test(expected=5, group_by_reps = T, group_by_type = T, directory = '/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/simulated_alignments/false_negatives/third_attempt/fn/independent/output/4_shift_2k', confidence_level = 0.95, p_value_correction = 'fdr')
  shift_four_10k <- parse_files_hmshj_CI_t_test(expected=5, group_by_reps = T, group_by_type = T, directory = '/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/simulated_alignments/false_negatives/third_attempt/fn/independent/output/4_shift_10k', confidence_level = 0.95, p_value_correction = 'fdr')
  shift_four_25k <- parse_files_hmshj_CI_t_test(expected=5, group_by_reps = T, group_by_type = T, directory = '/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/simulated_alignments/false_negatives/third_attempt/fn/independent/output/4_shift_25k', confidence_level = 0.95, p_value_correction = 'fdr')
  shift_four_50k <- parse_files_hmshj_CI_t_test(expected=5, group_by_reps = T, group_by_type = T, directory = '/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/simulated_alignments/false_negatives/third_attempt/fn/independent/output/4_shift_50k', confidence_level = 0.95, p_value_correction = 'fdr')
  
  ##nested
  
  #dirichlet all - eg every model is a dirichlet sample
  dirichlet_all_2k <- parse_files_hmshj_CI_t_test(expected=3, group_by_reps = T, group_by_type = T, directory = '/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/simulated_alignments/false_negatives/third_attempt/fn/nested/dirichlet_all/output/nested_dirichlet_all_2k_min4_buffer2_max150_HKY2', confidence_level = 0.95, p_value_correction = 'fdr')
  dirichlet_all_25k <- parse_files_hmshj_CI_t_test(expected=3, group_by_reps = T, group_by_type = T, directory = '/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/simulated_alignments/false_negatives/third_attempt/fn/nested/dirichlet_all/output/nested_dirichlet_all_25k_min4_buffer2_max150_HKY2', confidence_level = 0.95, p_value_correction = 'fdr')
  dirichlet_all_50k <- parse_files_hmshj_CI_t_test(expected=3, group_by_reps = T, group_by_type = T, directory = '/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/simulated_alignments/false_negatives/third_attempt/fn/nested/dirichlet_all/output/nested_dirichlet_all_50k_min4_buffer2_max150_HKY2', confidence_level = 0.95, p_value_correction = 'fdr')
  dirichlet_all_100k <- parse_files_hmshj_CI_t_test(expected=3, group_by_reps = T, group_by_type = T, directory = '/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/simulated_alignments/false_negatives/third_attempt/fn/nested/dirichlet_all/output/nested_dirichlet_all_100k_min4_buffer2_max150_HKY2', confidence_level = 0.95, p_value_correction = 'fdr')
  
  
  #fixing names
  shift_zero_2k$type <- shift_one_2k$type
  shift_zero_10k$type <- shift_one_2k$type
  shift_zero_50k$type <- shift_one_2k$type
}


#final plot
pdf(height=9.5*0.9, width=7.5*0.9, file='FP_FN_t.pdf')
plot_rate_comparison_CI_inset_fn(list(shift_zero_2k=shift_zero_2k, shift_zero_10k=shift_zero_10k, shift_zero_50k=shift_zero_50k), phylogram_list = list(exontree, introntree, mtdnatree, utrtree), y_lim = c(0,1.0), plot_type = 'se', global_title = "No shifts", nested=F, include_false_negatives = F, test_type = 't')
plot_rate_comparison_CI_inset(list(shift_one_2k=shift_one_2k, shift_one_10k=shift_one_10k, shift_one_25k=shift_one_25k, shift_one_50k=shift_one_50k), phylogram_list = list(exontree, introntree, mtdnatree, utrtree), y_lim = c(0,1.0), plot_type = 'se', global_title = "One derived shift", nested=F, test_type = 't')
plot_rate_comparison_CI_inset(list(shift_two_2k=shift_two_2k, shift_two_10k=shift_two_10k, shift_two_25k=shift_two_25k, shift_two_50k=shift_two_50k), phylogram_list = list(exontree, introntree, mtdnatree, utrtree), y_lim = c(0,1.0), plot_type = 'se', global_title = "Two derived shifts (independent)", nested=F, test_type = 't')
plot_rate_comparison_CI_inset(list(shift_three_2k=shift_three_2k, shift_three_10k=shift_three_10k, shift_three_25k=shift_three_25k, shift_three_50k=shift_three_50k), phylogram_list = list(exontree, introntree, mtdnatree, utrtree), y_lim = c(0,1.0), plot_type = 'se', global_title = "Three derived shifts (independent)", nested=F, test_type = 't')
plot_rate_comparison_CI_inset(list(shift_four_2k=shift_four_2k, shift_four_10k=shift_four_10k, shift_four_25k=shift_four_25k, shift_four_50k=shift_four_50k), phylogram_list = list(exontree, introntree, mtdnatree, utrtree), y_lim = c(0,1.0), plot_type = 'se', global_title = "Four derived shifts (independent)", nested=F, test_type = 't')

#nested shifts
plot_rate_comparison_CI_inset(list(dirichlet_all_2k=dirichlet_all_2k, dirichlet_all_25k=dirichlet_all_25k, dirichlet_all_50k=dirichlet_all_50k,dirichlet_all_100k=dirichlet_all_100k), phylogram_list = list(exontree, introntree, mtdnatree, utrtree), y_lim = c(0,1.0), plot_type = 'se', global_title = "Nested shifts", nested=T, test_type = 't')
dev.off()


#z test version
{
  ### dirichlet all root models, no shifts
  shift_zero_2k <- parse_files_hmshj_CI_z_test(expected=1, group_by_reps = T, group_by_type = T, directory = '/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/simulated_alignments/false_negatives/second_attempt/output/fp_random_2k/fp_random_2k', confidence_level = 0.95, p_value_correction = 'fdr')
  shift_zero_10k <- parse_files_hmshj_CI_z_test(expected=1, group_by_reps = T, group_by_type = T, directory = '/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/simulated_alignments/false_negatives/second_attempt/output/fp_random_10k/fp_attempt2_hmshj', confidence_level = 0.95, p_value_correction = 'fdr')
  shift_zero_50k <- parse_files_hmshj_CI_z_test(expected=1, group_by_reps = T, group_by_type = T, directory = '/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/simulated_alignments/false_negatives/second_attempt/output/fp_random_50k', confidence_level = 0.95, p_value_correction = 'fdr')
  
  ### dirichlet all shifts, independent
  shift_one_2k <- parse_files_hmshj_CI_z_test(expected=2, group_by_reps = T, group_by_type = T, directory = '/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/simulated_alignments/false_negatives/third_attempt/fn/independent/output/1_shift_2k', confidence_level = 0.95, p_value_correction = 'fdr')
  shift_one_10k <- parse_files_hmshj_CI_z_test(expected=2, group_by_reps = T, group_by_type = T, directory = '/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/simulated_alignments/false_negatives/third_attempt/fn/independent/output/1_shift_10k', confidence_level = 0.95, p_value_correction = 'fdr')
  shift_one_25k <- parse_files_hmshj_CI_z_test(expected=2, group_by_reps = T, group_by_type = T, directory = '/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/simulated_alignments/false_negatives/third_attempt/fn/independent/output/1_shift_25k', confidence_level = 0.95, p_value_correction = 'fdr')
  shift_one_50k <- parse_files_hmshj_CI_z_test(expected=2, group_by_reps = T, group_by_type = T, directory = '/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/simulated_alignments/false_negatives/third_attempt/fn/independent/output/1_shift_50k', confidence_level = 0.95, p_value_correction = 'fdr')
  
  shift_two_2k <- parse_files_hmshj_CI_z_test(expected=3, group_by_reps = T, group_by_type = T, directory = '/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/simulated_alignments/false_negatives/third_attempt/fn/independent/output/2_shift_2k', confidence_level = 0.95, p_value_correction = 'fdr')
  shift_two_10k <- parse_files_hmshj_CI_z_test(expected=3, group_by_reps = T, group_by_type = T, directory = '/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/simulated_alignments/false_negatives/third_attempt/fn/independent/output/2_shift_10k', confidence_level = 0.95, p_value_correction = 'fdr')
  shift_two_25k <- parse_files_hmshj_CI_z_test(expected=3, group_by_reps = T, group_by_type = T, directory = '/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/simulated_alignments/false_negatives/third_attempt/fn/independent/output/2_shift_25k', confidence_level = 0.95, p_value_correction = 'fdr')
  shift_two_50k <- parse_files_hmshj_CI_z_test(expected=3, group_by_reps = T, group_by_type = T, directory = '/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/simulated_alignments/false_negatives/third_attempt/fn/independent/output/2_shift_50k', confidence_level = 0.95, p_value_correction = 'fdr')
  
  shift_three_2k <- parse_files_hmshj_CI_z_test(expected=4, group_by_reps = T, group_by_type = T, directory = '/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/simulated_alignments/false_negatives/third_attempt/fn/independent/output/3_shift_2k', confidence_level = 0.95, p_value_correction = 'fdr')
  shift_three_10k <- parse_files_hmshj_CI_z_test(expected=4, group_by_reps = T, group_by_type = T, directory = '/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/simulated_alignments/false_negatives/third_attempt/fn/independent/output/3_shift_10k', confidence_level = 0.95, p_value_correction = 'fdr')
  shift_three_25k <- parse_files_hmshj_CI_z_test(expected=4, group_by_reps = T, group_by_type = T, directory = '/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/simulated_alignments/false_negatives/third_attempt/fn/independent/output/3_shift_25k', confidence_level = 0.95, p_value_correction = 'fdr')
  shift_three_50k <- parse_files_hmshj_CI_z_test(expected=4, group_by_reps = T, group_by_type = T, directory = '/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/simulated_alignments/false_negatives/third_attempt/fn/independent/output/3_shift_50k', confidence_level = 0.95, p_value_correction = 'fdr')
  
  shift_four_2k <- parse_files_hmshj_CI_z_test(expected=5, group_by_reps = T, group_by_type = T, directory = '/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/simulated_alignments/false_negatives/third_attempt/fn/independent/output/4_shift_2k', confidence_level = 0.95, p_value_correction = 'fdr')
  shift_four_10k <- parse_files_hmshj_CI_z_test(expected=5, group_by_reps = T, group_by_type = T, directory = '/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/simulated_alignments/false_negatives/third_attempt/fn/independent/output/4_shift_10k', confidence_level = 0.95, p_value_correction = 'fdr')
  shift_four_25k <- parse_files_hmshj_CI_z_test(expected=5, group_by_reps = T, group_by_type = T, directory = '/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/simulated_alignments/false_negatives/third_attempt/fn/independent/output/4_shift_25k', confidence_level = 0.95, p_value_correction = 'fdr')
  shift_four_50k <- parse_files_hmshj_CI_z_test(expected=5, group_by_reps = T, group_by_type = T, directory = '/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/simulated_alignments/false_negatives/third_attempt/fn/independent/output/4_shift_50k', confidence_level = 0.95, p_value_correction = 'fdr')
  
  ##nested
  
  #dirichlet all - eg every model is a dirichlet sample
  dirichlet_all_2k <- parse_files_hmshj_CI_z_test(expected=3, group_by_reps = T, group_by_type = T, directory = '/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/simulated_alignments/false_negatives/third_attempt/fn/nested/dirichlet_all/output/nested_dirichlet_all_2k_min4_buffer2_max150_HKY2', confidence_level = 0.95, p_value_correction = 'fdr')
  dirichlet_all_25k <- parse_files_hmshj_CI_z_test(expected=3, group_by_reps = T, group_by_type = T, directory = '/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/simulated_alignments/false_negatives/third_attempt/fn/nested/dirichlet_all/output/nested_dirichlet_all_25k_min4_buffer2_max150_HKY2', confidence_level = 0.95, p_value_correction = 'fdr')
  dirichlet_all_50k <- parse_files_hmshj_CI_z_test(expected=3, group_by_reps = T, group_by_type = T, directory = '/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/simulated_alignments/false_negatives/third_attempt/fn/nested/dirichlet_all/output/nested_dirichlet_all_50k_min4_buffer2_max150_HKY2', confidence_level = 0.95, p_value_correction = 'fdr')
  dirichlet_all_100k <- parse_files_hmshj_CI_z_test(expected=3, group_by_reps = T, group_by_type = T, directory = '/Users/cotinga/jsb439@cornell.edu/AnchoredEnrichment/bird2020/berv_alignments/simulated_alignments/false_negatives/third_attempt/fn/nested/dirichlet_all/output/nested_dirichlet_all_100k_min4_buffer2_max150_HKY2', confidence_level = 0.95, p_value_correction = 'fdr')
  
  #fixing names
  shift_zero_2k$type <- shift_one_2k$type
  shift_zero_10k$type <- shift_one_2k$type
  shift_zero_50k$type <- shift_one_2k$type
}

#final plot
pdf(height=9.5*0.9, width=7.5*0.9, file='FP_FN_z.pdf')
plot_rate_comparison_CI_inset_fn(list(shift_zero_2k=shift_zero_2k, shift_zero_10k=shift_zero_10k, shift_zero_50k=shift_zero_50k), phylogram_list = list(exontree, introntree, mtdnatree, utrtree), y_lim = c(0,1.0), plot_type = 'se', global_title = "No shifts", nested=F, include_false_negatives = F, test_type = 'z')
plot_rate_comparison_CI_inset(list(shift_one_2k=shift_one_2k, shift_one_10k=shift_one_10k, shift_one_25k=shift_one_25k, shift_one_50k=shift_one_50k), phylogram_list = list(exontree, introntree, mtdnatree, utrtree), y_lim = c(0,1.0), plot_type = 'se', global_title = "One derived shift", nested=F, test_type = 'z')
plot_rate_comparison_CI_inset(list(shift_two_2k=shift_two_2k, shift_two_10k=shift_two_10k, shift_two_25k=shift_two_25k, shift_two_50k=shift_two_50k), phylogram_list = list(exontree, introntree, mtdnatree, utrtree), y_lim = c(0,1.0), plot_type = 'se', global_title = "Two derived shifts (independent)", nested=F, test_type = 'z')
plot_rate_comparison_CI_inset(list(shift_three_2k=shift_three_2k, shift_three_10k=shift_three_10k, shift_three_25k=shift_three_25k, shift_three_50k=shift_three_50k), phylogram_list = list(exontree, introntree, mtdnatree, utrtree), y_lim = c(0,1.0), plot_type = 'se', global_title = "Three derived shifts (independent)", nested=F, test_type = 'z')
plot_rate_comparison_CI_inset(list(shift_four_2k=shift_four_2k, shift_four_10k=shift_four_10k, shift_four_25k=shift_four_25k, shift_four_50k=shift_four_50k), phylogram_list = list(exontree, introntree, mtdnatree, utrtree), y_lim = c(0,1.0), plot_type = 'se', global_title = "Four derived shifts (independent)", nested=F, test_type = 'z')

#nested shifts
plot_rate_comparison_CI_inset(list(dirichlet_all_2k=dirichlet_all_2k, dirichlet_all_25k=dirichlet_all_25k, dirichlet_all_50k=dirichlet_all_50k,dirichlet_all_100k=dirichlet_all_100k), phylogram_list = list(exontree, introntree, mtdnatree, utrtree), y_lim = c(0,1.0), plot_type = 'se', global_title = "Nested shifts", nested=T, test_type = 'z')
dev.off()


