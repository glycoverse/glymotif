#' Extract All Substructures (Motifs)
#'
#' Extract all unique connected subgraphs (motifs) from the input glycans up to a specified size.
#' This function can be useful combined with `count_motifs()` or `glydet::quantify_motifs()`.
#' If so, set `alignment` to "substructure" for these functions.
#'
#' @param glycans One of:
#'   - A [glyrepr::glycan_structure()] vector.
#'   - A glycan structure string vector. All formats supported by [glyparse::auto_parse()] are accepted.
#' @param max_size The maximum number of monosaccharides in the extracted motifs. Default is 3.
#'   Note that setting this value very large can be computationally expensive.
#'   Try the default value first, and increase it progressively if needed.
#'
#' @return A [glyrepr::glycan_structure()] vector containing the unique extracted motifs.
#'
#' @examples
#' glycan <- "Gal(b1-3)[GlcNAc(a1-6)]GalNAc(a1-"
#' extract_motif(glycan, max_size = 2)
#'
#' @importFrom igraph V E induced_subgraph incident edge_attr neighbors
#' @export
extract_motif <- function(glycans, max_size = 3) {
  glycans <- ensure_glycans_are_structures(glycans)
  glycans <- unique(glycans)
  structure_graphs <- glyrepr::get_structure_graphs(glycans, return_list = TRUE)

  extracted_subtrees <- list()

  for (g in structure_graphs) {
    # For each node in the graph, treat it as the root of a potential motif
    nodes <- as.numeric(igraph::V(g))

    for (root_id in nodes) {
      # Find all connected subgraphs rooted at 'root_id' with size <= max_size
      # We only consider "out" edges (descendants) to form the subgraph because
      # glycan structures are directed trees.

      # We use a recursive function to find all valid subsets of descendants
      subsets <- .find_connected_subsets(g, root_id, max_size)

      for (subset_nodes in subsets) {
        # Induce subgraph
        subtree <- igraph::induced_subgraph(g, subset_nodes)

        # Handle Anomer
        # The subtree needs an 'anomer' attribute.
        # If root_id is the root of g, use g$anomer.
        # Otherwise, use the linkage of the incoming edge to root_id in g.

        # Check if root_id is a root in g (indegree == 0)
        # Note: In glyrepr, usually there is only one root (reducing end).
        # But let's check properly.
        in_edges <- igraph::incident(g, root_id, mode = "in")

        if (length(in_edges) == 0) {
          # It is the root of the original glycan
          subtree$anomer <- g$anomer
        } else {
          # It is an internal node. Extract anomer from incoming linkage.
          linkage <- igraph::edge_attr(g, "linkage", in_edges)
          # Linkage format: "b1-4". We want "b1".
          anomer_part <- stringr::str_split_i(linkage, "-", 1)
          subtree$anomer <- anomer_part
        }

        extracted_subtrees[[length(extracted_subtrees) + 1]] <- subtree
      }
    }
  }

  if (length(extracted_subtrees) == 0) {
    return(glyrepr::glycan_structure())
  }

  res <- glyrepr::as_glycan_structure(extracted_subtrees)
  unique(res)
}

# Helper function to find connected subsets rooted at 'node'
# Returns a list of vectors, each vector containing node IDs.
.find_connected_subsets <- function(g, root, max_size) {
  # Current set is just {root}
  results <- list(c(root))

  if (max_size <= 1) {
    return(results)
  }

  # Find children of root
  children <- as.numeric(igraph::neighbors(g, root, mode = "out"))

  if (length(children) == 0) {
    return(results)
  }

  # We need to pick subsets of children to extend the current component.
  # This is a bit complex: we can pick child A and some of its descendants,
  # AND child B and some of its descendants.
  # total size <= max_size.

  # Recursive strategy:
  # For each child, get all its valid connected subsets (size <= max_size - 1).
  # Then combine these subsets such that total size <= max_size (including root).

  child_subsets_list <- list()
  for (child in children) {
    child_subsets_list[[as.character(child)]] <- .find_connected_subsets(g, child, max_size - 1)
  }

  # Now we need to combine these options.
  # Problem: "Knapsack"-like combination.
  # We have K children. Each child i offers a set of Subsets S_i.
  # We can choose one subset from S_i, OR choose nothing from child i.
  # Constraint: 1 (root) + sum(|chosen_subset_i|) <= max_size.

  # Since max_size is small (default 3), we can brute-force combination.

  # Simplify:
  # Create a list of "choices" for each child.
  # Choice 0: pick nothing (nodes = empty, size = 0).
  # Choice 1..M: pick one of the subsets from child_subsets_list.

  choices_per_child <- list()
  for (child in children) {
    subsets <- child_subsets_list[[as.character(child)]]
    # Add "empty" choice
    choices <- list(numeric(0))
    choices <- c(choices, subsets)
    choices_per_child[[length(choices_per_child) + 1]] <- choices
  }

  # Cartesian product of choices
  all_combinations <- expand.grid(purrr::map(choices_per_child, seq_along))

  # Iterate all combinations
  for (i in seq_len(nrow(all_combinations))) {
    combination_idx <- as.numeric(all_combinations[i, ])

    current_nodes <- c(root)
    current_size <- 1

    for (j in seq_along(combination_idx)) {
      choice_idx <- combination_idx[j]
      chosen_subset <- choices_per_child[[j]][[choice_idx]]

      size_incr <- length(chosen_subset)
      if (size_incr > 0) {
        current_nodes <- c(current_nodes, chosen_subset)
        current_size <- current_size + size_incr
      }
    }

    if (current_size <= max_size && current_size > 1) {
      # Don't add size 1 again (already added initially)
      results[[length(results) + 1]] <- current_nodes
    }
  }

  results
}
