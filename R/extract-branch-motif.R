#' Extract Branch Motifs
#'
#' @description
#' An N-glycan branching motif if the substructures linked to either the a3- or a6-core-mannose.
#'
#' For example:
#'
#' ```
#' Neu5Ac - Gal - GlcNAc - Man
#' ~~~~~~~~~~~~~~~~~~~~~      \
#'   A branching motif         Man - GlcNAc - GlcNAc -
#'                            /
#'                         Man
#' ```
#'
#' This function returns all the unique branching motifs found in the input glycans.
#' It is very useful for analyzing N-glycan antennary patterns combined with `glydet::quantify_motifs()`.
#'
#' @param glycans One of:
#'   - A [glyrepr::glycan_structure()] vector.
#'   - A glycan structure string vector. All formats supported by [glyparse::auto_parse()] are accepted.
#'
#' @details
#' The function works by:
#' 1. Converting the input to a set of unique `glycan_structure` objects.
#' 2. Searching for the N-glycan branch pattern: `HexNAc(??-?)Hex(??-?)Hex(??-?)HexNAc(??-?)HexNAc(??-`.
#' 3. For each match, identifying the root node of the branch (the leftmost HexNAc in the pattern).
#' 4. Extracting the full subtree rooted at that node.
#' 5. Preserving the correct anomeric configuration (e.g., "b1") by inspecting the linkage to the root node.
#'
#' @return A [glyrepr::glycan_structure()] vector containing the unique extracted branching motifs.
#'
#' @examples
#' glycans <- c(
#'   "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(a1-4)GlcNAc(b1-",
#'   "Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(a1-4)GlcNAc(b1-"
#' )
#' extract_branch_motif(glycans)
#'
#' @export
extract_branch_motif <- function(glycans) {
  # 1. Handle input types and deduplicate
  glycans <- ensure_glycans_are_structures(glycans)
  glycans <- unique(glycans)
  .assert_n_glycans(glycans)

  # 2. Define the motif
  # The motif represents an N-glycan branch, rooted at the GlcNAc attached to the Mannose core.
  motif <- glyrepr::as_glycan_structure("HexNAc(??-?)Hex(??-?)Hex(??-?)HexNAc(??-?)HexNAc(??-")

  # 3. Find matches
  # Use match_motif to find where the motif maps to the glycan
  matches <- match_motif(glycans, motif, alignment = "core")

  # 4. Extract subtrees
  # matches is a list of lists: [[glycan_1_matches], [glycan_2_matches], ...]
  # attributes(glycans, "structures") gives us the actual graphs

  structure_graphs <- glyrepr::get_structure_graphs(glycans, return_list = TRUE)

  extracted_subtrees <- list()

  for (i in seq_along(structure_graphs)) {
    g <- structure_graphs[[i]]
    g_matches <- matches[[i]]

    if (length(g_matches) == 0) next

    for (match_idx in g_matches) {
      # The match_idx is a vector of node indices in 'g' corresponding to nodes in 'motif'.
      # The first node of the motif (HexNAc) is the root of the branch.
      # We assume the motif string was parsed such that the first node is indeed the root HexNAc.
      # In glyrepr/igraph, nodes are 1-indexed.

      branch_root_id <- match_idx[1]

      # Extract the subtree rooted at branch_root_id
      # We use igraph::subcomponent with mode = "out" to get all descendants
      subtree_nodes <- igraph::subcomponent(g, branch_root_id, mode = "out")
      subtree <- igraph::induced_subgraph(g, subtree_nodes)

      # 5. Handle Anomer
      # The subtree needs an 'anomer' attribute (e.g., "b1").
      # This information is contained in the linkage of the INCOMING edge to branch_root_id in 'g'.

      # Find the incoming edge to the branch root in the original graph
      in_edges <- igraph::incident(g, branch_root_id, mode = "in")

      # There should be exactly one incoming edge for a tree structure (unless it's the absolute root)
      linkage <- in_edges$linkage
      # Linkage format is like "b1-4" or "a1-3" or "??-?".
      # We need to extract the anomer part (first 2 chars usually, e.g. "b1", "a1", "??")

      # Simple parsing logic (assuming standard format like "a1-..." or "b1-...")
      # We can split by "-"
      anomer_part <- stringr::str_split_i(linkage, "-", 1)
      subtree$anomer <- anomer_part

      extracted_subtrees[[length(extracted_subtrees) + 1]] <- subtree
    }
  }

  if (length(extracted_subtrees) == 0) {
    return(glyrepr::glycan_structure())
  }

  # Convert extracted graphs back to glycan_structure and return unique ones
  res <- glyrepr::as_glycan_structure(extracted_subtrees)
  unique(res)
}

.assert_n_glycans <- function(glycans) {
  n_motif <- glyrepr::as_glycan_structure("Hex(??-?)[Hex(??-?)]Hex(??-?)HexNAc(??-?)HexNAc(??-")
  have_n_motif <- have_motif(glycans, n_motif, alignment = "core")
  if (!all(have_n_motif)) {
    cli::cli_abort(c(
      "{.arg glycans} must be N-glycans.",
      "x" = "Some of {.arg glycans} do not have the N-glycan core motif."
    ))
  }
}