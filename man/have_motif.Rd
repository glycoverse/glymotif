% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/have-motif.R, R/vectorized.R
\name{have_motif}
\alias{have_motif}
\alias{have_motifs}
\title{Check if the Glycans have the Given Motif(s)}
\usage{
have_motif(glycans, motif, alignment = NULL, ignore_linkages = FALSE)

have_motifs(glycans, motifs, alignments = NULL, ignore_linkages = FALSE)
}
\arguments{
\item{glycans}{A 'glyrepr_structure' object, or an IUPAC-condensed structure string vector.}

\item{motif}{A 'glyrepr_structure' object, an IUPAC-condensed structure string,
or a known motif name (use \code{\link[=available_motifs]{available_motifs()}} to see all available motifs).}

\item{alignment}{A character string.
Possible values are "substructure", "core", "terminal" and "whole".
If not provided, the value will be decided based on the \code{motif} argument.
If \code{motif} is a motif name, the alignment in the database will be used.
Otherwise, "substructure" will be used.}

\item{ignore_linkages}{A logical value. If \code{TRUE}, linkages will be ignored in the comparison.}

\item{motifs}{A character vector of motif names, IUPAC-condensed structure strings,
or a 'glyrepr_structure' object.}

\item{alignments}{A character vector specifying alignment types for each motif.
Can be a single value (applied to all motifs) or a vector of the same length as motifs.}
}
\value{
\itemize{
\item \code{have_motif()}: A logical vector indicating if each \code{glycan} has the \code{motif}.
\item \code{have_motifs()}: A logical matrix where rows correspond to glycans and columns correspond to motifs.
Row names contain glycan identifiers and column names contain motif identifiers.
}
}
\description{
These functions check if the given \code{glycan}s have the given \code{motif}(s).
\itemize{
\item \code{have_motif()} checks a single motif against multiple glycans
\item \code{have_motifs()} checks multiple motifs against multiple glycans
}

Technically speaking, they perform subgraph isomorphism tests to
determine if the \code{motif}(s) are subgraphs of the \code{glycan}s.
Monosaccharides, linkages, and substituents are all considered.
}
\section{Graph mode and monosaccharide type}{
Both \code{glycan} and \code{motif} should be 'glyrepr_structure' objects
(see \code{\link[glyrepr:glycan_structure]{glyrepr::glycan_structure()}}).

Also, they can have different monosaccharide types
("concrete" or "generic", see \code{\link[glyrepr:get_mono_type]{glyrepr::get_mono_type()}}).
However, a "concrete" motif cannot be matched to a "generic" glycan.
}

\section{Linkages}{
Obscure linkages (e.g. "??-?") are allowed in the \code{motif} graph
(see \code{\link[glyrepr:possible_linkages]{glyrepr::possible_linkages()}}).
"?" in a motif graph means "anything could be OK",
so it will match any linkage in the \code{glycan} graph.
However, "?" in a \code{glycan} graph will only match "?" in the \code{motif} graph.
You can set \code{ignore_linkages = TRUE} to ignore linkages in the comparison.

Some examples:
\itemize{
\item "b1-?" in motif will match "b1-4" in glycan.
\item "b1-?" in motif will match "b1-?" in glycan.
\item "b1-4" in motif will NOT match "b1-?" in glycan.
\item "a1-?" in motif will NOT match "b1-4" in glycan.
\item "a1-?" in motif will NOT match "a?-4" in glycan.
}

Both motifs and glycans can have a "half-linkage" at the reducing end,
e.g. "GlcNAc(b1-".
The half linkage in the motif will be matched to any linkage in the glycan,
or the half linkage of the glycan.
e.g. Glycan "GlcNAc(b1-4)Gal(a1-" will have both "GlcNAc(b1-" and "Gal(a1-" motifs.
}

\section{Alignment}{
According to the \href{https://glycomotif.glyomics.org}{GlycoMotif} database,
a motif can be classified into four alignment types:
\itemize{
\item "substructure": The motif can be anywhere in the glycan. This is the default.
See \href{https://glycomotif.glyomics.org/glycomotif/Substructure_Alignment}{substructure}
for details.
\item "core": The motif must align with at least one connected substructure
(subtree) at the reducing end of the glycan.
See \href{https://glycomotif.glyomics.org/glycomotif/Glycan_Core_Alignment}{glycan core}
for details.
\item "terminal": The motif must align with at least one connected substructure
(subtree) at the nonreducing end of the glycan.
See \href{https://glycomotif.glyomics.org/glycomotif/Nonreducing-End_Alignment}{nonreducing end}
for details.
\item "whole": The motif must align with the entire glycan.
See \href{https://glycomotif.glyomics.org/glycomotif/Whole-Glycan_Alignment}{whole-glycan}
for details.
}

When using known motifs in the GlycoMotif GlyGen Collection,
the best practice is to not provide the \code{alignment} argument,
and let the function decide the alignment based on the motif name.
However, it is still possible to override the default alignments.
In this case, the user-provided alignments will be used,
but a warning will be issued.
}

\section{Substituents}{
Substituents (e.g. "Ac", "SO3") are matched in strict mode.
"Neu5Ac-9Ac" will only match "Neu5Ac-9Ac" but not "Neu5Ac",
and "Neu5Ac" will not match "Neu5Ac-9Ac".
Obscure linkage in the motif will match any linkage in the glycan.
e.g. Motif "Neu5Ac-?Ac" will match "Neu5Ac-9Ac" in the glycan.
}

\section{Implementation}{
Under the hood, the function uses \code{\link[igraph:subgraph_isomorphisms]{igraph::graph.get.subisomorphisms.vf2()}}
to get all possible subgraph isomorphisms between \code{glycan} and \code{motif}.
\code{color} vertex attributes are added to the graphs to distinguish monosaccharides.
For all possible matches, the function checks the following:
\itemize{
\item Alignment: using \code{alignment_check()}
\item Substituents: using \code{substituent_check()}
\item Linkages: using \code{linkage_check()}
\item Anomer: using \code{anomer_check()}
The function returns \code{TRUE} if any of the matches pass all checks.
}
}

\examples{
library(glyparse)
library(glyrepr)

(glycan <- o_glycan_core_2(mono_type = "concrete"))

# The glycan has the motif "Gal(b1-3)GalNAc"
have_motif(glycan, "Gal(b1-3)GalNAc")

# But not "Gal(b1-4)GalNAc" (wrong linkage)
have_motif(glycan, "Gal(b1-4)GalNAc")

# Set `ignore_linkages` to `TRUE` to ignore linkages
have_motif(glycan, "Gal(b1-4)GalNAc", ignore_linkages = TRUE)

# Different monosaccharide types are allowed
have_motif(glycan, "Hex(b1-3)HexNAc")

# Obscure linkages in the `motif` graph are allowed
have_motif(glycan, "Gal(b1-?)GalNAc")

# However, obscure linkages in `glycan` will only match "?" in the `motif` graph
glycan_2 <- parse_iupac_condensed("Gal(b1-?)[GlcNAc(b1-6)]GalNAc")
have_motif(glycan_2, "Gal(b1-3)GalNAc")
have_motif(glycan_2, "Gal(b1-?)GalNAc")

# The anomer of the motif will be matched to linkages in the glycan
have_motif(glycan_2, "GlcNAc(b1-")

# Alignment types
# The default type is "substructure", which means the motif can be anywhere in the glycan.
# Other options include "core", "terminal" and "whole".
glycan_3 <- parse_iupac_condensed("Gal(a1-3)Gal(a1-4)Gal(a1-6)Gal")
motifs <- c(
  "Gal(a1-3)Gal(a1-4)Gal(a1-6)Gal",
  "Gal(a1-3)Gal(a1-4)Gal",
           "Gal(a1-4)Gal(a1-6)Gal",
           "Gal(a1-4)Gal"
)

purrr::map_lgl(motifs, ~ have_motif(glycan_3, .x, alignment = "whole"))
purrr::map_lgl(motifs, ~ have_motif(glycan_3, .x, alignment = "core"))
purrr::map_lgl(motifs, ~ have_motif(glycan_3, .x, alignment = "terminal"))
purrr::map_lgl(motifs, ~ have_motif(glycan_3, .x, alignment = "substructure"))

# Substituents
glycan_4 <- "Neu5Ac9Ac(a2-3)Gal(b1-4)GlcNAc"
glycan_5 <- "Neu5Ac(a2-3)Gal(b1-4)GlcNAc"

have_motif(glycan_4, glycan_5)
have_motif(glycan_5, glycan_4)
have_motif(glycan_4, glycan_4)
have_motif(glycan_5, glycan_5)

# Vectorization with single motif
glycans <- c(glycan, glycan_2, glycan_3)
motif <- "Gal(b1-3)GalNAc"
have_motif(glycans, motif)

# Multiple motifs with have_motifs()
glycan1 <- o_glycan_core_2(mono_type = "concrete")
glycan2 <- parse_iupac_condensed("Gal(b1-?)[GlcNAc(b1-6)]GalNAc")
glycans <- c(glycan1, glycan2)

motifs <- c("Gal(b1-3)GalNAc", "Gal(b1-4)GalNAc", "GlcNAc(b1-6)GalNAc")
result <- have_motifs(glycans, motifs)
print(result)

}
\seealso{
\code{\link[=count_motif]{count_motif()}}, \code{\link[=count_motifs]{count_motifs()}}
}
