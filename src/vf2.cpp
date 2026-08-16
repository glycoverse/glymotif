#include <Rcpp.h>

#include <algorithm>
#include <cstddef>
#include <set>
#include <utility>
#include <vector>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/vf2_sub_graph_iso.hpp>

namespace {

using Graph = boost::adjacency_list<
  boost::vecS,
  boost::vecS,
  boost::bidirectionalS,
  boost::no_property,
  boost::property<boost::edge_index_t, int>
>;

struct InterruptState {
  std::size_t checks = 0;

  void tick() {
    ++checks;
    if ((checks & 65535U) == 0U) {
      Rcpp::checkUserInterrupt();
    }
  }
};

struct VertexCompatibility {
  const Rcpp::LogicalMatrix* compatibility;
  InterruptState* interrupt;

  bool operator()(
    Graph::vertex_descriptor motif_vertex,
    Graph::vertex_descriptor glycan_vertex
  ) const {
    interrupt->tick();
    return (*compatibility)(
      static_cast<int>(motif_vertex),
      static_cast<int>(glycan_vertex)
    ) == TRUE;
  }
};

struct EdgeCompatibility {
  const Graph* motif;
  const Graph* glycan;
  const Rcpp::LogicalMatrix* compatibility;

  bool operator()(
    Graph::edge_descriptor motif_edge,
    Graph::edge_descriptor glycan_edge
  ) const {
    const int motif_id = boost::get(
      boost::edge_index,
      *motif,
      motif_edge
    );
    const int glycan_id = boost::get(
      boost::edge_index,
      *glycan,
      glycan_edge
    );
    return (*compatibility)(motif_id, glycan_id) == TRUE;
  }
};

struct MappingCollector {
  const Graph* motif;
  std::vector<std::vector<int>>* mappings;
  bool first_only;

  template <typename MotifToGlycanMap, typename GlycanToMotifMap>
  bool operator()(
    MotifToGlycanMap motif_to_glycan,
    GlycanToMotifMap
  ) const {
    std::vector<int> mapping;
    mapping.reserve(boost::num_vertices(*motif));
    const auto vertices = boost::vertices(*motif);
    for (auto vertex = vertices.first; vertex != vertices.second; ++vertex) {
      mapping.push_back(
        static_cast<int>(boost::get(motif_to_glycan, *vertex)) + 1
      );
    }

    mappings->push_back(std::move(mapping));

    return !first_only;
  }
};

void sort_and_deduplicate_mappings(
  std::vector<std::vector<int>>* mappings
) {
  std::sort(mappings->begin(), mappings->end());
  std::set<std::vector<int>> mapped_vertex_sets;
  std::vector<std::vector<int>> unique_mappings;
  unique_mappings.reserve(mappings->size());

  for (const auto& mapping : *mappings) {
    std::vector<int> key = mapping;
    std::sort(key.begin(), key.end());
    if (mapped_vertex_sets.insert(std::move(key)).second) {
      unique_mappings.push_back(mapping);
    }
  }
  mappings->swap(unique_mappings);
}

Graph make_graph(
  int vertex_count,
  const Rcpp::IntegerMatrix& edges
) {
  Graph graph(vertex_count);
  auto edge_index = boost::get(boost::edge_index, graph);
  for (int edge_id = 0; edge_id < edges.nrow(); ++edge_id) {
    const int from = edges(edge_id, 0) - 1;
    const int to = edges(edge_id, 1) - 1;
    if (
      from < 0 || from >= vertex_count ||
      to < 0 || to >= vertex_count
    ) {
      Rcpp::stop("Edge endpoint is outside the graph vertex range.");
    }

    const auto added = boost::add_edge(from, to, graph);
    if (!added.second) {
      Rcpp::stop("Could not add graph edge.");
    }
    boost::put(edge_index, added.first, edge_id);
  }
  return graph;
}

} // namespace


// [[Rcpp::export]]
Rcpp::List cpp_vf2_subgraph_mono(
  int glycan_vertex_count,
  Rcpp::IntegerMatrix glycan_edges,
  int motif_vertex_count,
  Rcpp::IntegerMatrix motif_edges,
  Rcpp::LogicalMatrix vertex_compatibility,
  Rcpp::LogicalMatrix edge_compatibility,
  bool first_only,
  bool unique_vertex_sets
) {
  if (glycan_vertex_count < 0 || motif_vertex_count < 0) {
    Rcpp::stop("Vertex counts must be non-negative.");
  }
  if (glycan_edges.ncol() != 2 || motif_edges.ncol() != 2) {
    Rcpp::stop("Edge matrices must have two columns.");
  }
  if (
    vertex_compatibility.nrow() != motif_vertex_count ||
    vertex_compatibility.ncol() != glycan_vertex_count
  ) {
    Rcpp::stop("Vertex compatibility matrix has invalid dimensions.");
  }
  if (
    edge_compatibility.nrow() != motif_edges.nrow() ||
    edge_compatibility.ncol() != glycan_edges.nrow()
  ) {
    Rcpp::stop("Edge compatibility matrix has invalid dimensions.");
  }

  Graph glycan = make_graph(glycan_vertex_count, glycan_edges);
  Graph motif = make_graph(motif_vertex_count, motif_edges);
  std::vector<std::vector<int>> mappings;
  InterruptState interrupt;

  MappingCollector callback{
    &motif,
    &mappings,
    first_only
  };
  VertexCompatibility vertex_predicate{
    &vertex_compatibility,
    &interrupt
  };
  EdgeCompatibility edge_predicate{
    &motif,
    &glycan,
    &edge_compatibility
  };

  boost::vf2_subgraph_mono(
    motif,
    glycan,
    callback,
    boost::vertex_order_by_mult(motif),
    boost::edges_equivalent(edge_predicate)
      .vertices_equivalent(vertex_predicate)
  );

  if (unique_vertex_sets) {
    sort_and_deduplicate_mappings(&mappings);
  }

  Rcpp::List output(mappings.size());
  for (std::size_t i = 0; i < mappings.size(); ++i) {
    output[i] = mappings[i];
  }
  return output;
}
