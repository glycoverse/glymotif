#include <Rcpp.h>
#include <algorithm>
#include <map>
#include <string>
#include <unordered_set>
#include <vector>

namespace {
using Nodes = std::vector<int>;
using Subsets = std::vector<Nodes>;

std::string field(const std::string& value) {
  return std::to_string(value.size()) + ":" + value;
}

std::string field(SEXP value) {
  return value == NA_STRING ? "M" : field(std::string(Rf_translateCharUTF8(value)));
}

class MotifTree {
 public:
  std::vector<Nodes> children;
  std::vector<std::vector<std::string>> links;
  std::vector<std::string> labels;
  Rcpp::CharacterVector anomers;
  std::map<std::pair<int, int>, Subsets> cache;

  explicit MotifTree(const Rcpp::List& context) : anomers(context["anomers"]) {
    Rcpp::List child_list = context["children"];
    Rcpp::List link_list = context["child_linkages"];
    Rcpp::CharacterVector mono = context["mono"], sub = context["sub"];
    for (int i = 0; i < child_list.size(); ++i) {
      Rcpp::IntegerVector ids = Rcpp::as<Rcpp::IntegerVector>(child_list[i]);
      Nodes order(ids.begin(), ids.end());
      std::sort(order.begin(), order.end());
      children.push_back(Nodes());
      links.push_back(std::vector<std::string>());
      Rcpp::CharacterVector linkage = link_list[i];
      for (int id : order) {
        auto position = std::find(ids.begin(), ids.end(), id) - ids.begin();
        children.back().push_back(id - 1);
        links.back().push_back(field(linkage[position]));
      }
      labels.push_back("N" + field(mono[i]) + field(sub[i]));
    }
  }

  const Subsets& subsets(int root, int limit) {
    auto key = std::make_pair(root, limit);
    auto found = cache.find(key);
    if (found != cache.end()) return found->second;
    Subsets result{Nodes{root}};
    if (limit > 1 && !children[root].empty()) {
      Subsets combinations{Nodes{}};
      for (int child : children[root]) {
        const Subsets& choices = subsets(child, limit - 1);
        Subsets next = combinations; // Empty choice comes first.
        // expand.grid() varies the first child's choice fastest.
        for (const Nodes& choice : choices) {
          for (const Nodes& prefix : combinations) {
            Rcpp::checkUserInterrupt();
            if (prefix.size() + choice.size() >= static_cast<size_t>(limit)) continue;
            Nodes joined = prefix;
            joined.insert(joined.end(), choice.begin(), choice.end());
            next.push_back(std::move(joined));
          }
        }
        combinations = std::move(next);
      }
      for (const Nodes& combination : combinations) {
        if (combination.empty()) continue;
        Nodes joined{root};
        joined.insert(joined.end(), combination.begin(), combination.end());
        result.push_back(std::move(joined));
      }
    }
    return cache.emplace(key, std::move(result)).first->second;
  }

  std::string encode(int root, const std::vector<bool>& selected) const {
    std::vector<std::string> tokens;
    for (size_t i = 0; i < children[root].size(); ++i) {
      int child = children[root][i];
      if (selected[child]) tokens.push_back("E" + links[root][i] + field(encode(child, selected)));
    }
    std::sort(tokens.begin(), tokens.end());
    std::string joined;
    for (const auto& token : tokens) joined += token;
    return labels[root] + field(joined);
  }

  Nodes descendants(int root) const {
    Nodes result{root};
    for (size_t i = 0; i < result.size(); ++i) {
      Rcpp::checkUserInterrupt();
      const Nodes& next = children[result[i]];
      result.insert(result.end(), next.begin(), next.end());
    }
    return result;
  }
};
}

// Return first-occurrence node selections; R constructs only unique graphs.
// [[Rcpp::export]]
Rcpp::List cpp_extract_motif_candidates(Rcpp::List contexts, Rcpp::List roots,
                                       double max_size, bool branches) {
  std::unordered_set<std::string> seen;
  std::vector<Rcpp::List> result;
  for (int graph = 0; graph < contexts.size(); ++graph) {
    Rcpp::checkUserInterrupt();
    MotifTree tree(Rcpp::as<Rcpp::List>(contexts[graph]));
    Rcpp::IntegerVector graph_roots = roots[graph];
    int limit = max_size >= tree.children.size() ? tree.children.size() :
      static_cast<int>(std::max(1.0, max_size));
    for (int root_id : graph_roots) {
      int root = root_id - 1;
      Subsets full;
      if (branches) full.push_back(tree.descendants(root));
      const Subsets& candidates = branches ? full : tree.subsets(root, limit);
      for (const Nodes& nodes : candidates) {
        Rcpp::checkUserInterrupt();
        std::vector<bool> selected(tree.children.size(), false);
        for (int node : nodes) selected[node] = true;
        std::string key = "R" + field(tree.anomers[root]) + field(tree.encode(root, selected));
        if (!seen.insert(key).second) continue;
        Rcpp::IntegerVector ids(nodes.size());
        for (size_t j = 0; j < nodes.size(); ++j) ids[j] = nodes[j] + 1;
        result.push_back(Rcpp::List::create(Rcpp::_["graph"] = graph + 1,
          Rcpp::_["nodes"] = ids, Rcpp::_["anomer"] = Rcpp::CharacterVector::create(tree.anomers[root])));
      }
    }
  }
  return Rcpp::wrap(result);
}
