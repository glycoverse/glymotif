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

// Experimental object-to-result path. The VF2 code above is copied unchanged
// from glymotif/src/vf2.cpp (MIT). No R functions are called by the path below.
// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp14)]]
#include <unordered_map>
#include <string>
#include <functional>
using namespace Rcpp;
using std::string;
using std::vector;

namespace fused {
vector<string> split(const string& s, char sep) {
  vector<string> out; size_t start=0, pos;
  while ((pos=s.find(sep,start))!=string::npos) {
    if(pos>start) out.push_back(s.substr(start,pos-start));
    start=pos+1;
  }
  if(start<s.size()) out.push_back(s.substr(start));
  return out;
}
bool token(const string& g,const string& m,bool lenient) {
  return m=="?" || (lenient && g=="?") || g==m;
}
bool anomer(const string& g,const string& m,bool lenient) {
  return token(g.substr(0,1),m.substr(0,1),lenient) &&
    token(g.substr(1),m.substr(1),lenient);
}
bool linkage(const string& g,const string& m,bool lenient) {
  if (!anomer(g.substr(0,2),m.substr(0,2),lenient)) return false;
  auto gs=split(g.substr(3),'/'), ms=split(m.substr(3),'/');
  if(std::find(ms.begin(),ms.end(),"?")!=ms.end() ||
     (lenient && std::find(gs.begin(),gs.end(),"?")!=gs.end())) return true;
  bool any=false;
  for(auto& x:gs) {
    bool found=std::find(ms.begin(),ms.end(),x)!=ms.end();
    if(!lenient && !found) return false;
    any|=found;
  }
  return lenient?any:true;
}
bool subs(const string& g,const string& m,bool strict,bool lenient) {
  if(!strict && m.empty()) return true;
  if(g.empty() || m.empty()) return g.empty() && m.empty();
  auto gs=split(g,','),ms=split(m,',');
  if(gs.size()!=ms.size()) return false;
  vector<bool> used(ms.size(),false);
  std::function<bool(size_t)> assign=[&](size_t i) {
    if(i==gs.size()) return true;
    for(size_t j=0;j<ms.size();++j) {
      if(!used[j] && token(gs[i].substr(0,1),ms[j].substr(0,1),lenient) &&
         gs[i].substr(1)==ms[j].substr(1)) {
        used[j]=true; if(assign(i+1)) return true; used[j]=false;
      }
    }
    return false;
  };
  return assign(0);
}
struct Dictionary {
  std::unordered_map<string,string> generic;
  std::set<string> generic_names,known;
  Dictionary(DataFrame table) {
    CharacterVector c=table["concrete"], g=table["generic"];
    for(int i=0;i<c.size();++i) {
      string cs=as<string>(c[i]),gs=as<string>(g[i]);
      generic.emplace(cs,gs); generic_names.insert(gs);
      known.insert(cs);known.insert(gs);
    }
  }
  bool mono(const string& g,const string& m,bool lenient) const {
    if(g==m) return true;
    bool gg=generic_names.count(g),mg=generic_names.count(m);
    auto gi=generic.find(g),mi=generic.find(m);
    if(!gg && mg && gi!=generic.end()) return gi->second==m;
    return lenient && gg && !mg && mi!=generic.end() && mi->second==g;
  }
  bool residue(const string& g,const string& gs,const string& m,
               const string& ms,bool strict,bool lenient) const {
    if(mono(g,m,lenient) && subs(gs,ms,strict,lenient)) return true;
    if(!(ms.size() && ms[0]=='?') && ms.find(",?")==string::npos) return false;
    string base,sub;
    if(g=="Neu5Ac") {base="Neu";sub="5Ac";}
    else for(auto suffix:{string("NAc"),string("N")}) {
      if(g.size()>suffix.size() && g.compare(g.size()-suffix.size(),suffix.size(),suffix)==0) {
        string b=g.substr(0,g.size()-suffix.size());
        if(known.count(b)) {base=b;sub="?"+suffix;break;}
      }
    }
    return !base.empty() && mono(base,m,lenient) &&
      subs(sub+(gs.empty()?"":","+gs),ms,strict,lenient);
  }
};
struct Profile {
  int n; vector<string> mono,sub,links,anomers;
  vector<int> in,out; IntegerMatrix edges;
  bool informative=false;
  Profile(SEXP obj, bool forest=false):edges(0,2) {
    if(!Rf_inherits(obj,"igraph") || TYPEOF(obj)!=VECSXP || Rf_length(obj)!=10)
      stop("Unsupported igraph layout; prototype requires igraph 2.3.3 layout.");
    List graph(obj); n=as<int>(graph[0]);
    if(!as<bool>(graph[1])) stop("Expected directed graph.");
    List attrs=graph[8],ga=attrs[1],va=attrs[2],ea=attrs[3];
    for(auto key:{"floating_parts","floating_substituents"})
      if(!forest && ga.containsElementNamed(key) && Rf_length(ga[key])>0)
        stop("Floating localization is not implemented in this prototype.");
    mono=as<vector<string>>(va["mono"]);sub=as<vector<string>>(va["sub"]);
    if(n<1 || mono.size()!=size_t(n) || sub.size()!=size_t(n)) stop("Invalid vertices.");
    string root=as<string>(ga["anomer"]);
    in.assign(n,0);out.assign(n,0);anomers.assign(n,root);
    NumericVector from=as<NumericVector>(graph[2]),to=as<NumericVector>(graph[3]);
    if(from.size()!=to.size() || (!forest && from.size()!=n-1)) stop("Expected a rooted tree.");
    links=ea.containsElementNamed("linkage")?as<vector<string>>(ea["linkage"]):vector<string>();
    if(links.size()!=size_t(from.size())) stop("Invalid linkages.");
    edges=IntegerMatrix(from.size(),2); informative=root!="??";
    for(int e=0;e<from.size();++e) {
      int f=from[e],t=to[e];
      if(f<0||t<0||f>=n||t>=n) stop("Invalid endpoint.");
      edges(e,0)=f+1;edges(e,1)=t+1;out[f]++;in[t]++;
      anomers[t]=links[e].substr(0,links[e].find('-'));
      informative|=links[e]!="??-?";
    }
    if((!forest && std::count(in.begin(),in.end(),0)!=1) || *std::max_element(in.begin(),in.end())>1)
      stop("Expected single rooted tree.");
  }
};
struct Structures {
  vector<Profile> graphs;vector<SEXP> sources;vector<int> restore;CharacterVector codes;
  Structures(SEXP obj, bool forest=false):codes(Rf_xlength(obj)) {
    if(!Rf_inherits(obj,"glyrepr_structure") || TYPEOF(obj)!=VECSXP)
      stop("Expected glycan_structure object.");
    List x(obj);List pool=x.attr("graphs");CharacterVector keys=pool.names();
    std::unordered_map<string,int> pool_index,used;
    for(int i=0;i<keys.size();++i) pool_index[as<string>(keys[i])]=i;
    for(int i=0;i<x.size();++i) {
      CharacterVector value=x[i];
      if(value.size()!=1) stop("Invalid structure key.");
      codes[i]=value[0];
      if(CharacterVector::is_na(value[0])) {restore.push_back(-1);continue;}
      string key=as<string>(value[0]);auto it=used.find(key);
      if(it==used.end()) {
        auto p=pool_index.find(key);if(p==pool_index.end()) stop("Missing cached graph.");
        int id=graphs.size();graphs.emplace_back(pool[p->second],forest);sources.push_back(pool[p->second]);used[key]=id;restore.push_back(id);
      } else restore.push_back(it->second);
    }
  }
};
List match(const Profile& g,const Profile& m,const Dictionary& dict,
           const string& alignment,bool ignore,bool strict,bool lenient,
           SEXP degree,bool first) {
  if(g.n<m.n || (alignment=="whole" && g.n!=m.n)) return List();
  bool check=!ignore && m.informative;
  if(check && !lenient && !g.informative) return List();
  bool has_degree=!Rf_isNull(degree);LogicalVector deg;
  if(has_degree) {deg=LogicalVector(degree);if(deg.size()!=m.n) stop("Degree mask length.");}
  LogicalMatrix vc(m.n,g.n),ec(m.links.size(),g.links.size());
  for(int mi=0;mi<m.n;++mi) {
    bool any=false;
    for(int gi=0;gi<g.n;++gi) {
      bool ok=dict.residue(g.mono[gi],g.sub[gi],m.mono[mi],m.sub[mi],strict,lenient);
      if(has_degree && deg[mi]) ok=ok && m.in[mi]==g.in[gi] && m.out[mi]==g.out[gi];
      if(!has_degree && alignment=="core" && m.in[mi]==0) ok=ok && g.in[gi]==0;
      if(!has_degree && alignment=="terminal" && m.out[mi]==0) ok=ok && g.out[gi]==0;
      if(check && m.in[mi]==0) ok=ok && anomer(g.anomers[gi],m.anomers[mi],lenient);
      // Necessary degree lower bounds are valid for any subgraph monomorphism.
      ok=ok && g.in[gi]>=m.in[mi] && g.out[gi]>=m.out[mi];
      vc(mi,gi)=ok;any|=ok;
    }
    if(!any) return List();
  }
  for(int mi=0;mi<ec.nrow();++mi)for(int gi=0;gi<ec.ncol();++gi)
    ec(mi,gi)=!check || linkage(g.links[gi],m.links[mi],lenient);
  return cpp_vf2_subgraph_mono(g.n,g.edges,m.n,m.edges,vc,ec,first,!first);
}
} // namespace fused

// [[Rcpp::export]]
SEXP cpp_structure_match(SEXP glycans, SEXP motifs, DataFrame dictionary,
                         CharacterVector alignments, bool ignore_linkages,
                         bool strict_sub, bool lenient, List degrees,
                         string output="have") {
  fused::Structures gs(glycans),ms(motifs);fused::Dictionary dict(dictionary);
  int ng=gs.restore.size(),nm=ms.restore.size();
  if(nm==0 || alignments.size()!=nm || degrees.size()!=nm) stop("Invalid motif options.");
  if(output!="have" && output!="count" && output!="match") stop("Invalid output.");
  LogicalMatrix have(ng,nm);IntegerMatrix count(ng,nm);List matches(nm);
  for(int j=0;j<nm;++j) {
    checkUserInterrupt();string alignment=as<string>(alignments[j]);
    if(alignment!="substructure" && alignment!="core" && alignment!="terminal" && alignment!="whole") stop("Invalid alignment.");
    if(ms.restore[j]<0) stop("Missing motif unsupported.");
    vector<List> values;
    for(auto& g:gs.graphs) {
      checkUserInterrupt();
      values.push_back(fused::match(g,ms.graphs[ms.restore[j]],dict,alignment,
        ignore_linkages,strict_sub,lenient,degrees[j],output=="have"));
    }
    List column(ng);
    for(int i=0;i<ng;++i) {
      int id=gs.restore[i];
      if(id<0) {have(i,j)=NA_LOGICAL;count(i,j)=NA_INTEGER;column[i]=R_NilValue;}
      else {have(i,j)=values[id].size()>0;count(i,j)=values[id].size();column[i]=values[id];}
    }
    matches[j]=column;
  }
  SEXP gn=Rf_getAttrib(glycans,R_NamesSymbol),mn=Rf_getAttrib(motifs,R_NamesSymbol);

  if(output=="match") {
    // Public match_motifs() has motif-major lists and preserves glycan names.
    for(int j=0;j<nm;++j) {List col=matches[j];col.attr("names")=gn;}
    if(!Rf_isNull(mn)) matches.attr("names")=mn;
    return matches;
  }
  if(Rf_isNull(gn)) gn=gs.codes;
  List dims=List::create(gn,mn);
  if(output=="have") {have.attr("dimnames")=dims;return have;}
  count.attr("dimnames")=dims;return count;
}

namespace floating {
struct Attachment {
  bool part;
  int root=-1;
  vector<int> nodes, parents;
  string token;
};
string position(const string& token) {
  size_t i=0;
  while(i<token.size() && (std::isdigit(token[i]) || token[i]=='/' || token[i]=='?')) ++i;
  return token.substr(0,i);
}
vector<string> slots(int parent,const string& positions) {
  vector<string> result;
  if(positions=="?") return result;
  for(auto& p:fused::split(positions,'/')) result.push_back(std::to_string(parent)+"\r"+p);
  return result;
}
bool slot_assignment(vector<vector<string>> domains) {
  // All domains here are known; unknown positions never occupy a known slot.
  std::stable_sort(domains.begin(),domains.end(),[](const vector<string>& a,const vector<string>& b){return a.size()<b.size();});
  std::unordered_map<string,int> owner;
  std::function<bool(int,std::set<string>&)> augment=[&](int i,std::set<string>& seen) {
    for(auto& s:domains[i]) {
      if(!seen.insert(s).second) continue;
      auto it=owner.find(s);
      if(it==owner.end() || augment(it->second,seen)) {owner[s]=i;return true;}
    }
    return false;
  };
  for(int i=0;i<int(domains.size());++i) {std::set<string> seen;if(!augment(i,seen)) return false;}
  return true;
}
string append_sub(const string& existing,const string& extra) {
  auto tokens=fused::split(existing,',');tokens.push_back(extra);
  std::stable_sort(tokens.begin(),tokens.end(),[](const string& a,const string& b){
    int pa=a[0]=='?'?100:a[0]-'0',pb=b[0]=='?'?100:b[0]-'0';return pa<pb;
  });
  string result;for(auto& t:tokens) {if(!result.empty()) result+=",";result+=t;}return result;
}
struct Localizations {
  vector<fused::Profile> variants;
  vector<vector<int>> assignments;
  double raw=1;
  Localizations(SEXP obj,const fused::Profile& base,int maximum) {
    List graph(obj),attrs=graph[8],ga=attrs[1];
    vector<Attachment> attachments;
    int nparts=0;
    for(auto key:{"floating_parts","floating_substituents"}) {
      if(!ga.containsElementNamed(key)) continue;
      List meta=ga[key];
      for(SEXP item:meta) {
        List value(item);Attachment a;a.part=string(key)=="floating_parts";
        a.parents=as<vector<int>>(value["parents"]);
        for(auto& parent:a.parents) --parent;
        if(a.part) {
          ++nparts;a.root=as<int>(value["root"])-1;
          a.nodes=as<vector<int>>(value["nodes"]);for(auto& n:a.nodes)--n;
          a.token=as<string>(value["linkage"]);
        } else a.token=as<string>(value["substituent"]);
        if(a.parents.empty()) for(int i=0;i<base.n;++i)
          if(!a.part || std::find(a.nodes.begin(),a.nodes.end(),i)==a.nodes.end()) a.parents.push_back(i);
        for(int parent:a.parents) if(parent<0 || parent>=base.n) stop("Invalid floating parent.");
        raw*=a.parents.size();attachments.push_back(a);
      }
    }
    if(raw>maximum || !std::isfinite(raw)) stop("Floating localization count exceeds max_variants.");
    if(attachments.empty()) {variants.push_back(base);assignments.push_back({});return;}
    vector<int> membership(base.n,-1);
    for(int i=0;i<nparts;++i) for(int node:attachments[i].nodes) membership[node]=i;
    vector<vector<string>> fixed;
    for(int i=0;i<base.edges.nrow();++i) {
      auto s=slots(base.edges(i,0)-1,base.links[i].substr(base.links[i].find('-')+1));
      if(!s.empty()) fixed.push_back(s);
    }
    for(int i=0;i<base.n;++i) for(auto& t:fused::split(base.sub[i],',')) {
      auto s=slots(i,position(t));if(!s.empty()) fixed.push_back(s);
    }
    // expand.grid varies the first attachment fastest. Preserve this order for mappings.
    for(int combo=0;combo<int(raw);++combo) {
      checkUserInterrupt();int q=combo;vector<int> chosen;
      for(auto& a:attachments) {chosen.push_back(a.parents[q%a.parents.size()]);q/=a.parents.size();}
      bool valid=true;
      for(int i=0;i<nparts && valid;++i) {
        std::set<int> visited;int current=i;
        while(current>=0) {
          if(!visited.insert(current).second) {valid=false;break;}
          current=membership[chosen[current]];
        }
      }
      if(!valid) continue;
      auto domains=fixed;
      for(int i=0;i<int(attachments.size());++i) {
        auto& a=attachments[i];string pos=a.part?a.token.substr(a.token.find('-')+1):position(a.token);
        if(pos!="?") domains.push_back(slots(chosen[i],pos));
      }
      if(!slot_assignment(domains)) continue;
      fused::Profile localized=base;
      localized.edges=IntegerMatrix(base.edges.nrow()+nparts,2);
      for(int e=0;e<base.edges.nrow();++e)for(int j=0;j<2;++j) localized.edges(e,j)=base.edges(e,j);
      for(int i=0;i<int(attachments.size());++i) {
        auto& a=attachments[i];int parent=chosen[i];
        if(a.part) {
          int e=base.edges.nrow()+i;
          localized.edges(e,0)=parent+1;localized.edges(e,1)=a.root+1;
          localized.links.push_back(a.token);localized.out[parent]++;localized.in[a.root]++;
          localized.anomers[a.root]=a.token.substr(0,a.token.find('-'));
          localized.informative|=a.token!="??-?";
        } else localized.sub[parent]=append_sub(localized.sub[parent],a.token);
      }
      if(localized.edges.nrow()!=base.n-1 || std::count(localized.in.begin(),localized.in.end(),0)!=1)
        stop("Localized graph is not a rooted tree.");
      variants.push_back(localized);
      for(auto& x:chosen) ++x;assignments.push_back(chosen);
    }
    if(variants.empty()) stop("No conflict-free floating localization.");
  }
};
}

// Diagnostic export: compare localized graphs before testing final match outputs.
// [[Rcpp::export]]
List cpp_localizations(SEXP graph,int maximum=256) {
  fused::Profile base(graph,true);floating::Localizations loc(graph,base,maximum);
  List result(loc.variants.size());
  for(int i=0;i<result.size();++i) {
    auto& g=loc.variants[i];
    result[i]=List::create(_["edges"]=g.edges,_["linkage"]=g.links,_["sub"]=g.sub,_["parents"]=loc.assignments[i]);
  }
  return List::create(_["raw"]=loc.raw,_["variants"]=result);
}

// [[Rcpp::export]]
SEXP cpp_floating_match(SEXP glycans,SEXP motifs,DataFrame dictionary,
                       CharacterVector alignments,bool ignore_linkages,
                       bool strict_sub,bool lenient,List degrees,
                       string output="have",bool strict_floating=true,int maximum=256) {
  fused::Structures gs(glycans,true),ms(motifs);fused::Dictionary dict(dictionary);
  int ng=gs.restore.size(),nm=ms.restore.size();
  if(maximum<1 || nm==0 || alignments.size()!=nm || degrees.size()!=nm) stop("Invalid options.");
  if(output!="have" && output!="count" && output!="match") stop("Invalid output.");
  vector<floating::Localizations> localized;
  for(int i=0;i<int(gs.graphs.size());++i) localized.emplace_back(gs.sources[i],gs.graphs[i],maximum);
  LogicalMatrix have(ng,nm);IntegerMatrix count(ng,nm);List matches(nm);
  for(int j=0;j<nm;++j) {
    string alignment=as<string>(alignments[j]);
    if(alignment!="substructure" && alignment!="core" && alignment!="terminal" && alignment!="whole") stop("Invalid alignment.");
    if(ms.restore[j]<0) stop("Missing motif unsupported.");
    List column(ng);vector<int> counts;vector<List> mappings;
    for(auto& loc:localized) {
      int aggregate=-1;List merged;vector<SEXP> ordered;std::set<vector<int>> seen;
      // Keep R mapping objects protected while combining variants.
      List results(loc.variants.size());
      for(int v=0;v<int(loc.variants.size());++v) {
        List one=fused::match(loc.variants[v],ms.graphs[ms.restore[j]],dict,alignment,
          ignore_linkages,strict_sub,lenient,degrees[j],output=="have");
        results[v]=one;int n=one.size();
        if(aggregate<0) aggregate=n;
        else aggregate=strict_floating?std::min(aggregate,n):std::max(aggregate,n);
        if(output=="match") for(SEXP mapping:one) {
          auto key=as<vector<int>>(mapping);
          if(seen.insert(key).second) ordered.push_back(mapping);
        }
        if(output=="have" && ((strict_floating && aggregate==0) || (!strict_floating && aggregate>0))) break;
      }
      for(SEXP item:ordered) merged.push_back(item);
      counts.push_back(aggregate);mappings.push_back(merged);
    }
    for(int i=0;i<ng;++i) {
      int id=gs.restore[i];
      if(id<0) {have(i,j)=NA_LOGICAL;count(i,j)=NA_INTEGER;column[i]=R_NilValue;}
      else {have(i,j)=counts[id]>0;count(i,j)=counts[id];column[i]=mappings[id];}
    }
    matches[j]=column;
  }
  SEXP gn=Rf_getAttrib(glycans,R_NamesSymbol),mn=Rf_getAttrib(motifs,R_NamesSymbol);
  if(output=="match") {
    for(int j=0;j<nm;++j){List col=matches[j];col.attr("names")=gn;}
    if(!Rf_isNull(mn)) matches.attr("names")=mn;return matches;
  }
  if(Rf_isNull(gn)) gn=gs.codes;
  List dims=List::create(gn,mn);
  if(output=="have") {have.attr("dimnames")=dims;return have;}
  count.attr("dimnames")=dims;return count;
}
