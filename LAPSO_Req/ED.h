#include "LaPSO.hpp"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <string>
#include <map>

using namespace LaPSO;
using namespace std;
using namespace boost;

typedef vector<string> split_vector_type;
// we are only dealing with undirected graphs
typedef adjacency_list<listS, vecS, undirectedS, no_property, property<edge_index_t, size_t> > graph_t;
typedef property<edge_weight_t, double> WeightProperty;
typedef graph_traits<graph_t>::vertex_descriptor vertex_descriptor;
typedef graph_traits<graph_t>::edge_descriptor edge_descriptor;
typedef graph_traits<graph_t>::edge_iterator edge_iterator;
typedef pair<int, int> Edge;
typedef map<Edge,int> Edge_Int_Map;
typedef vector<Edge> EdgeVec;
typedef EdgeVec::const_iterator EdgeIter;

struct Commodity {
    int origin;
    int dest;
    vector<Edge> solution_edges;
    float solution_value;
    int comm_idx;
};

class ED : public LaPSO::UserHooks {

public:
    LaPSO::DblVec weights;
    Edge_Int_Map EIM;
    EdgeVec graph_edges;
    vector<EdgeVec> solution_edges; //edges used in each commodity SP
    ED(string graph_filename, string pairs_filename);
    //void solve_ED(EDParticle &p);
    int nEDsolves;		// number of times ED was solved
    int maxEDsolves;		// abort after this many
    ~ED();
    Status reducedCost(const Particle &p, DblVec &redCost);
    Status solveSubproblem(Particle& p);
    Status fixConstraint(const int constraint,
									 const Particle &p,
									 SparseVec &feas);
    Status heuristics(Particle &p);
    Status updateBest(Particle &p);
    const graph_t& getGraph() const {return g; }
	const vector<Commodity>& getComm() const {return commodities;}
  // map from edge in graph to index of lagrange vector
  // or of the primal vector (num.arcs * num.commodities)
  const size_t dualIdx(const edge_iterator &e) const {
    return get(edge_index,g,*e); }
  const size_t primalIdx(const edge_iterator &e,int c) const {
    return num_edges(g)*c+dualIdx(e); }
  const size_t primalIdx(int u,int v,int c) const {
    auto ei = EIM.find(Edge(u,v)); // Look up (u,v) = edge
    if(ei == EIM.end()) return 999999999; // error
    return num_edges(g)*c+ei->second; } // ei->second=edge index
  const size_t edgeIdx(size_t primalIdx) const {
    return primalIdx % num_edges(g); }
  const size_t edgeIdx(Edge a) const {
    auto ei = EIM.find(a);
    if(ei == EIM.end()) return 999999999; // error
    return ei->second; } // ei->second=edge index
  const size_t commIdx(size_t primalIdx) const {
    return primalIdx / num_edges(g); }
  const size_t primalIdx(size_t edge ,int c) const {
    return num_edges(g)*c+ edge; }

private:
    vector<Commodity> commodities;
    graph_t g;
    int num_nodes;
    void populate_graph(string filename);
    void populate_commodities(string filename);
};

class EDParticle : public LaPSO::Particle {
public:
    EDParticle(const EdgeVec &edges, const vector<Commodity> &comm, const int n, Edge_Int_Map em) : 
    LaPSO::Particle((int)edges.size() * (int)comm.size(),(int)edges.size())
    ,graph_edges(edges),commodities(comm),num_nodes(n){}
    EdgeVec graph_edges;		/// list of edges in graph
    vector<EdgeVec> solution_edges;		/// edges involved in ED solution
    vector<Commodity> commodities;
    int num_nodes;
    int max_lb = 0;
    int max_ub = 0;
};
