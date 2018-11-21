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
typedef boost::vector_property_map<double, boost::adj_list_edge_property_map<boost::undirected_tag, unsigned long, unsigned long&, unsigned long, boost::property<boost::edge_index_t, unsigned long, boost::no_property>, boost::edge_index_t> const> wM; 

struct Commodity {
    int origin;
    int dest;
    EdgeVec solution_edges;
    float solution_value;
    int comm_idx;
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

};

class ED : public LaPSO::UserHooks {

public:
    LaPSO::DblVec weights;
    Edge_Int_Map EIM;
    EdgeVec graph_edges;
    int solution_cost; // objective value of solution
    vector<EdgeVec> solution_edges; //edges used in each commodity SP
    ED(string graph_filename, string pairs_filename, bool _printing, bool _randComm, bool _fudge_factor);
    //void solve_ED(EDParticle &p);
    int nEDsolves;		// number of times ED was solved
    int maxEDsolves;		// abort after this many
    ~ED();
    Status reducedCost(const Particle &p, DblVec &redCost);
    void update_weightMap(bool fudge_factor, std::map<Edge, double>& added_weight, wM& weightMap,
    double& minWeight, int comm_idx, EDParticle& p);
    void update_comm_sol(EDParticle& p, DblVec distances,vector<vertex_descriptor> parents, double& total_paths_cost, int random_index, int comm_idx, vertex_descriptor& start, vertex_descriptor& end,bool printing);
    Status solveSubproblem(Particle& p);
    //void randomiseMethod(EDParticle &p);
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
  const int getCommSize() {return commodities.size();}
  void write_mip(vector<Particle*> &non_dom, double lb, double ub, string outfile_name);
  void setPrinting(bool p) {printing=p;}
  bool getPrinting() const {return printing;}
  bool getrandComm() {return randComm;}
  bool getfudgeFactor() {return fudge_factor;}


private:
    bool printing;
    vector<Commodity> commodities;
    graph_t g;
    int num_nodes;
    void populate_graph(string filename);
    void populate_commodities(string filename);
    bool randComm;
    bool fudge_factor;
};


