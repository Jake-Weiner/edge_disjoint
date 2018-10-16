#include "LaPSO.hpp"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <string>
#include <map>

using namespace LaPSO;
using namespace std;
using namespace boost;

typedef vector<string> split_vector_type;
typedef adjacency_list<listS, vecS, directedS, no_property, property<edge_weight_t, double>> graph_t;
typedef graph_traits<graph_t>::vertex_descriptor vertex_descriptor;
typedef graph_traits<graph_t>::edge_descriptor edge_descriptor;
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
    void solve_ED(EDParticle &p);
    int nEDsolves;		// number of times MST was solved
    int maxEDsolves;		// abort after this many
    ~ED();
    Status reducedCost(const Particle &p, DblVec &redCost);
    Status solveSubproblem(Particle& p);
    Status fixConstraint(const int constraint,
									 const Particle &p,
									 SparseVec &feas);
    Status heuristics(Particle &p);
    Status updateBest(Particle &p);
    

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
};
