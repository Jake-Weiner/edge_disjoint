#include "LaPSO.hpp"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <string>
#include <map>

using namespace LaPSO;
using namespace std;
using namespace boost;

typedef vector<string> split_vector_type;
typedef adjacency_list<listS, vecS, directedS, no_property, property<edge_weight_t, int>> graph_t;
typedef graph_traits<graph_t>::vertex_descriptor vertex_descriptor;
typedef graph_traits<graph_t>::edge_descriptor edge_descriptor;
typedef pair<int, int> Edge;
typedef map<Edge,int> Edge_Int_Map;



struct Commodity {
    int origin;
    int dest;
    vector<Edge> solution_edges;
    float solution_value;
};

class ED : public LaPSO::UserHooks {

public:
    LaPSO::DblVec weight;
    ED(string graph_filename, string pairs_filename);
    void solve_ED();
    int nMSTsolves;		// number of times MST was solved
    int maxMSTsolves;		// abort after this many
    vector<Commodity> get_commodities();
    ~ED(){}
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
    void populate_graph(string filename);
    void populate_commodities(string filename);
    void solve_shortest_path(Commodity& commodity);
};
