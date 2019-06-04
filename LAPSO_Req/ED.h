#include "LaPSO.hpp"
#include <ilcplex/ilocplex.h>
#include <string>
#include <map>


using namespace LaPSO;
using namespace std;

typedef vector<string> split_vector_type;
typedef pair<int, int> Edge;
typedef map<Edge,int> Edge_Int_Map;
typedef vector<Edge> EdgeVec;
typedef EdgeVec::const_iterator EdgeIter;

struct Commodity {
    int origin;
    int dest;
    EdgeVec solution_edges;
    float solution_value;
    int comm_idx;
};

struct Commodity_SP {
    int start;
    int end;
    vector<int> parents;
};



class EDParticle : public LaPSO::Particle {
public:
    EDParticle(const EdgeVec &edges, const vector<Commodity> &comm, const int n, Edge_Int_Map em) : 
    LaPSO::Particle((int)edges.size() * (int)comm.size(),(int)edges.size())
    ,graph_edges(edges),commodities(comm),num_nodes(n)
    {
        c.resize(commodities.size(),0.0);
        commodity_shortest_paths.resize(commodities.size());

        cut_set_size = commodities.size();
    }
    EdgeVec graph_edges;		/// list of edges in graph
    vector<EdgeVec> solution_edges;		/// edges involved in ED solution
    vector<Commodity> commodities;
    int num_nodes;
    int cut_set_size;
    vector<Commodity_SP> commodity_shortest_paths;
    vector<float> c; // shortest path costs in MIP
    vector<vector<int>> cutsets; // include commodities involved in cutset e.g {{0,1,2}, {1,2,4}, {3,6,9} etc...}
    vector<int> cut_set_sizes;
    int repair_iter = 0;
    vector<int> solve_mip(map<vector<int>,int>& global_constraint_map);
    
};

class ED : public LaPSO::UserHooks {

public:
    LaPSO::DblVec weights;
    Edge_Int_Map EIM;
    EdgeVec graph_edges;
    int solution_cost; // objective value of solution
    vector<EdgeVec> solution_edges; //edges used in each commodity SP
    ED(string graph_filename, string pairs_filename, bool _printing, bool _randComm, 
    bool _djikstras_naive, string repair_remove_edge,string repair_add_edge);
    //void solve_ED(EDParticle &p);
    int nEDsolves;		// number of times ED was solved
    int maxEDsolves;		// abort after this many
    ~ED();
    Status reducedCost(const Particle &p, DblVec &redCost);
    void update_comm_sol(EDParticle& p, double SP, vector<int> parents, double& total_paths_cost, int random_index,
    int start, int end, bool printing);
    Status solveSubproblem(Particle& p);
    //void randomiseMethod(EDParticle &p);
    Status fixConstraint(const int constraint,
									 const Particle &p,
									 SparseVec &feas);
    Status heuristics(Particle &p);
    Status updateBest(Particle &p);
    void localSearch(Particle& p_);
    void initialise_global_constraints();
    void add_constraints_mip(vector<pair<vector<int>, int>>& local_constraints);

	const vector<Commodity>& getComm() const {return commodities;}
  // map from edge in graph to index of lagrange vector
  // or of the primal vector (num.arcs * num.commodities)
  /*const size_t dualIdx(const edge_iterator &e) const {
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
  */
  const size_t primalIdx(const int edge_idx,const int comm_idx){
    return comm_idx*num_edges + edge_idx;
  }
  const size_t edgeIdx(const int primal_idx){
    return primal_idx % num_edges;
  }

  const size_t edgeIdx(Edge e){
    auto ei = EIM.find(e);
    if(ei == EIM.end()) return 999999999; // error
    return ei->second; 
  }

  const int getCommSize() {return commodities.size();}
  void write_mip(vector<Particle*> &non_dom, double lb, double ub, string outfile_name);
  void setPrinting(bool p) {printing=p;}
  bool getPrinting() const {return printing;}
  bool getrandComm() {return randComm;}
  int get_nodes() {return num_nodes;} 
  int get_edges() {return num_edges;} 



private:
    bool printing;
    vector<Commodity> commodities;
    int num_nodes;
    int num_edges;
    void populate_graph(string filename);
    void populate_commodities(string filename);
    bool randComm;
    bool dN; // djikstas_naive
    //IloRangeArray global_constraints;
    vector<pair<vector<int>,int>> global_constraints;
    map<vector<int>,int> constraint_map;
    string repair_remove_edge;
    string repair_add_edge;
    map<int,map<int,bool>> node_neighbours;
    void remove_commodity(EDParticle& p, IntVec& viol, int commodity_index);
    void add_commodity(EDParticle& p, IntVec& viol, vector<int>& parents, int start, int end, int commodity_index);
    vector<int> find_cutset_commodities(EDParticle& p, map<int, bool>& S_cutset, map<int, bool>& T_cutset);
    int find_cutset_edges(map<int, bool>& S_cutset, map<int, bool>& T_cutset);
    void initialise_NumVarArray();
    vector<int> solve_mip(EDParticle &p);
};


