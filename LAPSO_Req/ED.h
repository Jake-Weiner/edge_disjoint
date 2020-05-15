#ifndef ED_H
#define ED_H

#include "LaPSO.hpp"
#include "types.h"
#include <ilcplex/ilocplex.h>
#include <string>
#include <map>
#include <unordered_map>
#include <boost/functional/hash.hpp>


using namespace LaPSO;
using namespace std;

typedef vector<string> split_vector_type;
typedef pair<int, int> Edge;
typedef map<Edge,int> Edge_Int_Map;
typedef EdgeVec::const_iterator EdgeIter;
typedef std::pair<int,int> NodeEdgePair; // neighbouring node and edge index


struct Commodity {
    int origin;
    int dest;
    EdgeVec solution_edges;
    vector<Edge> solution_edges_nodes;
    float solution_value;
    int comm_idx;
};

struct Commodity_SP {
    int start;
    int end;
    vector<NodeEdgePair> parents;
};


template <typename Container> // we can make this generic for any container [1]
struct container_hash {
    std::size_t operator()(Container const& c) const {
        return boost::hash_range(c.begin(), c.end());
    }
};

struct pair_hash {
    template <class T1, class T2>
    std::size_t operator () (const std::pair<T1,T2> &p) const {
        auto h1 = std::hash<T1>{}(p.first);
        auto h2 = std::hash<T2>{}(p.second);

        // Mainly for demonstration purposes, i.e. works but is overly simple
        // In the real world, use sth. like boost.hash_combine
        return h1 ^ h2;  
    }
};

class EDParticle : public LaPSO::Particle {
public:
    EDParticle(const vector<Edge> &edges, const vector<Commodity> &comm, const int n, Edge_Int_Map em) : 
    LaPSO::Particle((int)edges.size() * (int)comm.size(),(int)edges.size())
    ,graph_edges(edges),commodities(comm),num_nodes(n)
    {
        c.resize(commodities.size(),0.0);
        commodity_shortest_paths.resize(commodities.size());

        cut_set_size = commodities.size();
    }
    vector<Edge> graph_edges;		/// list of edges in graph
    vector<EdgeVec> solution_edges;		/// edges involved in ED solution
    vector<Commodity> commodities;
    int num_nodes;
    int cut_set_size;
    vector<Commodity_SP> commodity_shortest_paths;
    vector<double> c; // shortest path costs in MIP
    //vector<vector<int>> cutsets; // include commodities involved in cutset e.g {{0,1,2}, {1,2,4}, {3,6,9} etc...}
    vector<int> cut_set_sizes;
    int repair_iter = 0;
    vector<int> solve_mip(map<vector<int>,int>& global_constraint_map);
    
};

class ED : public LaPSO::UserHooks {

public:
    LaPSO::DblVec weights;
    Edge_Int_Map EIM;
    unordered_map<int, Edge> IEM;  // int edge map - map edge index to edge
    vector<Edge> graph_edges;
    int solution_cost; // objective value of solution
    vector<EdgeVec> solution_edges; //edges used in each commodity SP
    vector<Edge> solution_edges_nodes; //edges used for MIP warmstart
    ED(string graph_filename, string pairs_filename, bool _printing, bool _randComm, 
    bool _djikstras_naive, MyTypes::repairRemoveEdgeMethod RREM, MyTypes::repairAddEdgeMethod RAEM);
    //void solve_ED(EDParticle &p);
    int nEDsolves;		// number of times ED was solved
    int maxEDsolves;		// abort after this many
    ~ED();
    Status reducedCost(const Particle &p, DblVec &redCost);
    void update_comm_sol(EDParticle& p, double SP,const vector<NodeEdgePair> &parents, double& total_paths_cost, int random_index,
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
  const size_t primalIdx(const int edge_idx,const int comm_idx) const{
    return comm_idx*num_edges + edge_idx;
  }
  const size_t edgeIdx(const int primal_idx) const{
    return primal_idx % num_edges;
  }

    const size_t edgeIdx(Edge e) const{
    std::cerr << "WARNING: edgeIdx may not be unique!!!!!!!!!!!!!!!!!!!!!\n";
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
  void storePath(EDParticle &p,int comm,int start,int end,const vector<NodeEdgePair> &parents,
		 vector<int> *viol=0,bool doPrint=false) const;
  vector<vector<NodeEdgePair>> get_node_neighbours(){return node_neighbours;}


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
    MyTypes::repairRemoveEdgeMethod RREM;
    MyTypes::repairAddEdgeMethod RAEM;
    void removeCommodity(EDParticle& p, IntVec& viol, int comm_idx_to_remove);
    // For each node_idx, node_neighbours contains NodeEdgePair --> {neighbour node_idx , edge_idx}
    vector<vector<NodeEdgePair>> node_neighbours;
    void remove_commodity(EDParticle& p, IntVec& viol, int commodity_index);
    void add_commodity(EDParticle& p, IntVec& viol, vector<int>& parents, int start, int end, int commodity_index);
        vector<int> find_cutset_commodities(const EDParticle& p, const vector<bool>& S_cutSet) const;
    int find_cutset_edges(const vector<bool>& S_cutSet) const;
    void initialise_NumVarArray();
    vector<int> solve_mip(EDParticle &p);

};


#endif 