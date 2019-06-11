#include "LaPSO.hpp"
#include <stdio.h>
#include <vector>
#include <queue>
#include <stack>
#include <cstdlib>
#include <map>

using namespace std;
using namespace LaPSO;

typedef pair<int, int> Edge;
typedef map<Edge,int> Edge_Int_Map;
typedef vector<string> split_vector_type;
typedef pair<int, int> NodeEdgePair;


double djikstras(double max_perturb, const Edge_Int_Map& EIM,
		 const vector<vector<NodeEdgePair>>& node_neighbours, int start, int end,
		 vector<NodeEdgePair>& parents, int num_nodes, int num_edges, 
		 DblVec& rc, IntVec& x, int comm_idx, int num_comm,double thresh=1e99);

double djikstras_naive(const Edge_Int_Map& EIM,
		       const vector<vector<NodeEdgePair>>& node_neighbours, int start, int end,
		       vector<NodeEdgePair>& parents, int num_nodes, int num_edges, 
		       DblVec& rc, IntVec& x, int comm_idx, int num_comm,double thresh=1e99);

double djikstras_naive_cutSet(const Edge_Int_Map& EIM,
    const vector<vector<NodeEdgePair>>& node_neighbours, int start, int end,
    vector<NodeEdgePair>& parents, int num_nodes, int num_edges,
    DblVec& rc, IntVec& x, int comm_idx, int num_comm,vector<bool>& S_cutset, double thresh);
