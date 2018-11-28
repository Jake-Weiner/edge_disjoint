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
typedef std::vector<Edge> EdgeVec;
typedef map<Edge,int> Edge_Int_Map;
typedef vector<int> vi;
typedef vector<double,int> di;  
typedef vector<Edge> EdgeVec;
typedef vector<string> split_vector_type;



double djikstras(Edge_Int_Map& EIM,
    map<int, map<int, bool>>& node_neighbours, int start, int end,
    vector<int>& parents, int num_nodes, int num_edges, 
    DblVec& rc, IntVec& x, int comm_idx, int num_comm);