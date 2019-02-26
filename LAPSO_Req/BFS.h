#include <stdio.h>
#include <vector>
#include <queue>
#include <stack>
#include <cstdlib>
#include <map>

using namespace std;

typedef pair<int, int> Edge;
typedef std::vector<Edge> EdgeVec;
typedef map<Edge,int> Edge_Int_Map;
typedef vector<int> vi;
typedef vector<double,int> di;  
typedef vector<Edge> EdgeVec;
typedef vector<string> split_vector_type;

void BFS(
    map<int, map<int, bool>>& node_neighbours, int start, map<int,bool>& set_S);

