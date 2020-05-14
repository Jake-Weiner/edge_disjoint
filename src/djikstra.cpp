#include "djikstra.h"
#include <iostream>

#define INF 0x3f3f3f3f

using namespace std;
using namespace LaPSO;

struct node {
    int idx;
    pair<double, int> properties;
    //double epsilon;
   
};

constexpr bool operator>(const pair<double,int> &A,const pair<double,int> &B){
  return (A.first > B.first) || (A.first==B.first && A.second > B.second);
}

struct node_naive {
    int idx;
    double dist;   
};

// STRUCT TEMPLATE used for comparison operator in priority_queue // what am i comparing? a node which will have
//distance and viol
struct mygreater { // functor for operator>
    constexpr bool operator()(const node& _Left, const node& _Right)
    { // apply operator> to operands
      // AE: should be consistent on use of epsilon (either use it also for equality or not for >
      return _Left.properties > _Right.properties;
      //return (_Left.properties.first > (_Right.properties.first /*+ _Right.epsilon*/) || ((_Left.properties.first == _Right.properties.first) && (_Left.properties.second > _Right.properties.second)));
    }
};
bool operator>(const node& Left,const node& Right){ return mygreater()(Left,Right); }

double djikstras(double max_perturb, const Edge_Int_Map& EIM,
    const vector<vector<NodeEdgePair>>&  node_neighbours, int start, int end,
    vector<NodeEdgePair>& parents, int num_nodes, int num_edges,
		 DblVec& rc, IntVec& x, int comm_idx, int num_comm,double thresh)
{
    //vector<pair<double, int>> dist; // for storing the distance of every other node from source.
    //vector<int> parents; // predecessor nodes
    priority_queue<node, vector<node>, mygreater> Q; // min heap
    vector<pair<double,int>> distances(num_nodes,{thresh,0});
    
    vector<bool> visited(num_nodes,false);
    if(parents.size() < num_nodes)
	parents.resize(num_nodes,NodeEdgePair(-1,-1));
    else
	std::fill(parents.begin(),parents.end(),NodeEdgePair(-1,-1));
    distances[start].first = 0;

    Q.push({start,{0.0,0}/*,max_perturb*num_edges*/});
    while (!Q.empty()) {
        node current_node = Q.top();

        if (current_node.idx == end || distances[current_node.idx].first >= thresh) {
            /*
            int last_node = end;
            while (parents[last_node] != -1) {
                // return the SP
                relaxed_solution(last_node + nodes * parents[last_node], 0) = 1;
                last_node = parents[last_node];
            }
            */
           return distances[end].first;
        }
        
        Q.pop();
	if(visited[current_node.idx]) {
        continue;
    }
    else{
         visited[current_node.idx] = true;
    }
        //cout << "node " << current_node.idx << " popped" << endl;
        //iterate through neighbours of current_node 
	for(const NodeEdgePair adj : node_neighbours[current_node.idx] ){
                
	    const int neighbour_node = adj.first;
	    if(visited[neighbour_node]) continue;
	    const int edge = adj.second;
	    const double edge_weight = rc[edge +  comm_idx*num_edges];
	    auto d = distances[current_node.idx];
	    d.first += edge_weight;
	    if(d.first >= thresh) continue; // too long to consider
	    int edge_used = 0;
	    for (int i=0; i< num_comm; i++){
		if(x[edge + i*num_edges] == 1)
		    edge_used += 1;
	    }
	    d.second += edge_used;
	    
	    // if distance of node is shorter, or if distances is same but violation is less
	    if ( distances[neighbour_node] > d ) {		
	      distances[neighbour_node] = d;
	      parents[neighbour_node] = NodeEdgePair(current_node.idx,edge);
	      Q.push({neighbour_node, {d.first,d.second}/*, max_perturb*num_edges*/});
	      //cout << "node " << neighbour_node << " pushed - dist " << distances[neighbour_node].first
	      //<< " viol - " << distances[neighbour_node].second << endl;
	    }
	}
    }
    
    return distances[end].first;
}


struct mygreater_naive
		{	// functor for operator>
			constexpr bool operator()(const node_naive& _Left, const node_naive& _Right) const
			{	// apply operator> to operands
				return (_Left.dist > _Right.dist);
			}
		};

double djikstras_naive(const Edge_Int_Map& EIM, // unused
		       const vector<vector<NodeEdgePair>> &node_neighbours, int start, int end,
		       vector<NodeEdgePair>& parents, int num_nodes, int num_edges,
		       DblVec& rc, IntVec& x, int comm_idx, int num_comm,double thresh){

    priority_queue<node_naive, vector<node_naive>, mygreater_naive> Q; // min heap
    vector<bool> visited(num_nodes,false);
    vector<double> distances(num_nodes,thresh);
    
    if(parents.size() < (size_t)num_nodes)
	parents.resize(num_nodes,NodeEdgePair(-1,-1));
    else
	std::fill(parents.begin(),parents.end(),NodeEdgePair(-1,-1));
    distances[start] = 0.0;
    
    Q.push({start,0.0});
    while (!Q.empty()) {
        node_naive current_node = Q.top();

        if (current_node.idx == end || distances[current_node.idx] >= thresh) {
            /*
            int last_node = end;
            while (parents[last_node] != -1) {
                // return the SP
                relaxed_solution(last_node + nodes * parents[last_node], 0) = 1;
                last_node = parents[last_node];
            }
            */
	    return distances[end]; 
        }
        
        Q.pop();
	if(visited[current_node.idx]) continue;
	visited[current_node.idx] = true;
        //cout << "node " << current_node.idx << " popped" << endl;
        //iterate through neighbours of current_node
	for(const NodeEdgePair adj : node_neighbours[current_node.idx] ){
                
	    const int neighbour_node = adj.first;
	    if(visited[neighbour_node]) continue;
	    const int edge = adj.second;
	    const double d = distances[current_node.idx]+rc[edge +  comm_idx*num_edges];
	    if(d >= thresh) continue;
	    // if distance of node is shorter, or if distances is same but violation is less
	    if ( distances[neighbour_node] > d){
		distances[neighbour_node] = d;
		parents[neighbour_node] = NodeEdgePair(current_node.idx,edge); 
		Q.push({neighbour_node, d});
		//cout << "node " << neighbour_node << " pushed - dist " << distances[neighbour_node].first
		//<< " viol - " << distances[neighbour_node].second << endl;
	    }
	}
    }
    
    return distances[end];
}

double djikstras_naive_cutSet(const Edge_Int_Map& EIM,
			      const vector<vector<NodeEdgePair>> &node_neighbours, int start, int end,
			      vector<NodeEdgePair>& parents, int num_nodes, int num_edges,
			      DblVec& rc, IntVec& x, int comm_idx, int num_comm,vector<bool>& S_cutSet, double thresh){

    priority_queue<node_naive, vector<node_naive>, mygreater_naive> Q; // min heap
    vector<double> distances(num_nodes,thresh);
    vector<bool> visited(num_nodes,false);
    
    if(parents.size() < num_nodes)
	parents.resize(num_nodes,NodeEdgePair(-1,-1));
    else
	std::fill(parents.begin(),parents.end(),NodeEdgePair(-1,-1));
    distances[start] = 0.0;

    Q.push({start,0.0});
    while (!Q.empty()) {
        node_naive current_node = Q.top();
	if(distances[current_node.idx] >= thresh) break;
        Q.pop();
	if(visited[current_node.idx]) continue;
        visited[current_node.idx] = true;
        //cout << "node " << current_node.idx << " popped" << endl;
        //iterate through neighbours of current_node
	for(const NodeEdgePair adj : node_neighbours[current_node.idx] ){
                
	    const int neighbour_node = adj.first;
	    if(visited[neighbour_node]) continue;
	    const int edge = adj.second;
	    const double d = distances[current_node.idx] + rc[edge +  comm_idx*num_edges];               
	    
	    // if distance of node is shorter, or if distances is same but violation is less
	    if ( distances[neighbour_node] > d){
		distances[neighbour_node] = d;
                      
		parents[neighbour_node] = NodeEdgePair(current_node.idx,edge);
		Q.push({neighbour_node, d});
		//cout << "node " << neighbour_node << " pushed - dist " << distances[neighbour_node].first
		//<< " viol - " << distances[neighbour_node].second << endl;
	    }
	}
    }
    // at this point, any nodes with distances greater than threshold mean they can't be reached
    if(S_cutSet.size() < num_nodes)
	S_cutSet.resize(num_nodes,false);
    for(int i=0;i<num_nodes;++i)
	S_cutSet[i] = (distances[i] < thresh);

    return distances[end];
}
