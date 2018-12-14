#include "djikstra.h"
#include <iostream>

#define INF 0x3f3f3f3f

using namespace std;
using namespace LaPSO;

struct node {
    int idx;
    pair<double, int> properties;
   
};

struct node_naive {
    int idx;
    double dist;
   
};

// STRUCT TEMPLATE used for comparison operator in priority_queue // what am i comparing? a node which will have
//distance and viol
struct mygreater { // functor for operator>
    constexpr bool operator()(const node& _Left, const node& _Right)
    { // apply operator> to operands
        return (_Left.properties.first > _Right.properties.first || ((_Left.properties.first == _Right.properties.first) && (_Left.properties.second > _Right.properties.second)));
    }
};

double djikstras(Edge_Int_Map& EIM,
    map<int, map<int, bool>>& node_neighbours, int start, int end,
    vector<int>& parents, int num_nodes, int num_edges,
      DblVec& rc, IntVec& x, int comm_idx, int num_comm)
{
    //vector<pair<double, int>> dist; // for storing the distance of every other node from source.
    //vector<int> parents; // predecessor nodes
    priority_queue<node, vector<node>, mygreater> Q; // min heap
    vector<pair<double,int>> distances;
    int current_edge_idx;
    double edge_weight;
    int highest_node_idx = node_neighbours.rbegin()->first;
    for (int i=0; i<highest_node_idx+1; i++){
        distances.push_back({INF,0});
    }
    map <int,bool> visited;
    parents.assign(highest_node_idx+1, -1);
    distances[start].first = 0;

    Q.push({start,{0.0,0}});
    while (!Q.empty()) {
        node current_node = Q.top();

        if (current_node.idx == end) {
            /*
            int last_node = end;
            while (parents[last_node] != -1) {
                // return the SP
                relaxed_solution(last_node + nodes * parents[last_node], 0) = 1;
                last_node = parents[last_node];
            }
            */
           return current_node.properties.first;
        }
        
        Q.pop();
        visited[current_node.idx] = 1;
        //cout << "node " << current_node.idx << " popped" << endl;
        //iterate through neighbours of current_node
        if (node_neighbours.find(current_node.idx) != node_neighbours.end()) {
            //it refers to each <int,bool> pair
            for (map<int, bool>::iterator it = node_neighbours[current_node.idx].begin();
                 it != node_neighbours[current_node.idx].end(); it++) {
                
                int neighbour_node = it->first; 
                // has node been visited
                if (visited.find(neighbour_node) != visited.end()){
                    continue;
                }
                int edge_used = 0;
                /*
                cout << "current node " << current_node.idx << endl;
                cout << "neighbour node " << neighbour_node << endl;
                */
                Edge current_edge = Edge(current_node.idx, neighbour_node);
                if (EIM.find(current_edge) != EIM.end()) {
                    current_edge_idx = EIM[current_edge];
                }
                else{
                    cerr << "edge not found in EIM --> djikstras" << endl;
                    return -1;
                }
                //cout << "current edge idx " << current_edge_idx << endl;
                edge_weight = rc[current_edge_idx +  comm_idx*num_edges];
                //cout << "edge weight is " << edge_weight << endl;
                for (int i=0; i< num_comm; i++){
                    if(x[current_edge_idx + i*num_edges] == 1){
                        edge_used += 1;
                    }
                }

                // if distance of node is shorter, or if distances is same but violation is less
                    if ( (distances[neighbour_node].first > (distances[current_node.idx].first + edge_weight))
                    || ( (distances[neighbour_node].first == (distances[current_node.idx].first + edge_weight)) &&
                       (distances[neighbour_node].second > (distances[current_node.idx].second + edge_used) ) ) ) {

                           distances[neighbour_node].first = distances[current_node.idx].first + edge_weight;
                           distances[neighbour_node].second = distances[current_node.idx].second + edge_used;
                           parents[neighbour_node] = current_node.idx;
                           Q.push({neighbour_node, {distances[neighbour_node].first,distances[neighbour_node].second}});
                           //cout << "node " << neighbour_node << " pushed - dist " << distances[neighbour_node].first
                           //<< " viol - " << distances[neighbour_node].second << endl;
                    }
            }
        } else {
            cerr << "node not found in neighbours list --> djikstras" << endl;
            return -1;
        }
	}
    
    return -1;
}


struct mygreater_naive
		{	// functor for operator>
			constexpr bool operator()(const node_naive& _Left, const node_naive& _Right) const
			{	// apply operator> to operands
				return (_Left.dist > _Right.dist);
			}
		};

double djikstras_naive(Edge_Int_Map& EIM,
    map<int, map<int, bool>>& node_neighbours, int start, int end,
    vector<int>& parents, int num_nodes, int num_edges,
      DblVec& rc, IntVec& x, int comm_idx, int num_comm){

    priority_queue<node_naive, vector<node_naive>, mygreater_naive> Q; // min heap
    vector<double> distances;
    int current_edge_idx;
    double edge_weight;
    int highest_node_idx = node_neighbours.rbegin()->first;
    
    map <int,bool> visited;
    parents.assign(highest_node_idx+1, -1);
    distances.assign(highest_node_idx+1,INF);
    distances[start] = 0.0;

    Q.push({start,0.0});
    while (!Q.empty()) {
        node_naive current_node = Q.top();

        if (current_node.idx == end) {
            /*
            int last_node = end;
            while (parents[last_node] != -1) {
                // return the SP
                relaxed_solution(last_node + nodes * parents[last_node], 0) = 1;
                last_node = parents[last_node];
            }
            */
           return current_node.dist;
        }
        
        Q.pop();
        visited[current_node.idx] = 1;
        //cout << "node " << current_node.idx << " popped" << endl;
        //iterate through neighbours of current_node
        if (node_neighbours.find(current_node.idx) != node_neighbours.end()) {
            //it refers to each <int,bool> pair
            for (map<int, bool>::iterator it = node_neighbours[current_node.idx].begin();
                 it != node_neighbours[current_node.idx].end(); it++) {
                
                int neighbour_node = it->first; 
                // has node been visited
                if (visited.find(neighbour_node) != visited.end()){
                    continue;
                }
            
               
                Edge current_edge = Edge(current_node.idx, neighbour_node);
                if (EIM.find(current_edge) != EIM.end()) {
                    current_edge_idx = EIM[current_edge];
                }
                else{
                    cerr << "edge not found in EIM --> djikstras" << endl;
                    return -1;
                }
                edge_weight = rc[current_edge_idx +  comm_idx*num_edges];

               

                // if distance of node is shorter, or if distances is same but violation is less
                    if ( distances[neighbour_node] > (distances[current_node.idx] + edge_weight)){
                        distances[neighbour_node] = distances[current_node.idx] + edge_weight;
                      
                        parents[neighbour_node] = current_node.idx;
                        Q.push({neighbour_node, distances[neighbour_node]});
                        //cout << "node " << neighbour_node << " pushed - dist " << distances[neighbour_node].first
                        //<< " viol - " << distances[neighbour_node].second << endl;
                    }
            }
        } else {
            cerr << "node not found in neighbours list --> djikstras" << endl;
            return -1;
        }
	}
    
    return -1;
}