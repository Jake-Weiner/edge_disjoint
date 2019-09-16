#include "MIPsolver_test.h"

void generate_nodeNeighbours_test(vector<vector<NodeEdgePair>>& node_neighbours){

    MIP_Solver ms;
    ms.generate_nodeNeighbours(node_neighbours);
    unordered_map<int, vector<int>> nN = ms.get_nodeNeighbours();
    //node_neighbours (from ED class) output
    for (int node_num = 0; node_num<node_neighbours.size(); node_num++ ){
        cout << "neighbours for node " << node_num << " are - ";
        for (auto& neighbours_info : node_neighbours[node_num]){
            cout << neighbours_info.first << " ";
        }
        cout << endl;
    }

    // nodeNeighbours output (from MIP_Solver class)
    for (auto& node_info : nN){
        cout << "neighbours for node " << node_info.first << " are - ";
        for (auto& neighbours : nN[node_info.first]){
            cout << neighbours << " ";
        }
        cout << endl;
    }

}