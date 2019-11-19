#ifndef MIPSOLVER_H
#define MIPSOLVER_H


#include <ilcplex/ilocplex.h>
#include <fstream>
#include <iostream>
#include <list>
#include <vector>
#include <map>
#include "ED.h"

typedef IloArray<IloNumVarArray> NumVarMatrix;
typedef IloArray<IloArray<IloNumVarArray>> NumVar3Matrix;



class EDP_MIPsolver{
    public:
        EDP_MIPsolver(){};
        ~EDP_MIPsolver(){};
        
        //generate_model()
        void generate_nodeNeighbours(std::vector<std::vector<NodeEdgePair>>& node_neighbours);
        void solve_EDP(int edges, const vector<Commodity>& commodities, vector<vector<NodeEdgePair>>& node_neighbours);
        std::unordered_map<int, std::vector<int>> get_nodeNeighbours(){return nodeNeighbours;}

    private:

        IloRangeArray generate_EDP_constraints(IloEnv& env, const vector<Commodity>& commodities, std::unordered_map<int, std::vector<int>>& nodeNeighbours, 
            IloNumVarArray& z_var, NumVar3Matrix& x_var);
        std::unordered_map<int, std::vector<int>> nodeNeighbours;
        // x[3][2] --> x[3][0] , x[3][4] --> x[3][1] etc..
        std::unordered_map<std::pair<int, int>, int, pair_hash> neighbour_index_mapping; 
        //x x[3][0] --> x[3][2] , x[3][1] --> x[3][4]
        std::unordered_map<std::pair<int, int>, int, pair_hash> index_neighbour_mapping;
};

#endif

