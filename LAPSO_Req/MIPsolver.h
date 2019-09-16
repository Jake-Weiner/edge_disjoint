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

class MIP_Solver{
    public:
        MIP_Solver(){};
        ~MIP_Solver(){};
        IloRangeArray generate_EDP_constraints(IloEnv& env, IloNumVarArray& z_var, NumVar3Matrix& x_var);
        //generate_model()
        void generate_nodeNeighbours(std::vector<std::vector<NodeEdgePair>>& node_neighbours);
        std::unordered_map<int, std::vector<int>> get_nodeNeighbours(){return nodeNeighbours;}

    private:
        std::unordered_map<int, std::vector<int>> nodeNeighbours;

};

#endif

