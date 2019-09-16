#include <ilcplex/ilocplex.h>
#include <fstream>
#include <iostream>
#include <list>
#include <vector>


#include <map>


class MIP_Solver{
    public:
    MIP_Solver();
    ~MIP_Solver();
    
    IloRangeArray generate_EDP_constraints(IloEnv& env, IloNumVarArray& z_var, NumVar3Matrix& x_var);
    //generate_model()
    void generate_nodeNeighbours(vector<vector<NodeEdgePair>>& node_neighbours);
    void test_generate_nodeNeighbours();
    unordered_map<int, vector<int>> nodeNeighbours;

}