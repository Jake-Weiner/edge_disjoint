#include "MIPsolver.h"
#include "ED.h"



using namespace std;



void MIP_Solver::generate_nodeNeighbours(vector<vector<NodeEdgePair>>& node_neighbours){
    for (int node_idx = 0; node_idx < node_neighbours.size(); ++node_idx) {
        vector<int> neighbours;
        for (auto& neighbour_node_info : node_neighbours[node_idx]){
            int neighbour_node = neighbour_node_info.first;
            neighbours.push_back(neighbour_node);
        }
        nodeNeighbours[node_idx] = neighbours;
    }
}



// void MIP_Solver::test_generate_nodeNeighbours(){

// }

// IloRangeArray MIP_Solver::generate_EDP_constraints(IloEnv& env, IloNumVarArray& z_var, NumVar3Matrix& x_var){

    
//       // flow constraints
//     IloRangeArray constraints_to_add(env);
//     for (int k = 0; k < num_commodities; k++) {
//         for (int vertex = 0; vertex < node_neighbours.size(); ++vertex) {
//             IloExpr out_flow(env);
//             IloExpr in_flow(env);
//             IloExpr flow_constr(env);
//             for (auto& neighbours : node_neighbours[vertex]) {
//                 int neighbour_node = neighbours.first;
//                 out_flow += x_var[vertex][neighbour_node][k];
//                 in_flow += x_var[neighbour_node][vertex][k];
//             }
//             flow_constr = out_flow - in_flow;
//             if (vertex != commodities[k].origin && vertex != commodities[k].dest) {
//                 IloRange r1(env, flow_constr == 0);
//                 constraints_to_add.add(r1);
//             } else if (vertex == commodities[k].origin) {
//                 IloRange r1(env, flow_constr == z_var[k]);
//                 constraints_to_add.add(r1);
//             } else if (vertex == commodities[k].dest) {
//                 IloRange r1(env, flow_constr == -z_var[k]);
//                 constraints_to_add.add(r1);
//             }
//         }
       
//     }

// }

// MIP_Solver::generate_model(int edges, vector<Commodity>& commodities, vector<vector<NodeEdgePair>>& node_neighbours)
// {
//     int num_commodities = commodities.size();
//     int num_vertices = node_neighbours.size();
//     vector<bool> z;
//     z.resize(num_commodities, 0);
//     vector<bool> x;
//     x.resize(edges * num_commodities, 0);

//     IloEnv env;
//     IloModel model(env);
//     IloNumVarArray z_var(env);

//     for (int i = 0; i < num_commodities; i++) {
//         z_var.add(IloBoolVar(env));
//     }

//     NumVar3Matrix x_var(env, num_vertices);
//     for (int i = 0; i < num_vertices; i++) {
//         x_var[i] = NumVarMatrix(env, node_neighbours[i].size());
//         for (int j = 0; j < node_neighbours[i].size(); j++) {
//             x_var[i][j] = IloNumVarArray(env, num_commodities);
//             for (int k = 0; k < num_commodities; k++) {
//                 x_var[i][j][k] = IloBoolVar(env);
//             }
//         }
//     }

   
//     IloExpr obj_exp(env);
//     for (int i = 0; i < num_commodities; ++i) {
//         obj_exp += z_var[i];
//     }

//     IloObjective obj_fn = IloMaximize(env, obj_exp);
//     model.add(obj_fn);
   
   


//     // K=Set( k for k in keys(od_pairs) if od_pairs[k][source] in vertices && od_pairs[k][sink] in vertices)
//     // @variable(mip,x[i=vertices,j=neighbours[i],k=K],Bin, start=0)
//     // @variable(mip,z[k=K],Bin, start = 0) # is this commodity used?

//     // @constraint(mip,flow[i=vertices,k=K], sum(x[i,j,k] for j in neighbours[i]) - sum( x[j,i,k] for j in neighbours[i] ) ==
//     //                                               ((i==od_pairs[k][source] ? z[k] : 0) + (i==od_pairs[k][sink] ? -z[k] : 0)))
//     // # @constraint(mip,test_flow[k=K], )

//     // @constraint(mip,cap[i=vertices,j=neighbours[i];i<j], sum(x[i,j,k]+x[j,i,k] for k=K) <= 1)


//     model.add(constraints_to_add);
        
    
//     IloCplex cplex(model);
//     cplex.setParam(IloCplex::Threads, 1); // each particle gets 1 thread only
//     cplex.setOut(env.getNullStream());
        

//     // Optimize the problem and obtain solution.
//     if (!cplex.solve()) {
//         env.error() << "Failed to optimize LP" << endl;
//         throw(-1);
//     }

//     IloNumArray z_vals(env);
//         //populate y and return it

//     cout << "Solution status = " << cplex.getStatus() << endl;
//     cout << "Solution value  = " << cplex.getObjValue() << endl;
//     cplex.getValues(z_vals, z_var);

//         //cout << "Values        = " << vals << endl;

//         // double improvement = 0;
//         // for (int i = 0; i < vals.getSize(); i++) {
//         //     y[i] = floor(vals[i] + 0.1); // make sure we don't get rounding error
//         //     if (c[i] < 1)
//         //         improvement += (1 - vals[i]) * (1 - c[i]);
//         //     //cout << "i = " << i << " vals[i] = " << vals[i] << endl;
//         // }

//     model.remove(obj_fn);
//     model.remove(constraints_to_add);

//     // } catch (IloException& e) {
//     //     cout << e << endl;
//     // } catch (...) {
//     //     cout << "Unknown exception caught" << endl;
//     // }
//     env.end();

//     //return y;
// }
