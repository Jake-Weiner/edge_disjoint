#include "EDP_MIPsolver.h"
#include "ED.h"

using namespace std;

void EDP_MIPsolver::generate_nodeNeighbours(vector<vector<NodeEdgePair>>& node_neighbours)
{
    for (int node_idx = 0; node_idx < node_neighbours.size(); ++node_idx) {
        vector<int> neighbours;
        for (auto& neighbour_node_info : node_neighbours[node_idx]) {
            int neighbour_node = neighbour_node_info.first;
            neighbours.push_back(neighbour_node);
        }
        nodeNeighbours[node_idx] = neighbours;
    }
}

IloRangeArray EDP_MIPsolver::generate_EDP_constraints(IloEnv& env, const vector<Commodity>& commodities, std::unordered_map<int, std::vector<int>>& nodeNeighbours,
    IloNumVarArray& z_var, NumVar3Matrix& x_var)
{

    int num_commodities = commodities.size();


    // flow constraints
    IloRangeArray constraints_to_add(env);
    for (auto& vertex_info : nodeNeighbours) {
        int vertex = vertex_info.first;
        for (int k = 0; k < num_commodities; k++) {
            IloExpr out_flow(env);
            IloExpr in_flow(env);
            IloExpr flow_constr(env);
            for (auto& neighbour : nodeNeighbours[vertex]) {
                int neighbour_node = neighbour;
                // identify the mapping key for neighbours
                pair<int, int> out_neighbour_idx_mapping_key = {vertex, neighbour_node};
                int out_neighbour_idx_mapping_value = neighbour_index_mapping[out_neighbour_idx_mapping_key];
        
                pair<int, int> in_neighbour_idx_mapping_key = { neighbour_node, vertex };
                int in_neighbour_idx_mapping_value = neighbour_index_mapping[in_neighbour_idx_mapping_key];
                
                out_flow += x_var[vertex][out_neighbour_idx_mapping_value][k];
                cout << "variable x[" << vertex << "][" << out_neighbour_idx_mapping_value << "][" << k << "] is used" << endl;
                in_flow += x_var[neighbour_node][in_neighbour_idx_mapping_value][k];
            }
            flow_constr = out_flow - in_flow;
            cout << "flow constraint " << flow_constr;
            cout << "vertex is " << vertex <<  "origin node is " << commodities[k].origin << " dest node is " << commodities[k].dest << endl;
            if (vertex != commodities[k].origin && vertex != commodities[k].dest) {
                IloRange r1(env, flow_constr == 0);
                cout << "= 0" << endl;
                constraints_to_add.add(r1);
            } else if (vertex == commodities[k].origin) {
                cout << "origin flow const = " << flow_constr << endl;
                IloRange r1(env, flow_constr == z_var[k]);
                cout << "= " << z_var[k] << endl;
                constraints_to_add.add(r1);
            } else if (vertex == commodities[k].dest) {
                IloRange r1(env, flow_constr == -z_var[k]);
                cout << "= -" << z_var[k] << endl;
                constraints_to_add.add(r1);
            }
        }
    }

    cout << "creating capacity constraints " << endl;
    //capacity_constraints
    for (int vertex = 0; vertex < nodeNeighbours.size(); ++vertex) {
        cout << "vertex checked -" << vertex << endl;
        for (auto& neighbour : nodeNeighbours[vertex]) {
            int neighbour_node = neighbour;
            if (vertex < neighbour_node) {
                pair<int, int> out_neighbour_idx_mapping_key = { vertex, neighbour_node };
                int out_neighbour_idx_mapping_value = neighbour_index_mapping[out_neighbour_idx_mapping_key];
                pair<int, int> in_neighbour_idx_mapping_key = { neighbour_node, vertex };
                int in_neighbour_idx_mapping_value = neighbour_index_mapping[in_neighbour_idx_mapping_key];
                IloExpr edge_out_flow(env);
                IloExpr edge_in_flow(env);
                IloExpr edge_cap(env);
                for (int k = 0; k < num_commodities; k++) {
                    edge_out_flow += x_var[vertex][out_neighbour_idx_mapping_value][k];
                    edge_in_flow += x_var[neighbour_node][in_neighbour_idx_mapping_value][k];
                }
                edge_cap = edge_in_flow + edge_out_flow;
                IloRange r1(env, edge_cap <= 1);
                constraints_to_add.add(r1);
            }
            else{
                cout << "larger neighbour node found" << endl;
            }
        }
    }

    return constraints_to_add;
}

void EDP_MIPsolver::solve_EDP(int edges, const vector<Commodity>& commodities, vector<vector<NodeEdgePair>>& node_neighbours)
{
    // generate node neighbours in a hashmap format
    generate_nodeNeighbours(node_neighbours);
    int num_commodities = commodities.size();
    int num_vertices = nodeNeighbours.size();
    vector<bool> z;
    z.resize(num_commodities, 0);
    vector<bool> x;
    x.resize(edges * num_commodities, 0);

    IloEnv env;
    IloModel model(env);
    IloNumVarArray z_var(env);
    IloRangeArray EDP_constraints(env);

    // z variables used as decision variables to send/not send commodities
    for (int i = 0; i < num_commodities; i++) {
        z_var.add(IloBoolVar(env));
    }

    // x variables used to represent whether an edge is used or not.
    NumVar3Matrix x_var(env, num_vertices);
    for (int vertex_index = 0; vertex_index < num_vertices; vertex_index++) {
        x_var[vertex_index] = NumVarMatrix(env, nodeNeighbours[vertex_index].size());
        for (int neighbour_index = 0; neighbour_index < nodeNeighbours[vertex_index].size(); neighbour_index++) {
            //create indexes for edge mapping to save space on problem size
            int neighbour_node = nodeNeighbours[vertex_index][neighbour_index];

            // map edge pair to a neighbour index x[1][2] --> x[1][0] x[1][4] --> x[1][1] etc.. 
            pair<int, int> neighbour_index_pair = {vertex_index, neighbour_node};
            neighbour_index_mapping[neighbour_index_pair] = neighbour_index;
            // mapping back from neighbour index to actual neighbour node
            pair<int, int> index_neighbour_pair = {vertex_index, neighbour_index};
            index_neighbour_mapping[index_neighbour_pair] = neighbour_node;
            x_var[vertex_index][neighbour_index] = IloNumVarArray(env, num_commodities);
            for (int commodity_index = 0; commodity_index < num_commodities; commodity_index++) {
                x_var[vertex_index][neighbour_index][commodity_index] = IloBoolVar(env);
            }
        }
    }

    IloExpr obj_exp(env);
    for (int i = 0; i < num_commodities; ++i) {
        obj_exp += z_var[i];
    }
    IloObjective obj_fn = IloMaximize(env, obj_exp);
    model.add(obj_fn);
    //generate flow and capacity constraints for the EDP problem
    EDP_constraints = generate_EDP_constraints(env, commodities, nodeNeighbours, z_var, x_var);
    model.add(EDP_constraints);

    // K=Set( k for k in keys(od_pairs) if od_pairs[k][source] in vertices && od_pairs[k][sink] in vertices)
    // @variable(mip,x[i=vertices,j=neighbours[i],k=K],Bin, start=0)
    // @variable(mip,z[k=K],Bin, start = 0) # is this commodity used?

    // @constraint(mip,flow[i=vertices,k=K], sum(x[i,j,k] for j in neighbours[i]) - sum( x[j,i,k] for j in neighbours[i] ) ==
    //                                               ((i==od_pairs[k][source] ? z[k] : 0) + (i==od_pairs[k][sink] ? -z[k] : 0)))
    // # @constraint(mip,test_flow[k=K], )

    // @constraint(mip,cap[i=vertices,j=neighbours[i];i<j], sum(x[i,j,k]+x[j,i,k] for k=K) <= 1)

    //model.add(constraints_to_add);

    IloCplex cplex(model);
    cplex.setParam(IloCplex::Threads, 1); // each particle gets 1 thread only
    cplex.setOut(env.getNullStream());

    // Optimize the problem and obtain solution.
    if (!cplex.solve()) {
        env.error() << "Failed to optimize LP" << endl;
        throw(-1);
    }

    IloNumArray z_vals(env);
    IloNumArray x_vals(env);
    //populate y and return it

    cout << "Solution status = " << cplex.getStatus() << endl;
    cout << "Solution value  = " << cplex.getObjValue() << endl;
    cplex.getValues(z_vals, z_var);

    for (int i = 0; i < x_var.getSize(); ++i) {
        for (int j = 0; j < x_var[i].getSize(); ++j) {
            for (int k = 0; k < num_commodities; ++k) {
                IloNum x_val = cplex.getValue(x_var[i][j][k]);
                // cout << x_val << endl;
                double x_val_rounded = floor(x_val + 0.1);
                int neighbour_node_idx = j;
                pair<int, int> edge_idx = { i, neighbour_node_idx };
                int neighbour_node = index_neighbour_mapping[edge_idx];
                cout << "x[" << i << "][" << neighbour_node << "][" << k << "] = " << x_val_rounded << endl;
              
            }
        }
    }

    cout << "Z Values        = " << z_vals << endl;

    // double improvement = 0;
    // for (int i = 0; i < vals.getSize(); i++) {
    //     y[i] = floor(vals[i] + 0.1); // make sure we don't get rounding error
    //     if (c[i] < 1)
    //         improvement += (1 - vals[i]) * (1 - c[i]);
    //     //cout << "i = " << i << " vals[i] = " << vals[i] << endl;
    // }

    // model.remove(obj_fn);
    // model.remove(constraints_to_add);

    // } catch (IloException& e) {
    //     cout << e << endl;
    // } catch (...) {
    //     cout << "Unknown exception caught" << endl;
    // }
    env.end();

    //return y;
}
