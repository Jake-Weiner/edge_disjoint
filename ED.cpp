#include "ED.h"
#include "Random.h"
#include "VolVolume.hpp"
#include "djikstra.h"
#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>
#include <fstream>
#include <iostream>
#include <list>
#include <map>

#include <type_traits>
#include <typeinfo>
#ifndef _MSC_VER
#include <cxxabi.h>
#endif
#include <cstdlib>
#include <memory>
#include <string>

template <class T>
std::string
type_name()
{
    typedef typename std::remove_reference<T>::type TR;
    std::unique_ptr<char, void (*)(void*)> own(
#ifndef _MSC_VER
        abi::__cxa_demangle(typeid(TR).name(), nullptr,
            nullptr, nullptr),
#else
        nullptr,
#endif
        std::free);
    std::string r = own != nullptr ? own.get() : typeid(TR).name();
    if (std::is_const<TR>::value)
        r += " const";
    if (std::is_volatile<TR>::value)
        r += " volatile";
    if (std::is_lvalue_reference<T>::value)
        r += "&";
    else if (std::is_rvalue_reference<T>::value)
        r += "&&";
    return r;
}

using namespace boost;
using namespace std;

ED::ED(string graph_filename, string pairs_filename, bool _printing, bool _randComm, bool _djikstras_naive,
    string _repair_remove_edge, string _repair_add_edge)
    : nEDsolves(0)
    , maxEDsolves(10)
{
    ED::populate_commodities(pairs_filename);
    ED::populate_graph(graph_filename);
    solution_edges.resize(commodities.size());
    printing = _printing;
    randComm = _randComm;
    dN = _djikstras_naive;
    repair_remove_edge = _repair_remove_edge;
    repair_add_edge = _repair_add_edge;
}

ED::~ED() {}

void ED::populate_commodities(string filename)
{
    int idx = 0;
    ifstream input(filename);
    split_vector_type SplitVec;
    if (input.is_open()) {
        while (!input.eof()) {
            string line;
            getline(input, line);
            split(SplitVec, line, is_any_of(" \t"), token_compress_on);
            if (SplitVec.size() > 1) {
                commodities.push_back(Commodity{ stoi(SplitVec[0]), stoi(SplitVec[1]), {}, 0.0, idx });
                idx++;
            } else
                break;
        }
    }
}

void ED::populate_graph(string filename)
{
    //populate the graph
    ifstream input(filename);
    split_vector_type SplitVec;
    size_t edge_number = 0;
    //double weight_init = 0;
    Edge current_edge;
    int start_node;
    int end_node;
    if (input.is_open()) {
        int line_count = 0;
        while (!input.eof()) {
            string line;
            getline(input, line); //read number
            split(SplitVec, line, is_any_of(" \t"), token_compress_on); // SplitVec == { "hello abc","ABC","aBc goodbye" }
            if (line_count == 0) {
                num_nodes = stoi(SplitVec[0]);
                //g = graph_t(num_nodes);
                line_count++;
                continue;
            }

            if (SplitVec.size() >= 3) {

                start_node = stoi(SplitVec[0]);
                end_node = stoi(SplitVec[1]);
                current_edge = Edge(start_node, end_node);
                if (EIM.find(current_edge) != EIM.end()) {
                    std::cout << "WARNING: " << filename << " contains 2 edges between "
                              << current_edge.first << " & " << current_edge.second
                              << " - ignored\n"; // this may give us
                    continue;
                }
                //add_edge(current_edge.first, current_edge.second, edge_number, g);
                graph_edges.push_back(current_edge);
                EIM[current_edge] = edge_number;
                // store which nodes are connected to each other

                if (node_neighbours[start_node].find(end_node) == node_neighbours[start_node].end()) {
                    node_neighbours[start_node][end_node] = 1;
                    //cout << start_node << " " << end_node << endl;
                }

                if (node_neighbours[end_node].find(start_node) == node_neighbours[end_node].end()) {
                    node_neighbours[end_node][start_node] = 1;
                    //cout << "added_node" << endl;
                }

                current_edge = Edge(current_edge.second, current_edge.first);
                // same edge number in both directions
                EIM[current_edge] = edge_number;

                edge_number += 1;
            }
        }
    }
    // he's removed this part?
    // add the dummy edge between orig and dest nodes
    // for (vector<Commodity>::iterator itr = commodities.begin(); itr < commodities.end(); itr++) {
    // current_edge = Edge((*itr).origin, (*itr).dest);
    // EIM[current_edge] = edge_number;
    // edge_number += 1;
    // weights.push_back(1);
    // }
    num_edges = graph_edges.size();
}

Status ED::reducedCost(const Particle& p, DblVec& redCost)
{
    // update reduced costs for all edges except for the origin->node dummy edge weight
    for (size_t c = 0; c < commodities.size(); ++c) {
        for (EdgeIter ei = graph_edges.begin(); ei != graph_edges.end(); ei++) {
            if (EIM.find(*ei) != EIM.end()) {
                redCost[primalIdx(EIM[*ei], c)] = -p.dual[EIM[*ei]];
            }
        }
    }
    return OK;
}

void ED::update_comm_sol(EDParticle& p, double SP, vector<int> parents, double& total_paths_cost, int index,
    int start, int end, bool printing)
{
    int current;
    if (SP >= 1) {
        if (printing)
            std::cout << "\tCost for " << p.commodities[index].comm_idx << " (" << (p.commodities[index]).origin
                      << "," << (p.commodities[index]).dest << ") " << SP << std::endl;
        p.ub += 1;
        total_paths_cost += 1;
        p.commodities[index].solution_value = 1; // penalty for not including path

    } else {

        p.commodities[index].solution_value = 0;
        // Iterate over path and add to primal solution

        if (printing == true) {
            std::cout << "\tPath for " << p.commodities[index].comm_idx << " (" << p.commodities[index].origin
                      << "," << p.commodities[index].dest << ") " << SP << ": ";
        }

        for (current = end; current != start; current = parents[current]) {
            if (printing == true)
                std::cout << " " << current;
            if (parents[current] == -1) {
                cout << "issue with update_comm_sol - parents array incorrect" << endl;
                // exit;
            }
            Edge current_edge = Edge(parents[current], current);
            p.x[primalIdx(EIM[current_edge], p.commodities[index].comm_idx)] += 1;
            // solution is stored in reverse at here
            p.commodities[index].solution_edges.push_back(current_edge);
        }
        if (printing == true)
            std::cout << " " << start << std::endl;
        // reversing each time is not really necessary but nice
        reverse(p.commodities[index].solution_edges.begin(), p.commodities[index].solution_edges.end());
        total_paths_cost += SP;
    }
}

Status ED::solveSubproblem(Particle& p_)
{

    ++nEDsolves;
    EDParticle& p(static_cast<EDParticle&>(p_));
    //clear previous primal sol
    p.x = 0; // clear would resize the array

    // create/update the edge_weights in the graph using reduced costs
    // solve each commodity pair independantly

    double total_paths_cost = 0.0;
    vector<int> parents;

    int start;
    int end;
    int current;
    p.ub = p.commodities.size();
    p.lb = 0;
    double max_perturb = 0.0;

    for (size_t i = 0; i < p.perturb.size(); i++) {
        if (p.perturb[i] > max_perturb) {
            max_perturb = p.perturb[i];
            continue;
        }
        if (-1 * p.perturb[i] > max_perturb)
            max_perturb = -1 * p.perturb[i];
    }

    //solve subproblem selecting commodity iteration order randomly
    if (randComm) {
        // randomise order of commodities
        vector<int> random_indices;
        for (int i = 0; i < p.commodities.size(); i++) {
            random_indices.push_back(i);
        }
        random_shuffle(random_indices.begin(), random_indices.end());

        // loop through commodities
        for (int i = 0; i < random_indices.size(); i++) {
            //int random_index = random_indices[i];
            int random_index = i;

            //solve SP
            p.commodities[random_index].solution_edges.clear();
            start = p.commodities[random_index].origin;
            end = p.commodities[random_index].dest;
            double SP;

            //ensure no negative edge weights
            for (int i = 0; i < p.rc.size(); i++) {
                if (p.rc[i] < 0) {
                    p.rc[i] = 0;
                }
            }
            if (dN) { // djikstras without violation consideration
                SP = djikstras_naive(EIM, node_neighbours, start, end, parents, num_nodes, num_edges,
                    p.rc, p.x, p.commodities[random_index].comm_idx, commodities.size());
            }

            else {
                SP = djikstras(max_perturb, EIM, node_neighbours, start, end, parents, num_nodes, num_edges,
                    p.rc, p.x, p.commodities[random_index].comm_idx, commodities.size());
            }
            if (SP == -1) {
                cerr << "error with djikstras" << endl;
                exit;
            }

            Commodity_SP sp = {
                start, end, parents
            };
            //update solution for current commodity
            //shortest path info required to update solutions
            p.commodity_shortest_paths[p.commodities[random_index].comm_idx] = sp;
            // p.c is cost of path -- used in MIP solver
            p.c[p.commodities[random_index].comm_idx] = SP;
            //update_comm_sol(p, SP, parents, total_paths_cost, random_index, start, end, printing);
        }
    }

    //iterate through commodities in standard order
    else {
        int loop_idx = 0;
        for (vector<Commodity>::iterator itr = p.commodities.begin(); itr < p.commodities.end(); ++itr) {
            const int comm_idx = itr->comm_idx;

            itr->solution_edges.clear();

            //solve SP
            p.commodities[loop_idx].solution_edges.clear();
            start = p.commodities[loop_idx].origin;
            end = p.commodities[loop_idx].dest;
            double SP;
            if (dN) { // djikstras without violation consideration
                SP = djikstras_naive(EIM, node_neighbours, start, end, parents, num_nodes, num_edges,
                    p.rc, p.x, p.commodities[loop_idx].comm_idx, commodities.size());
            } else {
                SP = djikstras(max_perturb, EIM, node_neighbours, start, end, parents, num_nodes, num_edges,
                    p.rc, p.x, p.commodities[loop_idx].comm_idx, commodities.size());
            }

            if (SP == -1) {
                cerr << "error with djikstras" << endl;
                exit;
            }
            //update solution for current commodity
            update_comm_sol(p, SP, parents, total_paths_cost, loop_idx, start, end, printing);
            loop_idx++;
        }
    }

    vector<int> y = solve_mip(p);
    int path_saved = 0;

    for (int i = 0; i < y.size(); i++) {
        if (y[i] == 1) {
            vector<int> parents = p.commodity_shortest_paths[i].parents;
            int start = p.commodity_shortest_paths[i].start;
            int end = p.commodity_shortest_paths[i].end;
            double SP = p.c[i];
            update_comm_sol(p, SP, parents, total_paths_cost, i, start, end, printing);
        } else {
            if (p.c[i] < 1) {
                path_saved += 1;
            }
            update_comm_sol(p, 1.5, parents, total_paths_cost, i, start, end, printing);
        }
    }
    p.path_saved = path_saved;
    p.commodity_shortest_paths.clear();
    p.commodity_shortest_paths.resize(commodities.size());

    // edges violated in the ED solution
    p.viol = 1; // sets every violation to 1
    for (size_t i = 0; i < p.x.size(); i++) {
        p.viol[edgeIdx(i)] -= p.x[i]; // sum over all commodites
    }

    // is feasible if no edges are violated
    p.isFeasible = (p.viol.min() >= 0);
    p.lb = (total_paths_cost + p.dual.sum()) - num_edges * max_perturb;
    //update particles best local solution
    if (p.lb > p_.best_lb) {
        p_.best_lb = p.lb;
        double sum_viol = 0;
        for (int i = 0; i < p.viol.size(); i++) {
            // count every violation
            if (p.viol[i] < 0) {
                sum_viol += p.viol[i];
            }
        }
        p_.best_lb_viol = sum_viol;
        p_.best_lb_sol.clear();
        for (vector<Commodity>::iterator itr = p.commodities.begin(); itr < p.commodities.end(); ++itr) {
            for (EdgeIter E = itr->solution_edges.begin(); E != itr->solution_edges.end(); E++) {
                p_.best_lb_sol.push_back(*E);
            }
        }
        if (p_.best_lb_sol.empty()) {
            std::cout << "best sol is empty" << endl;
        }
    }

    if (p.isFeasible && (p.ub < p_.best_ub)) {
        for (vector<Commodity>::iterator itr = p.commodities.begin(); itr < p.commodities.end(); ++itr) {
            if (p_.best_ub_sol.find(itr->comm_idx) != p_.best_ub_sol.end()) {
                p_.best_ub_sol[itr->comm_idx].clear();
            }
            for (EdgeIter E = itr->solution_edges.begin(); E != itr->solution_edges.end(); E++) {
                p_.best_ub_sol[itr->comm_idx].push_back(*E);
            }
        }

        p_.best_ub = p.ub;
    }

    if (printing == true) {
        std::cout << "Subproblem solve " << nEDsolves << "/" << maxEDsolves << ": "
                  << " lb=" << p.lb << " ub=" << p.ub
                  << "\trange of dual = " << p.dual.min() << " to " << p.dual.max() << std::endl
                  << "\trange of viol = " << p.viol.min() << " to " << p.viol.max() << std::endl;
        if (!p.perturb.empty())
            std::cout << "\trange of perturb = " << p.perturb.min() << " to " << p.perturb.max() << std::endl;
        std::cout << "\tpath edge counts:";
        for (auto c = p.commodities.begin(); c != p.commodities.end(); ++c)
            std::cout << " " << c->solution_edges.size();
        std::cout << std::endl;
    }

    return (nEDsolves < maxEDsolves) ? OK : ABORT;
}

Status ED::fixConstraint(const int constraint,
    const Particle& p,
    SparseVec& feas) { return OK; }

Status ED::heuristics(Particle& p_)
{

    EDParticle& p(static_cast<EDParticle&>(p_));
    if (p.isFeasible) {
        if (printing)
            std::cout << "Heuristic called with feasible solution "
                      << p.ub << std::endl;
        return OK; // already feasible
    }
    int largest_com_idx = 0;
    int largest_viol = 0;
    int current_viol = 0;

    vector<int> cut_set_commodities;
    int cut_set_edges;
    vector<int> potential_cutset;
    map<int, bool> S_cutSet;
    map<int, bool> T_cutSet;
    IntVec viol(p.viol.size(), 1);
    p.x = 0;
    int viol_sum = 0;
    // iterate over solution edges --> locally set viol intvec and make sure p.x is correct
    for (auto it = p.commodities.begin(); it != p.commodities.end(); it++)
        for (const Edge& e : it->solution_edges) {
            p.x[primalIdx(EIM[e], it->comm_idx)] += 1;
            viol[edgeIdx(e)] -= 1;
        }

    for (int i = 0; i < viol.size(); i++) {
        if (viol[i] < 0) {
            viol_sum += 1;
        }
    }
    p.viol_sum = viol_sum;

    while (viol.min() < 0) {
        if (printing == true)
            cout << "\tviol sum is " << viol.sum() << endl;

        if (repair_remove_edge.compare("random") == 0) {
            cout << "using random repair" << endl;
            vector<int> random_indices;
            for (int i = 0; i < p.commodities.size(); i++) {
                random_indices.push_back(i);
            }
            random_shuffle(random_indices.begin(), random_indices.end());

            // loop through commodities randomly
            for (int i = 0; i < random_indices.size(); i++) {
                //int random_index = random_indices[i];
                int random_index = i;
                bool remove_path = false;
                if (!(p.commodities[random_index].solution_edges.empty())) {
                    for (const Edge& e : p.commodities[random_index].solution_edges) {
                        //contains an edge which is violated
                        if (viol[edgeIdx(e)] < 0) {
                            remove_path = true;
                            break;
                        }
                    }
                }
                if (remove_path == true) {
                    for (const Edge& e : p.commodities[random_index].solution_edges) {
                        p.x[primalIdx(edgeIdx(e), p.commodities[random_index].comm_idx)] -= 1;
                        viol[edgeIdx(e)] += 1;
                    }
                    p.commodities[random_index].solution_edges.clear();

                    //try and re-add in path

                    DblVec temp_rc;
                    temp_rc.resize(p.rc.size());
                    for (int i = 0; i < p.rc.size(); i++) {
                        temp_rc[i] = (viol[edgeIdx(i)] == 1) ? 1 : num_edges;
                    }

                    vector<int> parents;
                    int start = p.commodities[random_index].origin;
                    int end = p.commodities[random_index].dest;
                    int current;
                    double SP = djikstras_naive(EIM, node_neighbours, start, end, parents, num_nodes, num_edges,
                        temp_rc, p.x, p.commodities[random_index].comm_idx, commodities.size());
                    // if no feasible sp exists between orig and dest nodes

                    double thresh = num_edges;
                    if (SP < thresh) {
                        //if (SP < 1) {

                        // Iterate over path and add to primal solution
                        for (current = end; current != start; current = parents[current]) {
                            if (parents[current] == -1) {
                                cerr << "issue with repair - parents array incorrect" << endl;
                                exit;
                            }
                            const Edge e(Edge(parents[current], current));
                            p.x[primalIdx(edgeIdx(e), random_index)] += 1;
                            // solution is stored in reverse at here

                            p.commodities[random_index].solution_edges.push_back(e);
                            viol[edgeIdx(e)] -= 1;
                        }
                        // reversing each time is not really necessary but nice
                        reverse(p.commodities[random_index].solution_edges.begin(), p.commodities[random_index].solution_edges.end());
                        if (printing == true)
                            cout << "\t\trepaired path" << endl;
                    }
                }
            }

        }

        else {

            auto largest_edge_violation = std::min_element(viol.begin(), viol.end());
            // identify edge_idx with the largest violation
            int largest_viol_idx = distance(viol.begin(), largest_edge_violation);
            if (printing == true) {
                cout << "\tmin viol is at index " << largest_viol_idx
                     << " with value " << viol[largest_viol_idx] << endl;
            }

            largest_viol = 0;
            double largest_perturb = -1;

            if (repair_remove_edge.compare("largest_viol") == 0) {
                //find commodity with the largest violation, that contains this edge
                for (auto it = p.commodities.begin(); it != p.commodities.end(); it++) {
                    if (p.x[primalIdx(largest_viol_idx, it->comm_idx)] > 0) { // the commodity contains this edge
                        current_viol = 0;
                        for (size_t i = 0; i < it->solution_edges.size(); i++) { // sum up violations this commodity has
                            if (viol[edgeIdx(it->solution_edges[i])] < 1.0e-6) {
                                current_viol += viol[edgeIdx(it->solution_edges[i])];
                            }
                            if (current_viol < largest_viol) {
                                largest_com_idx = it->comm_idx;
                                largest_viol = current_viol;
                            }
                        }
                    }
                }
            } else if (repair_remove_edge.compare("perturb") == 0) {
                //remove with highest perturbation value (this commodity is least likely to use this edge)
                for (auto it = p.commodities.begin(); it != p.commodities.end(); it++) {
                    if (p.x[primalIdx(largest_viol_idx, it->comm_idx)] > 1.0e-6) { // the commodity contains this edge

                        if (p.perturb[primalIdx(largest_viol_idx, it->comm_idx)] > largest_perturb) {
                            largest_perturb = p.perturb[primalIdx(largest_viol_idx, it->comm_idx)];
                            largest_com_idx = it->comm_idx;
                        }
                    }
                }
            }

            else {
                cout << "repair method not set properly" << endl;
                exit(1);
            }

            if (printing == true)
                cout << "\tlargest comm idx " << largest_com_idx << " with " << -largest_viol << endl;
            // remove this path from both x AND from the solution_edges
            for (const Edge& e : p.commodities[largest_com_idx].solution_edges) {
                p.x[primalIdx(edgeIdx(e), largest_com_idx)] -= 1;
                viol[edgeIdx(e)] += 1;
            }
            p.commodities[largest_com_idx].solution_edges.clear();

            // try find a new path for this OD-PAIR, using perturbations

            IntVec distances_heur(num_nodes); // integer distances
            DblVec temp_rc;
            temp_rc.resize(p.rc.size());
            double min_perturb = p.perturb.min();

            if (min_perturb > 0) {
                min_perturb = 0;
            }

            //cout << min_perturb << endl;

            for (int i = 0; i < p.rc.size(); i++) {
                if (repair_add_edge.compare("pert_repair_0") == 0) {
                    temp_rc[i] = (viol[edgeIdx(i)] == 1) ? p.perturb[i] : 1;
                } else if (repair_add_edge.compare("pert_repair_min") == 0) {
                    temp_rc[i] = (viol[edgeIdx(i)] == 1) ? p.perturb[i] - min_perturb : 1;
                    if (temp_rc[i] < 0)
                        cout << temp_rc[i] << endl;
                } else if (repair_add_edge.compare("rc_repair") == 0) {
                    temp_rc[i] = (viol[edgeIdx(i)] == 1) ? p.rc[i] : 1;
                } else if (repair_add_edge.compare("arb_repair") == 0) {
                    temp_rc[i] = (viol[edgeIdx(i)] == 1) ? 1 : num_edges;
                } else {
                    cout << "repair method - add edge not set properly" << endl;
                }

                if (temp_rc[i] < 0) {
                    temp_rc[i] = 0;
                }
                //temp_rc[i] = (viol[edgeIdx(i)] == 1) ? p.rc[i] : 1;
            }

            vector<int> parents;
            int start = p.commodities[largest_com_idx].origin;
            int end = p.commodities[largest_com_idx].dest;
            int current;

            // if no feasible sp exists between orig and dest nodes

            double thresh;
            if (repair_add_edge.compare("arb_repair") == 0) {
                thresh = num_edges;
            } else {
                thresh = 1;
            }
            double SP = djikstras_naive_cutSet(EIM, node_neighbours, start, end, parents, num_nodes, num_edges,
                temp_rc, p.x, p.commodities[largest_com_idx].comm_idx, commodities.size(), S_cutSet, T_cutSet, thresh);

            if (SP < thresh) {
                //if (SP < 1) {

                // Iterate over path and add to primal solution
                for (current = end; current != start; current = parents[current]) {
                    if (parents[current] == -1) {
                        cerr << "issue with repair - parents array incorrect" << endl;
                        exit;
                    }
                    const Edge e(Edge(parents[current], current));
                    p.x[primalIdx(edgeIdx(e), largest_com_idx)] += 1;
                    // solution is stored in reverse at here

                    p.commodities[largest_com_idx].solution_edges.push_back(e);
                    viol[edgeIdx(e)] -= 1;
                }
                // reversing each time is not really necessary but nice
                reverse(p.commodities[largest_com_idx].solution_edges.begin(), p.commodities[largest_com_idx].solution_edges.end());
                if (printing == true)
                    cout << "\t\trepaired path" << endl;
            }

            else {
                cut_set_commodities = find_cutset_commodities(p, S_cutSet, T_cutSet);
                cut_set_edges = find_cutset_edges(S_cutSet, T_cutSet);
                p.cutsets.push_back(cut_set_commodities); // include commodities involved in cutset e.g {{0,1,2}, {1,2,4}, {3,6,9} etc...}
                p.cut_set_sizes.push_back(cut_set_edges); // number of edges connecting the 2 sets S and T
                if (cut_set_commodities.size() > cut_set_edges) {
                    p_.local_constraints.push_back({ cut_set_commodities, cut_set_edges });
                }

                //reset cut_set_edges and cut_set_commodities
                //reset_vector<int>(cut_set_commodities);
                cut_set_edges = 0;
                S_cutSet.clear();
                T_cutSet.clear();
            }
        }
    }

    vector<Commodity> commodities_to_add;
    for (auto it = p.commodities.begin(); it != p.commodities.end(); it++) {
        // commodity checks and reset p.x
        if (it->solution_edges.empty()) {
            commodities_to_add.push_back(*it);
        }
    }
    DblVec temp_rc;
    temp_rc.resize(p.rc.size());
    for (auto it = commodities_to_add.begin(); it != commodities_to_add.end(); it++) {
        int comm_idx = it->comm_idx;
        double min_perturb = p.perturb.min();
        for (int i = 0; i < p.rc.size(); i++) {
            if (repair_add_edge.compare("pert_repair_0") == 0) {
                temp_rc[i] = (viol[edgeIdx(i)] == 1) ? p.perturb[i] : 1;
            } else if (repair_add_edge.compare("pert_repair_min") == 0) {
                temp_rc[i] = (viol[edgeIdx(i)] == 1) ? p.perturb[i] - min_perturb : 1;
            } else if (repair_add_edge.compare("rc_repair") == 0) {
                temp_rc[i] = (viol[edgeIdx(i)] == 1) ? p.rc[i] : 1;
            } else if (repair_add_edge.compare("arb_repair") == 0) {
                temp_rc[i] = (viol[edgeIdx(i)] == 1) ? 1 : num_edges;
            }

            if (temp_rc[i] < 0) {
                temp_rc[i] = 0 + 1e-16;
            }
        }

        // calculate theshold that paths require for successful path connection
        double thresh;
        if (repair_add_edge.compare("arb_repair") == 0) {
            thresh = num_edges;
        } else {
            thresh = 1;
        }

        vector<int> parents;
        int start = it->origin;
        int end = it->dest;
        int current;
        vector<int> cutset;
        double SP = djikstras_naive_cutSet(EIM, node_neighbours, start, end, parents, num_nodes, num_edges,
            temp_rc, p.x, comm_idx, commodities.size(), S_cutSet, T_cutSet, thresh);

        if (SP < thresh) {

            // Iterate over path and add to primal solution
            for (current = end; current != start; current = parents[current]) {
                if (parents[current] == -1) {
                    cerr << "issue with repair - parents array incorrect" << endl;
                    exit;
                }
                const Edge e(Edge(parents[current], current));
                p.x[primalIdx(edgeIdx(e), largest_com_idx)] += 1;
                // solution is stored in reverse at here
                p.commodities[largest_com_idx].solution_edges.push_back(e);
                viol[edgeIdx(e)] -= 1;
            }
            // reversing each time is not really necessary but nice
            reverse(p.commodities[largest_com_idx].solution_edges.begin(), p.commodities[largest_com_idx].solution_edges.end());
            if (printing == true)
                cout << "\t\trepaired path" << endl;

        } else {
            cut_set_commodities = find_cutset_commodities(p, S_cutSet, T_cutSet);
            cut_set_edges = find_cutset_edges(S_cutSet, T_cutSet);
            p.cutsets.push_back(cut_set_commodities); // include commodities involved in cutset e.g {{0,1,2}, {1,2,4}, {3,6,9} etc...}
            p.cut_set_sizes.push_back(cut_set_edges);

            // add this information to the MIP structures
            if (cut_set_commodities.size() > cut_set_edges) {
                p_.local_constraints.push_back({ cut_set_commodities, cut_set_edges });
            }

            //reset cut_set_edges and cut_set_commodities
            //reset_vector<int>(cut_set_commodities);
            cut_set_edges = 0;
            S_cutSet.clear();
            T_cutSet.clear();
        }
    }
    p.ub = 0;

    for (auto it = p.commodities.begin(); it != p.commodities.end(); it++) {
        if (it->solution_edges.empty())
            p.ub += 1;
    }

    p.isFeasible = true;
    if (printing == true)
        cout << "Heuristic found solution with " << p.ub << " missing paths\n";

    if (p.isFeasible && (p.ub < p_.best_ub)) {
        for (vector<Commodity>::iterator itr = p.commodities.begin(); itr < p.commodities.end(); ++itr) {
            if (p_.best_ub_sol.find(itr->comm_idx) != p_.best_ub_sol.end()) {
                p_.best_ub_sol[itr->comm_idx].clear();
            }
            for (EdgeIter E = itr->solution_edges.begin(); E != itr->solution_edges.end(); E++) {
                p_.best_ub_sol[itr->comm_idx].push_back(*E);
            }
        }
        p_.best_ub = p.ub;
    }
    return (nEDsolves < maxEDsolves) ? OK : ABORT;
}

Status ED::updateBest(Particle& p_)
{
    EDParticle& p(static_cast<EDParticle&>(p_));
    if (solution_edges.size() != commodities.size())
        solution_edges.resize(commodities.size()); // should not be necessary
    solution_cost = 0;
    for (size_t c = 0; c < commodities.size(); ++c) {
        solution_edges[c] = p.commodities[c].solution_edges; // copy path
        if (p.commodities[c].solution_edges.empty())
            ++solution_cost;
    }
    if (printing == true)
        cout << "################## " << solution_cost << " ############\n";
    return OK;
}

void ED::localSearch(Particle& p_)
{

    EDParticle& p(static_cast<EDParticle&>(p_));

    if (!p.isFeasible) {
        return;
    }

    IntVec viol(p.viol.size(), 1);
    DblVec temp_rc;
    temp_rc.resize(p.rc.size());

    // need to double check why this needs to be reset
    p.x = 0;
    vector<Commodity> commodities_to_add;
    for (auto it = p.commodities.begin(); it != p.commodities.end(); it++) {
        // commodity checks and reset p.x
        if (it->solution_edges.empty()) {
            commodities_to_add.push_back(*it);
        } else {
            for (const Edge& e : it->solution_edges) {
                p.x[primalIdx(EIM[e], it->comm_idx)] += 1;
                viol[edgeIdx(e)] -= 1;
            }
        }
    }
    //randomise commodity iteration order
    random_shuffle(commodities_to_add.begin(), commodities_to_add.end());

    double perturb_min = p.perturb.min();
    if (perturb_min > 0) {
        perturb_min = 0;
    }
    //try and add in commodities 1 at a time using previously_unused edges
    for (auto it = commodities_to_add.begin(); it != commodities_to_add.end(); it++) {
        int comm_idx = it->comm_idx;
        double min_perturb = p.perturb.min();
        for (int i = 0; i < p.rc.size(); i++) {
            if (repair_add_edge.compare("pert_repair_0") == 0) {
                temp_rc[i] = (viol[edgeIdx(i)] == 1) ? p.perturb[i] : 1;
            } else if (repair_add_edge.compare("pert_repair_min") == 0) {
                temp_rc[i] = (viol[edgeIdx(i)] == 1) ? p.perturb[i] - perturb_min : 1;
            } else if (repair_add_edge.compare("rc_repair") == 0) {
                temp_rc[i] = (viol[edgeIdx(i)] == 1) ? p.rc[i] : 1;
            } else if (repair_add_edge.compare("arb_repair") == 0) {
                temp_rc[i] = (viol[edgeIdx(i)] == 1) ? 1 : num_edges;
            }

            if (temp_rc[i] < 0) {
                temp_rc[i] = 0 + 1e-16;
            }
        }
        vector<int> parents;
        int start = it->origin;
        int end = it->dest;
        int current;
        double SP = djikstras_naive(EIM, node_neighbours, start, end, parents, num_nodes, num_edges,
            temp_rc, p.x, comm_idx, commodities.size());
        // if no feasible sp exists between orig and dest nodes

        double thresh;
        if (repair_add_edge.compare("arb_repair") == 0) {
            thresh = num_edges;
        } else {
            thresh = 1;
        }
        if (SP < thresh) {
            //if (SP < 1) {
            //cout << "added commodity in local_search" << endl;
            // Iterate over path and add to primal solution
            for (current = end; current != start; current = parents[current]) {
                if (parents[current] == -1) {
                    cerr << "issue with local search - parents array incorrect" << endl;
                    exit;
                }
                const Edge e(Edge(parents[current], current));
                p.x[primalIdx(edgeIdx(e), comm_idx)] += 1;
                // solution is stored in reverse at here
                p.commodities[comm_idx].solution_edges.push_back(e);
                viol[edgeIdx(e)] -= 1;
            }
            // reversing each time is not really necessary but nice
            reverse(p.commodities[comm_idx].solution_edges.begin(), p.commodities[comm_idx].solution_edges.end());
        }
    }

    //make sure solution is feasible
    if (!(viol.min() < 0)) {
        return;
    }

    p.ub = 0;

    //recalculate ub
    for (auto it = p.commodities.begin(); it != p.commodities.end(); it++) {
        if (it->solution_edges.empty())
            p.ub += 1;
    }

    //update best_ub solution
    if (p.ub < p_.best_ub) {
        for (vector<Commodity>::iterator itr = p.commodities.begin(); itr < p.commodities.end(); ++itr) {
            if (p_.best_ub_sol.find(itr->comm_idx) != p_.best_ub_sol.end()) {
                p_.best_ub_sol[itr->comm_idx].clear();
            }
            for (EdgeIter E = itr->solution_edges.begin(); E != itr->solution_edges.end(); E++) {
                p_.best_ub_sol[itr->comm_idx].push_back(*E);
            }
        }
        p_.best_ub = p.ub;
    }
}

void ED::write_mip(vector<Particle*>& non_dom, double lb, double ub, string outfile_name)
{
    /*
    std::ofstream outfile;
    try {
        outfile.open(outfile_name, std::ios_base::app);
    } catch (std::ofstream::failure e) {
        std::cerr << "Exception opening output file\n";
        return;
    }

    outfile << num_edges(g) << endl;
    outfile << num_vertices(g) << endl;
    outfile << lb << endl;
    outfile << ub << endl;
    */
}
/* Stuff for emacs/xemacs to conform to the "Visual Studio formatting standard"
 * Local Variables:
 * tab-width: 4
 * eval: (c-set-style "stroustrup")
 * End:
*/

vector<int> ED::find_cutset_commodities(EDParticle& p, map<int, bool>& S_cutSet, map<int, bool>& T_cutSet)
{
    vector<int> cut_set;
    for (vector<Commodity>::iterator comm_it = p.commodities.begin(); comm_it < p.commodities.end(); comm_it++) {
        int commodity_orig = comm_it->origin;
        int commodity_dest = comm_it->dest;
        if ((S_cutSet.find(commodity_orig) != S_cutSet.end()) && (T_cutSet.find(commodity_dest) != T_cutSet.end())) {
            cut_set.push_back(comm_it->comm_idx);
        }
        if ((T_cutSet.find(commodity_orig) != T_cutSet.end()) && (S_cutSet.find(commodity_dest) != S_cutSet.end())) {
            cut_set.push_back(comm_it->comm_idx);
        }
    }
    return cut_set;
}

int ED::find_cutset_edges(map<int, bool>& S_cutSet, map<int, bool>& T_cutSet)
{

    int cut_set_edges = 0;
    for (map<int, bool>::iterator S_cut_it = S_cutSet.begin(); S_cut_it != S_cutSet.end(); S_cut_it++) {
        for (map<int, bool>::iterator T_cut_it = T_cutSet.begin(); T_cut_it != T_cutSet.end(); T_cut_it++) {
            if (node_neighbours[S_cut_it->first].find(T_cut_it->first) != node_neighbours[S_cut_it->first].end()) {
                cut_set_edges += 1;
            }
        }
    }

    return cut_set_edges;
}

/*
void ED::initialise_global_constraints(){
    global_constraints = IloRangeArray(env);
}
*/

void ED::add_constraints_mip(vector<pair<vector<int>, int>>& local_constraints)
{

    vector<int> constraint_vars;
    int cut_set_edges;
    
    for (vector<pair<vector<int>, int>>::iterator cs_it = local_constraints.begin(); cs_it != local_constraints.end(); cs_it++) {
        
        constraint_vars = cs_it->first;
        cut_set_edges = cs_it->second;

        if (constraint_map.find(constraint_vars) == constraint_map.end()) {
            constraint_map[constraint_vars] = cut_set_edges;
        } else {
            if (constraint_map[constraint_vars] > cut_set_edges) {
                constraint_map[constraint_vars] = cut_set_edges;
            }
        }
    }
}

vector<int> EDParticle::solve_mip(map<vector<int>, int>& constraint_map)
{
    vector<int> y;
    y.resize(c.size(), 0);
    IloEnv env;
    IloModel model(env);
    IloNumVarArray var(env);
    IloRangeArray con(env);

    for (int i = 0; i < c.size(); i++) {
        var.add(IloBoolVar(env));
    }

    try {
        IloExpr con_exp(env);
        IloRangeArray constraints_to_add(env);
        IloExpr con_exp_temp(env);
        for (map<vector<int>, int>::iterator it = constraint_map.begin(); it != constraint_map.end(); it++) {
            IloExpr con_exp(env);

            //constraint pair <variables, |Cutset Edges|>
            vector<int> constraint_vars = it->first;
            int constraint_bound = it->second;
            for (vector<int>::iterator cs_it = constraint_vars.begin(); cs_it != constraint_vars.end(); cs_it++) {
                con_exp += var[*cs_it];
            }
            IloRange r1(env, 0, con_exp, constraint_bound);
            constraints_to_add.add(r1);
        }

        model.add(constraints_to_add);
        double cost;
        IloExpr obj_exp(env);
        for (int i = 0; i < c.size(); i++) {
            if (c[i] == 1) {
                cost = 1.01;
            } else if (c[i] < 0) {
                cost = 0;
            } else {
                cost = c[i];
            }
            obj_exp += (1 - cost) * var[i];
        }
        IloObjective obj_fn = IloMaximize(env, obj_exp);
        model.add(obj_fn);

        IloCplex cplex(model);

        //cplex.exportModel("lpex1.lp");

        // Optimize the problem and obtain solution.
        if (!cplex.solve()) {
            env.error() << "Failed to optimize LP" << endl;
            throw(-1);
        }

        IloNumArray vals(env);
        //populate y and return it

        //cout << "Solution status = " << cplex.getStatus() << endl;
        //cout << "Solution value  = " << cplex.getObjValue() << endl;
        cplex.getValues(vals, var);

        //cout << "Values        = " << vals << endl;

        cout << "vals size is " << vals.getSize() << endl;
        for (int i = 0; i < vals.getSize(); i++) {
            y[i] = vals[i];
            //cout << "i = " << i << " vals[i] = " << vals[i] << endl;
        }

        model.remove(obj_fn);
        cout << "number of constraints are " << cplex.getNrows() << endl;
        model.remove(constraints_to_add);

    } catch (IloException& e) {
        cout << e << endl;
    } catch (...) {
        cout << "Unknown exception caught" << endl;
    }
    env.end();
    return y;
}

vector<int> ED::solve_mip(EDParticle& p)
{
    return (p.solve_mip(constraint_map));
    /*vector<int> y;
    y.resize(p.c.size(), 0);
    IloEnv env;
    IloModel model(env);
    IloNumVarArray var(env);
    IloRangeArray con(env);

    for (int i = 0; i < p.c.size(); i++) {
        var.add(IloBoolVar(env));
    }


    try {
        IloExpr con_exp(env);
        IloRangeArray constraints_to_add(env);
        IloExpr con_exp_temp(env);
        for (map<vector<int>, int>::iterator it = constraint_map.begin(); it != constraint_map.end(); it++) {
            IloExpr con_exp(env);

            //constraint pair <variables, |Cutset Edges|>
            vector<int> constraint_vars = it->first;
            int constraint_bound = it->second;
            for (vector<int>::iterator cs_it = constraint_vars.begin(); cs_it != constraint_vars.end(); cs_it++) {
                con_exp += var[*cs_it];
            }
            IloRange r1(env, 0, con_exp, constraint_bound);
            constraints_to_add.add(r1);
        }

        model.add(constraints_to_add);
        double cost;
        IloExpr obj_exp(env);
        for (int i = 0; i < p.c.size(); i++) {
            if (p.c[i] == 1) {
                cost = 1.01;
            } else if (p.c[i] < 0) {
                cost = 0;
            } else {
                cost = p.c[i];
            }
            obj_exp += (1 - cost) * var[i];
        }
        IloObjective obj_fn = IloMaximize(env, obj_exp);
        model.add(obj_fn);

        IloCplex cplex(model);

        //cplex.exportModel("lpex1.lp");

        // Optimize the problem and obtain solution.
        if (!cplex.solve()) {
            env.error() << "Failed to optimize LP" << endl;
            throw(-1);
        }

        IloNumArray vals(env);
        //populate y and return it

        //cout << "Solution status = " << cplex.getStatus() << endl;
        //cout << "Solution value  = " << cplex.getObjValue() << endl;
        cplex.getValues(vals, var);

        //cout << "Values        = " << vals << endl;

        cout << "vals size is " << vals.getSize() << endl;
        for (int i = 0; i < vals.getSize(); i++) {
            y[i] = vals[i];
            //cout << "i = " << i << " vals[i] = " << vals[i] << endl;
        }

        model.remove(obj_fn);
        cout << "number of constraints are " << cplex.getNrows() << endl;
        model.remove(constraints_to_add);

    } catch (IloException& e) {
        cout << e << endl;
    } catch (...) {
        cout << "Unknown exception caught" << endl;
    }
    env.end();
    return y;
    */
}
