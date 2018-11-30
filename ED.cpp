#include "Random.h"
#include "VolVolume.hpp"
#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>
#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include "ED.h"
#include "djikstra.h"



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

ED::ED(string graph_filename, string pairs_filename, bool _printing, bool _randComm)
    : nEDsolves(0)
    , maxEDsolves(10)
{
    ED::populate_commodities(pairs_filename);
    ED::populate_graph(graph_filename);
    solution_edges.resize(commodities.size());
    printing = _printing;
    randComm = _randComm;
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
            if (parents[current] == -1){
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
    p.ub = p.lb = 0;

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
            int random_index = random_indices[i];
            

            //solve SP
            p.commodities[random_index].solution_edges.clear();
            start = p.commodities[random_index].origin;
            end = p.commodities[random_index].dest;
            double SP = djikstras(EIM, node_neighbours, start, end,parents, num_nodes,num_edges, 
            p.rc, p.x, p.commodities[random_index].comm_idx,commodities.size());

            if (SP == -1) {
                cerr << "error with djikstras" << endl;
                exit;
            }
            //update solution for current commodity
            update_comm_sol(p, SP, parents, total_paths_cost, random_index, start, end, printing);
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
            double SP = djikstras(EIM, node_neighbours, start, end,parents, num_nodes,num_edges, 
            p.rc, p.x, p.commodities[loop_idx].comm_idx,commodities.size());

            if (SP == -1) {
                cerr << "error with djikstras" << endl;
                exit;
            }
            //update solution for current commodity
            update_comm_sol(p, SP, parents, total_paths_cost, loop_idx, start, end, printing);
            loop_idx++;
        }
    }

    // edges violated in the ED solution
    // violation = b - Ax = 1 - sum_c x_ic
    p.viol = 1; // sets every violation to 1
    for (size_t i = 0; i < p.x.size(); i++) {
        p.viol[edgeIdx(i)] -= p.x[i]; // sum over all commodites
    }
    // is feasible if no edges are violated
    p.isFeasible = (p.viol.min() >= 0);
    p.lb = (total_paths_cost + p.dual.sum());
    //update particles best local solution
    if (p.lb > p_.best_lb) {
        p_.best_lb = p.lb;
        int sum = 0;
        for (int i = 0; i < p.viol.size(); i++) {
            sum += p.viol[i];
        }
        p_.best_lb_viol = sum;
        p_.best_lb_sol.clear();
        for (vector<Commodity>::iterator itr = p.commodities.begin(); itr < p.commodities.end(); ++itr) {
            for (EdgeIter E = itr->solution_edges.begin(); E != itr->solution_edges.end(); E++) {
                p_.best_lb_sol.push_back(*E);
            }
        }
        if (p_.best_lb_sol.empty()) {
            cout << "best sol is empty" << endl;
        }
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

    IntVec viol(p.viol.size(), 1);
    p.x = 0;
    
    // iterate over solution edges
    for (auto it = p.commodities.begin(); it != p.commodities.end(); it++)
        for (const Edge& e : it->solution_edges) {
            p.x[primalIdx(EIM[e], it->comm_idx)] += 1;
            viol[edgeIdx(e)] -= 1;
        }
    while (viol.min() < 0) {
        if (printing == true)
            cout << "\tviol sum is " << viol.sum() << endl;
        auto result = std::min_element(viol.begin(), viol.end());
        // identify edge_idx with the largest violation
        int largest_viol_idx = distance(viol.begin(), result);
        if (printing == true) {
            cout << "\tmin viol is at index " << largest_viol_idx
                 << " with value " << viol[largest_viol_idx] << endl;
        }
        largest_viol = 0;
        //find commodity with the largest violation, that contains this edge
        for (auto it = p.commodities.begin(); it != p.commodities.end(); it++) {
            if (p.x[primalIdx(largest_viol_idx, it->comm_idx)] > 0) { // the commodity contains this edge
                current_viol = 0;
                for (size_t i = 0; i < it->solution_edges.size(); i++) // sum up violations this commodity has
                    if (viol[edgeIdx(it->solution_edges[i])] < 1.0e-6)
                        current_viol += viol[edgeIdx(it->solution_edges[i])];
                if (current_viol < largest_viol) {
                    largest_com_idx = it->comm_idx;
                    largest_viol = current_viol;
                }
            }
        }

        if (printing == true)
            cout << "\tlargest comm idx " << largest_com_idx << " with " << -largest_viol << endl;
        // remove this path from both x AND from the solution_edges
        for (const Edge& e : p.commodities[largest_com_idx].solution_edges) {
            p.x[primalIdx(edgeIdx(e), largest_com_idx)] -= 1;
            viol[edgeIdx(e)] += 1;
        }
        p.commodities[largest_com_idx].solution_edges.clear();

        // try find a new path for this OD-PAIR


        IntVec distances_heur(num_nodes); // integer distances
        DblVec temp_rc;
        temp_rc.resize(p.rc.size());

        for (int i=0; i<p.rc.size(); i++){
            temp_rc[i] = (viol[edgeIdx(i)] == 1) ? 1 : num_edges;
        }


        vector<int> parents; 
        int start = p.commodities[largest_com_idx].origin;
        int end = p.commodities[largest_com_idx].dest;
        int current;
        double SP = djikstras(EIM, node_neighbours, start, end, parents, num_nodes,num_edges, 
            temp_rc, p.x, p.commodities[largest_com_idx].comm_idx,commodities.size());
        // if no feasible sp exists between orig and dest nodes

        
        if (SP < num_edges) {

            // Iterate over path and add to primal solution
            for (current = end; current != start; current = parents[current]) {
                if (parents[current] == -1){
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
        
    }
    p.ub = 0;

    for (auto it = p.commodities.begin(); it != p.commodities.end(); it++) {
        if (it->solution_edges.empty())
            p.ub += 1;
    }

    p.isFeasible = true;
    if (printing == true)
        cout << "Heuristic found solution with " << p.ub << " missing paths\n";
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
