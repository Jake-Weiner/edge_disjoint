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

ED::ED(string graph_filename, string pairs_filename, bool _printing, bool _randComm, bool _djikstras_naive, MyTypes::repairRemoveEdgeMethod RREM,
    MyTypes::repairAddEdgeMethod RAEM)
    : nEDsolves(0)
    , maxEDsolves(10)
{
    ED::populate_commodities(pairs_filename);
    ED::populate_graph(graph_filename);
    solution_edges.resize(commodities.size());
    solution_edges_nodes.resize(commodities.size());
    printing = _printing;
    randComm = _randComm;
    dN = _djikstras_naive;
    this->RREM = RREM;
    this->RAEM = RAEM;
}

ED::~ED() {}

void ED::populate_commodities(string filename)
{
    int comm_idx = 0;
    ifstream input(filename);
    split_vector_type SplitVec;
    if (input.is_open()) {
        cout << "pairs file found" << endl;
        string line;
        getline(input, line);
        while (!input.eof()) {
            split(SplitVec, line, is_any_of(" \t"), token_compress_on);
            if (SplitVec.size() > 1) {
                cout << line << endl;
                int origin_node = stoi(SplitVec[0]);
                int dest_node = stoi(SplitVec[1]);
                commodities.push_back(Commodity{ origin_node, dest_node, {}, {}, 0.0, comm_idx });
                comm_idx++;
                getline(input, line);
            }
        }
    }
    cout << "finished reading in pairs" << endl;
}

void ED::populate_graph(string filename)
{
    //populate the graph
    ifstream input(filename);
    split_vector_type SplitVec;
    size_t edge_number = 0;
    Edge current_edge;
    const int start_node_idx = 0;
    const int end_node_idx = 1;
    int start_node;
    int end_node;
    if (input.is_open()) {
        int line_count = 0;
        while (!input.eof()) {
            string line;
            getline(input, line); //read number
            // split the line based on \t and whitespace
            split(SplitVec, line, is_any_of(" \t"), token_compress_on);
            // first line contains the number of nodes
            if (line_count == 0) {
                num_nodes = stoi(SplitVec[0]) + 1; // in case the data is indexed from 1..n
                node_neighbours.resize(num_nodes);
                ++line_count;
                continue;
            }
            // second line contains the number of edges
            else if(line_count == 1){
                ++line_count;
                continue;
            }
            // edges are written as Node1 Node2 Weight
            if (SplitVec.size() >= 3) {
                try{
                    start_node = stoi(SplitVec[start_node_idx]);
                    end_node = stoi(SplitVec[end_node_idx]);
                }
                catch(...){
                    cout << "error in converting string to int in populate_graph()" << endl;
                    exit(1);
                }

                // ask andreas about duplicate edge check?
                current_edge = Edge(start_node, end_node);
                if (EIM.find(current_edge) != EIM.end()) {
                    std::cout << "WARNING: " << filename << " contains 2 edges between "
                              << current_edge.first << " & " << current_edge.second
                              << " - duplicates accepted\n"; // this may give us
                }
                graph_edges.push_back(current_edge);
                EIM[current_edge] = edge_number;
                IEM[edge_number] = current_edge;
                // store which nodes are connected to each other
                node_neighbours[start_node].push_back(NodeEdgePair(end_node, edge_number));
                node_neighbours[end_node].push_back(NodeEdgePair(start_node, edge_number));
                ++edge_number;
            }
        }
    }
    // check for 0/1 indexing...
    if (node_neighbours[num_nodes - 1].empty()){
        node_neighbours.resize(--num_nodes);
    }
        

    num_edges = graph_edges.size();

    cout << "finished reading in graph" << endl;
}

Status ED::reducedCost(const Particle& p, DblVec& redCost)
{
    // update reduced costs for all edges except for the origin->node dummy edge weight
    for (size_t c = 0; c < commodities.size(); ++c) {
        for (int e = 0; e < get_edges(); ++e)
            redCost[primalIdx(e, c)] = -p.dual[e];
    }
    return OK;
}

void ED::update_comm_sol(EDParticle& p, double SP, const vector<NodeEdgePair>& parents, double& total_paths_cost, int index,
    int start, int end, bool printing)
{
    if (SP >= 1 || parents[end].first == -1) { // no (short enough) path to end
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
        storePath(p, index, start, end, parents, 0, printing);
        if (printing == true)
            std::cout << " " << start << std::endl;
        total_paths_cost += SP;
    }
}

void ED::storePath(EDParticle& p, int comm, int start, int end, const vector<NodeEdgePair>& parents,
    vector<int>* viol, bool doPrint) const
{
    for (int current = end; current != start; current = parents[current].first) {
        if (doPrint)
            std::cout << " " << current;
        if (parents[current].first == -1) {
            cerr << "issue storingPath - parents array incorrect" << endl;
        }
        const int eidx = parents[current].second;

        p.x[primalIdx(eidx, comm)] += 1;

       
        // solution is stored in reverse at here

        p.commodities[comm].solution_edges.push_back(eidx);

        int parent_node = parents[current].first;

        // ensure correct edge directions for mip solver
        // p.commodities[comm].solution_edges_nodes.

        if (p.commodities[comm].solution_edges_nodes.empty()) {
            p.commodities[comm].solution_edges_nodes.push_back({ parent_node, current });
        }


        // else{
        //     Edge last_edge = p.commodities[comm].solution_edges_nodes.back();

        //     if(current == last_edge.second){
        //         p.commodities[comm].solution_edges_nodes.push_back({current,parent_node});
        //     }
        //     else if(parent_node ==  last_edge.second){
        //         p.commodities[comm].solution_edges_nodes.push_back({parent_node,current});
        //     }
        //     else if(current == last_edge.first){
        //         p.commodities[comm].solution_edges_nodes.push_back({parent_node,current});
        //     }
        //     else{
        //         p.commodities[comm].solution_edges_nodes.push_back({current,parent_node});
        //     }
        // }

        // for repair method, keep track of violations
        if (viol != nullptr)
            (*viol)[eidx] -= 1;


        // else if (current == p.commodities[comm].origin){
        //     p.commodities[comm].solution_edges_nodes.push_back({current,parent_node});
        // }
        // else if (current == p.commodities[comm].dest){
        //     p.commodities[comm].solution_edges_nodes.push_back({parent_node,current});
        // }
        // order does not matter for these edges
        // else {
        //     p.commodities[comm].solution_edges_nodes.push_back({current,parent_node});
        // }

        //cout << parents[current].first << " " << current << endl;
    }
    // reversing each time is not really necessary but nice
    //reverse(p.commodities[comm].solution_edges_nodes.begin(), p.commodities[comm].solution_edges_nodes.end());
}

Status ED::solveSubproblem(Particle& p_)
{

    // updated randomshuffle
    ++nEDsolves;
    EDParticle& p(static_cast<EDParticle&>(p_));
    //clear previous primal sol
    p.x = 0; // clear would resize the array

    // create/update the edge_weights in the graph using reduced costs
    // solve each commodity pair independantly

    double total_paths_cost = 0.0;
    vector<NodeEdgePair> parents;

    int start;
    int end;
    p.ub = p.commodities.size();
    p.lb = 0;
    double max_perturb = 0.0;

    // find max perturbation
    for (size_t i = 0; i < p.perturb.size(); i++) {
        if (p.perturb[i] > max_perturb) {
            max_perturb = p.perturb[i];
        } else if (-1 * p.perturb[i] > max_perturb)
            max_perturb = -1 * p.perturb[i];
    }

    vector<int> commodity_order;
    for (int i = 0; i < p.commodities.size(); i++) {
        commodity_order.push_back(i);
    }
    // shuffle order if random selection selected
    if (randComm) {
        random_shuffle(commodity_order.begin(), commodity_order.end());
    }

    // loop through commodities to solve shortest path problems
    for (int loop_idx = 0; loop_idx < commodity_order.size(); loop_idx++) {
        int comm_index = commodity_order[loop_idx];
        //solve SP
        p.commodities[comm_index].solution_edges.clear();
        p.commodities[comm_index].solution_edges_nodes.clear();
        start = p.commodities[comm_index].origin;
        end = p.commodities[comm_index].dest;
        double SP;

        //ensure no negative edge weights
        for (int i = 0; i < p.rc.size(); i++) {
            if (p.rc[i] < 0) {
                // cout << "edge weight is " << p.rc[i] << endl;
                // cout << "error, negative edge weights detected " << endl;
                // exit(1);
                p.rc[i] = 0;
            }
        }

        if (dN) { // djikstras without violation consideration
            SP = djikstras_naive(EIM, node_neighbours, start, end, parents, num_nodes, num_edges,
                p.rc, p.x, p.commodities[comm_index].comm_idx, commodities.size(), 1.0);
        } else { // Threshold = 1.0 means don't return paths longer than 1.0
            SP = djikstras(max_perturb, EIM, node_neighbours, start, end, parents, num_nodes, num_edges,
                p.rc, p.x, p.commodities[comm_index].comm_idx, commodities.size(), 1.0);
        }
        if (SP == -1) {
            cerr << "error with djikstras" << endl;
            exit(2);
        }

        const Commodity_SP sp = {
            start, end, parents
        };
        //update solution for current commodity
        //shortest path info required to update solutions
        p.commodity_shortest_paths[p.commodities[comm_index].comm_idx] = sp;
        // p.c is cost of path -- used in MIP solver
        p.c[p.commodities[comm_index].comm_idx] = SP;
        //update_comm_sol(p, SP, parents, total_paths_cost, random_index, start, end, printing);
    }

    for (int i = 0; i < p.c.size(); i++) {
        const vector<NodeEdgePair>& parents = p.commodity_shortest_paths[i].parents;
        int start = p.commodity_shortest_paths[i].start;
        int end = p.commodity_shortest_paths[i].end;

        double SP = p.c[i];
        update_comm_sol(p, SP, parents, total_paths_cost, i, start, end, printing);
    }

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
        p_.best_lb_sol.resize(p.commodities.size());
        for (vector<Commodity>::iterator itr = p.commodities.begin(); itr < p.commodities.end(); ++itr) {
            p_.best_lb_sol[itr->comm_idx] = itr->solution_edges_nodes;
        }
    }

    if (p.isFeasible && (p.ub < p_.best_ub)) {
        p_.best_ub_sol.resize(p.commodities.size());
        for (vector<Commodity>::iterator itr = p.commodities.begin(); itr < p.commodities.end(); ++itr) {
            p_.best_ub_sol[itr->comm_idx] = itr->solution_edges_nodes;
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

    IntVec viol(p.viol.size(), 1);
    p.x = 0;
    int viol_sum = 0;
    // iterate over solution edges --> locally set viol intvec and make sure p.x is correct
    for (auto it = p.commodities.begin(); it != p.commodities.end(); it++)
        for (const int e : it->solution_edges) {
            p.x[primalIdx(e, it->comm_idx)] += 1;
            viol[edgeIdx(e)] -= 1;
        }

    for (int i = 0; i < viol.size(); i++) {
        if (viol[i] < 0) {
            viol_sum += 1;
        }
    }
    p.viol_sum = viol_sum;

    // keep removing paths until feasible solution found
    while (viol.min() < 0) {
        if (printing == true)
            cout << "\tviol sum is " << viol.sum() << endl;

        
        int comm_idx_to_remove = -1;

        // remove commodities based on random selection
        if (RREM == MyTypes::random) {
            cout << "using random repair" << endl;
            vector<int> random_indices;
            for (int i = 0; i < p.commodities.size(); i++) {
                random_indices.push_back(i);
            }
            random_shuffle(random_indices.begin(), random_indices.end());
            bool remove_path = false;
            int loop_idx = 0;
            while (remove_path == false) {
                // loop through commodities randomly until a commodity with a violation is found
                int random_index = random_indices[loop_idx];
                if (!(p.commodities[random_index].solution_edges.empty())) {
                    for (const int e : p.commodities[random_index].solution_edges) {
                        //contains an edge which is violated
                        if (viol[e] < 0) {
                            remove_path = true;
                            comm_idx_to_remove = random_index;
                        }
                    }
                }
                loop_idx++;
            }
            // remove paths based on the which paths have the most amount of violations
        } else if (RREM == MyTypes::largest_viol) {
            auto largest_edge_violation = std::min_element(viol.begin(), viol.end());
            // identify edge_idx with the largest violation
            int largest_viol_idx = distance(viol.begin(), largest_edge_violation);
            if (printing == true) {
                cout << "\tmin viol is at index " << largest_viol_idx
                     << " with value " << viol[largest_viol_idx] << endl;
            }

            int largest_viol = 0;
            // identify which commodity uses this edge idx and remove the commodity with the largest number of violations
            for (auto it = p.commodities.begin(); it != p.commodities.end(); it++) {
                if (p.x[primalIdx(largest_viol_idx, it->comm_idx)] > 0) { // the commodity contains this edge
                    current_viol = 0;
                    for (size_t i = 0; i < it->solution_edges.size(); i++) { // sum up violations this commodity has
                        if (viol[edgeIdx(it->solution_edges[i])] < 1.0e-6) {
                            current_viol += viol[edgeIdx(it->solution_edges[i])];
                        }
                        if (current_viol < largest_viol) {
                            comm_idx_to_remove = it->comm_idx;
                            largest_viol = current_viol;
                        }
                    }
                }
            }
            // remove paths based on the perturbation values which the paths have for the edge with the largest viol
        } else if (RREM == MyTypes::perturb) {
            auto largest_edge_violation = std::min_element(viol.begin(), viol.end());
            // identify edge_idx with the largest violation
            int largest_viol_idx = distance(viol.begin(), largest_edge_violation);
            if (printing == true) {
                cout << "\tmin viol is at index " << largest_viol_idx
                     << " with value " << viol[largest_viol_idx] << endl;
            }
            double largest_perturb = -1;
            //remove with highest perturbation value (this commodity is least likely to use this edge)
            for (auto it = p.commodities.begin(); it != p.commodities.end(); it++) {
                if (p.x[primalIdx(largest_viol_idx, it->comm_idx)] > 1.0e-6) { // the commodity contains this edge
                    if (p.perturb[primalIdx(largest_viol_idx, it->comm_idx)] > largest_perturb) {
                        largest_perturb = p.perturb[primalIdx(largest_viol_idx, it->comm_idx)];
                        comm_idx_to_remove = it->comm_idx;
                    }
                }
            }
        }

        // remove commodity chosen
        removeCommodity(p, viol, comm_idx_to_remove);

        // try find a new path for this removed commodity using either perturbations, reduced costs or arbitrary values
        DblVec temp_rc;
        temp_rc.resize(p.rc.size());
        const double min_perturb = std::min(0.0, p.perturb.min());
        const double max_perturb = std::max(0.0, p.perturb.max());
        const double scale = 1.0 / num_nodes * 1.0 / (max_perturb - min_perturb + 1e-5);
        for (int i = 0; i < p.rc.size(); i++) {
            // using perturbations where negative values are shifted to 0
            if (RAEM == MyTypes::pert_repair_0) {
                temp_rc[i] = (viol[edgeIdx(i)] == 1) ? scale * p.perturb[i] : 1;
            } // using perturbations where negative values are up by min perturb value
            else if (RAEM == MyTypes::pert_repair_min) {
                temp_rc[i] = (viol[edgeIdx(i)] == 1) ? scale * (p.perturb[i] - min_perturb) : 1;
            } // using the reduced costs as repair guide
            else if (RAEM == MyTypes::rc_repair) {
                temp_rc[i] = (viol[edgeIdx(i)] == 1) ? std::min(1.0 / num_nodes, p.rc[i] / num_nodes) : 1;
            } // every available edge has the same edge weight
            else if (RAEM == MyTypes::arb_repair) {
                temp_rc[i] = (viol[edgeIdx(i)] == 1) ? 1.0 / num_nodes : 1.0;
            } else {
                cout << "repair method - add edge not set properly" << endl;
            }

            if (temp_rc[i] < 0) { // rounding error
                temp_rc[i] = 0.0;
            }
            //temp_rc[i] = (viol[edgeIdx(i)] == 1) ? p.rc[i] : 1;
        }

        vector<NodeEdgePair> parents;
        int start = p.commodities[comm_idx_to_remove].origin;
        int end = p.commodities[comm_idx_to_remove].dest;

        // // if no feasible sp exists between orig and dest nodes

        const double thresh = 1.0;
        double SP = djikstras_naive(EIM, node_neighbours, start, end, parents, num_nodes, num_edges,
            temp_rc, p.x, p.commodities[comm_idx_to_remove].comm_idx, commodities.size(), thresh);

        if (SP < thresh) {
            storePath(p, comm_idx_to_remove, start, end, parents, &viol);
            if (printing == true)
                cout << "\t\trepaired path" << endl;
        }
    }

    // Final feasibility check
    for (int i = 0; i < p.x.size(); i++) {
        if (p.x[i] > 1) {
            cout << "error, infeasible solution found" << endl;
            exit(1);
        }
    }

    p.isFeasible = true;
    p.ub = 0;

    // missing paths, ub +=1. For true LB to original problem, this should be num_commodities - ub
    for (auto it = p.commodities.begin(); it != p.commodities.end(); it++) {
        if (it->solution_edges.empty())
            p.ub += 1;
    }

    
    if (printing == true)
        cout << "Heuristic found solution with " << p.ub << " missing paths\n";

    if (p.isFeasible && (p.ub < p_.best_ub)) {
        p.best_ub_sol.resize(p.commodities.size());
        for (vector<Commodity>::iterator itr = p.commodities.begin(); itr < p.commodities.end(); ++itr) {
            p.best_ub_sol[itr->comm_idx] = itr->solution_edges_nodes;
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

void ED::removeCommodity(EDParticle& p, IntVec& viol, int comm_idx_to_remove)
{

    for (const int e : p.commodities[comm_idx_to_remove].solution_edges) {
        p.x[primalIdx(e, p.commodities[comm_idx_to_remove].comm_idx)] -= 1;
        viol[e] += 1;
    }
    p.commodities[comm_idx_to_remove].solution_edges.clear();
    p.commodities[comm_idx_to_remove].solution_edges_nodes.clear();
}

void ED::localSearch(Particle& p_)
{

    return;
    /*
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
            for (const int& e : it->solution_edges) {
                p.x[primalIdx(e, it->comm_idx)] += 1;
                viol[e] -= 1;
            }
        }
    }
    //randomise commodity iteration order
    random_shuffle(commodities_to_add.begin(), commodities_to_add.end());
    const double min_perturb = std::min(0.0,p.perturb.min());
	const double max_perturb = std::max(0.0,p.perturb.max());
	const double scale =1.0/num_nodes * 1.0/(max_perturb-min_perturb+1e-5);

    //try and add in commodities 1 at a time using previously_unused edges
    for (auto it = commodities_to_add.begin(); it != commodities_to_add.end(); it++) {
        int comm_idx = it->comm_idx;
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

            if (temp_rc[i] < 0) temp_rc[i] = 0;
        }
        vector<NodeEdgePair> parents;
        int start = it->origin;
        int end = it->dest;
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
            storePath(p, comm_idx, start, end, parents, &viol);
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
        p_.best_ub_sol.resize(p.commodities.size());
        for (vector<Commodity>::iterator itr = p.commodities.begin(); itr < p.commodities.end(); ++itr) {
            p_.best_ub_sol[itr->comm_idx] = itr->solution_edges;
        }
        p_.best_ub = p.ub;
    }
    */
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

vector<int> ED::find_cutset_commodities(const EDParticle& p, const vector<bool>& S_cutSet) const
{
    vector<int> cut_set;
    for (const auto& comm : p.commodities) {
        if (S_cutSet[comm.origin] != S_cutSet[comm.dest])
            cut_set.push_back(comm.comm_idx);
    }
    return cut_set;
}

int ED::find_cutset_edges(const vector<bool>& S_cutSet) const
{

    int cut_set_edges = 0;
    for (const Edge& edge : graph_edges)
        if (S_cutSet[edge.first] != S_cutSet[edge.second])
            ++cut_set_edges;

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

    // K=Set( k for k in keys(od_pairs) if od_pairs[k][source] in vertices && od_pairs[k][sink] in vertices)
    // @variable(mip,x[i=vertices,j=neighbours[i],k=K],Bin, start=0)
    // @variable(mip,z[k=K],Bin, start = 0) # is this commodity used?
    // # for k in keys(od_pairs)
    // # 	setvalue(z[k],0)
    // # for key in neighbours
    // # 	setvalue(x[key,neighbours[key]])

    // open(reduced_edges_file) do reduced_edge_input
    // 	for (index, data) in enumerate(eachline(reduced_edge_input))
    // 		a = split(data)
    // 		first_node = parse(Int16,a[1])
    // 		second_node = parse(Int16,a[2])
    // 		commodity_number = parse(Int16,a[3]) + 1

    // 		setlowerbound(x[first_node,second_node ,commodity_number],0)
    // 		setlowerbound(x[second_node,first_node ,commodity_number],0)
    // 		setupperbound(x[first_node,second_node ,commodity_number],0)
    // 		setupperbound(x[second_node,first_node ,commodity_number],0)
    // 		# x[first_node,second_node ,commodity_number] == 0
    // 		# x[second_node,first_node ,commodity_number] == 0

    // 	end
    // end

    // for obj in warm_start_edges
    // 	setupperbound(x[obj[1],obj[2],obj[3]],1)
    //  	setvalue(x[obj[1],obj[2],obj[3]], 1)
    //  	setvalue(z[obj[3]], 1)
    // end

    // @objective(mip,Max,sum( z[k] for k=K ))

    // @constraint(mip,flow[i=vertices,k=K], sum(x[i,j,k] for j in neighbours[i]) - sum( x[j,i,k] for j in neighbours[i] ) ==
    //                                               ((i==od_pairs[k][source] ? z[k] : 0) + (i==od_pairs[k][sink] ? -z[k] : 0)))
    // # @constraint(mip,test_flow[k=K], )

    // @constraint(mip,cap[i=vertices,j=neighbours[i];i<j], sum(x[i,j,k]+x[j,i,k] for k=K) <= 1)

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

        /*
        if (constraint_map.empty()){
            IloExpr con_exp(env);
            for (int i = 0 ; i<c.size() ; i++){
                con_exp += var[i];
            }
            IloRange r1(env, 0, con_exp, c.size());
            constraints_to_add.add(r1);
        }
        */

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
        cplex.setParam(IloCplex::Threads, 1); // each particle gets 1 thread only
        cplex.setOut(env.getNullStream());
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

        double improvement = 0;
        for (int i = 0; i < vals.getSize(); i++) {
            y[i] = floor(vals[i] + 0.1); // make sure we don't get rounding error
            if (c[i] < 1)
                improvement += (1 - vals[i]) * (1 - c[i]);
            //cout << "i = " << i << " vals[i] = " << vals[i] << endl;
        }

        model.remove(obj_fn);
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
