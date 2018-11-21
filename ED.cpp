#include "Random.h"
#include "VolVolume.hpp"
#include <algorithm>
#include <boost/config.hpp>
#include <fstream>
#include <iostream>
#include <list>
#include <map>

#include <boost/algorithm/string.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/tokenizer.hpp>

#include "ED.h"

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

ED::ED(string graph_filename, string pairs_filename, bool _printing, bool _randComm, bool _fudge_factor)
    : nEDsolves(0)
    , maxEDsolves(10)
{
    ED::populate_commodities(pairs_filename);
    ED::populate_graph(graph_filename);
    solution_edges.resize(commodities.size());
    printing = _printing;
    randComm = _randComm;
    fudge_factor = _fudge_factor;
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
    if (input.is_open()) {
        int line_count = 0;
        while (!input.eof()) {
            string line;
            getline(input, line); //read number
            split(SplitVec, line, is_any_of(" \t"), token_compress_on); // SplitVec == { "hello abc","ABC","aBc goodbye" }
            if (line_count == 0) {
                num_nodes = stoi(SplitVec[0]);
                g = graph_t(num_nodes);
                line_count++;
                continue;
            }
            if (SplitVec.size() >= 3) {
                current_edge = Edge(stoi(SplitVec[0]), stoi(SplitVec[1]));
                if (EIM.find(current_edge) != EIM.end()) {
                    std::cout << "WARNING: " << filename << " contains 2 edges between "
                              << current_edge.first << " & " << current_edge.second
                              << " - ignored\n"; // this may give us
                    continue;
                }
                add_edge(current_edge.first, current_edge.second, edge_number, g);
                graph_edges.push_back(current_edge);
                EIM[current_edge] = edge_number;
                // same edge number in both directions
                EIM[Edge(current_edge.second, current_edge.first)] = edge_number;
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
}

Status ED::reducedCost(const Particle& p, DblVec& redCost)
{
    // update reduced costs for all edges except for the origin->node dummy edge weight
    for (size_t c = 0; c < commodities.size(); ++c) {
        edge_iterator ei, ei_end;
        for (tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
            //  const int e = get(edge_index_t(),g,*ei);
            redCost[primalIdx(ei, c)] = -p.dual[dualIdx(ei)];
    }
    return OK;
}


void ED::update_weightMap(bool fudge_factor, std::map<Edge, double>& added_weight, wM& weightMap,
    double& minWeight, int comm_idx, EDParticle& p)
{
    edge_iterator ei, ei_end;
    const double max_rand = 0;
    double random_val = 0;
    double total_random_vals;
    //std::cout << "decltype(eim) is " << type_name<decltype(eim)>() << '\n';
    if (fudge_factor == true) {
        for (tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {
            // if edge is not used, use shortest path using standard costs
            if (p.x[primalIdx(ei, comm_idx)] != 1) {
                weightMap[*ei] = p.rc[primalIdx(ei, comm_idx)];
                minWeight = std::min(minWeight, weightMap[*ei]);
            }

            // if edges are used, add a random perturbation during the sub-problem solve to bias towards edges not used
            else {
                double min_rand = p.dual.min();
                std::uniform_real_distribution<double> unif(min_rand, max_rand);
                std::default_random_engine re;
                random_val = unif(re);
                weightMap[*ei] = p.rc[primalIdx(ei, comm_idx)] - random_val;
                added_weight[Edge(source(*ei, g), target(*ei, g))] = random_val;
                added_weight[Edge(target(*ei, g), source(*ei, g))] = random_val;
                // - because edge weights are negative of the dual
                minWeight = std::min(minWeight, weightMap[*ei]);
            }
        }
    } else {
        for (tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {
            weightMap[*ei] = p.rc[primalIdx(ei, comm_idx)];
            minWeight = std::min(minWeight, weightMap[*ei]);
        }
    }
}

void ED::update_comm_sol(EDParticle& p, DblVec distances, vector<vertex_descriptor> parents, double& total_paths_cost, int random_index, int comm_idx,
vertex_descriptor& start, vertex_descriptor& end, bool printing)
{
    vertex_descriptor current;

    if (distances[end] >= 1) {
        if (printing)
            std::cout << "\tCost for " << comm_idx << " (" << (p.commodities[random_index]).origin
                      << "," << (p.commodities[random_index]).dest << ") " << distances[end] << std::endl;
        p.ub += 1;
        total_paths_cost += 1;
        p.commodities[random_index].solution_value = 1; // penalty for not including path

    } else {
        p.commodities[random_index].solution_value = 0;
        // Iterate over path and add to primal solution
        std::vector<vertex_descriptor>::iterator it;
        if (printing == true) {
            std::cout << "\tPath for " << comm_idx << " (" << p.commodities[random_index].origin
                      << "," << p.commodities[random_index].dest << ") " << distances[end] << ": ";
        }
        for (current = end; current != start; current = parents[current]) {
            if (printing == true)
                std::cout << " " << current;
            p.x[primalIdx(parents[current], current, comm_idx)] += 1;
            // solution is stored in reverse at here
            p.commodities[comm_idx].solution_edges.push_back(Edge(parents[current], current));
        }
        if (printing == true)
            std::cout << " " << start << std::endl;
        // reversing each time is not really necessary but nice
        reverse(p.commodities[random_index].solution_edges.begin(), p.commodities[random_index].solution_edges.end());
        total_paths_cost += distances[end];
    }
}

Status ED::solveSubproblem(Particle& p_)
{

    ++nEDsolves;
    EDParticle& p(static_cast<EDParticle&>(p_));
    //clear previous primal sol
    p.x = 0; // clear would resize the array
    double total_paths_cost = 0;
    // create/update the edge_weights in the graph using reduced costs
    // solve each commodity pair independantly

    const auto& eim = get(edge_index, g); // map edge_descriptor to index
   
    vector_property_map<double, __typeof__(eim)> weightMap(num_edges(g), eim);
    //std::cout << "decltype(weightMap) is " << type_name<decltype(weightMap)>() << '\n';
    std::vector<vertex_descriptor> parents(num_vertices(g));
    DblVec distances(num_vertices(g)); // floating point distances
    edge_iterator ei, ei_end;
    vertex_descriptor start;
    vertex_descriptor end;
    vertex_descriptor current;
    p.ub = p.lb = 0;

    //solve subproblem selecting commodity iteration order randomly
    if (randComm) {
        // randomise order of commodities
        vector<int> random_indices;
        for (int i = 0; i < p.commodities.size(); i++) {
            random_indices.push_back(i);
        }
        random_shuffle(random_indices.begin(), random_indices.end());
        for (int i = 0; i < random_indices.size(); i++) {
            int random_index = random_indices[i];
            const int comm_idx = p.commodities[random_index].comm_idx;
            double minWeight = 1e99;
            std::map<Edge, double> added_weight;

            //void ED::update_weightMap(bool fudge_factor, std::map<Edge, double>& added_weight, wM& weightMap,
            //  double& minWeight, int comm_idx, EDParticle& p)
            // update the edge weights in the graph
            update_weightMap(fudge_factor, added_weight, weightMap, minWeight, comm_idx, p);

            if (minWeight < 0.0) { // this may happen due to perturbations
                if (minWeight < -1e-2)
                    std::cerr << "WARNING: negative edge weight " << minWeight << std::endl;
                for (tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
                    weightMap[*ei] = std::max(0.0, weightMap[*ei]);
            }
            //solve SP
            p.commodities[random_index].solution_edges.clear();
            start = vertex(p.commodities[random_index].origin, g);
            end = vertex(p.commodities[random_index].dest, g);
            dijkstra_shortest_paths(g, start,
                predecessor_map(parents.data())
                    .distance_map(distances.data())
                    .weight_map(weightMap));

            //update distance vector if fudge_factor technique is used
            if (fudge_factor) {
                for (current = end; current != start; current = parents[current]) {
                    if (added_weight.find(Edge(current, parents[current])) != added_weight.end() || added_weight.find(Edge(parents[current], current)) != added_weight.end()) {
                        distances[end] += added_weight[Edge(current, parents[current])];
                    }
                }
            }
               //update solution for current commodity
            update_comm_sol(p, distances, parents, total_paths_cost, random_index, comm_idx, start, end, printing);
        }
    }
     
        //iterate through commodities in standard order
        else
        {
            int loop_idx = 0;
            for (vector<Commodity>::iterator itr = p.commodities.begin(); itr < p.commodities.end(); ++itr) {
                const int comm_idx = itr->comm_idx;
                double minWeight = 1e99;
                std::map<Edge, double> added_weight;

                // update the edge weights in the graph
                update_weightMap(fudge_factor, added_weight, weightMap, minWeight, comm_idx, p);

                if (minWeight < 0.0) { // this may happen due to perturbations
                    if (minWeight < -1e-2)
                        std::cerr << "WARNING: negative edge weight " << minWeight << std::endl;
                    for (tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
                        weightMap[*ei] = std::max(0.0, weightMap[*ei]);
                }
                itr->solution_edges.clear();
                start = vertex(itr->origin, g);
                end = vertex(itr->dest, g);
                dijkstra_shortest_paths(g, start,
                    predecessor_map(parents.data())
                        .distance_map(distances.data())
                        .weight_map(weightMap));

                //update distances[end] to account for added_weight
                //cout << "distances before = " << distances[end] << endl;
                if (fudge_factor) {
                    for (current = end; current != start; current = parents[current]) {
                        if (added_weight.find(Edge(current, parents[current])) != added_weight.end() || added_weight.find(Edge(parents[current], current)) != added_weight.end()) {
                            distances[end] += added_weight[Edge(current, parents[current])];
                        }
                    }
                }
                update_comm_sol(p, distances, parents, total_paths_cost, loop_idx, comm_idx, start, end, printing);
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
    double max_perturb = 0;
    for (size_t i = 0; i < p.perturb.size(); i++) {
        if (p.perturb[i] > max_perturb) {
            max_perturb = p.perturb[i];
            continue;
        }
        if (-1 * p.perturb[i] > max_perturb)
            max_perturb = -1 * p.perturb[i];
    }
    p.lb = (total_paths_cost + p.dual.sum()) - ((num_nodes - 1) * max_perturb);
    //update particles best local solution
    if (p.lb > p_.best_lb) {
        p_.best_lb = p.lb;
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
                  << " max_perturb=" << max_perturb << std::endl
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

    edge_iterator ei, ei_end;
    IntVec viol(p.viol.size(), 1);
    p.x = 0;
    for (auto it = p.commodities.begin(); it != p.commodities.end(); it++)
        for (const Edge& e : it->solution_edges) {
            p.x[primalIdx(e.first, e.second, it->comm_idx)] += 1;
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
            if (p.x[primalIdx(largest_viol_idx, it->comm_idx)] > 1.0e-6) { // the commodity contains this edge
                current_viol = 0;
                for (size_t i = 0; i < it->solution_edges.size(); i++) // sum up violations this commodity has
                    if (viol[edgeIdx(it->solution_edges[i])] < 0)
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
            p.x[primalIdx(e.first, e.second, largest_com_idx)] -= 1;
            viol[edgeIdx(e)] += 1;
        }
        p.commodities[largest_com_idx].solution_edges.clear();

        // try find a new path for this OD-PAIR
        vertex_descriptor start;
        vertex_descriptor end;
        vertex_descriptor current;

        const auto& eim = get(edge_index, g); // map edge_descriptor to index
        // should be consistent - here we are doing all integer weights
        vector_property_map<int, __typeof__(eim)> weightMap_heur(num_edges(g), eim);
        std::vector<vertex_descriptor> parents(num_vertices(g));
        IntVec distances_heur(num_vertices(g)); // integer distances
        const int num_edges_g = num_edges(g);
        for (tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
            weightMap_heur[*ei] = (viol[dualIdx(ei)] == 1) ? 1 : num_edges_g;

        start = vertex(p.commodities[largest_com_idx].origin, g);
        end = vertex(p.commodities[largest_com_idx].dest, g);

        dijkstra_shortest_paths(g, start, predecessor_map(parents.data()).distance_map(distances_heur.data()).weight_map(weightMap_heur));
        // if no feasible sp exists between orig and dest nodes
        if (distances_heur[end] < num_edges_g) {

            // Iterate over path and add to primal solution
            std::vector<vertex_descriptor>::iterator it;
            for (current = end; current != start; current = parents[current]) {
                p.x[primalIdx(parents[current], current, largest_com_idx)] += 1;
                // solution is stored in reverse at here
                const Edge e(Edge(parents[current], current));
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
