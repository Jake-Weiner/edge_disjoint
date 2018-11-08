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

using namespace boost;
using namespace std;

ED::ED(string graph_filename, string pairs_filename, bool _printing)
    : nEDsolves(0)
    , maxEDsolves(10)
{
    ED::populate_commodities(pairs_filename);
    ED::populate_graph(graph_filename);
    solution_edges.resize(commodities.size());
    printing = _printing;
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

Status ED::solveSubproblem(Particle& p_)
{
    ++nEDsolves;
    EDParticle& p(static_cast<EDParticle&>(p_));
    //clear previous primal sol
    p.x = 0; // clear would resize the array
    double total_paths_cost = 0;
    // create/update the edge_weights in the graph using reduced costs
    //g = graph_t(graph_edges.data(), graph_edges.data() + graph_edges.size(), p.rc.data(), num_nodes);
    // solve each commodity pair independantly
    //int comm_idx = 0;
    vertex_descriptor start;
    vertex_descriptor end;
    vertex_descriptor current;
    const auto& eim = get(edge_index, g); // map edge_descriptor to index
    vector_property_map<double, __typeof__(eim)> weightMap(num_edges(g), eim);
    std::vector<vertex_descriptor> parents(num_vertices(g));
    DblVec distances(num_vertices(g)); // floating point distances

    p.ub = p.lb = 0;
    for (vector<Commodity>::iterator itr = p.commodities.begin(); itr < p.commodities.end(); ++itr) {
        const int comm_idx = itr->comm_idx;
        edge_iterator ei, ei_end;
        double minWeight = 1e99;
        // update the edge weights in the graph -> rc
        for (tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {
            weightMap[*ei] = p.rc[primalIdx(ei, comm_idx)];
            minWeight = std::min(minWeight, weightMap[*ei]);
        }
        if (minWeight < 0.0) { // this should never happen!
            if (minWeight < -1e-4)
                std::cerr << "WARNING: negative edge weight " << minWeight << std::endl;
            for (tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
                weightMap[*ei] = std::max(0.0, weightMap[*ei]);
        }
        itr->solution_edges.clear();
        start = vertex((*itr).origin, g);
        end = vertex((*itr).dest, g);
        dijkstra_shortest_paths(g, start,
            predecessor_map(parents.data())
                .distance_map(distances.data())
                .weight_map(weightMap));

        // if no feasible sp exists between orig and dest nodes
        if (distances[end] >= 1) {
            std::cout << "\tCost for " << comm_idx << " (" << (*itr).origin
                      << "," << (*itr).dest << ") " << distances[end] << std::endl;  
            p.ub += 1;
            total_paths_cost += 1;
            itr->solution_value = 1; // penalty for not including path
            continue;
        }
        itr->solution_value = 0;
        // Iterate over path and add to primal solution
        std::vector<vertex_descriptor>::iterator it;
        if (printing == true) {
            std::cout << "\tPath for " << comm_idx << " (" << (*itr).origin
                      << "," << (*itr).dest << ") " << distances[end] << ": ";
        }
        for (current = end; current != start; current = parents[current]) {
            if (printing == true)
                std::cout << " " << current;
            p.x[primalIdx(parents[current], current, comm_idx)] += 1;
            // solution is stored in reverse at here
            itr->solution_edges.push_back(Edge(parents[current], current));
        }
        if (printing == true)
            std::cout << " " << start << std::endl;
        // reversing each time is not really necessary but nice
        reverse(itr->solution_edges.begin(), itr->solution_edges.end());
        total_paths_cost += distances[end];
        //commodity.solution_value = distances[commodity.dest];
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
    if (printing == true) {
        std::cout << "Subproblem solve " << nEDsolves << "/" << maxEDsolves << ": "
                  << " lb=" << p.lb << " ub=" << p.ub
                  << " max_perturb=" << max_perturb << std::endl
                  << "\trange of dual = " << p.dual.min() << " to " << p.dual.max() << std::endl
                  << "\trange of viol = " << p.viol.min() << " to " << p.viol.max() << std::endl;
	if(! p.perturb.empty() )
	  std::cout << "\trange of perturb = " << p.perturb.min() << " to " << p.perturb.max() << std::endl;
        std::cout << "\tpath edge counts:";
        for (auto c = p.commodities.begin(); c != p.commodities.end(); ++c)
            std::cout << " " << c->solution_edges.size();
        std::cout << std::endl;
    }
    if (p.lb > p.max_lb) 
        p.max_lb = p.lb;

    if (p.isFeasible && (p.commodities.size() - p.ub) > p.max_ub) {
        p.max_ub = p.commodities.size() - p.ub;
    }

    return (nEDsolves < maxEDsolves) ? OK : ABORT;
}

Status ED::fixConstraint(const int constraint,
    const Particle& p,
    SparseVec& feas) { return OK; }

Status ED::heuristics(Particle& p_)
{

    EDParticle& p(static_cast<EDParticle&>(p_));
    if (p.isFeasible){
      if(printing) std::cout << "Heuristic called with feasible solution "
			     << p.ub << std::endl;      
      return OK; // already feasible
    }
    int largest_com_idx = 0;
    int largest_viol = 0;
    int current_viol = 0;

    edge_iterator ei, ei_end;
    DblVec viol(p.viol.size(), 1.0);
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
            if (p.x[primalIdx(largest_viol_idx, it->comm_idx)] != 0) { // the commodity contains this edge
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
        for (tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {
            viol[dualIdx(ei)] += p.x[primalIdx(ei, largest_com_idx)];
            p.x[primalIdx(ei, largest_com_idx)] = 0;
        }
        p.commodities[largest_com_idx].solution_edges.clear();

        // try find a new path for this OD-PAIR
        vertex_descriptor start;
        vertex_descriptor end;
        vertex_descriptor current;

        const auto& eim = get(edge_index, g); // map edge_descriptor to index
        vector_property_map<double, __typeof__(eim)> weightMap_heur(num_edges(g), eim);
        std::vector<vertex_descriptor> parents(num_vertices(g));
        IntVec distances_heur(num_vertices(g)); // integer distances
        int num_edges_g = num_edges(g);
        for (tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {
            if (viol[dualIdx(ei)] == 1) {
                weightMap_heur[*ei] = 1;
            } else {
                weightMap_heur[*ei] = num_edges_g;
            }
        }

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
    if(printing == true)
      cout<< "################## " << solution_cost  << " ############\n";
    return OK;
}

/* Stuff for emacs/xemacs to conform to the "Visual Studio formatting standard"
 * Local Variables:
 * tab-width: 4
 * eval: (c-set-style "stroustrup")
 * End:
*/
