#include "Random.h"
#include "VolVolume.hpp"
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

ED::ED(string graph_filename, string pairs_filename)
{
    ED::populate_commodities(pairs_filename);
    ED::populate_graph(graph_filename);
}

ED::~ED() {}

void ED::populate_commodities(string filename)
{
    int idx = 0;
    ifstream input(filename);
    split_vector_type SplitVec; // #2: Search for tokens
    if (input.is_open()) {
        while (!input.eof()) {
            string line;
            getline(input, line); //read number
            split(SplitVec, line, is_any_of(" \t"), token_compress_on); // SplitVec == { "hello abc","ABC","aBc goodbye" }
            if (SplitVec.size() > 1) {
                commodities.push_back(Commodity{ stoi(SplitVec[0]), stoi(SplitVec[1]), {}, idx });
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
    split_vector_type SplitVec; // #2: Search for tokens
    int edge_number = 0;
    double weight_init = 0;
    Edge current_edge;
    if (input.is_open()) {
        int line_count = 0;
        while (!input.eof()) {
            string line;
            getline(input, line); //read number
            split(SplitVec, line, is_any_of(" \t"), token_compress_on); // SplitVec == { "hello abc","ABC","aBc goodbye" }
            if (line_count == 0) {
                num_nodes = stoi(SplitVec[0]);
            }
            line_count++;
            if (SplitVec.size() != 4)
                continue;
            current_edge = Edge(stoi(SplitVec[0]), stoi(SplitVec[1]));
            graph_edges.push_back(current_edge);
            EIM[current_edge] = edge_number;
            edge_number += 1;
            weights.push_back(weight_init);
            current_edge = Edge(stoi(SplitVec[1]), stoi(SplitVec[0]));
            graph_edges.push_back(current_edge);
            EIM[current_edge] = edge_number;
            edge_number += 1;
            weights.push_back(weight_init);
        }
    }
    // add the dummy edge between orig and dest nodes
    for (vector<Commodity>::iterator itr = commodities.begin(); itr < commodities.end(); itr++) {
        current_edge = Edge((*itr).origin, (*itr).dest);
        EIM[current_edge] = edge_number;
        edge_number += 1;
        weights.push_back(1);
    }
}

Status ED::reducedCost(const Particle& p, DblVec& redCost)
{
    int edge_count = 0;
    // update reduced costs for all edges except for the origin->node dummy edge weight
    for (EdgeIter e = graph_edges.begin(); e != graph_edges.end(); ++e)
        if (edge_count < graph_edges.size() - commodities.size()) {
            redCost[EIM[*e]] = weights[EIM[*e]] - p.dual[EIM[*e]];
        }
    return OK;
}

Status ED::solveSubproblem(Particle& p_)
{
    ++maxEDsolves;
    EDParticle& p(static_cast<EDParticle&>(p_));
    //clear previous primal sol
    p.x.clear();
    g.clear();
    double total_paths_cost = 0;
    // create/update the edge_weights in the graph using reduced costs
    g = graph_t(graph_edges.data(), graph_edges.data() + graph_edges.size(), p.rc.data(), num_nodes);
    // solve each commodity pair independantly
    int comm_idx = 0;
    vertex_descriptor start;
    vertex_descriptor end;
    vertex_descriptor current;
    vector<vertex_descriptor> path;
    property_map<graph_t, edge_weight_t>::type weightmap = get(edge_weight, g);
    std::vector<vertex_descriptor> parents(num_vertices(g));
    std::vector<int> distances(num_vertices(g));

    p.ub = p.lb = 0;

    for (vector<Commodity>::iterator itr = p.commodities.begin(); itr < p.commodities.end(); itr++) {
        path.clear();
        start = vertex((*itr).origin, g);
        end = vertex((*itr).dest, g);
        dijkstra_shortest_paths(g, start, predecessor_map(make_iterator_property_map(parents.begin(), get(vertex_index, g))).distance_map(make_iterator_property_map(distances.begin(), get(vertex_index, g))));

        // if no feasible sp exists between orig and dest nodes
        if (distances[end] == 1) {
            p.ub += 1;
            comm_idx++;
            total_paths_cost += 1;
            continue;
        }
        current = end;
        while (current != start) {
            path.push_back(current);
            current = parents[current];
        }
        path.push_back(start);
        //This prints the path reversed use reverse_iterator and rbegin/rend
        std::vector<vertex_descriptor>::iterator it;
        for (int i = 0; i < path.size() - 1; i++) {
            //solution_edges[comm_idx].push_back(Edge(path[i + 1], path[i]));
            p.x[EIM[Edge(path[i + 1], path[i])]] += 1;
        }
        total_paths_cost += distances[end];
        //commodity.solution_value = distances[commodity.dest];
        comm_idx++;
    }
    // edges violated in the ED solution
    for (int i = 0; i < p.x.size(); i++) {
        p.viol[i] = p.x[i] - 1;
    }
    // is feasible if no edges are violated
    p.isFeasible = (p.viol.max() <= 0);
    double max_perturb = 0;
    for (int i = 0; i < p.perturb.size(); i++) {
        if (p.perturb[i] > max_perturb) {
            max_perturb = p.perturb[i];
            continue;
        }
        if (-1 * p.perturb[i] > max_perturb)
            max_perturb = -1 * p.perturb[i];
    }
    p.lb = (total_paths_cost - p.dual.sum()) - ((num_nodes - 1) * max_perturb);
    return (nEDsolves < maxEDsolves) ? OK : ABORT;
}

Status ED::fixConstraint(const int constraint,
    const Particle& p,
    SparseVec& feas) {}
Status ED::heuristics(Particle& p) {}
Status ED::updateBest(Particle& p) {}
