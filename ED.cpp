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

ED::ED(string graph_filename, string pairs_filename){
    ED::populate_graph(graph_filename);
    ED::populate_commodities(pairs_filename);
}

ED::~ED(){}

void ED::solve_ED(){
    if (ED::commodities.empty()) {
        cout << "there are no commodity pairs" << endl;
        return;
    }
    for (vector<Commodity>::iterator itr = commodities.begin(); itr < commodities.end(); itr++) {
        solve_shortest_path(*itr);
    }
}

vector<Commodity> ED::get_commodities(){
    return commodities;
}

void ED::populate_graph(string filename){
    //populate the graph
    ifstream input(filename);
    int number_nodes = 0;
    split_vector_type SplitVec; // #2: Search for tokens
    map<Edge, int> edge_index;
    vector<Edge> edges;
    vector<int> weights;
    vector<vector<Edge>> Paths;
    int edge_number = 0;
    if (input.is_open()) {
        int line_count = 0;
        while (!input.eof()) {
            string line;
            getline(input, line); //read number
            split(SplitVec, line, is_any_of(" \t"), token_compress_on); // SplitVec == { "hello abc","ABC","aBc goodbye" }
            if (line_count == 0) {
                number_nodes = stoi(SplitVec[0]);
            }
            line_count++;
            if (SplitVec.size() != 4)
                continue;
            Edge current_edge = Edge(stoi(SplitVec[0]), stoi(SplitVec[1]));
            edges.push_back(current_edge);
            edge_index[current_edge] = edge_number;
            edge_number += 1;
            weights.push_back(1);
            current_edge = Edge(stoi(SplitVec[1]), stoi(SplitVec[0]));
            edges.push_back(current_edge);
            edge_index[current_edge] = edge_number;
            edge_number += 1;
            weights.push_back(1);
        }
    }
    g = graph_t(edges.data(), edges.data() + edges.size(), weights.data(), number_nodes);
}

Status ED::reducedCost(const Particle &p, DblVec &redCost){}
Status ED::solveSubproblem(Particle& p){}
Status ED::fixConstraint(const int constraint,
									 const Particle &p,
									 SparseVec &feas){}
Status ED::heuristics(Particle &p){}
Status ED::updateBest(Particle &p){}

void ED::populate_commodities(string filename){
    ifstream input(filename);
    split_vector_type SplitVec; // #2: Search for tokens
    if (input.is_open()) {
        while (!input.eof()) {
            string line;
            getline(input, line); //read number
            split(SplitVec, line, is_any_of(" \t"), token_compress_on); // SplitVec == { "hello abc","ABC","aBc goodbye" }
            if (SplitVec.size() > 1)
                commodities.push_back(Commodity{ stoi(SplitVec[0]), stoi(SplitVec[1]), {} });
            else
                break;
        }
    }
}

void ED::solve_shortest_path(Commodity& commodity){
    commodity.solution_edges.clear();
    commodity.solution_value = -1;
    property_map<graph_t, edge_weight_t>::type weightmap = get(edge_weight, g);
    std::vector<vertex_descriptor> parents(num_vertices(g));
    std::vector<int> distances(num_vertices(g));
    vertex_descriptor start = vertex(commodity.origin, g);
    vertex_descriptor end = vertex(commodity.dest, g);
    dijkstra_shortest_paths(g, start, predecessor_map(boost::make_iterator_property_map(parents.begin(), 
            get(boost::vertex_index, g))).distance_map(boost::make_iterator_property_map(distances.begin(), 
            get(boost::vertex_index, g))));
    vector<vertex_descriptor> path;
    vertex_descriptor current = end;
    while (current != start) {
        path.push_back(current);
        current = parents[current];
    }
    path.push_back(start);

    //This prints the path reversed use reverse_iterator and rbegin/rend
    std::vector<vertex_descriptor>::iterator it;
    for (int i = 0; i < path.size() - 1; i++) {
        commodity.solution_edges.push_back(Edge(path[i + 1], path[i]));
    }
    commodity.solution_value = distances[commodity.dest];
}
