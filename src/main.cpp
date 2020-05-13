#include <boost/config.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <map>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>

using namespace boost;
using namespace std;

typedef vector<string> split_vector_type;
typedef adjacency_list<listS, vecS, directedS, no_property, property<edge_weight_t, int>> graph_t;
typedef graph_traits<graph_t>::vertex_descriptor vertex_descriptor;
typedef graph_traits<graph_t>::edge_descriptor edge_descriptor;
typedef pair<int, int> Edge;

struct Commodity{
  int origin;
  int dest;
  vector<Edge> solution_edges;
  float solution_value;
};  



/*how to change weight in the graph, will need for later 
std::pair<edge_descriptor, bool> ed = boost::edge(v1,v2,g);
int weight = get(boost::edge_weight_t(), g, ed.first);
int weightToAdd = 10;
boost::put(boost::edge_weight_t(), g, ed.first, weight+weightToAdd);
*/

//global variables

graph_t g;

vector<Commodity> initialise_commodities(string filename){ 
  vector<Commodity> commodities;
  ifstream input(filename);
  split_vector_type SplitVec; // #2: Search for tokens
  if (input.is_open()){
    while (!input.eof()){
      string line;
      getline(input, line);                                       //read number
      split(SplitVec, line, is_any_of(" \t"), token_compress_on); // SplitVec == { "hello abc","ABC","aBc goodbye" }
      if (SplitVec.size() >1)
        commodities.push_back(Commodity{stoi(SplitVec[0]), stoi(SplitVec[1]), {}});
      else
        break;
    }
  }
  return commodities;
}

void populate_graph(string filename){
    //populate the graph
  ifstream input(filename);
  int number_nodes = 0;
  split_vector_type SplitVec; // #2: Search for tokens
  map<Edge,int> edge_index; 
  vector<Edge> edges;
  vector<int> weights;
  vector< vector<Edge> > Paths;
  int edge_number = 0;
  if (input.is_open()){
    int line_count = 0;
    while (!input.eof())
    {
      string line;
      getline(input, line);                                       //read number
      split(SplitVec, line, is_any_of(" \t"), token_compress_on); // SplitVec == { "hello abc","ABC","aBc goodbye" }
      if (line_count==0){
          number_nodes = stoi(SplitVec[0]);
      }
      line_count++;
      if (SplitVec.size() != 4)
        continue;
      Edge current_edge = Edge(stoi(SplitVec[0]),stoi(SplitVec[1]));
      edges.push_back(current_edge);
      edge_index[current_edge] = edge_number;
      edge_number += 1;
      weights.push_back(1);
      add_edge(current_edge.first,current_edge.second,1,g);
      current_edge = Edge(stoi(SplitVec[1]),stoi(SplitVec[0]));
      edges.push_back(current_edge);
      edge_index[current_edge] = edge_number;
      edge_number += 1;
      weights.push_back(1);
      add_edge(current_edge.first,current_edge.second,1,g)
    }
  }

  g = graph_t(edges.data(), edges.data() + edges.size(), weights.data(), number_nodes);
}

void solve_shortest_path(Commodity& commodity){
  commodity.solution_edges.clear();
  commodity.solution_value = -1;
  property_map<graph_t, edge_weight_t>::type weightmap = get(edge_weight, ::g);
  std::vector<vertex_descriptor> parents(num_vertices(::g));
  std::vector<int> distances(num_vertices(::g));
  vertex_descriptor start = vertex(commodity.origin, ::g);
  vertex_descriptor end = vertex(commodity.dest, ::g);
  dijkstra_shortest_paths(::g, start, predecessor_map(boost::make_iterator_property_map(parents.begin(), 
      get(boost::vertex_index, ::g))).distance_map(boost::make_iterator_property_map(distances.begin(), 
      get(boost::vertex_index, ::g))));
  
  vector<vertex_descriptor> path;
  vertex_descriptor current=end;

  while(current!=start) {
    path.push_back(current);
    current=parents[current];
  }
  path.push_back(start);

  //This prints the path reversed use reverse_iterator and rbegin/rend
  std::vector<vertex_descriptor>::iterator it;
  for (int i=0;i<path.size()-1;i++){
    commodity.solution_edges.push_back(Edge(path[i+1],path[i]));
  }
  commodity.solution_value = distances[commodity.dest];
  /*
  for (it=path.end(); it != path.begin(); --it) {

    std::cout << *it << " ";
  }
  */
}

int main(int, char *[])
{
  //initialise the commodities
  string pairs_filename = "/home/jake/PhD/Edge_Disjoint/LP/Data/pairs/AS-BA.R-Wax.v100e190.rpairs.10.1";
  vector<Commodity> commodities = initialise_commodities(pairs_filename);

  //change the output to a file due to problem with visual studio code debug
  ofstream out("/home/jake/PhD/Edge_Disjoint/c++/Lagrangian_Relax/output");
  streambuf *coutbuf = std::cout.rdbuf(); //save old buf
  std::cout.rdbuf(out.rdbuf());                //redirect std::cout to out.txt!

  //initialise graph
  string graph_file = "/home/jake/PhD/Edge_Disjoint/LP/Data/AS-BA.R-Wax.v100e190.bb";
  populate_graph(graph_file);

  for (vector<Commodity>::iterator itr = commodities.begin(); itr < commodities.end(); itr++){
    solve_shortest_path(*itr);
  }

 int a = 0;
  /*
  property_map<graph_t, edge_weight_t>::type weightmap = get(edge_weight, g);
  std::vector<vertex_descriptor> parents(num_vertices(g));
  std::vector<int> distances(num_vertices(g));
  vertex_descriptor start = vertex(0, g);
  vertex_descriptor end = vertex(20, g);

  dijkstra_shortest_paths(g, start,
                          predecessor_map(boost::make_iterator_property_map(parents.begin(), get(boost::vertex_index, g))).distance_map(boost::make_iterator_property_map(distances.begin(), get(boost::vertex_index, g))));

  /*std::cout << "distances and parents:" << std::endl;
  graph_traits<graph_t>::vertex_iterator vi, vend;
  for (boost::tie(vi, vend) = vertices(g); vi != vend; ++vi)
  {
    std::cout << "distance(" << *vi << ") = " << distances[*vi] << ", ";
    std::cout << "parent(" << *vi << ") = " << parents[*vi] << std::endl;
  }
  std::cout << std::endl;

vector<vertex_descriptor> path;
vertex_descriptor current=end;

while(current!=start) {
    path.push_back(current);
    current=parents[current];
}
path.push_back(start);


//This prints the path reversed use reverse_iterator and rbegin/rend
std::vector< graph_traits< graph_t >::vertex_descriptor >::iterator it;
for (it=path.begin(); it != path.end(); ++it) {

    std::cout << *it << " ";
}
std::cout << endl;
std::cout << distances[20] << endl;
std::cout << std::endl;

*/
return EXIT_SUCCESS;
}