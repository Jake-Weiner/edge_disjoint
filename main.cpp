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
typedef pair<int, int> Edge;

struct Path{
  int origin;
  int dest;
  vector<Edge> solution;
};  

void initialise_paths(string filename, vector<Path>& paths){
  ifstream input(filename);
  split_vector_type SplitVec; // #2: Search for tokens

  if (input.is_open()){
    while (!input.eof()){
      string line;
      getline(input, line);                                       //read number
      split(SplitVec, line, is_any_of(" \t"), token_compress_on); // SplitVec == { "hello abc","ABC","aBc goodbye" }
      paths.push_back(Path{stoi(SplitVec[0]), stoi(SplitVec[1]), {}});
    }
  }
}

int main(int, char *[])
{
  vector<Path> paths;
  initialise_paths("/home/jake/PhD/Edge_Disjoint/LP/Data/pairs/AS-BA.R-Wax.v100e190.rpairs.10.1",paths);
  ofstream out("/home/jake/PhD/Edge_Disjoint/c++/Lagrangian_Relax/output");
  streambuf *coutbuf = std::cout.rdbuf(); //save old buf
  std::cout.rdbuf(out.rdbuf());                //redirect std::cout to out.txt!
  ifstream input("/home/jake/PhD/Edge_Disjoint/LP/Data/AS-BA.R-Wax.v100e190.bb");
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
      current_edge = Edge(stoi(SplitVec[1]),stoi(SplitVec[0]));
      edges.push_back(current_edge);
      edge_index[current_edge] = edge_number;
      edge_number += 1;
      weights.push_back(1);

    }
  }

  int num_nodes = edges.size();

  graph_t g(edges.data(), edges.data() + edges.size(), weights.data(), num_nodes);
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
*/
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
return EXIT_SUCCESS;
}