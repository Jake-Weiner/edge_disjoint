#include "LaPSO.hpp"
#include "ED.h"

using namespace std;
using namespace LaPSO;

class EDParticle : public LaPSO::Particle {
    public:
        EdParticle(const EdgeVec &graph_edges);
        EdgeVec graph_edges;		/// List of edges sorted by reduced cost
        EdgeVec solution_edges;		/// minimum spanning tree edges
};