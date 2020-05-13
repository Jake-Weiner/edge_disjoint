#include "EDParticle.h"

using namespace LaPSO;

EDParticle::EDParticle(const EdgeVec &edges){
    LaPSO::Particle((int)graph_edges.size(), (int)graph_edges.size())
    ,graph_edges(edges);
}

