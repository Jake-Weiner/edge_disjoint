#include "LaPSO.hpp"
#include "VolVolume.hpp"
#include "Random.h"
#include <iostream>
#include <fstream>
#include <list>
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
using namespace LaPSO;

class ReducedCostCmp {
    const LaPSO::DblVec &rc;
public:
    ReducedCostCmp(const LaPSO::DblVec &reducedCost) : rc(reducedCost) {}
    bool operator () (const Edge &a,const Edge &b) const
	{ return rc[a.idx] < rc[b.idx]; }
};


class DcmstParticle : public LaPSO::Particle {
public:
    DcmstParticle(const EdgeVec &edges,int nNodes) :
	LaPSO::Particle((int)edges.size(),nNodes), edge(edges) {}
    EdgeVec edge;		/// List of edges sorted by reduced cost
    EdgeVec mst;		/// minimum spanning tree edges
};

class DCMST : public LaPSO::UserHooks {
public:    
    LaPSO::DblVec cost;
    EdgeVec edge;		///< List of edges 
    EdgeVec mst;		///< Minimum spanning tree edges
    int n;			// number of nodes
    const int degree;			// degree constraint
    int nMSTsolves;		// number of times MST was solved
    int maxMSTsolves;		// abort after this many
    DCMST(const char *filename,int degree_) : degree(degree_), nMSTsolves(0) {
	std::ifstream in(filename);
	double d;
	for(n=1; ! in.eof() ;++n){
	    for(int j=0;j < n;++j){
		in >> d;
		if( in.fail() ) return; // failed to read another number
		edge.push_back(Edge(j,n,cost.size()));
		cost.push_back(d);
	    }
	}
    }
    ~DCMST() {}			// nothing to do in destructor
    Status reducedCost(const LaPSO::Particle &p, DblVec &redCost) {
	for(EdgeIter e=edge.begin();e!=edge.end();++e)
	    redCost[e->idx] = cost[e->idx] - p.dual[e->i] - p.dual[e->j];
	return OK;
    }
    /// kruskal's spanning tree implementation
    /// @param p (IN) reduced costs for edges
    /// @param ensureFeas (IN) if true creates a heuristic feasible solution
    /// @param treeEdge (OUT) set of edges in the tree (exactly n)
    void kruskal(DcmstParticle &p,bool ensureFeas=false) {
	++nMSTsolves;
	sort(p.edge.begin(),p.edge.end(),ReducedCostCmp(p.rc));
	std::vector<std::list<int> > tree(n); // nodes in each subtree
	IntVec treeIdx(n); // tree that each node belongs to
	IntVec deg(n,0);
	for(int i=0;i < n;++i){
	    tree[i].push_back(i);
	    treeIdx[i]=i;
	}	
	p.mst.clear();
	for(EdgeIter e=p.edge.begin();e!=p.edge.end();++e){
	    if( treeIdx[e->i] == treeIdx[e->j] ) continue; // part of same tree
	    if(ensureFeas){
		if( deg[e->i] < degree && deg[e->j] < degree)
		    ++deg[e->i], ++deg[e->j]; // update degree counts
		else
		    continue;	// not feasible
	    }
	    p.mst.push_back(*e);
	    int i=treeIdx[e->i],j=treeIdx[e->j];
	    if(tree[i].size()+tree[j].size() >= (size_t)n)
		return; // all done ------- exit function
	    if(tree[i].size() < tree[j].size()) std::swap(i,j);
	    for(std::list<int>::iterator nd=tree[j].begin();
		nd!=tree[j].end(); ++nd)
		treeIdx[*nd] = i;
	    tree[i].splice(tree[i].end(),tree[j]);
	}
	if(ensureFeas)
	    for(int i=0;i < n;++i)
		if(deg[i] > degree)
		    fprintf(stderr,"ERROR: degree of %d violated by %d\n",
			    i,deg[i]-degree);
	
    }
    Status solveSubproblem(Particle &p_) {
	DcmstParticle &p(static_cast<DcmstParticle &>(p_));
	kruskal(p);	// solve minimum spanning tree problem
	// calculate violation & store tree solution
	p.viol = degree;
	p.x = 0;
	p.ub = p.lb = 0;
	for(EdgeIter e=p.mst.begin();e!=p.mst.end();++e){
	    p.viol[e->i] -= 1.0;
	    p.viol[e->j] -= 1.0;
	    p.x[e->idx] = 1;
	    p.ub += cost[e->idx]; // primal solution value
	    p.lb += p.rc[e->idx]; // lagrangean/reduced cost value
	}
        // 	printf("MST: size=%d, cost=%g, maxDegree=%d\n",
        // 	       p.mst.size(),p.ub,degree - (int)p.viol.min());
	p.isFeasible = (p.viol.min() >= 0);
	// add lagrangian constant
	for(int i=0; i<n ;++i)
	    p.lb += degree*p.dual[i];
	// add perturbation constant
	// What is the maximum perturbation that could have been included
	// in any feasible DCMST ?
	// as we have a tree we have at least one edge incident to each node
	// so take the n maximum edges (one per node).
	// Could also try (n-1)*max_edge cost.
	// Alternatively could look for d * max edge cost + d * 2nd max edge
	// cost + ... until we reached n-1 edge costs (but this is harder to do)
	DblVec maxEdgeCost(n,-INF);
	for(EdgeIter e=edge.begin();e!=edge.end();++e){
	    maxEdgeCost[e->i] = std::max(maxEdgeCost[e->i],
					 p.perturb[e->idx]);
	    maxEdgeCost[e->j] = std::max(maxEdgeCost[e->j],
					 p.perturb[e->idx]);
	}
	double pConst1 = 0,pConst2=-1e99;
	for(int i=0;i<n;++i){
	  pConst1 += maxEdgeCost[i];
	  if(maxEdgeCost[i] > pConst2) pConst2 = maxEdgeCost[i];
	}
	pConst2 = (n-1) * pConst2;
	p.lb -= std::min(pConst1,pConst2);
	return (nMSTsolves < maxMSTsolves) ? OK : ABORT;
    }

    Status fixConstraint(const int node,const Particle &p_,SparseVec &feas) {
	const DcmstParticle &p(static_cast<const DcmstParticle &>(p_));
	int d=0;
	feas.clear();
	feas.reserve(n);
	// go through in cost order, the "degree" cheapest edges incident
	// to node should go to 1, the others to zero.
	for(EdgeIter e=p.edge.begin(); e!=p.edge.end();++e){
	    if(e->i == node || e->j == node){
		if( ++d <= degree )
		    feas.add(e->idx,1);
		else
		    feas.add(e->idx,0);
	    }
	}
	return OK;
    }

    Status heuristics(Particle &p_) {
	DcmstParticle &p(static_cast<DcmstParticle &>(p_));
	if(p.isFeasible) return OK; // save 1 kruskal() call
	kruskal(p,true);	// run min.spanning tree only accepting feas.
	p.isFeasible = true;
	p.x = 0;		// copy & evaluate primal solution
	p.ub = 0;
	for(EdgeIter e=p.mst.begin();e!=p.mst.end();++e){
	    p.x[e->idx] = 1;
	    p.ub += cost[e->idx]; // primal solution value
	}
	return (nMSTsolves < maxMSTsolves) ? OK : ABORT;
    }
    Status updateBest(Particle &p_) {
	DcmstParticle &p(static_cast<DcmstParticle &>(p_));
	mst = p.mst;		// store best solution
	return OK;
    }
};
