#include "ED.h"
#include "Random.h"
#include "LaPSO.hpp"
#include "VolVolume.hpp"

using namespace LaPSO;

int main(int argc,const char** argv){
    string graph_file = "AS-BA.R-Wax.v100e190.bb";
    string pairs_filename = "AS-BA.R-Wax.v100e190.rpairs.10.1";
    std::cout << "Running LaPSO for disjoint paths problem with "
	      << graph_file  << " & " << pairs_filename << std::endl;
    ED ed(graph_file,pairs_filename);
    const int nnode = num_vertices(ed.getGraph());
    const int nedge = num_edges(ed.getGraph());
    const int ncomm = (int) ed.getComm().size();
    std::cout << "Read data " << nnode << " nodes, "
	      << nedge << " edges, " << ncomm << " ODs\n";
    LaPSO::Problem solver(nedge*ncomm,nedge); // number of variables & (relaxed) constraints
    //-------- set default parameter values
    solver.param.absGap = 0.999; // any two solutions must differ by at least 1
    solver.param.printFreq = 1;
    solver.param.printLevel = 1;
    solver.param.maxIter = 200;
    solver.param.subgradFactor=0.01; // start with a small step size
    solver.param.subgradFmin=0.0001; // allow very small steps
    // override defaults with command line arguments
    solver.param.parse(argc,argv); 
    ed.maxEDsolves = solver.param.maxIter;
    solver.best.ub = ncomm; // every path excluded
    solver.best.lb = 0; // no path left out
    solver.dualUB = 0; // all constraints are <= so lagrange multipliers are <= 0
    Uniform rand;
    if(solver.param.randomSeed == 0) rand.seedTime();
    else rand.seed(solver.param.randomSeed);
    //if(useVol){ const EdgeVec &edges, const vector<Commodity> &comm, const int n, Edge_Int_Map em) 
        EDParticle *p = new EDParticle(ed.graph_edges,ed.getComm(),nnode,ed.EIM);
        for(int j=0;j<solver.dsize;++j)
	  p->dual[j] = -1.0/nnode; //-rand(0,0.5); // random initial point
        VOL_LaPSO_adaptor vol(solver,ed,p);
        vol.prob.parm.lambdainit = solver.param.subgradFactor;
        vol.prob.parm.greentestinvl = 3;
        vol.prob.parm.redtestinvl = 5;
        std::cout << "Using volume algorithm with random initial point\n";
        vol.solve(true);
        solver.best = vol.best;

}
