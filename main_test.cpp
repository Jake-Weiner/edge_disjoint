#include "ED.h"
#include "LaPSO.hpp"
#include "Random.h"
#include "VolVolume.hpp"
#include <fstream>
#include <string>

using namespace LaPSO;

int main(int argc, const char** argv)
{
    bool printing = false;
    string graph_file = "";
    string pairs_filename = "";
    bool useVol = true;
    bool particle_param = false;
    bool pert_param = false;
    bool globalFact_param = false;
    bool velocity_param = false;
    for (int i = 1; i < argc; ++i) {
        if (argv[i][0] == '-')
            continue;
        const size_t arglen = strlen(argv[i]);
        if (arglen > 3 && string(&argv[i][arglen - 3]) == ".bb")
            graph_file = argv[i];
        else if (arglen > 10 && string(argv[i]).find(".rpairs") != std::string::npos)
            pairs_filename = argv[i];
        else if (argv[i][0] == 'L')
            useVol = false;
        else if (string(argv[i]) == "P")
            particle_param = true;
        else if (string(argv[i]) == "Pe")
            pert_param = true;
        else if (string(argv[i]) == "gF")
            globalFact_param = true;
        else if (string(argv[i]) == "v")
            velocity_param = true;
    }

    std::cout << "Running " << (useVol ? "Volume" : "LaPSO") << " for disjoint paths problem with "
              << graph_file << " & " << pairs_filename << std::endl;
    ED ed(graph_file, pairs_filename, printing);
    const int nnode = num_vertices(ed.getGraph());
    const int nedge = num_edges(ed.getGraph());
    const int ncomm = (int)ed.getComm().size();
    std::cout << "Read data " << nnode << " nodes, "
              << nedge << " edges, " << ncomm << " ODs\n";
    LaPSO::Problem solver(nedge * ncomm, nedge); // number of variables & (relaxed) constraints
    //-------- set default parameter values
    solver.param.absGap = 0.999; // any two solutions must differ by at least 1
    solver.param.printFreq = 1;
    solver.param.printLevel = 1;
    solver.param.heurFreq = 1;
    solver.param.maxIter = 200;
    solver.param.subgradFactor = 0.01; // start with a small step size
    solver.param.subgradFmin = 0.0001; // allow very small steps
    // override defaults with command line arguments
    solver.param.parse(argc, argv);
    ed.maxEDsolves = solver.param.maxIter;
    ed.solution_cost = solver.best.ub = ncomm; // every path excluded
    solver.best.lb = 0; // no path left out
    solver.dualUB = 0; // all constraints are <= so lagrange multipliers are <= 0

    Uniform rand;
    if (solver.param.randomSeed == 0)
        rand.seedTime();
    else
        rand.seed(solver.param.randomSeed);
    if (useVol) { //const EdgeVec &edges, const vector<Commodity> &comm, const int n, Edge_Int_Map em)
        EDParticle* p = new EDParticle(ed.graph_edges, ed.getComm(), nnode, ed.EIM);
        for (int j = 0; j < solver.dsize; ++j)
            p->dual[j] = -1.0 / nnode; //-rand(0,0.5); // random initial point
        VOL_LaPSO_adaptor vol(solver, ed, p);
        vol.prob.parm.lambdainit = solver.param.subgradFactor;
        vol.prob.parm.greentestinvl = 3;
        vol.prob.parm.redtestinvl = 5;
        std::cout << "Using volume algorithm with random initial point\n";
        vol.solve(true);
        solver.best = vol.best;
    } else {
        for (int i = 0; i < solver.param.nParticles; ++i) {
            EDParticle* p = new EDParticle(ed.graph_edges, ed.getComm(), nnode, ed.EIM);
            for (int j = 0; j < solver.dsize; ++j)
                p->dual[j] = -rand(0, 0.1) / nnode; // random initial point
            solver.swarm.push_back(p);
        }
        std::cout << "set up solver with " << solver.param.nParticles
                  << " particles\n";
        solver.solve(ed);
    }

    if (printing = true) {
        std::cout << "Best solution missing " << solver.best.ub
                  << " paths, lower bound " << solver.best.lb
                  << std::endl
                  << "CPU time = " << solver.cpuTime()
                  << " elapsed = " << solver.wallTime() << " sec"
                  << std::endl;
    }
    std::ofstream outfile;

    string output_file = "";
    if (particle_param == true) {
        output_file = "/home/jake/PhD/Edge_Disjoint/c++/Outputs/particles_param.csv";
        std::cout << "read in P" << endl;
        std::cout << "running particle param test with " << solver.param.nParticles << endl;
    } else if (pert_param == true) {
        output_file = "/home/jake/PhD/Edge_Disjoint/c++/Outputs/pert_param.csv";
        std::cout << "read in Pe" << endl;
        std::cout << "running pert param test with " << solver.param.perturbFactor << endl;
    } else if (globalFact_param == true) {
        output_file = "/home/jake/PhD/Edge_Disjoint/c++/Outputs/globalFactor_param.csv";
        std::cout << "read in gF" << endl;
        std::cout << "running gF param test with " << solver.param.globalFactor << endl;
    }
      else if (velocity_param == true) {
        output_file = "/home/jake/PhD/Edge_Disjoint/c++/Outputs/velocity_param.csv";
        std::cout << "read in v" << endl;
        std::cout << "running v param test with " << solver.param.velocityFactor << endl;
    }

    try {
        outfile.open(output_file, std::ios_base::app);
        outfile << graph_file << "," << pairs_filename << "," << solver.param.nCPU << "," << solver.param.nParticles << "," << solver.param.absGap << "," << solver.param.maxIter << "," << solver.param.perturbFactor << "," << solver.param.subgradFactor
                << "," << solver.param.subgradFmult << "," << solver.param.subgradFmin << "," << solver.param.velocityFactor
                << "," << solver.param.globalFactor << ","
                << solver.cpuTime() << "," << solver.best.lb << ","
                << ed.getCommSize() - solver.best.ub << endl;
        std::cout << "writing to " << output_file << endl;
        outfile.close();
    } catch (std::ofstream::failure e) {
        std::cerr << "Exception opening/reading/closing output file\n";
    }

      vector<Particle *> population = solver.swarm;
    return 0;
}
