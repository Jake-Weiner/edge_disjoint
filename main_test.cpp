#include "ED.h"
#include "LaPSO.hpp"
#include "Random.h"
#include "VolVolume.hpp"
#include "prep_mip.h"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
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
    bool subgrad_param = false;
    bool mult_update = false;
    bool mult_random_update = false;
    bool write_out_edges = false;
    bool randComm = false;

    bool test_rand_fudge = false;

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
        else if (string(argv[i]) == "S")
            subgrad_param = true;
        else if (string(argv[i]) == "M")
            mult_update = true;
        else if (string(argv[i]) == "MR")
            mult_random_update = true;
        else if (string(argv[i]) == "rC")
            randComm = true;
        else if (string(argv[i]) == "trf") {
            test_rand_fudge = true;
        } else if (string(argv[i]) == "e") {
            write_out_edges = true;
        }
    }

    std::cout << "Running " << (useVol ? "Volume" : "LaPSO") << " for disjoint paths problem with "
              << graph_file << " & " << pairs_filename << std::endl;
    ED ed(graph_file, pairs_filename, printing, randComm);
    const int nnode = ed.get_nodes();
    const int nedge = ed.get_edges();
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
    printing = solver.param.printLevel > 0;
    ed.setPrinting(printing);
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

    //if (printing == true) { // always show the final result
    std::cout << "Best solution missing " << solver.best.ub
              << " paths, lower bound " << solver.best.lb
              << std::endl
              << "CPU time = " << solver.cpuTime()
              << " elapsed = " << solver.wallTime() << " sec"
              << std::endl;
    //}
    std::ofstream outfile;

    string output_file = "out.csv";
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
    } else if (velocity_param == true) {
        output_file = "/home/jake/PhD/Edge_Disjoint/c++/Outputs/velocity_param.csv";
        std::cout << "read in v" << endl;
        std::cout << "running v param test with " << solver.param.velocityFactor << endl;
    } else if (subgrad_param == true) {
        output_file = "/home/jake/PhD/Edge_Disjoint/c++/Outputs/subgrad_param.csv";
        std::cout << "read in S" << endl;
        std::cout << "running subgrad param test" << endl;
    } else if (mult_update == true) {
        output_file = "/home/jake/PhD/Edge_Disjoint/c++/Outputs/mult_update.csv";
        std::cout << "read in M" << endl;
        std::cout << "running mult param test" << endl;
    } else if (mult_random_update == true) {
        output_file = "/home/jake/PhD/Edge_Disjoint/c++/Outputs/mult_update_randomised.csv";
        std::cout << "read in MR" << endl;
        std::cout << "running mult_rand param test" << endl;
    } 

    /*
    output_file = "/home/jake/PhD/Edge_Disjoint/c++/Outputs/new_djikstras.csv";

    try {
        outfile.open(output_file, std::ios_base::app);
        outfile << graph_file << "," << pairs_filename << "," << solver.param.nCPU << "," << solver.param.nParticles << "," << solver.param.absGap << "," << solver.param.maxIter << "," << solver.param.perturbFactor << "," << solver.param.subgradFactor
                << "," << solver.param.subgradFmult << "," << solver.param.subgradFmin << "," << solver.param.velocityFactor
                << "," << solver.param.globalFactor << ","
                << solver.cpuTime() << "," << solver.best.lb << ","
                << ed.getCommSize() - solver.best.ub << "," << ed.getrandComm() << endl;
        std::cout << "writing to " << output_file << endl;
        outfile.close();
    } catch (std::ofstream::failure e) {
        std::cerr << "Exception opening/reading/closing output file\n";
    }
    */
    vector<Particle *> non_dom_set = sort_non_dom(solver.swarm);

    
    if (write_out_edges) {
        map<Edge, bool> edges_used;
        int edge_count = 0;

        //write best feasible / ub solution for MIP solver
        for (Edge_Int_Map::iterator it = ed.EIM.begin(); it != ed.EIM.end(); ++it) {
            for (int comm = 0; comm < ed.getCommSize(); comm++) {
                // contains edge idx it->second
                if (solver.best.x[ed.primalIdx(it->second, comm)] == 1) {
                    // if edge is already accounted for
                    if (edges_used.find(it->first) != edges_used.end()
                        || edges_used.find(Edge((it->first).second, (it->first).first)) != edges_used.end()) {
                        continue;
                    } else {
                        edges_used[it->first] = true;
                        edge_count++;
                        cout << (it->first).first << " " << (it->first).second << " " << comm << endl;
                    }
                }
            }
        }
        //cout << "total edges used is " << edge_count << endl;
        //iterate through swarm
        for (int idx = 0; idx < non_dom_set.size(); ++idx) {
            Problem::ParticleIter p(non_dom_set, idx);


            if (p->best_lb_sol.empty()) {
                //cout << "empty" << endl;
            }
            for (vector<Edge>::iterator sol_it = p->best_lb_sol.begin(); sol_it != p->best_lb_sol.end(); sol_it++) {
                if (edges_used.find(*sol_it) != edges_used.end() || edges_used.find(Edge(sol_it->second, sol_it->first)) != edges_used.end()) {
                    continue;
                } else {
                    edges_used[*sol_it] = true;
                    cout << sol_it->first << " " << sol_it->second << endl;
                    edge_count++;
                }
            }
        }

        //cout << "best edges are" << endl;
        // include edges from best feasible solution

        //cout << "non_zero =" << not_zero << endl;
        //cout << "total edges used is " << edge_count << endl;
        /*
        string edge_filename = "/home/jake/PhD/Edge_Disjoint/c++/Outputs/edge_reduction2.csv";
        try {
            outfile.open(edge_filename, std::ios_base::app);
            outfile << graph_file << "," << ed.get_edges() << "," << edge_count << endl;
            outfile.close();
        } catch (std::ofstream::failure e) {
            std::cerr << "Exception opening/reading/closing output file\n";
        }
        */
        
    }

 
           
           //sort particles into non-dominated set

            //ed.write_mip(solver.swarm, solver.best.lb, ed.getCommSize() - solver.best.ub, outfile_name);

    return 0;
}
    
