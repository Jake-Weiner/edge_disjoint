#include "ArgParser.h"
#include "ED.h"
#include "EDP_MIPsolver.h"
#include "LaPSO.hpp"
#include "MIPsolver_test.h"
#include "Random.h"
#include "VolVolume.hpp"
#include "prep_mip.h"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/regex.hpp>
#include <fstream>
#include <iostream>
#include <string>

using std::cout;
using std::endl;

template <typename T>
void write_iteration_checks(string filename, vector<T>& input_vec);
void initLaPSOParameters(MainArg::ArgParser& ap, LaPSO::Problem& solver);

int main(int argc,const char** argv)
{
    MainArg::ArgParser ap;
    ap.parse(argc, argv);

    cout << "ap.printing is" << ap.printing << endl;
    cout << "ap.RAEM is" << ap.RAEM << endl;
    cout << "iterations checks is " << ap.iteration_checks << endl;
    
    std::cout << "Running " << (ap.useVol ? "Volume" : "LaPSO") << " for disjoint paths problem with "
              << ap.graph_filename << " & " << ap.pairs_filename << std::endl;
    ED ed(ap.graph_filename, ap.pairs_filename, ap.printing, ap.randComm, ap.djikstras_naive, ap.RREM, ap.RAEM);
    //ed.initialise_global_constraints();
    const int nnode = ed.get_nodes();
    const int nedge = ed.get_edges();
    const int ncomm = (int)ed.getComm().size();
    std::cout << "Read data " << nnode << " nodes, "
              << nedge << " edges, " << ncomm << " ODs\n";
    LaPSO::Problem solver(nedge * ncomm, nedge); // number of variables & (relaxed) constraints
    //-------- set default parameter values
    initLaPSOParameters(ap, solver);
    // override defaults with command line arguments

    //TODO - put initLaPSO Parameters as default values in solver.param.parse
    solver.param.parse(argc, argv);
    //encapsulate in function...
    cout << "cpu Time limit is" << solver.param.maxCPU << endl;
    solver.best.lb = 0; // no path left out
    solver.dualUB = 0; // all constraints are <= so lagrange multipliers are <= 0
    Uniform rand;
    if (solver.param.randomSeed == 0)
        rand.seedTime();
    else
        rand.seed(solver.param.randomSeed);

    ed.maxEDsolves = solver.param.maxIter;
    ed.solution_cost = solver.best.ub = ncomm; // every path excluded

    if (ap.useVol) { //const EdgeVec &edges, const vector<Commodity> &comm, const int n, Edge_Int_Map em)
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

        double initial_dual = 0.0;
        for (int i = 0; i < solver.param.nParticles; ++i) {
            EDParticle* p = new EDParticle(ed.graph_edges, ed.getComm(), nnode, ed.EIM);
            for (int j = 0; j < solver.dsize; ++j) {
                initial_dual = -rand(0, ap.initial_dual_max);
                p->dual[j] = initial_dual; // random initial point
            }
            solver.swarm.push_back(p);
        }
        std::cout << "set up solver with " << solver.param.nParticles
                  << " particles"
                  << " using veloctiy factor of " << solver.param.velocityFactor << endl;
        solver.solve(ed);
    }

    // always show the final result
    std::cout << "Best solution missing " << solver.best.ub
              << " paths, lower bound " << solver.best.lb
              << std::endl
              << "CPU time = " << solver.cpuTime()
              << " elapsed = " << solver.wallTime() << " sec"
              << " primal cpu time " << solver.primal_cpu_time
              << " dual cpu time " << solver.dual_cpu_time
              << std::endl;

    // if output_file given
    if (ap.output_filename.compare("") != 0) {
        // if (!ap.parameter_output_file.empty())
        //     ap.output_filename = ap.parameter_output_file;
        try {
            std::ofstream outfile;
            outfile.open(ap.output_filename, std::ios_base::app);
            outfile << ap.graph_filename << "," << ap.pairs_filename << "," << solver.param.nCPU << "," << solver.param.nParticles << "," << solver.param.absGap << "," << solver.param.maxIter << "," << solver.param.perturbFactor << "," << solver.param.subgradFactor
                    << "," << solver.param.subgradFmult << "," << solver.param.subgradFmin << "," << solver.param.velocityFactor
                    << "," << solver.param.globalFactor << ","
                    << solver.cpuTime() << "," << ed.getCommSize() - solver.best.lb << ","
                    << ed.getCommSize() - solver.best.ub << "," << ed.getCommSize() - solver.lb_primal_cpu_time << ","
                    << solver.primal_cpu_time << ","
                    << solver.best_nIter
                    << endl;
            std::cout << "writing to " << ap.output_filename << endl;
            outfile.close();
        } catch (std::ofstream::failure e) {
            std::cerr << "Exception opening/reading/closing output file\n";
        }
    }

    // Output iteration tracking

    if (ap.iteration_checks){
    // dual euclid
    write_iteration_checks(ap.dual_euclid_filename, solver.dual_euclid);
    // perturb euclid
    write_iteration_checks(ap.perturb_euclid_filename, solver.perturb_euclid);
    // best lb
    write_iteration_checks(ap.best_lb_filename, solver.best_lb_tracking);
    // best ub
    write_iteration_checks(ap.best_ub_filename, solver.best_ub_tracking);
    // average lb
    write_iteration_checks(ap.average_lb_filename, solver.average_lb_tracking);
    //viol tracking
    write_iteration_checks(ap.average_viol_filename, solver.average_viol_tracking);
    //path_saved tracking
    write_iteration_checks(ap.average_path_saved_filename, solver.average_path_saved_tracking);
    // average ub
    write_iteration_checks(ap.average_ub_filename, solver.average_ub_tracking);
    //dual_0
    write_iteration_checks(ap.dual_0_filename, solver.dual_0_tracking);
    }
    
    return 0;
}

void initLaPSOParameters(MainArg::ArgParser& ap, LaPSO::Problem& solver)
{
    solver.param.absGap = 0.999; // any two solutions must differ by at least 1
    solver.param.printFreq = 1;
    solver.param.printLevel = 1;
    solver.param.heurFreq = 1;
    solver.param.localSearchFreq = 3;
    solver.param.maxIter = 200;
    solver.param.subgradFactor = 0.01; // start with a small step size
    solver.param.subgradFmin = 0.0001; // allow very small steps

    solver.param.zeroInitial = ap._zeroInitial;
    solver.param.particle_tracking = ap.particle_tracking;
    solver.param.localSearch = ap._localSearch;
    solver.param.iteration_checks = ap.iteration_checks;
    cout << "iteration_checks is " << ap.iteration_checks << endl;
    solver.param.time_limit_checks = ap.time_limit_checks;
    solver.param.convergence_output = ap.convergence_filename;
}


// plots iteration number against quantity recorded
template <typename T>
void write_iteration_checks(string filename, vector<T>& input_vec)
{
    std::ofstream outfile;
    outfile.open(filename, std::ios_base::app);
    if (outfile) {
        for (typename vector<T>::iterator it = input_vec.begin(); it != input_vec.end(); it++) {
            outfile << distance(input_vec.begin(), it) << "," << *(it) << endl;
        }
    }
    else
    {
        cout << "output file- " << filename << " not found for iteration checks" << endl;
    }
    
    outfile.close();
    
}


    // for (int i = 1; i < argc; ++i) {
    //     if (argv[i][0] == '-')
    //         continue;
    //     const size_t arglen = strlen(argv[i]);
    //     if (arglen > 3 && string(&argv[i][arglen - 3]) == ".bb")
    //         graph_file = argv[i];
    //     else if (arglen > 10 && string(argv[i]).find(".rpairs") != std::string::npos)
    //         pairs_filename = argv[i];
    //     else if (arglen > 3 && string(&argv[i][arglen - 3]) == ".csv")
    //         output_filename = argv[i];

    // else if (argv[i][0] == 'L')
    //     useVol = false;
    // //output files
    // else if (string(argv[i]) == "wes") {
    //     write_edges_stats = true;
    //     if (string(argv[i + 1]).find(".csv") != std::string::npos)
    //         edgestats_filename = string(argv[i + 1]);
    // } else if (string(argv[i]) == "wme") {
    //     write_mip_edges = true;
    //     mip_edges_folder = string(argv[i + 1]);
    // } else if (string(argv[i]) == "wre") {
    //     write_reduced_edges = true;
    //     reduced_edges_folder = string(argv[i + 1]);
    // } else if (string(argv[i]) == "wo") {
    //     write_outputs = true;
    //     if (string(argv[i + 1]).find(".csv") != std::string::npos)
    //         output_filename = string(argv[i + 1]);
    // } else if (string(argv[i]) == "bc") {
    //     lb_time_comparisons_filename = string(argv[i + 1]);
    // } else if (string(argv[i]) == "cT") {
    //     convergence_test = true;
    //     char* end;
    //     int cT_int = strtol(argv[i + 1], &end, 10);
    //     if (*end != '\0') {
    //         std::cout << "invalid cT integer.\n";
    //         return 1;
    //     }

    //         int current_index = i;
    //         for (int i = current_index + 2; i < current_index + 2 + cT_int; i++) {
    //             std::cout << argv[i] << endl;
    //             if (string(argv[i]).find("dual") != std::string::npos)
    //                 dual_euclid_filename = string(argv[i]);
    //             if (string(argv[i]).find("perturb") != std::string::npos)
    //                 perturb_euclid_filename = string(argv[i]);
    //             if (string(argv[i]).find("best_lb") != std::string::npos)
    //                 best_lb_filename = string(argv[i]);
    //             if (string(argv[i]).find("best_noncutset_lb") != std::string::npos)
    //                 best_noncutset_lb_filename = string(argv[i]);
    //             if (string(argv[i]).find("best_ub") != std::string::npos)
    //                 best_ub_filename = string(argv[i]);
    //             if (string(argv[i]).find("average_lb") != std::string::npos)
    //                 average_lb_filename = string(argv[i]);
    //             if (string(argv[i]).find("average_noncutset_lb") != std::string::npos)
    //                 average_noncutset_lb_filename = string(argv[i]);
    //             if (string(argv[i]).find("average_viol") != std::string::npos)
    //                 average_viol_filename = string(argv[i]);
    //             if (string(argv[i]).find("average_path_saved") != std::string::npos)
    //                 average_path_saved_filename = string(argv[i]);
    //             if (string(argv[i]).find("average_ub") != std::string::npos)
    //                 average_ub_filename = string(argv[i]);
    //             if (string(argv[i]).find("dual_0") != std::string::npos)
    //                 dual_0_filename = string(argv[i]);
    //         }
    //     }
    // }

    /*
    //Cutset related informatiom
    vector<vector<NodeEdgePair>> node_neighbours = ed.get_node_neighbours();
    solve_EDP_MIP(ed.get_edges(),ed.getComm(),node_neighbours)
    */

 //x variables test
    // if (write_reduced_edges) {
    //     size_t pos = graph_file.find("/graphs"); //find location of word
    //     graph_file.erase(0, pos + 8); //delete everything before /Graphs
    //     string instance = pairs_filename;
    //     pos = pairs_filename.find("rpairs.");
    //     instance.erase(0, pos + 7);
    //     ofstream reduced_edges_outfile;
    //     string mip_edges_filename = reduced_edges_folder + "/" + graph_file + "." + instance;
    //     reduced_edges_outfile.open(mip_edges_filename);
    //     for (int primal_idx = 0; primal_idx < solver.x_total.size(); primal_idx++) {
    //         if (solver.x_total[primal_idx] < 4) {
    //             int edge_idx = ed.edgeIdx(primal_idx);
    //             Edge e = ed.IEM[edge_idx];
    //             int commodity_idx = (primal_idx - edge_idx) / ed.get_edges();
    //             reduced_edges_outfile << e.first << " " << e.second << " " << commodity_idx << endl;
    //         }
    //     }
    // }

    //  lb comparisons write out # NEED TO CHECK WHAT THIS IS?
    // try {
    //     outfile.open(lb_time_comparisons_filename, std::ios_base::app);
    //     for (vector<double>::iterator it = solver.lb_comparisons.begin(); it != solver.lb_comparisons.end();
    //          it++) {
    //         outfile << distance(solver.lb_comparisons.begin(), it) << "," << *(it) << endl;
    //         cout << distance(solver.lb_comparisons.begin(), it) << "," << *(it) << endl;
    //     }
    //     outfile.close();
    // } catch (std::ofstream::failure& writeErr) {
    //     cout << "error writing to lb comparison file" << endl;
    // }

    // if (write_mip_edges) {
    //     size_t pos = graph_file.find("/graphs"); //find location of word
    //     graph_file.erase(0, pos + 8); //delete everything before /Graphs
    //     string instance = pairs_filename;
    //     pos = pairs_filename.find("rpairs.");
    //     instance.erase(0, pos + 7);
    //     int edge_count = 0;

    //     ofstream mip_edges_outfile;
    //     string mip_edges_filename = mip_edges_folder + "/" + graph_file + "." + instance;
    //     mip_edges_outfile.open(mip_edges_filename);
    //     //best feasible soln

    //     for (int i = 0; i < solver.best.best_ub_sol.size(); ++i) {
    //         for (auto& edge_pair : solver.best.best_ub_sol[i]) {
    //             mip_edges_outfile << edge_pair.first << " " << edge_pair.second << " " << i << endl;
    //         }
    //     }
    // }

    // parameter testing suite
    // string parameter_output_file = "";
    // if (ap.particle_param == true) {
    //     parameter_output_file = "/home/jake/PhD/Edge_Disjoint/c++/Outputs/particles_param.csv";
    //     std::cout << "read in P" << endl;
    //     std::cout << "running particle param test with " << solver.param.nParticles << endl;
    // } else if (ap.pert_param == true) {
    //     parameter_output_file = "/home/jake/PhD/Edge_Disjoint/c++/Outputs/pert_param.csv";
    //     std::cout << "read in Pe" << endl;
    //     std::cout << "running pert param test with " << solver.param.perturbFactor << endl;
    // } else if (ap.globalFact_param == true) {
    //     parameter_output_file = "/home/jake/PhD/Edge_Disjoint/c++/Outputs/globalFactor_param.csv";
    //     std::cout << "read in gF" << endl;
    //     std::cout << "running gF param test with " << solver.param.globalFactor << endl;
    // } else if (ap.velocity_param == true) {
    //     parameter_output_file = "/home/jake/PhD/Edge_Disjoint/c++/Outputs/velocity_param.csv";
    //     std::cout << "read in v" << endl;
    //     std::cout << "running v param test with " << solver.param.velocityFactor << endl;
    // } else if (ap.subgrad_param == true) {
    //     parameter_output_file = "/home/jake/PhD/Edge_Disjoint/c++/Outputs/subgrad_param.csv";
    //     std::cout << "read in S" << endl;
    //     std::cout << "running subgrad param test" << endl;
    // } else if (ap.mult_update == true) {
    //     parameter_output_file = "/home/jake/PhD/Edge_Disjoint/c++/Outputs/mult_update.csv";
    //     std::cout << "read in M" << endl;
    //     std::cout << "running mult param test" << endl;
    // } else if (ap.mult_random_update == true) {
    //     parameter_output_file = "/home/jake/PhD/Edge_Disjoint/c++/Outputs/mult_update_randomised.csv";
    //     std::cout << "read in MR" << endl;
    //     std::cout << "running mult_rand param test" << endl;
    // }