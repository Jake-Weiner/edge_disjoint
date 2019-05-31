#include "ED.h"
#include "LaPSO.hpp"
#include "Random.h"
#include "VolVolume.hpp"
#include "prep_mip.h"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/regex.hpp>
#include <fstream>
#include <iostream>
#include <string>

using namespace LaPSO;

int main(int argc, const char** argv)
{

    bool printing = false;
    string graph_file = "";
    string pairs_filename = "";
    string output_filename = "";
    string mip_edges_folder = "";
    string mip_edges_map = "";
    string edgestats_filename = "";
    string dual_euclid_filename = "";
    string perturb_euclid_filename = "";
    string best_lb_filename = "";
    string best_ub_filename = "";
    string average_lb_filename = "";
    string average_viol_filename = "";
    string average_path_saved_filename = "";
    string average_ub_filename = "";
    string dual_0_filename = "";
    string best_bounds_tracking = "";
    string convergence_filename = "";
    string repair_add_edge = "";
    string repair_edge_removal = "";

    bool useVol = false;
    bool particle_param = false;
    bool pert_param = false;
    bool globalFact_param = false;
    bool velocity_param = false;
    bool subgrad_param = false;
    bool mult_update = false;
    bool mult_random_update = false;
    bool write_edges_stats = false;
    bool randComm = false;
    bool write_outputs = false;
    bool write_mip_edges = false;
    bool djikstras_naive = false;
    bool _zeroInitial = false;
    bool particle_tracking = false;
    bool convergence_test = false;
    bool _localSearch = false;
    bool print_initial_costs = false;

    double set_initial = 0.1;
    for (int i = 1; i < argc; ++i) {
        if (argv[i][0] == '-')
            continue;
        const size_t arglen = strlen(argv[i]);
        if (arglen > 3 && string(&argv[i][arglen - 3]) == ".bb")
            graph_file = argv[i];
        else if (arglen > 10 && string(argv[i]).find(".rpairs") != std::string::npos)
            pairs_filename = argv[i];
        else if (arglen > 3 && string(&argv[i][arglen - 3]) == ".csv")
            output_filename = argv[i];
        else if (argv[i][0] == 'L')
            useVol = false;

        // Parameter controls
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
        else if (string(argv[i]) == "dN") {
            djikstras_naive = true;
        } else if (string(argv[i]) == "zI") {
            _zeroInitial = true;
        } else if (string(argv[i]) == "lS") {
            _localSearch = true;
        } else if (string(argv[i]) == "sI") {
            set_initial = atof(argv[i + 1]);
        } else if (string(argv[i]) == "remove_largest_viol") {
            repair_edge_removal = "largest_viol";
        } else if (string(argv[i]) == "remove_perturb") {
            repair_edge_removal = "perturb";
        } else if (string(argv[i]) == "remove_random") {
            repair_edge_removal = "random";
        }

        else if (string(argv[i]) == "add_pert_0") {
            repair_add_edge = "pert_repair_0";
        } else if (string(argv[i]) == "add_pert_min") {
            repair_add_edge = "pert_repair_min";
        } else if (string(argv[i]) == "add_rc") {
            repair_add_edge = "rc_repair";
        } else if (string(argv[i]) == "add_arb") {
            repair_add_edge = "arb_repair";
        }

        //output files
        else if (string(argv[i]) == "wes") {
            write_edges_stats = true;
            if (string(argv[i + 1]).find(".csv") != std::string::npos)
                edgestats_filename = string(argv[i + 1]);
        } else if (string(argv[i]) == "wme") {
            write_mip_edges = true;
            if (string(argv[i + 1]).find("Reduced_Graph") != std::string::npos)
                mip_edges_folder = string(argv[i + 1]);
            if (string(argv[i + 2]).find("maps") != std::string::npos)
                mip_edges_map = string(argv[i + 2]);
        } else if (string(argv[i]) == "wo") {
            write_outputs = true;
            if (string(argv[i + 1]).find(".csv") != std::string::npos)
                output_filename = string(argv[i + 1]);
        } else if (string(argv[i]) == "cT") {
            convergence_test = true;
            if (string(argv[i + 1]).find(".csv") != std::string::npos)
                dual_euclid_filename = string(argv[i + 1]);
            if (string(argv[i + 2]).find(".csv") != std::string::npos)
                perturb_euclid_filename = string(argv[i + 2]);
            if (string(argv[i + 3]).find(".csv") != std::string::npos)
                best_lb_filename = string(argv[i + 3]);
            if (string(argv[i + 4]).find(".csv") != std::string::npos)
                best_ub_filename = string(argv[i + 4]);
            if (string(argv[i + 5]).find(".csv") != std::string::npos)
                average_lb_filename = string(argv[i + 5]);
            if (string(argv[i+6]).find(".csv") != std::string::npos)
                average_viol_filename = string(argv[i+6]);
            if (string(argv[i+7]).find(".csv") != std::string::npos)
                average_path_saved_filename = string(argv[i+7]);
            if (string(argv[i+8]).find(".csv") != std::string::npos)
                average_ub_filename = string(argv[i+8]);
            if (string(argv[i+9]).find(".csv") != std::string::npos)
                dual_0_filename = string(argv[i+9]);
            
        }
    }

    ED ed(graph_file, pairs_filename, printing, randComm, djikstras_naive, repair_edge_removal, repair_add_edge);
    //ed.initialise_global_constraints();
    const int nnode = ed.get_nodes();
    const int nedge = ed.get_edges();
    const int ncomm = (int)ed.getComm().size();

    LaPSO::Problem solver(nedge * ncomm, nedge); // number of variables & (relaxed) constraints
    //-------- set default parameter values
    solver.param.absGap = 0.999; // any two solutions must differ by at least 1
    solver.param.printFreq = 1;
    solver.param.printLevel = 1;
    solver.param.heurFreq = 1;
    solver.param.localSearchFreq = 3;
    solver.param.maxIter = 200;
    solver.param.subgradFactor = 0.01; // start with a small step size
    solver.param.subgradFmin = 0.0001; // allow very small steps
    // override defaults with command line arguments
    solver.param.parse(argc, argv);
    solver.param.zeroInitial = _zeroInitial;
    solver.param.particle_tracking = particle_tracking;
    solver.param.localSearch = _localSearch;
    solver.param.convergence_test = convergence_test;
    solver.param.convergence_output = convergence_filename;

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
        double max_initial_dual = 0.0;
        double initial_dual = 0.0;
        for (int i = 0; i < solver.param.nParticles; ++i) {
            EDParticle* p = new EDParticle(ed.graph_edges, ed.getComm(), nnode, ed.EIM);
            for (int j = 0; j < solver.dsize; ++j) {
                initial_dual = -rand(0, set_initial);

                p->dual[j] = initial_dual; // random initial point
                if (initial_dual < max_initial_dual) {
                    max_initial_dual = initial_dual;
                }
            }
            solver.swarm.push_back(p);
        }
    

        solver.solve(ed);
    }

    //if (printing == true) { // always show the final result


    string parameter_output_file = "";
    if (particle_param == true) {
        parameter_output_file = "/home/jake/PhD/Edge_Disjoint/c++/Outputs/particles_param.csv";
        std::cout << "read in P" << endl;
        std::cout << "running particle param test with " << solver.param.nParticles << endl;
    } else if (pert_param == true) {
        parameter_output_file = "/home/jake/PhD/Edge_Disjoint/c++/Outputs/pert_param.csv";
        std::cout << "read in Pe" << endl;
        std::cout << "running pert param test with " << solver.param.perturbFactor << endl;
    } else if (globalFact_param == true) {
        parameter_output_file = "/home/jake/PhD/Edge_Disjoint/c++/Outputs/globalFactor_param.csv";
        std::cout << "read in gF" << endl;
        std::cout << "running gF param test with " << solver.param.globalFactor << endl;
    } else if (velocity_param == true) {
        parameter_output_file = "/home/jake/PhD/Edge_Disjoint/c++/Outputs/velocity_param.csv";
        std::cout << "read in v" << endl;
        std::cout << "running v param test with " << solver.param.velocityFactor << endl;
    } else if (subgrad_param == true) {
        parameter_output_file = "/home/jake/PhD/Edge_Disjoint/c++/Outputs/subgrad_param.csv";
        std::cout << "read in S" << endl;
        std::cout << "running subgrad param test" << endl;
    } else if (mult_update == true) {
        parameter_output_file = "/home/jake/PhD/Edge_Disjoint/c++/Outputs/mult_update.csv";
        std::cout << "read in M" << endl;
        std::cout << "running mult param test" << endl;
    } else if (mult_random_update == true) {
        parameter_output_file = "/home/jake/PhD/Edge_Disjoint/c++/Outputs/mult_update_randomised.csv";
        std::cout << "read in MR" << endl;
        std::cout << "running mult_rand param test" << endl;
    }

    cout << solver.best.ub << " " << solver.cpuTime() <<  endl;	
    return 0;
}
