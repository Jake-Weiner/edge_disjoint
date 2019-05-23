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

    std::cout << "Running " << (useVol ? "Volume" : "LaPSO") << " for disjoint paths problem with "
              << graph_file << " & " << pairs_filename << std::endl;
    ED ed(graph_file, pairs_filename, printing, randComm, djikstras_naive, repair_edge_removal, repair_add_edge);
    //ed.initialise_global_constraints();
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
    rand.seed(1);
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
        std::cout << "set up solver with " << solver.param.nParticles
                  << " particles"
                  << " using veloctiy factor of " << solver.param.velocityFactor << endl;

        solver.solve(ed);
    }

    //if (printing == true) { // always show the final result
    std::cout << "Best solution missing " << solver.best.ub
              << " paths, lower bound " << solver.best.lb
              << std::endl
              << "CPU time = " << solver.cpuTime()
              << " elapsed = " << solver.wallTime() << " sec"
              << " primal cpu time " << solver.primal_cpu_time
              << " dual cpu time " << solver.dual_cpu_time
              << std::endl;
    //}
    std::ofstream outfile;

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

    if (write_outputs) {
        if (!parameter_output_file.empty())
            output_filename = parameter_output_file;

        try {
            outfile.open(output_filename, std::ios_base::app);
            outfile << graph_file << "," << pairs_filename << "," << solver.param.nCPU << "," << solver.param.nParticles << "," << solver.param.absGap << "," << solver.param.maxIter << "," << solver.param.perturbFactor << "," << solver.param.subgradFactor
                    << "," << solver.param.subgradFmult << "," << solver.param.subgradFmin << "," << solver.param.velocityFactor
                    << "," << solver.param.globalFactor << ","
                    << solver.cpuTime() << "," << ed.getCommSize() - solver.best.lb << ","
                    << ed.getCommSize() - solver.best.ub << "," << ed.getCommSize() - solver.lb_primal_cpu_time << ","
                    << solver.primal_cpu_time << ","
                    << solver.best_nIter
                    << endl;
            std::cout << "writing to " << output_filename << endl;
            outfile.close();
        } catch (std::ofstream::failure e) {
            std::cerr << "Exception opening/reading/closing output file\n";
        }
    }

    if (write_mip_edges) {
        vector<Particle*> swarm_unsorted;
        if (solver.param.nParticles == 1) {
            swarm_unsorted = solver.best_particles_primal_time;
        } else {
            swarm_unsorted = solver.swarm_primal_time;
        }
        vector<Particle*> non_dom_set = sort_non_dom(swarm_unsorted);
        size_t pos = graph_file.find("/Graphs"); //find location of word
        graph_file.erase(0, pos + 8); //delete everything before /Graphs
        string instance = pairs_filename;
        pos = pairs_filename.find("rpairs.");
        instance.erase(0, pos + 7);
        int edge_count = 0;

        ofstream mip_edges_outfile;
        string mip_edges_filename = mip_edges_folder + "/" + graph_file + "." + instance;
        mip_edges_outfile.open(mip_edges_filename);
        // this file is used for the mip solver
        bool solved_optimality = false;
        if (solver.best.lb + solver.param.absGap <= solver.best.ub) {
            solved_optimality = true;
        }
        ofstream mip_edges_map_outfile;
        mip_edges_map_outfile.open(mip_edges_map, std::ios_base::app);
        mip_edges_map_outfile << mip_edges_filename << "," << pairs_filename << ","
                              << solver.param.maxCPU - solver.primal_cpu_time
                              << "," << solved_optimality << endl;

        map<Edge, bool> edges_used;

        //best feasible soln
        for (map<int, EdgeVec>::iterator it = solver.best.best_ub_sol.begin(); it != solver.best.best_ub_sol.end(); it++) {
            for (vector<Edge>::iterator edge_it = solver.best.best_ub_sol[it->first].begin();
                 edge_it != solver.best.best_ub_sol[it->first].end(); edge_it++) {
                if (edges_used.find(*edge_it) != edges_used.end()
                    || edges_used.find(Edge(edge_it->second, edge_it->first)) != edges_used.end()) {
                    continue;
                } else {
                    edges_used[*edge_it] = true;
                    edge_count++;
                    if (write_mip_edges) {
                        mip_edges_outfile << edge_it->first << " " << edge_it->second << " " << it->first << endl;
                    }
                }
            }
        }

        //write out edges from non_dom set
        for (int idx = 0; idx < non_dom_set.size(); ++idx) {
            Problem::ParticleIter p(non_dom_set, idx);
            for (vector<Edge>::iterator sol_it = p->best_lb_sol.begin(); sol_it != p->best_lb_sol.end(); sol_it++) {
                if (edges_used.find(*sol_it) != edges_used.end() || edges_used.find(Edge(sol_it->second, sol_it->first)) != edges_used.end()) {
                    continue;
                } else {
                    edges_used[*sol_it] = true;
                    if (write_mip_edges) {
                        mip_edges_outfile << sol_it->first << " " << sol_it->second << endl;
                    }
                    edge_count++;
                }
            }
        }
        if (write_edges_stats) {
            try {
                outfile.open(edgestats_filename, std::ios_base::app);
                outfile << graph_file << "," << ed.get_edges() << "," << edge_count << endl;
                outfile.close();
            } catch (std::ofstream::failure e) {
                std::cerr << "Exception opening/reading/closing output file\n";
            }
        }
    }

    if (convergence_test) {

        // dual euclid
        std::ofstream outfile;
        outfile.open(dual_euclid_filename, std::ios_base::app);
        for (vector<double>::iterator it = solver.dual_euclid.begin(); it != solver.dual_euclid.end();
             it++) {
            outfile << distance(solver.dual_euclid.begin(), it) << "," << *(it) << endl;
        }
        outfile.close();

        // perturb euclid
        outfile.open(perturb_euclid_filename, std::ios_base::app);
        for (vector<double>::iterator it = solver.perturb_euclid.begin(); it != solver.perturb_euclid.end();
             it++) {
            outfile << distance(solver.perturb_euclid.begin(), it) << "," << *(it) << endl;
        }
        outfile.close();

        // best lb
        outfile.open(best_lb_filename);
        for (vector<double>::iterator it = solver.best_lb_tracking.begin(); it != solver.best_lb_tracking.end();
             it++) {
            outfile << distance(solver.best_lb_tracking.begin(), it) << "," << *(it) << endl;
        }
        outfile.close();

        // best ub
        outfile.open(best_ub_filename);
        for (vector<double>::iterator it = solver.best_ub_tracking.begin(); it != solver.best_ub_tracking.end();
             it++) {
            outfile << distance(solver.best_ub_tracking.begin(), it) << "," << *(it) << endl;
        }
        outfile.close();

        // sum lb
        outfile.open(average_lb_filename);
        for (vector<double>::iterator it = solver.average_lb_tracking.begin(); it != solver.average_lb_tracking.end();
             it++) {
            outfile << distance(solver.average_lb_tracking.begin(), it) << "," << *(it) << endl;
        }
        outfile.close();

        //viol tracking
        outfile.open(average_viol_filename);
        for (vector<double>::iterator it = solver.average_viol_tracking.begin(); it != solver.average_viol_tracking.end();
             it++) {
            outfile << distance(solver.average_viol_tracking.begin(), it) << "," << *(it) << endl;
        }
        outfile.close();

        //path_saved tracking
        outfile.open(average_path_saved_filename);
        for (vector<double>::iterator it = solver.average_path_saved_tracking.begin(); it != solver.average_path_saved_tracking.end();
             it++) {
            outfile << distance(solver.average_path_saved_tracking.begin(), it) << "," << *(it) << endl;
        }
        outfile.close();


        // sum ub
        outfile.open(average_ub_filename);
        for (vector<double>::iterator it = solver.average_ub_tracking.begin(); it != solver.average_ub_tracking.end();
             it++) {
            outfile << distance(solver.average_ub_tracking.begin(), it) << "," << *(it) << endl;
        }
        outfile.close();

         // sum ub
        outfile.open(dual_0_filename);
        for (vector<double>::iterator it = solver.dual_0_tracking.begin(); it != solver.dual_0_tracking.end();
             it++) {
            outfile << distance(solver.dual_0_tracking.begin(), it) << "," << *(it) << endl;
        }
        outfile.close();
    }

    /*
        std::ofstream outfile; 
        outfile.open(LB_tracking_filename, std::ios_base::app);
        for (int i=0; i<solver.param.nParticles; i++){
            for (int j=0; j<solver.lower_bounds_tracking[i].size(); j++){
                outfile << i << ","  // this is the particle num
                << solver.lower_bounds_tracking[i][j] << "," 
                << j // iter num
                << endl;
            }
        }
        outfile.close();

        
        outfile.open(UB_tracking_filename, std::ios_base::app);
        for (int i=0; i<solver.param.nParticles; i++){
            for (int j=0; j<solver.upper_bounds_tracking[i].size(); j++){
                outfile << i << ","  // this is the particle num
                << solver.upper_bounds_tracking[i][j] << "," 
                << j // iter num
                << endl;
            }
        }
        outfile.close();
    

/*
    if (convergence_test == true){
        std::ofstream outfile; 
        outfile.open(convergence_filename, std::ios_base::app);
        //writeoutputs here
        for (int i =0; i < solver.convergence_info.size(); i++){
        outfile << solver.convergence_info[i].euclid_dist << "," << solver.convergence_info[i].nIter << "," 
        << solver.convergence_info[i].best_lb << "," << solver.convergence_info[i].best_ub
        << endl;
        }
        outfile.close();
    }
*/

    return 0;
}
