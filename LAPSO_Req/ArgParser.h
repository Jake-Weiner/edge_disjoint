//Parameters used for experiments

#ifndef __ARGPARSER__
#define __ARGPARSER__


#include <map>
#include <limits>
#include "types.h"
#include <string>
#include <string.h>
#include "anyoption.h"

using std::string;

namespace MainArg {


class ArgParser {

    /// load commandline arguments to set parameter values.
    /// This function will take any argument pair of the form
    /// --<name> <value> (where name is one of the attributes like
    /// maxIter) and set the corresponding parameter to the value

    public:

        void parse(int argc, const char** argv);
        void parseStrings(AnyOption& parser);
        void parseBools(AnyOption& parser);
        void parseDoubles(AnyOption& parser);
        void parseRepairMethod(AnyOption& parser);
        //Parameters

        // Boolean Arguments
        
        bool printing = false;
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
        bool write_reduced_edges = false;
        bool djikstras_naive = false;
        bool _zeroInitial = false;
        bool particle_tracking = false;
        bool iteration_checks = false;
        bool time_limit_checks = false;
        bool _localSearch = false;
        bool print_initial_costs = false;

        // Double
        double initial_dual_max = 0.1;

        // Strings
        string graph_filename = "";
        string pairs_filename = "";
        string output_filename = "";
        string edgestats_filename = "";
        string dual_euclid_filename = "";
        string perturb_euclid_filename = "";
        string lb_time_comparisons_filename = "";
        string best_lb_filename = "";
        string best_ub_filename = "";
        string average_lb_filename = "";
        string average_ub_filename = "";
        string best_bounds_tracking = "";
        string average_viol_filename = "";
        string average_path_saved_filename = "";
        string dual_0_filename = "";
        string convergence_filename = "";
        string repair_add_edge_method = "";
        string repair_edge_removal_method = "";

        MyTypes::repairAddEdgeMethod RAEM;
        MyTypes::repairRemoveEdgeMethod RREM;

        


    private:
       
    };
};
#endif