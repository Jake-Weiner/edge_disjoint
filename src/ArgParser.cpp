#include "ArgParser.h"

namespace MainArg{

bool charToBool(const char* inputString);



    void ArgParser::parse(int argc, const char** argv){
        AnyOption parser(100);

        // Bool Arguments
        parser.setOption("printing");
        parser.setOption("useVol");
        parser.setOption("randComm");
        parser.setOption("write_edges_stats");
        parser.setOption("_zeroInitial");
        parser.setOption("_localSearch");
        parser.setOption("particle_param");
        parser.setOption("pert_param");
        parser.setOption("globalFact_param");
        parser.setOption("velocity_param");
        parser.setOption("subgrad_param");
        parser.setOption("mult_update");
        parser.setOption("mult_random_update");
        parser.setOption("write_outputs");
        parser.setOption("write_mip_edges");
        parser.setOption("write_reduced_edges");
        parser.setOption("djikstras_naive");
        parser.setOption("particle_tracking");
        parser.setOption("convergence_test");
        parser.setOption("iteration_checks");
        parser.setOption("time_limit_checks");
        parser.setOption("print_initial_costs");
        
        // String Arguments
        parser.setOption("graph_filename");
        parser.setOption("pairs_filename");
        parser.setOption("output_filename");
        parser.setOption("edgestats_filename");
        parser.setOption("dual_euclid_filename");
        parser.setOption("perturb_euclid_filename");
        parser.setOption("best_lb_filename");
        parser.setOption("lb_time_comparisons_filename");
        parser.setOption("best_lb_filename");
        parser.setOption("best_ub_filename");
        parser.setOption("average_lb_filename");
        parser.setOption("average_ub_filename");
        parser.setOption("best_bounds_tracking");
        parser.setOption("average_viol_filename");
        parser.setOption("average_path_saved_filename");
        parser.setOption("dual_0_filename");
        parser.setOption("convergence_filename");
        parser.setOption("repair_add_edge_method");
        parser.setOption("repair_edge_removal_method");
        
        // doubles 
        parser.setOption("initial_dual_max");

        if (argc == -1) // abuse of function
            parser.processFile(*argv);
        else
            parser.processCommandArgs(argc, argv);

        // Parameter Parsing
        parseStrings(parser);
        parseRepairMethod(parser);
        parseBools(parser);
        parseDoubles(parser);
    }

    void ArgParser::parseStrings(AnyOption& parser){

        // Input/output files
        if (parser.getValue("graph_filename"))  
            graph_filename = string(parser.getValue("graph_filename"));
        if (parser.getValue("pairs_filename"))  
            pairs_filename = string(parser.getValue("pairs_filename"));
        if (parser.getValue("output_filename"))  
            output_filename = string(parser.getValue("output_filename"));
        if (parser.getValue("edgestats_filename"))  
            edgestats_filename = string(parser.getValue("edgestats_filename"));
        if (parser.getValue("dual_euclid_filename"))  
            dual_euclid_filename = string(parser.getValue("dual_euclid_filename"));
        if (parser.getValue("perturb_euclid_filename"))  
            perturb_euclid_filename = string(parser.getValue("perturb_euclid_filename"));
        if (parser.getValue("best_lb_filename"))  
            best_lb_filename = string(parser.getValue("best_lb_filename"));
        if (parser.getValue("lb_time_comparisons_filename"))  
            lb_time_comparisons_filename = string(parser.getValue("lb_time_comparisons_filename"));
        if (parser.getValue("best_lb_filename"))  
            best_lb_filename = string(parser.getValue("best_lb_filename"));
        if (parser.getValue("best_ub_filename"))  
            best_ub_filename = string(parser.getValue("best_ub_filename"));
        if (parser.getValue("average_lb_filename"))  
            average_lb_filename = string(parser.getValue("average_lb_filename"));
        if (parser.getValue("average_ub_filename"))  
            average_ub_filename = string(parser.getValue("average_ub_filename"));
        if (parser.getValue("best_bounds_tracking"))  
            best_bounds_tracking = string(parser.getValue("best_bounds_tracking"));
        if (parser.getValue("average_viol_filename"))  
            average_viol_filename = string(parser.getValue("average_viol_filename"));
        if (parser.getValue("average_path_saved_filename"))  
            average_path_saved_filename = string(parser.getValue("average_path_saved_filename"));
        if (parser.getValue("dual_0_filename"))  
            dual_0_filename = string(parser.getValue("dual_0_filename"));
        if (parser.getValue("convergence_filename"))  
            convergence_filename = string(parser.getValue("convergence_filename"));
        
    }

    // parse in repair methods

    void ArgParser::parseBools(AnyOption& parser){

        if (parser.getValue("printing"))
                printing  = charToBool(parser.getValue("printing"));
        if (parser.getValue("useVol"))
                useVol  = charToBool(parser.getValue("useVol"));
        if (parser.getValue("particle_param"))
                particle_param = charToBool(parser.getValue("particle_param"));
        if (parser.getValue("pert_param"))
                pert_param = parser.getValue("pert_param");
        if (parser.getValue("globalFact_param"))
                globalFact_param = charToBool(parser.getValue("globalFact_param"));
        if (parser.getValue("velocity_param"))
                velocity_param = charToBool(parser.getValue("velocity_param"));
        if (parser.getValue("subgrad_param"))
                subgrad_param = charToBool(parser.getValue("subgrad_param"));
        if (parser.getValue("mult_update"))
                mult_update = charToBool(parser.getValue("mult_update"));
        if (parser.getValue("mult_random_update"))
                mult_random_update = charToBool(parser.getValue("mult_random_update"));
        if (parser.getValue("write_edges_stats"))
                write_edges_stats = charToBool(parser.getValue("write_edges_stats"));
        if (parser.getValue("randComm"))
                randComm = charToBool(parser.getValue("randComm"));
        if (parser.getValue("write_outputs"))
                write_outputs = charToBool(parser.getValue("write_outputs"));
        if (parser.getValue("write_mip_edges"))
                write_mip_edges = charToBool(parser.getValue("write_mip_edges"));
        if (parser.getValue("write_reduced_edges"))
                write_reduced_edges = charToBool(parser.getValue("write_reduced_edges"));
        if (parser.getValue("djikstras_naive"))
                djikstras_naive = charToBool(parser.getValue("djikstras_naive"));
        if (parser.getValue("_zeroInitial"))
                _zeroInitial = charToBool(parser.getValue("_zeroInitial"));
        if (parser.getValue("particle_tracking"))
                particle_tracking = charToBool(parser.getValue("particle_tracking"));
        if (parser.getValue("iteration_checks"))
            iteration_checks = charToBool(parser.getValue("iteration_checks"));  
        if (parser.getValue("time_limit_checks"))
                time_limit_checks = charToBool(parser.getValue("time_limit_checks"));


        if (parser.getValue("_localSearch"))
                _localSearch = charToBool(parser.getValue("_localSearch"));
        if (parser.getValue("print_initial_costs"))
                print_initial_costs = charToBool(parser.getValue("print_initial_costs"));
        
    }

    void ArgParser::parseDoubles(AnyOption& parser){

        if (parser.getValue("initial_dual_max"))
            initial_dual_max = atof(parser.getValue("initial_dual_max"));
    }


    // parse in the repair method chosen
    void ArgParser::parseRepairMethod(AnyOption& parser){

        // edge addition method
        if (parser.getValue("repair_add_edge_method")){
            if (strcmp(parser.getValue("repair_add_edge_method"), "pert_repair_0") == 0){
                RAEM = MyTypes::pert_repair_0;
                cout << "pert_repair_0 selected" << endl;
            }
            else if(strcmp(parser.getValue("repair_add_edge_method"),"pert_repair_min") == 0){
                RAEM = MyTypes::pert_repair_min;
                cout << "pert_repair_min selected" << endl;
            }
            else if(strcmp(parser.getValue("repair_add_edge_method"), "rc_repair")  == 0){
                RAEM = MyTypes::rc_repair;
            }
            else if(strcmp(parser.getValue("repair_add_edge_method"),"arb_repair") == 0){
                RAEM = MyTypes::arb_repair;
            }
        }

        
        // edge removal method
        if (parser.getValue("repair_edge_removal_method")){
            if (parser.getValue("repair_edge_removal_method") == "largest_viol"){
                RREM = MyTypes::largest_viol;
            }
            else if(parser.getValue("repair_edge_removal_method") == "perturb"){
                RREM = MyTypes::perturb;
            }
            else if(parser.getValue("repair_edge_removal_method") == " random"){
                RREM =  MyTypes::random;
            }
           
        }

    }


    //Parameters

    // Convert char to bool
    bool charToBool(const char* inputString){

        bool return_val = false;
        if (inputString!=nullptr){
            if (strcmp(inputString, "true") == 0){
                return_val = true;
            }
        }

        return return_val;
        
    }


};