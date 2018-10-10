#include "ED.h"

int main(){
    string graph_file = "/home/jake/PhD/Edge_Disjoint/LP/Data/AS-BA.R-Wax.v100e190.bb";
    string pairs_filename = "/home/jake/PhD/Edge_Disjoint/LP/Data/pairs/AS-BA.R-Wax.v100e190.rpairs.10.1";
    ED ed(graph_file,pairs_filename);
    ed.solve_ED();
    vector<Commodity> test = ed.get_commodities();
    int a = 0;
}