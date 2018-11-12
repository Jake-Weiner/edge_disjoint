#include <iostream>
#include <list>
#include <vector>
#include <string>
#include "LaPSO.hpp"
#include <fstream>
using namespace std;
using namespace LaPSO;

//does a dominate b
bool particle_comp(Particle a,Particle b);
vector<Particle *> sort_non_dom(vector<Particle *> swarm);