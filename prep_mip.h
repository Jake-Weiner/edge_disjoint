#include <iostream>
#include <list>
#include <vector>
#include "LaPSO.hpp"

using namespace std;
using namespace LaPSO;

//does a dominate b
bool particle_comp(Particle a,Particle b);
vector<Particle> sort_non_dom(vector<Particle>);