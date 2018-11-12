#include "prep_mip.h"
#include "ED.h"

using namespace std;
using namespace LaPSO;

typedef pair<int, int> u_point;
typedef vector<u_point>::iterator points_it;

int sum(DblVec x)
{
    int sum = 0;
    for (int i = 0; i < x.size(); i++) {
        sum += x[i];
    }
    return sum;
}
//does a dominate b
bool point_comp(Particle a, Particle b)
{

    if ((a.lb > b.lb) && (sum(a.viol) <= sum(b.viol))
        || (a.lb >= b.lb) && (sum(a.viol) < sum(b.viol))) {
        return true;
    } else
        return false;
}

vector<Particle*> sort_non_dom(vector<Particle*> swarm)
{
    // need to automatically put in best
    vector<Particle*> non_dom;
    non_dom.push_back(swarm[0]);
    for (auto swarm_it = swarm.begin(); swarm_it != swarm.end(); ++swarm_it) {
        bool add_dom = true;
        for (auto dom_it = non_dom.begin(); dom_it != non_dom.end(); ++dom_it) {

            // if potential point dominates
            if (point_comp(**swarm_it, **dom_it)) {
                non_dom.erase(dom_it);
                continue;
            }

            if (point_comp(**dom_it, **swarm_it)) {
                add_dom = false;
                break;
            }
        }
        if (add_dom) {
            non_dom.push_back(*swarm_it);
        }
    }

    for (auto dom_it = non_dom.begin(); dom_it != non_dom.end(); ++dom_it) {
        cout << (*dom_it)->lb << sum((*dom_it)->viol) << endl;
    }

    return non_dom;
}




    /* 
    test code 

    Particle a = Particle();
    Particle b = Particle();
    Particle c = Particle();
    Particle d = Particle();

    a.lb = 3;
    b.lb = 4;
    c.lb = 2;
    d.lb = 1;

    int temp[4] = {1,0,0,0};
    a.viol.assign(temp,temp+4);
    temp[1] = 1;
    b.viol.assign(temp,temp+4);
    temp[2] = 1;
    temp[3] = 1;
    c.viol.assign(temp,temp+4);
    temp[3] = 0;
    d.viol.assign(temp,temp+4);

    vector<Particle *> swarm;
    swarm.push_back(&a);
    swarm.push_back(&b);
    swarm.push_back(&c);
    swarm.push_back(&d);


    vector<Particle *> non_dom_set = sort_non_dom(swarm);

    */