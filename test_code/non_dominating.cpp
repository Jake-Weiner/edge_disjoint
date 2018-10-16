
#include <iostream>
#include <list>
#include <vector>

using namespace std;

typedef pair<int, int> u_point;
typedef vector<u_point>::iterator points_it;

//does a dominate b
bool point_comp(u_point a, u_point b)
{

    if ((a.first > b.first) && (a.second <= b.second)
        || (a.first >= b.first) && (a.second < b.second)) {
        return true;
    } else
        return false;
}

int main()
{

    vector<u_point> points;

    points.push_back(u_point(8, 5));
    points.push_back(u_point(9, 2));
    points.push_back(u_point(11, 3));
    points.push_back(u_point(12, 1));
    points.push_back(u_point(16, 2));

    vector<u_point> dominating_set;

    dominating_set.push_back(points[0]);

    for (points_it pot_point = points.begin(); pot_point != points.end(); ++pot_point) {
        bool add_dom = true;
        for (points_it dom_point = dominating_set.begin(); dom_point != dominating_set.end(); ++dom_point) {

            // if potential point dominates
            if (point_comp(*pot_point, *dom_point)) {
                dominating_set.erase(dom_point);
                continue;
            }

            if (point_comp(*dom_point, *pot_point)) {
                add_dom = false;
                break;
            }
        }
        if (add_dom) {
            dominating_set.push_back(*pot_point);
        }

        for (points_it dom_point = dominating_set.begin(); dom_point != dominating_set.end(); ++dom_point) {
            cout << (*dom_point).first << " " << (*dom_point).second << endl;
        }
    }
    return 0;
}