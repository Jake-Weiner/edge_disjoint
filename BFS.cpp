#include "BFS.h"


/*
Algorithm
Place the starting node s on the queue.
If the queue is empty, return failure and stop.
If the first element on the queue is a goal node g, return success and stop Otherwise,
Remove and expand the first element from the queue and place all the children at the end of the queue in any order.
Return to step 2.
*/

using namespace std;

void BFS(
    map<int, map<int, bool>>& node_neighbours, int start, map<int,bool>& cutSet)
    {
        queue<int> Q;
        Q.push(start);
        while (!Q.empty()){
            int head = Q.front();
            Q.pop();
            if(cutSet.find(head) == cutSet.end()){
                cutSet[head] = 1;
                // add in children of head node to queue
                if (node_neighbours.find(head) != node_neighbours.end()) {
                    for (map<int, bool>::iterator it = node_neighbours[head].begin();
                    it != node_neighbours[head].end(); it++) {
                        int child = it->first; 
                        Q.push(child);
                    }
                }
            }
        }
    }