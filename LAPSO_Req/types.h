#ifndef __TYPES__
#define __TYPES__

#include <vector>

namespace MyTypes{

// repair methods availabe - TODO Include descriptions of each
enum repairAddEdgeMethod{
  pert_repair_0, pert_repair_min, rc_repair, arb_repair
};

enum repairRemoveEdgeMethod{
  largest_viol, perturb, random
};

};

#endif
