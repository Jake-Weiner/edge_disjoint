#include "LaPSO.hpp"
#include "ED.h"

class ReducedCostCmp {
    const LaPSO::DblVec &rc;
public:
    ReducedCostCmp(const LaPSO::DblVec &reducedCost) : rc(reducedCost) {}
    bool operator () (const Edge &a,const Edge &b) const
	{ return rc[a.idx] < rc[b.idx]; }
};