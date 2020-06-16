/// Generic solver for combinatorial optimisation problems using a combination
/// of Wedelin's lagrangian relaxation heuristic and particle swarm
/// optimisation. This header file follows the general pattern of the Volume
/// Algorithm implementation in the COIN-OR library to make it easy to swap
/// algorithms and make comparisons between the two methods.
///
/// The current implementation only allows for binary variables.
///
/// Autor: Andreas Ernst   $ Date: January 2010 $  $ Revision: 0.0 $
#include "LaPSO.hpp"
#include "Random.h"
#include "anyoption.h"
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <stdlib.h>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#else
#include <time.h>
#define omp_get_wtime() (double)time(NULL)
#define omp_set_num_threads(n) // do nothing
#endif

namespace LaPSO {
lagged_fibonacci1279 _randGen;
void Param::parse(int argc, const char** argv)
{
    AnyOption parser(100);
    parser.setFlag("help");
    parser.setFlag("randomSeed");
    parser.setOption("subgradFactor");
    parser.setOption("subgradFmult");
    parser.setOption("subgradFmin");
    parser.setOption("perturbFactor");
    parser.setOption("velocityFactor");
    parser.setOption("globalFactor");
    parser.setOption("maxIter");
    parser.setOption("maxCPU");
    parser.setOption("maxWallTime");
    parser.setOption("nParticles");
    parser.setOption("printLevel");
    parser.setOption("printFreq");
    parser.setOption("heurFreq");
    parser.setOption("eps");
    parser.setOption("absGap");
    parser.setOption("relGap");
    parser.setOption("nCPU");
    if (argc == -1) // abuse of function
        parser.processFile(*argv);
    else
        parser.processCommandArgs(argc, argv);
    if (parser.getFlag("help"))
        parser.printAutoUsage();
    if (parser.getFlag("randomSeed"))
        randomSeed = 0; //
    if (parser.getValue("subgradFactor"))
        subgradFactor = atof(parser.getValue("subgradFactor"));
    if (parser.getValue("subgradFmult"))
        subgradFmult = atof(parser.getValue("subgradFmult"));
    if (parser.getValue("subgradFmin"))
        subgradFmin = atof(parser.getValue("subgradFmin"));
    if (parser.getValue("perturbFactor"))
        perturbFactor = atof(parser.getValue("perturbFactor"));
    if (parser.getValue("velocityFactor"))
        velocityFactor = atof(parser.getValue("velocityFactor"));
    if (parser.getValue("globalFactor"))
        globalFactor = atof(parser.getValue("globalFactor"));
    if (parser.getValue("maxIter"))
        maxIter = atoi(parser.getValue("maxIter"));
    if (parser.getValue("maxCPU"))
        maxCPU = atof(parser.getValue("maxCPU"));
    if (parser.getValue("maxWallTime"))
        maxWallTime = atof(parser.getValue("maxWallTime"));
    if (parser.getValue("nParticles"))
        nParticles = atoi(parser.getValue("nParticles"));
    if (parser.getValue("printLevel"))
        printLevel = atoi(parser.getValue("printLevel"));
    if (parser.getValue("printFreq"))
        printFreq = atoi(parser.getValue("printFreq"));
    if (parser.getValue("heurFreq"))
        heurFreq = atoi(parser.getValue("heurFreq"));
    if (parser.getValue("eps"))
        eps = atof(parser.getValue("eps"));
    if (parser.getValue("absGap"))
        absGap = atof(parser.getValue("absGap"));
    if (parser.getValue("relGap"))
        relGap = atof(parser.getValue("relGap"));
    if (parser.getValue("nCPU"))
        nCPU = std::max(1, atoi(parser.getValue("nCPU")));
}
void Param::parse(const char* filename) { parse(-1, &filename); }
void Particle::resize(size_t numVar, size_t numConstr)
{
    dual.resize(numConstr, 0.0);
    perturb.resize(numVar, 0.0);
    x.resize(numVar, 0);
    viol.resize(numConstr, 0.0);
    rc.resize(numVar, 0.0);
    isFeasible = false;
    dVel.resize(numConstr, 0.0);
    ub = INF;
    lb = -INF;
    pVel.resize(numVar, 0.0);
}
// set UB for all lagrangian multipliers to be 0
Problem::Problem(int nVar, int nConstr)
    : psize(nVar)
    , dsize(nConstr)
    , dualLB(dsize, -INF)
    , dualUB(dsize, 0)
    , nIter(-1)
{
    _wallTime = omp_get_wtime();
    best.ub = INF;
    best.lb = -INF;
    best.isFeasible = false;
    best.perturb.resize(nVar, 0.0);
    best.dual.resize(nConstr, 0.0);
}

void Problem::iteration_updates(const int& commodities){
     double sum_lb = 0;
    int sum_ub = 0;
    int sum_viol = 0;
    int path_saved = 0;
    for (int idx = 0; idx < param.nParticles; ++idx) {
        ParticleIter p(swarm, idx);
        sum_lb += p->lb;
        sum_ub += p->ub;
        sum_viol += p->viol_sum;
        path_saved += p->path_saved;
    }
    average_lb_tracking.push_back(commodities - (sum_lb / param.nParticles));
    average_viol_tracking.push_back(sum_viol / param.nParticles);
    average_path_saved_tracking.push_back(path_saved / param.nParticles);
    average_ub_tracking.push_back(commodities - (sum_ub / param.nParticles));
    dual_euclid.push_back(euclideanDistance(swarm, "dual"));
    perturb_euclid.push_back(euclideanDistance(swarm, "perturb"));
    best_lb_tracking.push_back(commodities - best.lb);
    best_ub_tracking.push_back(commodities - best.ub);
}

// Main method
void Problem::solve(UserHooks& hooks)
{

    // commodities size
    double commodities = best.perturb.size() / best.dual.size();

    int non_improv = 0;
    int ub_non_improv = 0;
    double improv_thresh = 1.01;
    const int lbCheckFreq = 10;
    initialise(hooks); // set up the swarm
    //------------- initial solution for swarm --------------
    if (param.printLevel > 1)
        printf("Initialisation:\n");
    DblVec bestLB(param.nParticles);
    IntVec bestIter(param.nParticles, 0);
    std::vector<Uniform> rand(param.nParticles); //

    if (param.randomSeed != 0)
        rand[0].seed(param.randomSeed);
    else
        rand[0].seedTime();

    Status status = OK;
#pragma omp parallel for schedule(static, std::max(1, param.nParticles / param.nCPU))
    for (int idx = 0; idx < param.nParticles; ++idx) {
        ParticleIter p(swarm, idx);
        if (p.idx > 0) { // create random see from previous generator
            rand[p.idx].seed((uint32_t)rand[p.idx - 1](
                0, std::numeric_limits<uint32_t>::max()));
        }
        for (int i = 0; i < dsize; ++i) { // check initial point is valid
            if (p->dual[i] < dualLB[i]) {
                if (param.printLevel)
                    printf("ERROR: initial dual[%d] %.2f below LB %.2f\n",
                        i, p->dual[i], dualLB[i]);
                p->dual[i] = dualLB[i];
            }
            if (p->dual[i] > dualUB[i]) {
                if (param.printLevel)
                    printf("ERROR: initial dual[%d] %.2f below UB %.2f\n",
                        i, p->dual[i], dualUB[i]);
                p->dual[i] = dualUB[i];
            }
        }

        if (hooks.reducedCost(*p, p->rc) == ABORT) {
            status = ABORT;
        }
        /*
        for (int i = 0; i < dsize; ++i) { // check initial point is vali
            printf("dual value is %f\n",p->dual[i]);

        }*/
        p->rc += p->perturb;
        if (hooks.solveSubproblem(*p) == ABORT) {
            status = ABORT;
        }
        /*
         printf("solved\n");
        for (int i = 0; i < dsize; ++i) { // check initial point is valid
 
            printf("dual value is %f\n",p->dual[i]);
        }
        */
        bestLB[p.idx] = p->lb;
        if (param.printLevel > 1) {
            printf("\tp%02d: LB=%g UB=%g feas=%d viol=%g -- %g\n",
                p.idx, p->lb, p->ub, p->isFeasible,
                p->viol.min(), p->viol.max());
            if (param.printLevel > 2) {
                printf("\t\tRedCst %g - %g\n", p->rc.min(), p->rc.max());
                printf("\t\tdVel %g - %g\n", p->dVel.min(), p->dVel.max());
                printf("\t\tpVel %g - %g\n", p->pVel.min(), p->pVel.max());
                printf("\t\tdual %g - %g\n", p->dual.min(), p->dual.max());
                printf("\t\tpert %g - %g\n",
                    p->perturb.min(), p->perturb.max());
            }
        }
        p->pVel = 0;
        p->dVel = 0;
    }
    best.x.resize(psize, 0);
    updateBest(hooks, 0);
    if (status == ABORT)
        return;
    DblVec subgradFactor(param.nParticles, param.subgradFactor);
    DblVec perturbFactor(param.nParticles, param.perturbFactor);
    for (int i = 0; i < param.nParticles; ++i)
        perturbFactor[i] = i * param.perturbFactor / param.nParticles;
    if (param.printLevel > 1)
        printf("Initial LB=%g UB=%g\n", best.lb, best.ub);
    int noImproveIter = 0, nReset = 0, maxNoImprove = 1000;

    // lb time limits for CPLEX/Volume algorithm comparisons.
    vector<double> lb_time_limits = { 90.0, 180.0, 270.0, 360.0, 450.0, 540.0, 630.0, 720.0, 810.0, 900.0 };
    int current_lb_time_limits_idx = 0;
    double last_lb_comparison;

    for (nIter = 1; nIter < param.maxIter && cpuTime() < param.maxCPU && wallTime() < param.maxWallTime && status == OK && (best.lb + param.absGap <= best.ub) && fabs(best.ub - best.lb / best.ub) > param.relGap;
         ++nIter) {

        printf("iter num = %d\n", nIter);
        printf("best upper bound so far is %f \n", commodities - best.lb);
        // cpu time is within the limit
        // if()
        if (param.time_limit_checks == true) {
            if (cpuTime() < lb_time_limits[current_lb_time_limits_idx]) {
                last_lb_comparison = commodities - best.lb;
            } else { // update next limit, use the last found best solution that was within the limit
                printf("Current Limit is %f with bound of %f \n", lb_time_limits[current_lb_time_limits_idx], last_lb_comparison);

                current_lb_time_limits_idx++;
                lb_comparisons.push_back(last_lb_comparison);
                if (current_lb_time_limits_idx > lb_time_limits.size() - 1) {
                    break;
                }
            }
        }

        if (param.iteration_checks == true) {
            iteration_updates(commodities);
        }

        if (param.printLevel > 1 && nIter % param.printFreq == 0)
            printf("Iteration %d --------------------\n", nIter);
#pragma omp parallel for schedule(static, std::max(1, param.nParticles / param.nCPU))
        for (int idx = 0; idx < param.nParticles; ++idx) {
            ParticleIter p(swarm, idx);

            // --------- calculate step size/direction ----------
            double norm = 0;
            for (int i = 0; i < dsize; ++i)
                norm += p->viol[i] * p->viol[i];
            DblVec perturbDir(psize, 0.0);
            perturbationDirection(hooks, *p, perturbDir);

            const double stepSize = rand[p.idx](
                                        subgradFactor[p.idx] / 2.0, subgradFactor[p.idx])
                * (best.ub - bestLB[p.idx]) / norm;
            // 		const double stepSize = randAscent * param.subgradFactor *
            // 					(best.ub-best.lb)  / norm;
            //-------------- update velocity (step) --------------
            double randGlobal = rand[p.idx]() * param.globalFactor;

            for (int i = 0; i < dsize; ++i)
                p->dVel[i] = param.velocityFactor * p->dVel[i] + stepSize * p->viol[i] + randGlobal * (best.dual[i] - p->dual[i]);

            for (int i = 0; i < psize; ++i) {
                double randAscent = rand[p.idx](); // 0,1 uniform
                randGlobal = rand[p.idx]() * param.globalFactor;

                p->pVel[i] = param.velocityFactor * p->pVel[i] + randAscent * perturbFactor[p.idx] * perturbDir[i]
                    + perturbFactor[p.idx] * randGlobal * (1 - 2 * best.x[i]) + randGlobal * (best.perturb[i] - p->perturb[i]);
            }
            //---------- make a step ------------------------------
            for (int i = 0; i < dsize; ++i) {
                p->dual[i] += p->dVel[i];
                p->dual[i] = std::max(dualLB[i], std::min(dualUB[i], p->dual[i]));
            }

            p->perturb *= 0.5;
            p->perturb += p->pVel; // add step in perturbation velocity
            //---------- solve subproblems -------------------------
            if (hooks.reducedCost(*p, p->rc) == ABORT) {
                status = ABORT;
                continue;
            }

            if (nIter % lbCheckFreq != 0) {
                p->rc += p->perturb; // if we care more about ub than lb
                if (hooks.solveSubproblem(*p) == ABORT)
                    status = ABORT;
            } else { // temporarily set pertubation to zero
                DblVec zero(p->perturb.size(), 0.0);
                std::swap(p->perturb, zero);
                if (hooks.solveSubproblem(*p) == ABORT)
                    status = ABORT;
                std::swap(p->perturb, zero); // swap back
            }

            if (status == ABORT)
                continue;
            //if( p->isFeasible ){
            //    // make sure our next step reduces perturbation
            //    p->pVel = p->perturb;
            //    p->pVel *= -0.5/param.velocityFactor;
            //}
            //if( p->lb < bestLB[p.idx]) // downplay current velocity
            //    p->dVel *= 0.5;
            if (nIter % lbCheckFreq == 0) {
                if (p->lb > bestLB[p.idx]) {
                    bestLB[p.idx] = p->lb;
                    bestIter[p.idx] = nIter;
                    //subgradFactor[p.idx] *= 1.1;
                    // 			printf("\tsubgradFactor[%d] increased to %.4f (bestLB = %.2f)\n",
                    // 			       p.idx,subgradFactor[p.idx],bestLB[p.idx]);
                } else { //if( (nIter - bestIter[p.idx]) % 2 == 0 ){
                    subgradFactor[p.idx] = // slow down
                        std::max(param.subgradFmin,
                            param.subgradFmult * subgradFactor[p.idx]);
                    // 			printf("\tsubgradFactor[%d] decreased to %.4f\n",
                    // 			       p.idx,subgradFactor[p.idx]);
                }
            }
            if (param.heurFreq > 0 && nIter % param.heurFreq == 0
                && !p->isFeasible) { // only run heuristic if not alreay feas.
                if (hooks.heuristics(*p) == ABORT) {
                    status = ABORT;
                    continue;
                } else {
                    if (((p->ub <= p->localSearch_thresh) || nIter % param.localSearchFreq == 0) && param.localSearch) {
                        p->localSearch_thresh = p->ub;
                        hooks.localSearch(*p);
                    }
                }
            }
            if (param.printLevel > 1 && nIter % param.printFreq == 0) {
                printf("\tp%02d: LB=%g UB=%g feas=%d minViol=%g\n",
                    p.idx, p->lb, p->ub, p->isFeasible,
                    p->viol.min());
                if (param.printLevel > 2) {
                    printf("\t\tstepSize=%g randGlobal=%g\n",
                        stepSize, randGlobal);
                    printf("\t\tRedCst %g - %g\n", p->rc.min(), p->rc.max());
                    printf("\t\tdVel %g - %g\n", p->dVel.min(), p->dVel.max());
                    printf("\t\tpDir %g - %g\n", perturbDir.min(), perturbDir.max());
                    printf("\t\tpVel %g - %g\n", p->pVel.min(), p->pVel.max());
                    printf("\t\tdual %g - %g\n", p->dual.min(), p->dual.max());
                    printf("\t\tpert %g - %g\n",
                        p->perturb.min(), p->perturb.max());
                }
            }
            if (nIter == 1 && param.nParticles == 1) {
                best_particles.push_back(&(*p)); //initialise best_particles_sols if only 1 particle is used
            }

            // write out particle info to see how it behaves after each iteration
        }
        for (int idx = 0; idx < param.nParticles; ++idx) {
            ParticleIter p(swarm, idx);
            for (int x_idx = 0; x_idx < p->x.size(); ++x_idx) {
                x_total[x_idx] += p->x[x_idx];
            }
        }

        if (updateBest(hooks, nIter)) {
            noImproveIter = 0;
        } else {
            if (++noImproveIter > maxNoImprove) { // reset duals
                ++nReset;
                swarm[0]->dual = 0;
                double minLB = INF, maxLB = -INF;
                for (int idx = 0; idx < param.nParticles; ++idx) {
                    if (swarm[idx]->lb < minLB)
                        minLB = swarm[idx]->lb;
                    if (swarm[idx]->lb > maxLB)
                        maxLB = swarm[idx]->lb;
                }
                hooks.reducedCost(*swarm[0], swarm[0]->rc);
                double maxdual = std::max(0.0, swarm[0]->rc.max()) - std::min(0.0, swarm[0]->rc.min());
                if (maxdual < 1e-9)
                    maxdual = 1;
#pragma omp parallel for schedule(static, std::max(1, param.nParticles / param.nCPU))
                for (int idx = 0; idx < param.nParticles; ++idx) {
                    ParticleIter p(swarm, idx);
                    for (int j = 0; j < dsize; ++j) { // generate pertubation about
                        if (p->dual[j] == 0 && rand[idx]() < 2.0 / dsize) {
                            p->dual[j] = rand[idx](
                                std::max(dualLB[j], -maxdual),
                                std::min(dualUB[j], maxdual));
                        } else
                            p->dual[j] = best.dual[j] * rand[idx](1 - 0.2 / nReset, 1 + 0.2 / nReset);
                        p->dual[j] = std::max(
                            dualLB[j], std::min(dualUB[j], p->dual[j]));
                    }
                    p->perturb = 0.0;
                    p->dVel = 0.0;
                    p->pVel = 0.0;

                    for (int i = 0; i < psize; ++i) // discourage current best
                        if (rand[idx](0, 1) < 0.1) // try to get diversity
                            p->pVel[i] = rand[idx](0, 1) * perturbFactor[idx]
                                * maxdual * (2 * best.x[i] - 1);

                    subgradFactor[p.idx] = param.subgradFactor / nReset;
                    perturbFactor[p.idx] *= 1.1;
                    if (hooks.reducedCost(*p, p->rc) == ABORT) {
                        status = ABORT;
                    }
                    if (hooks.solveSubproblem(*p) == ABORT) {
                        status = ABORT;
                    }
                }
                param.heurFreq = std::max(param.heurFreq / 2, 1);
                param.subgradFmin *= 0.5;
                updateBest(hooks, nIter);
                noImproveIter = 0;
                maxNoImprove += 1000;
                if (param.printLevel > 0) {
                    printf("%2d: Reset swarm LBs were %.2f - %.2f now",
                        nIter, minLB, maxLB);
                    minLB = INF;
                    maxLB = -INF;
                    for (int idx = 0; idx < param.nParticles; ++idx) {
                        if (swarm[idx]->lb < minLB)
                            minLB = swarm[idx]->lb;
                        if (swarm[idx]->lb > maxLB)
                            maxLB = swarm[idx]->lb;
                    }
                    printf(" %.2f - %.2f\n", minLB, maxLB);
                }
            }
        }
        if (param.printLevel && nIter % param.printFreq == 0)
            printf("%2d: LB=%.2f UB=%.2f\n", nIter, best.lb, best.ub);
    }


    if (param.iteration_checks == true) {
        iteration_updates(commodities);
    }
        // //set duals to 0
        // for (int idx = 0; idx < param.nParticles; ++idx) {
        //     ParticleIter p(swarm, idx);
        //     p->dual = 0;
        //     p->rc = 0;
        //     p->perturb = 0;
        //     hooks.solveSubproblem(*p);
        // }

        // sum_lb = 0;
        // for (int idx = 0; idx < param.nParticles; ++idx) {
        //     ParticleIter p(swarm, idx);
        //     sum_lb += p->lb;
        // }
        // dual_0_tracking.push_back(commodities - (sum_lb / param.nParticles));
    // }
}

double Problem::swarmRadius() const
{
    return 0;
} // not yet implemented

double Problem::wallTime() const
{
    return omp_get_wtime() - _wallTime;
}
// set up an initial swarm
void Problem::initialise(UserHooks& hooks)
{
    nIter = 0;
    x_total.resize(psize, 0);
    omp_set_num_threads(param.nCPU);
    _wallTime = omp_get_wtime();
    timer.reset();
    if (swarm.size() >= (size_t)param.nParticles)
        return; // all done
    // figure out what range of duals makes sense
    if (best.dual.size() < (size_t)dsize)
        best.dual.resize(dsize, 0.0);
    if (best.perturb.size() < (size_t)psize)
        best.perturb.resize(psize, 0.0);
    if (best.rc.size() < (size_t)psize)
        best.rc.resize(psize, 0.0);
    hooks.reducedCost(best, best.rc);
    // find biggest absolute value reduced cost
    const double maxCost = std::max(
        *std::max_element(best.rc.begin(), best.rc.end()),
        -*std::min_element(best.rc.begin(), best.rc.end()));
    swarm.reserve(param.nParticles);
    Uniform rand;
    rand.seed(1);
    while (swarm.size() < (size_t)param.nParticles) {
        Particle* p = new Particle(psize, dsize);
        swarm.push_back(p);
        for (int j = 0; j < dsize; ++j) { // generate values up to maxCost
            if (param.zeroInitial == false) {
                p->dual[j] = rand(std::max(dualLB[j], -maxCost),
                    std::min(dualUB[j], maxCost));
            } else {
                p->dual[j] = 0.0;
            }
        }
        p->perturb = 0.0;
        p->dVel = 0.0;
        p->pVel = 0.0;
    }
}
bool Problem::updateBest(UserHooks& hooks, int nIter)
{
    Particle* bestP = 0;
    bool improved = false;
    for (ParticleIter p(swarm); !p.done(); ++p) {
        if (p->isFeasible && p->ub < best.ub) {
            best.ub = p->ub;
            best.isFeasible = true;
            best.x = p->x;
            best.best_ub_sol = p->best_ub_sol;
            best.perturb = p->perturb; // perturbation gives good feasible
            bestP = &(*p);
            if (param.nParticles == 1) {
                best_particles_primal_time = best_particles;
            } else {
                swarm_primal_time = swarm;
            }
            primal_cpu_time = cpuTime();
            best_nIter = nIter;
            if (p->lb > best.lb) {
                lb_primal_cpu_time = p->lb;
            } else {
                lb_primal_cpu_time = best.lb;
            }
        }

        if (p->lb > best.lb) {
            best.lb = p->lb;
            best.dual = p->dual;
            best.viol = p->viol;
            improved = true;
            dual_cpu_time = cpuTime();
        }
        if (p->lb >= best.lb) {
            Particle p1 = *(p);
            if (param.nParticles == 1) {
                best_particles.push_back(&(p1));
            }
        }
    }

    if (bestP) // only call once if current swarm improved optimum
        hooks.updateBest(*bestP);
    return improved || (bestP != 0);
}

/** Calculate a direction for the perturbation that is likely
		to lead to more feasible solutions (using hooks.fixConstraint) */
void Problem::perturbationDirection(UserHooks& hooks, const Particle& p,
    DblVec& dir) const
{
    const double eps = param.eps;
    dir = 0.0; // set all to zero
    for (int i = 0; i < dsize; ++i) {
        if ((p.viol[i] > eps && dualLB[i] < -eps) || (p.viol[i] < -eps && dualUB[i] > eps)) { // constraint violated
            SparseVec feas;
            hooks.fixConstraint(i, p, feas);
            for (SparseVec::iterator x = feas.begin(); x != feas.end(); ++x) {
                if (x->second >= 1)
                    dir[x->first] -= 1;
                if (x->second <= 0)
                    dir[x->first] += 1;
            }
        }
    }
}

double Problem::euclideanDistance(vector<Particle*> swarm, std::string component)
{

    double euclidean_dist = 0.0;

    if (component.compare("dual") == 0) {
        for (int idx = 0; idx < param.nParticles; ++idx) {
            ParticleIter p1(swarm, idx);
            for (int idx_2 = idx + 1; idx_2 < param.nParticles; ++idx_2) {
                double pair_euclid = 0.0;
                ParticleIter p2(swarm, idx_2);
                for (int i = 0; i < p1->dual.size(); i++) {
                    pair_euclid += (p1->dual[i] - p2->dual[i]) * (p1->dual[i] - p2->dual[i]);
                    //printf("pair euclid is %f\n",pair_euclid);
                }
                euclidean_dist += std::sqrt(pair_euclid);
            }
        }

    }

    else if (component.compare("perturb") == 0) {
        for (int idx = 0; idx < param.nParticles; ++idx) {
            ParticleIter p1(swarm, idx);
            for (int idx_2 = idx + 1; idx_2 < param.nParticles; ++idx_2) {
                double pair_euclid = 0.0;
                ParticleIter p2(swarm, idx_2);
                for (int i = 0; i < p1->perturb.size(); i++) {
                    pair_euclid += (p1->perturb[i] - p2->perturb[i]) * (p1->perturb[i] - p2->perturb[i]);
                }
                euclidean_dist += std::sqrt(pair_euclid);
            }
        }
    }

    else {
        printf("error in eucld dist - method not set properly");
        exit(1);
    }

    if (component.compare("dual") == 0) {
        euclidean_dist = (euclidean_dist / dsize) / ((param.nParticles * (param.nParticles - 1) / 2));
    } else if (component.compare("perturb") == 0) {
        euclidean_dist = (euclidean_dist / psize) / ((param.nParticles * (param.nParticles - 1) / 2));
    }

    return euclidean_dist;
}

}; // end namespace

/* Stuff for emacs/xemacs to conform to the "Visual Studio formatting standard"
 * Local Variables:
 * tab-width: 4
 * eval: (c-set-style "stroustrup")
 * End:
*/
