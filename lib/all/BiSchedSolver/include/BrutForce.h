//  Copyright (C) 2024
//  Laboratoire d'Informatique Fondamentale et Appliquée de Tours, Tours, France
//
//  DIGEP, Politecnico di Torino, Corso Duca degli Abruzzi 24, Torino, Italy
//  This file is part of bilevel-scheduling.
//
//  bilevel-scheduling is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published
//  by the Free Software Foundation, either version 3 of the License,
//  or (at your option) any later version.
//
//  bilevel-scheduling is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty
//  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//  See the GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with bilevel-scheduling. If not, see <https://www.gnu.org/licenses/>.

#ifndef BILEVEL_SCHEDULING_BRUTFORCE_H
#define BILEVEL_SCHEDULING_BRUTFORCE_H

// uncomment this line to debug the method
//#define DEBUG_BRUTFORCE

#include "ISolver.h"
#include "Math.h"
#include <algorithm>
#include <map>
#include <limits>

/**œ
 * Class that implement a brut force method to solve the bilevel scheduling problem
 */
class BrutForce : public ISolver {
public:
    /************************/
    /*      CONSTRUCTOR     */
    /************************/


    BrutForce() : ISolver() {}

    /**
     * Constructor of a solver from a instance.
     * @param instance The instance to solve
     */
    explicit BrutForce(Instance *instance) : ISolver(instance) {}


    BrutForce(Instance *instance, nlohmann::json &object);

    /********************/
    /*      METHODS     */
    /********************/

    /**
     * Method that solve by brut force the bilevel problem. The instance MUST HAVE BEEN SORTED according LPT rule
     */
    virtual void solve() override;

    /**
     * Method that save the result of the instance in a file
     * @param fileOutputName The name of the file
     * @param outputFile The stream of the file
     */
    void printOutput(std::string &fileOutputName, std::ofstream &outputFile);


    /**
     * Method that from a list of jobs computes all lists where jobs with same processing times are permuted.
     * @param listPermutedJobs The lists where each permutation is added.
     * @param listJobs The list of index of jobs where all permutation of jobs with same processing time are computed.
     */
    void
    listOfPermutedJobWithSamePi(std::vector<std::vector<unsigned int>> &listPermutedJobs, const std::vector<unsigned int> &listJobs) {
        // first we create a list of jobs group by processing time
        std::vector<std::vector<unsigned int>> listJobsGroupByPi;
        std::map<unsigned int, std::vector<unsigned int>> map;
        for (const unsigned int &index: listJobs) {
            map[instance->getListJobs()[index].getPi()].push_back(index);
        }
        // transform the map on vector using move to not copy
        listJobsGroupByPi.reserve(map.size());
        // we add jobs in list according LPT rule (iterate by the end)
        for (auto it = map.rbegin(); it != map.rend(); ++it) {
            listJobsGroupByPi.push_back(std::move(it->second));
        }
        generateConcatenatedPermutations(listPermutedJobs, listJobsGroupByPi);
    }


    /**
     * Method that solves recursively for all sub-vector of the given size in the list of all jobs. This set the Solution pointer of Brutforce class. This method don't manage the
     * pointer life. So you must delete after use.
     * @param subVector The sub-vector to compute that contains jobs index. You must give an empty at first call
     * @param listAllJobs The list of all jobs given by the instance
     * @param listOfLocationForSumCj The list of all location for selected jobs. This list must contains all permuted position to explore all search space.
     * @param nbToSelect The number of job that must be selected.
     * @param nbSelectedJob The number of selected Job
     * @param depth The depth in the call stack of recursion
     */
    void solveForAllSubVector(std::vector<unsigned int> &subVector, const std::vector<Job> &listAllJobs, const std::vector<std::vector<std::pair<unsigned int, unsigned int>>> &listOfLocationForSumCj
                              , unsigned int nbToSelect, unsigned int nbSelectedJob, unsigned depth) {
        // if we have already an optimal solution (with Sum wjUj == 0)
        if (solution->getSumWjUj() == 0) return;
        else if (depth == nbToSelect) {
            #ifdef DEBUG_BRUTFORCE
            std::cout << "sub-jobs = [ " ;
            for (auto index  : subVector)  std::cout << index  << " ";
            std::cout << "]" << std::endl;
            #endif
            // base case : we have a sub-vector with the good size. We search the best solution
            // start time to measure performance
            const auto start = std::chrono::steady_clock::now();

            // if we have spent more time than the limit
            if (time_elapsed.count() > time_limits.count()) return;

            Solution *currentSol = findBestSolFromListOfJobsAndListPos(subVector, listOfLocationForSumCj);
            if (currentSol->getSumWjUj() < solution->getSumWjUj()) {
                solution->reset();
                *solution = *currentSol;
            }
            // stop time to measure performance
            const auto end{std::chrono::steady_clock::now()};
            time_elapsed = std::chrono::duration<double>{end - start};
            delete currentSol;


        } else {
            // recursion : loop over the list of all jobs and recurse for adding in sub-vector
            for (unsigned int j = nbSelectedJob + 1; j < listAllJobs.size(); ++j) {
                // add the jth job to the subVector
                subVector.push_back(j);
                // compute sub-vector with nbSelected + 1 elements
                solveForAllSubVector(subVector, listAllJobs, listOfLocationForSumCj, nbToSelect, j, depth + 1);
                // remove the jth added job for the next iteration
                subVector.pop_back();
            }
        }
    }

    /**
     * Method that find the best solution given a list of all positions on machines and list of Jobs. This method compute all permutations of same jobs' processing times for explore all
     * solution space. The method don't manage the pointer that be returned. You must delete it after use.
     * @param listJobs The list of jobs to be scheduled
     * @param listOfAllPosition The list of all positions allowed on machines
     * @return A pointer on the optimal solution.
     */
    Solution *findBestSolFromListOfJobsAndListPos(const std::vector<unsigned int> &listJobs, const std::vector<std::vector<std::pair<unsigned int, unsigned int>>> &listOfAllPosition) {
        std::vector<std::vector<unsigned int>> listAllPermutedJob;
        // compute the list of permuted identical jobs
        listOfPermutedJobWithSamePi(listAllPermutedJob, listJobs);

        Solution *bestSol = new Solution(instance);
        bestSol->setSumWjUj(std::numeric_limits<double>::infinity());
        Solution *currentSol = new Solution(instance);

        //loop over each list of jobs (where identical jobs were permuted)
        for (const std::vector<unsigned int> &listPermutedJob: listAllPermutedJob) {
            //loop over all feasible locations on machines
            for (const auto &listLocations: listOfAllPosition) {
                // add jobs to the right location. We need to loop from the end of both list of permuted
                // job and location (because we create it to schedule from the RIGHT to the LEFT)
                for (unsigned int i = instance->getNbToSelectJob(); i-- > 0;) {
                    auto position = listLocations[i];
                    (*currentSol).add_job(position.first, instance->getListJobs()[listPermutedJob[i]]);
                }
                // evaluate the solution
                currentSol->evaluate();
                if (currentSol->getSumWjUj() < bestSol->getSumWjUj()) {
                    bestSol->reset();
                    *bestSol = *currentSol;
                }
                currentSol->reset();
            }
        }
        delete currentSol;
        return bestSol;
    }

};

#endif //BILEVEL_SCHEDULING_BRUTFORCE_H
