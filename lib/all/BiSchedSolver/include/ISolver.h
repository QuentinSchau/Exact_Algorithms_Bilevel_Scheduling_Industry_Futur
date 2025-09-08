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

#ifndef BILEVEL_SCHEDULING_ISOLVER_H
#define BILEVEL_SCHEDULING_ISOLVER_H

#include <iostream>
#include <fstream>
#include "Solution.h"
#include "Node.h"
#include <chrono>
#include <nlohmann/json.hpp>

// epsilon for equality between two double : a - b < EPSILON <=> a~=b
#define EPSILON 1E-6

/**
 * Interface for the Solver. Each new method or algorithm to solve the bilevel scheduling problem must implement with
 * its own class that use ISolver interface
 */
class ISolver {
protected:
    Solution *solution; // the solution of the problem
    Instance *instance; // a reference on instance
    char verbose;
    std::chrono::steady_clock::time_point start;
    std::chrono::duration<double> time_elapsed;
    std::chrono::duration<double> time_limits = std::chrono::seconds(60);

public:

    /************************/
    /*      CONSTRUCTOR     */
    /************************/

    explicit ISolver() : solution(new Solution()), instance(nullptr), verbose(0) {}

    /**
     * Constructor of a solver from a instance.
     * @param instance The instance to solve
     */
    explicit ISolver(Instance *instance) :
            solution(new Solution(instance)), instance(instance), verbose(0) {}

    /**
     * Constructor of a solver from an json object that solve the instance.
     * @param instance The instance to solve
     * @param object The json object with the parameters of the solver
     */
    explicit ISolver(Instance *instance, [[maybe_unused]] nlohmann::json &object) :
            solution(new Solution(instance)), instance(instance), verbose(0) {}

    /***********************/
    /*      DESTRUCTOR     */
    /***********************/

    virtual ~ISolver() {
        // save the solution in a "solutions" directory at same location as the instance file
        std::string strPath = instance->getInstancePath().parent_path().string();
        strPath.append("/solutions/");
        strPath.append(instance->getInstanceName()).append(".sol");

        std::filesystem::path fileSolutionPath = std::filesystem::path(strPath);
        std::filesystem::create_directories(fileSolutionPath.lexically_normal().parent_path());
        std::ofstream outputFileStream;
        outputFileStream.open(fileSolutionPath, std::ios::out | std::ios::trunc);
        outputFileStream << *solution;
        outputFileStream.close();
        delete solution;
    }


    /********************/
    /*      METHODS     */
    /********************/

    /**
     * Pure Virtual Method that implement a algorithm to solve the current instance.
     */
    virtual void solve() = 0;

    /**
     * Method that solves all sub-problems with identical jobs, by rearranging them in an optimal way so that the weighted sum of tardy equal-size jobs is minimal
     * @param sol The solution to be improved.
     */
    void solveSubProblemIdenticalJobs(Solution &sol) {
        auto blockStruct = sol.toBlockStruct(instance);
        // declare the binary tree of completion times that is not assigned with the corresponding position,
        // i.e. key is the completion time and the value is (index of machine, index of block). This is sorted in increasing order
        typedef std::multimap<double, std::pair<unsigned int, unsigned int>> TreeCj;
        // first we determine all identical jobs in the schedule, and we minimize the weight optimally
        // use a vector of (index in list grouped jobs, tree with completion time) we will solve optimally all elements
        // of this vector
        std::vector<std::pair<unsigned int, TreeCj>> listIdenticalJobAndTheirTree;
        unsigned int indexBlock = instance->getE().size(); // start from this end of the schedule
        auto &mapJobToItsGroupIdenticalJobs = instance->getMapListJobToListGroupedJobs();
        while (indexBlock) {
            --indexBlock;
            for (auto [indexMachine, indexBlockInStruct]: instance->getE()[indexBlock]) {
                auto job = blockStruct[indexMachine][indexBlockInStruct].first;
                // check if a job is scheduled
                if (job != nullptr) {
                    // check if there exist identical jobs
                    unsigned int indexInListOfGroupedJobs = mapJobToItsGroupIdenticalJobs[job->getIndex()];
                    assert(indexInListOfGroupedJobs < instance->getListGrpJobs().size());
                    auto &groupJobs = instance->getListGrpJobs()[indexInListOfGroupedJobs];
                    if (groupJobs.size() > 1) {
                        double completionTime = blockStruct[indexMachine][indexBlockInStruct].second; // get the completion time of the job
                        auto predIndexIdenticalGrp = [indexInListOfGroupedJobs](const std::pair<unsigned int, TreeCj> &pair) { return pair.first == indexInListOfGroupedJobs; };
                        // if we have already a tree with Cj and location for this group of jobs
                        auto itGroupAndItsTree = std::find_if(listIdenticalJobAndTheirTree.begin(), listIdenticalJobAndTheirTree.end(), predIndexIdenticalGrp);
                        if (itGroupAndItsTree != listIdenticalJobAndTheirTree.end()) {
                            itGroupAndItsTree->second.insert({completionTime, {indexMachine, indexBlockInStruct}});
                        }//else we add the index and its group
                        else {
                            TreeCj treeCj;
                            treeCj.insert({completionTime, {indexMachine, indexBlockInStruct}});
                            listIdenticalJobAndTheirTree.emplace_back(indexInListOfGroupedJobs, treeCj);
                        }
                    }

                }
            }
        }
        for (unsigned int indexLoopTreeCj = 0; indexLoopTreeCj < listIdenticalJobAndTheirTree.size(); ++indexLoopTreeCj) {
            auto &[indexListOfIdenticalJobs, listCjAndAvailablePosition] = listIdenticalJobAndTheirTree[indexLoopTreeCj];
            // use a copy of the list of grouped jobs, because it's modifying by the method
            std::vector<Job> listOfIdenticalJobs(instance->getListGrpJobs()[indexListOfIdenticalJobs]);
            solveProblemWithIdenticalJobs(nullptr, &blockStruct, listCjAndAvailablePosition, listOfIdenticalJobs);
        }
        sol.fromBlockStruct(blockStruct);
    }

    /**
     * Method that solves all sub-problems with identical jobs, by rearranging them in an optimal way so that the weighted sum of tardy equal-size jobs is minimal
     * @param blockStruct The solution to be improved.
     */
    void solveSubProblemIdenticalJobs(Solution::BlockStructure &blockStruct) {

        // declare the binary tree of completion times that is not assigned with the corresponding position,
        // i.e. key is the completion time and the value is (index of machine, index of block). This is sorted in increasing order
        typedef std::multimap<double, std::pair<unsigned int, unsigned int>> TreeCj;
        // first we determine all identical jobs in the schedule, and we minimize the weight optimally
        // use a vector of (index in list grouped jobs, tree with completion time) we will solve optimally all elements
        // of this vector
        std::vector<std::pair<unsigned int, TreeCj>> listIdenticalJobAndTheirTree;
        unsigned int indexBlock = instance->getE().size(); // start from this end of the schedule
        auto &mapJobToItsGroupIdenticalJobs = instance->getMapListJobToListGroupedJobs();
        while (indexBlock) {
            --indexBlock;
            for (auto [indexMachine, indexBlockInStruct]: instance->getE()[indexBlock]) {
                auto job = blockStruct[indexMachine][indexBlockInStruct].first;
                // check if a job is scheduled
                if (job != nullptr) {
                    // check if there exist identical jobs
                    unsigned int indexInListOfGroupedJobs = mapJobToItsGroupIdenticalJobs[job->getIndex()];
                    assert(indexInListOfGroupedJobs < instance->getListGrpJobs().size());
                    auto &groupJobs = instance->getListGrpJobs()[indexInListOfGroupedJobs];
                    if (groupJobs.size() > 1) {
                        double completionTime = blockStruct[indexMachine][indexBlockInStruct].second; // get the completion time of the job
                        auto predIndexIdenticalGrp = [indexInListOfGroupedJobs](const std::pair<unsigned int, TreeCj> &pair) { return pair.first == indexInListOfGroupedJobs; };
                        // if we have already a tree with Cj and location for this group of jobs
                        auto itGroupAndItsTree = std::find_if(listIdenticalJobAndTheirTree.begin(), listIdenticalJobAndTheirTree.end(), predIndexIdenticalGrp);
                        if (itGroupAndItsTree != listIdenticalJobAndTheirTree.end()) {
                            itGroupAndItsTree->second.insert({completionTime, {indexMachine, indexBlockInStruct}});
                        }//else we add the index and its group
                        else {
                            TreeCj treeCj;
                            treeCj.insert({completionTime, {indexMachine, indexBlockInStruct}});
                            listIdenticalJobAndTheirTree.emplace_back(indexInListOfGroupedJobs, treeCj);
                        }
                    }

                }
            }
        }
        for (unsigned int indexLoopTreeCj = 0; indexLoopTreeCj < listIdenticalJobAndTheirTree.size(); ++indexLoopTreeCj) {
            auto &[indexListOfIdenticalJobs, listCjAndAvailablePosition] = listIdenticalJobAndTheirTree[indexLoopTreeCj];
            // use a copy of the list of grouped jobs, because it's modifying by the method
            std::vector<Job> listOfIdenticalJobs(instance->getListGrpJobs()[indexListOfIdenticalJobs]);
            solveProblemWithIdenticalJobs(nullptr, &blockStruct, listCjAndAvailablePosition, listOfIdenticalJobs);
        }
    }

    /**
     * Solves the problem with identical jobs.
     *
     * @param node The current node to solve for.
     * @param blockStruct A pointer to the block structure that will be updated.
     * @param pj The processing time of identical jobs.
     * @param listCjAndAvailablePosition A multimap of completion times and their corresponding positions in the block.
     * @param listOfIdenticalJobs A vector of jobs that are identical.
     */
    inline void solveProblemWithIdenticalJobs(Node *node, Solution::BlockStructure *blockStruct, std::multimap<double, std::pair<unsigned int, unsigned int>> &listCjAndAvailablePosition
                                              , std::vector<Job> &listOfIdenticalJobs) {
        unsigned int numberElementToAdd = listCjAndAvailablePosition.size();
        // use iterator to iterate over completion time from the back (in order to get High Cj first)
        auto itCompletionTime = listCjAndAvailablePosition.rbegin();
        auto itJobs = listOfIdenticalJobs.rbegin();

        // declare the queue for max weighted early job
        auto queueEarlyJob = std::priority_queue<std::pair<double, std::vector<Job>::reverse_iterator>, std::vector<std::pair<double, std::vector<Job>::reverse_iterator >>>();

        // we use this list to delete the jobs after determine on with completion time the job will be scheduled.
        std::vector<unsigned int> indicesOfJobsToDelete;
        indicesOfJobsToDelete.reserve(listOfIdenticalJobs.size());
        // compute the set of all early jobs for each completion times, we start from the end because we remove element
        while (itCompletionTime != listCjAndAvailablePosition.rend() && itJobs != listOfIdenticalJobs.rend() && queueEarlyJob.size() < numberElementToAdd) {
            // If the job can be schedule on time, i.e., dj >= Ci, we add it
            if (isSmallerOrEqual(itCompletionTime->first, itJobs->getDi())) {
                queueEarlyJob.emplace(itJobs->getWi(), itJobs);
                ++itJobs;
            } else if (!queueEarlyJob.empty()) {
                // else if we have some jobs that early (means there are in the queue), so we affect the one with max weight
                auto itAssignedJob = queueEarlyJob.top().second;
                queueEarlyJob.pop();
                // add the job to the block structure
                if (node) {
                    node->scheduleOneJob(itCompletionTime->second.first, itCompletionTime->second.second, &instance->getListJobs()[itAssignedJob->getIndex()], itCompletionTime->first);
                }
                if (blockStruct) {
                    (*blockStruct)[itCompletionTime->second.first][itCompletionTime->second.second] = {&instance->getListJobs()[itAssignedJob->getIndex()], itCompletionTime->first};
                }

                // remove the completion time and the job
                itCompletionTime = std::reverse_iterator(
                        listCjAndAvailablePosition.erase(std::next(itCompletionTime).base()));
                //add the index to the list of indices of job to delete
                indicesOfJobsToDelete.emplace_back(itAssignedJob->getIndex());

                --numberElementToAdd;
            } else {
                // we pass to the next Cj
                ++itCompletionTime;
            }
        }

        // We clear out the queue if it not empty
        while (!queueEarlyJob.empty() && itCompletionTime != listCjAndAvailablePosition.rend()) {
            auto itAssignedJob = queueEarlyJob.top().second;
            queueEarlyJob.pop();
            if (node) {
                node->scheduleOneJob(itCompletionTime->second.first, itCompletionTime->second.second, &instance->getListJobs()[itAssignedJob->getIndex()], itCompletionTime->first);
            }
            if (blockStruct) {
                (*blockStruct)[itCompletionTime->second.first][itCompletionTime->second.second] = {&instance->getListJobs()[itAssignedJob->getIndex()], itCompletionTime->first};
            }
            // remove the completion time and the job
            itCompletionTime = std::reverse_iterator(
                    listCjAndAvailablePosition.erase(std::next(itCompletionTime).base()));
            //add the index to the list of indices of job to delete
            indicesOfJobsToDelete.emplace_back(itAssignedJob->getIndex());
        }

        //remove all jobs that are already scheduled
        auto pred_remove_job = [&indicesOfJobsToDelete](const Job &job) {
            return std::find(indicesOfJobsToDelete.begin(), indicesOfJobsToDelete.end(), job.getIndex()) != indicesOfJobsToDelete.end();
        };
        auto itRemovedJobs = std::remove_if(listOfIdenticalJobs.begin(), listOfIdenticalJobs.end(), pred_remove_job);
        listOfIdenticalJobs.erase(itRemovedJobs, listOfIdenticalJobs.end());

        // now, we schedule not assigned jobs, i.e., the late ones, so iterate through the list of jobs
        itCompletionTime = listCjAndAvailablePosition.rbegin();
        // sort the remaining job according non-increasing weight, and loop over it from the end (to get the smallest weight first)
        std::sort(listOfIdenticalJobs.begin(), listOfIdenticalJobs.end(), [](const Job &lhs, const Job &rhs) { return lhs.getWi() >= rhs.getWi(); });
        itJobs = listOfIdenticalJobs.rbegin();
        while (itCompletionTime != listCjAndAvailablePosition.rend() && itJobs != listOfIdenticalJobs.rend()) {
            if (node) {
                node->scheduleOneJob(itCompletionTime->second.first, itCompletionTime->second.second, &instance->getListJobs()[itJobs->getIndex()], itCompletionTime->first);
            }
            if (blockStruct) {
                (*blockStruct)[itCompletionTime->second.first][itCompletionTime->second.second] = {&instance->getListJobs()[itJobs->getIndex()], itCompletionTime->first};
            }
            itJobs = std::reverse_iterator(
                    listOfIdenticalJobs.erase(std::next(itJobs).base()));
            itCompletionTime = std::reverse_iterator(
                    listCjAndAvailablePosition.erase(std::next(itCompletionTime).base()));
        }
    }

    /********************/
    /*      GETTER      */
    /********************/

    bool isWithinTimeLimit(bool raiseException=false,bool updateElapsedTime=false) {
        // measure the time
        const auto measureTime{std::chrono::steady_clock::now()};
        auto timeElapsed = std::chrono::duration<double>{measureTime - start};
        if ((time_elapsed.count() + timeElapsed.count()) > time_limits.count()) {
            if (verbose >= 2) std::cout << "Stop not enough time " << (time_elapsed.count() + timeElapsed.count()) << "/" << time_limits.count() << "s" << std::endl;
            if (raiseException) throw BiSchTimeOutException();
        }
        if (updateElapsedTime) time_elapsed = time_elapsed + timeElapsed;
        return false;
    }

    [[nodiscard]] Solution *getSolution() const { return solution; }

    [[nodiscard]] const Instance *getInstance() const { return instance; }

    [[nodiscard]] char levelVerbose() const { return verbose; }

    [[nodiscard]] const std::chrono::duration<double> &getTimeElapsed() const { return time_elapsed; }

    /********************/
    /*      SETTER      */
    /********************/

    void setSolution(Solution *solution) { ISolver::solution = solution; }

    void setInstance(Instance *instance) { ISolver::instance = instance; }

    void setVerbose(char verbose) { ISolver::verbose = verbose; }

    void setTimeLimit(unsigned int seconds) { time_limits = std::chrono::seconds(seconds); };

    void setTimeLimitInMilliSecond(unsigned int newMilli_seconds) { time_limits = std::chrono::milliseconds(newMilli_seconds); };

    void setTimeElapsed(unsigned int seconds) { time_elapsed = std::chrono::seconds(seconds);; }
};

template<typename T>
inline std::ostream &operator<<(std::ostream &os, const std::vector<T> &vector) {

    os << "[";
    if (vector.empty()) os << "]";
    else {
        for (size_t indexLoopVector = 0; indexLoopVector < vector.size() - 1; ++indexLoopVector) {
            const T &element = vector[indexLoopVector];
            os << element << ",";
        }
        os << vector.back() << "]";
    }
    return os;
}


#endif //BILEVEL_SCHEDULING_ISOLVER_H
