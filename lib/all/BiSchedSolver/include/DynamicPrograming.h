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

#ifndef BILEVEL_SCHEDULING_DYNAMICPROGRAMING_H
#define BILEVEL_SCHEDULING_DYNAMICPROGRAMING_H

// uncomment this line to debug the method
//#define DEBUG_DP

#include "ISolver.h"
#include "Position.h"
#include <map>
#include <iterator>

class DynamicPrograming : public ISolver {

public:

    struct State_recursion {
        std::vector<bool> A;
        std::vector<int> C;
        double minSumWj;

        //defined move constructor
        inline State_recursion(std::vector<bool> &&A, std::vector<int> &&C, double minSumWj) noexcept {
            State_recursion::A = std::move(A);
            State_recursion::C = std::move(C);
            State_recursion::minSumWj = minSumWj;
        }
    };

    /************************/
    /*      CONSTRUCTOR     */
    /************************/

    DynamicPrograming() = default;


    /**
     * Constructor of a solver from a instance. It set the memory size at 500 Mo by default
     * @param instance The instance to solve, sorted according to the LPT rule.
     */
    explicit DynamicPrograming(Instance *instance);

    explicit DynamicPrograming(Instance *instance, nlohmann::json &object);

    /********************/
    /*      METHODS     */
    /********************/

    size_t computeMemorySize(size_t nbMegabytesAllowedByUser);


    /**
     * Method that looks at if the job is late from completion time. It's return the weight if it is late, otherwise zero.
     * @param completionTime The completion time at the position on the machine
     * @param speed The speed of the machine
     * @param job The job to add at the position corresponding to the completion time
     * @return The weight of the job if it's late, zero otherwise
     */
    static double weightedTobeLate(int completionTime, double speed, const Job &job);

    /**
     * Method that computes the minimum and maximum for the sum of completion time that is realizable for bilevel problem.
     * We solve for the sum Cj objective function, and compute the min and max at each block. By adding them, we can obtain
     * the min and max of the sum of completion time, if we swap job inside a block.
     * @param listSelectedJobs The list of index for selected jobs
     * @return (min Sum Cj for Vmax, max Sum Cj for Vmax, min Sum Cj for V0, max Sum Cj for V0)
     */
    std::array<unsigned int, 4> computeMinMaxSumCompletionTimes(std::vector<size_t> &listIndexSelectedJobs);

    /**
     * Recursive method, that increments the value of C[m] from min Sum Cj to max Sum Cj for the corresponding machine. So, this method allow the exploration
     * of all vector C starting from min Sum Cj to max Sum Cj.
     * @param C The vector of completion times
     * @param m The index of the machine
     * @param minValue The minValue find by the method
     * @param minMaxCj The tuple (min Sum Cj for Vmax, max Sum Cj for Vmax, min Sum Cj for V0, max Sum Cj for V0)
     * @return The optimal value from the bilevel problem by try all completion times
     */
    double tryAllSumCj(std::vector<int> &C, unsigned int m, double minValue, std::array<unsigned int, 4> &minMaxCj);

    /**
     * Method that saves the result of the instance in a file
     * @param fileOutputName The name of the file
     * @param outputFile The stream of the file
     */
    void printOutput(std::string &fileOutputName, std::ofstream &outputFile);

    /**
     * Method that solves by dynamic programing the bilevel problem. The instance MUST HAVE BEEN SORTED according LPT rule
     */
    virtual void solve() override;

private:

    // Unordered_map to keep the state that be already explored. For each state already explored, we keep in mind the corresponding optimal value
    std::unordered_map<std::vector<int>, double, HashFunction> memoization;

    size_t absMemoizationSize;

    unsigned int nbCleaning;
    /*******************************/
    /*     SUB-PROBLEM METHODS     */
    /*******************************/

    /**
     * Method that takes the max weighted Job from a priorityQueue and the current position, i.e. the machine and the block, and add to a list of affectations.
     * @param earlyJobs The priority queue that contains early jobs that match to the position given by indexMachine and indexBlock
     * @param listOfJobs The list of jobs, in the form of multimap, where the added job is removed
     * i.e. each element is a tuple (index of Machine, index of block, index of Job)
     */
    void affectEarlyJobsToPosition(
            std::priority_queue<std::pair<double, unsigned int>, std::vector<std::pair<double, unsigned int>>> &earlyJobs, std::multimap<unsigned int, unsigned int, std::greater<>> &listOfJobs);

    /**
     * Method that computes the number of jobs to schedule on last block from an index of job j and the number of jobs
     * with processing times jobs
     * @param n The number of job that have been selected by the leader
     * @param n_j A pair that contains : (the number of identical jobs, the processing time of this jobs)
     * @return mu_j the number of identical jobs in the last block
     */
    unsigned int computeNumberOfJobsInLastBlock(unsigned int n, unsigned int n_j);

    /**
     * Method that updates the vector A according to given position. I.e. if the job J and the job J+1 are in two different block on two type of machine,
     * then we make available the good ones.
     * @param positionJobJ The position of the job J
     * @param positionNextJobAfterJ The position of the job J+1
     * @param A The vector of boolean, which characterize if a machine is available.
     */
    void updateAvailableMachine(Position &positionJobJ, Position &positionNextJobAfterJ, std::vector<bool> &A);

    /**
     * Method that computes the list of completion times and corresponding available position from a combination of the
     * affection of the job to the last block. It update the vector A and C according the affectation.
     * @param pj The processing time of the identical jobs
     * @param n_j A pair that contains : (the number of identical jobs, the processing time of this jobs)
     * @param n The number of job that have been selected by the leader
     * @param mu_j The number of identical Jobs on the last block
     * @param A The vector of available machines. The index is in range of 0 to the number of machines.
     * The first index is for high-speed machines, the last one for the low-speed machines. It will use a copy of vector.
     * @param C The vector of completion times for each machines. The index is in range of 0 to the number of machines.
     * The first index is for high-speed machines, the last one for the low-speed machines. It will use a copy of vector.
     * @param combinationAffectationLastBlock A combination of affectation of the first job on the last Block
     * @param completionTimes The list of completion times to compute
     * @param listAvailablePosition The map of available position to compute. Each key is the completion time whereas the value is a pair form of (index machine, index block)
     */
    void computeListCompletionTimeAndAvailablePosition(int pj, unsigned int n_j, unsigned int n, unsigned int mu_j, std::vector<bool> &A, std::vector<int> &C
                                                       , std::vector<bool> &combinationAffectationLastBlock, std::multimap<double, std::pair<unsigned int, unsigned int>> &listCjAndAvailablePosition);

    /**
     * This method solves the sub-problem by assign the job with same processing time. The affection is optimal for our problem.
     * This method uses a function f which returns the machine and the block position for an index 'n' of a job.
     * Its implementation is in Position structure.
     * @param A The vector of available machines. The index is in range of 0 to the number of machines.
     * The first index is for high-speed machines, the last one for the low-speed machines.
     * @param C The vector of completion times for each machines. The index is in range of 0 to the number of machines.
     * The first index is for high-speed machines, the last one for the low-speed machines.
     * @param indexIdenticalJobs The vector of identical jobs that will be scheduled.
     * @param n The number of jobs that the leader have already selected.
     * @param listStateRecursion A list of recursion states that is generate by all assignment on last block.
     */
    void solveSubProblem(std::vector<bool> &A, std::vector<int> &C, std::vector<unsigned int> &indexIdenticalJobs, unsigned int n, std::vector<DynamicPrograming::State_recursion> &listStateRecursion);

    /**************************************/
    /*     DYNAMIC PROGRAMING METHODS     */
    /**************************************/

    /**
     * Recursive method that solves the dynamic programing.
     * @param A The vector of size : number of machine. Each element is a boolean, if it is true, then the corresponding machine is available, otherwise is not.
     * @param C The vector of size : number of machine. Each element corresponds to the completion time on the machine
     * @param j The index of the considering job
     * @param indexIdenticalJobs The vector of identical jobs that will be scheduled.
     * @param n The number of job that the leader have already selected
     * @return The minimum of the objective value
     */
    double solveDynamicPrograming(std::vector<bool> A, std::vector<int> C, unsigned int j, std::vector<unsigned int> indexIdenticalJobs, unsigned int n);

};

/*******************/
/*     METHODS     */
/*******************/

inline double DynamicPrograming::weightedTobeLate(int completionTime, double speed, const Job &job) {
    if ((double(completionTime) / speed) > job.getDi()) return job.getWi();
    else return 0.0;
}

inline std::array<unsigned int, 4>
DynamicPrograming::computeMinMaxSumCompletionTimes(std::vector<size_t> &listIndexSelectedJobs) {

    int lastBlockHighSpeed = 0;
    int lastBlockLowSpeed = 0;

    unsigned int minPjHigh = 0;
    unsigned int maxPjHigh = 0;
    unsigned int minPjLow = 0;
    unsigned int maxPjLow = 0;

    std::array<unsigned int, 4> result{0, 0, 0, 0};

    for (size_t i = 0; i < listIndexSelectedJobs.size(); ++i) {
        Position position = Position(instance, i + 1);
        // if we are in new block for high speed machines
        if (position.machines[0] == 0) {
            if (position.blocks[0] != lastBlockHighSpeed) {
                lastBlockHighSpeed = position.blocks[0];
                result[0] += minPjHigh;
                result[1] += maxPjHigh;
                minPjHigh = std::numeric_limits<unsigned int>::max();
                maxPjHigh = 0;
            }
            minPjHigh = std::min(static_cast<unsigned int>(instance->getListJobs()[listIndexSelectedJobs[i]].getPi()), minPjHigh);
            maxPjHigh = std::max(static_cast<unsigned int>(instance->getListJobs()[listIndexSelectedJobs[i]].getPi()), maxPjHigh);
        }
        // if we are in new block for low speed machines
        if (position.machines[0] == 1) {
            if (position.blocks[0] != lastBlockLowSpeed) {
                lastBlockLowSpeed = position.blocks[0];
                result[2] += minPjLow;
                result[3] += maxPjLow;
                minPjLow = std::numeric_limits<unsigned int>::max();
                maxPjLow = 0;
            }
            minPjLow = std::min(static_cast<unsigned int>(instance->getListJobs()[listIndexSelectedJobs[i]].getPi()), minPjLow);
            maxPjLow = std::max(static_cast<unsigned int>(instance->getListJobs()[listIndexSelectedJobs[i]].getPi()), maxPjLow);
        }
        // if we are in new block for low speed machines with high speed machines
        if (position.machines[1] != -1) {
            if (position.blocks[1] != lastBlockLowSpeed) {
                lastBlockLowSpeed = position.blocks[1];
                result[2] += minPjLow;
                result[3] += maxPjLow;
                minPjLow = std::numeric_limits<unsigned int>::max();
                maxPjLow = 0;
            }
            minPjLow = std::min(static_cast<unsigned int>(instance->getListJobs()[listIndexSelectedJobs[i]].getPi()), minPjLow);
            maxPjLow = std::max(static_cast<unsigned int>(instance->getListJobs()[listIndexSelectedJobs[i]].getPi()), maxPjLow);
        }
    }

    // add the maximum of the last block
    result[1] += maxPjHigh;
    result[3] += maxPjLow;

    return result;
}


inline double DynamicPrograming::tryAllSumCj(std::vector<int> &C, unsigned int m, double minValue, std::array<unsigned int, 4> &minMaxCj) {


    // stop time to measure performance
    const auto end{std::chrono::steady_clock::now()};
    auto time = std::chrono::duration<double>{end - start};
    // check if we have enough time
    if (time.count() > time_limits.count()) return minValue;

    // if we have already an optimal solution (with Sum wjUj == 0)
    if (solution->getSumWjUj() == 0) minValue = 0.0;
        // terminal condition of the recursion. We can solve for the current value of C
    else if (m >= C.size()) {

        //create empty solution
        Solution sol = Solution(instance);
        std::vector<bool> A(instance->getNbMachines(), true);
        // the list of identical jobs, at least n jobs
        std::vector<unsigned int> indexIdenticalJobs;
        indexIdenticalJobs.reserve(instance->getNbToSelectJob());
        // the list of affected jobs, there is exactly n jobs
        std::vector<std::pair<unsigned int, unsigned int>> affectedPosition;
        affectedPosition.reserve(instance->getNbToSelectJob());

        #ifdef DEBUG_DP
        std::cout << (double(memoization.size()) / double(absMemoizationSize)) * 100 << " Memory use | " << "C = [ ";
        for (auto cj: C) std::cout << cj << " ";
        std::cout << "]" << std::endl;
        #endif

        // check if we have fill all memoization
        if (memoization.size() >= absMemoizationSize) {
            if (verbose >= 2) std::cout << "MEMORY FULL : clearing it ..." << std::endl;
            memoization.clear();
            ++nbCleaning;
        }

        minValue = std::min(minValue, solveDynamicPrograming(A, C, 0, indexIdenticalJobs, 0));

    } else {
        // try all sum Cj value for the m^th machine
        unsigned int startSumCj = (m < instance->getNbOfHighSpeedMachines()) ? std::get<0>(minMaxCj) : std::get<2>(
                minMaxCj);
        unsigned int endSumCj = (m < instance->getNbOfHighSpeedMachines()) ? std::get<1>(minMaxCj) : std::get<3>(
                minMaxCj);

        for (unsigned int i = startSumCj; i <= endSumCj; ++i) {
            // change completion time value
            C[m] = int(i);
            minValue = std::min(minValue, tryAllSumCj(C, m + 1, minValue, minMaxCj));
        }
    }
    return minValue;
}

/*******************************/
/*     SUB-PROBLEM METHODS     */
/*******************************/

inline void DynamicPrograming::affectEarlyJobsToPosition(
        std::priority_queue<std::pair<double, unsigned int>, std::vector<std::pair<double, unsigned int>>> &earlyJobs, std::multimap<unsigned int, unsigned int, std::greater<>> &listOfJobs) {
    // take the job with the max weight
    auto maxIndexJob = earlyJobs.top();
    // remove this from the queue
    earlyJobs.pop();
    // remove the affected job on the list of jobs
    auto range = listOfJobs.equal_range(
            instance->getListJobs()[maxIndexJob.second].getDi()); // get the range of jobs with same dj
    // iterator to the next element after removed job
    auto newItJobs = range.first;
    bool haveRemoved = false; // boolean to know if we have removed the maxIndexJob
    // remove the affect job
    while (newItJobs != range.second && !haveRemoved) {
        if (newItJobs->second == maxIndexJob.second) {
            listOfJobs.erase(newItJobs);
            haveRemoved = true;
        } else {
            newItJobs++;
        }
    }
}

inline void DynamicPrograming::updateAvailableMachine(Position &positionJobJ, Position &positionNextJobAfterJ, std::vector<bool> &A) {
    // defined the index of the first and last machines in the block of job J
    auto rangeIndexBlock = positionJobJ.range(instance);
    // defined the index of the first and last machines in the block of job J+1
    auto rangeIndexNextBlock = positionNextJobAfterJ.range(instance);

    // if we change the type of machine
    if (positionJobJ.machines != positionNextJobAfterJ.machines) {
        auto itANextBlockBegin = std::next(A.begin(), static_cast<long>(rangeIndexNextBlock.first));
        auto itANextBlockEnd = std::next(A.begin(), static_cast<long>(rangeIndexNextBlock.second));
        std::replace(itANextBlockBegin, itANextBlockEnd, false, true);
    }
        // else we stay in the same type of machine, but go on new block, so we make all machines available
    else if (positionJobJ.blocks != positionNextJobAfterJ.blocks) {
        auto itAEnd = std::next(A.begin(), static_cast<long>(rangeIndexBlock.second));
        std::replace(A.begin(), itAEnd, false, true);
    }
}

inline void DynamicPrograming::computeListCompletionTimeAndAvailablePosition(int pj, unsigned int n_j, unsigned int n, unsigned int mu_j, std::vector<bool> &A, std::vector<int> &C
                                                                             , std::vector<bool> &combinationAffectationLastBlock
                                                                             , std::multimap<double, std::pair<unsigned int, unsigned int>> &listCjAndAvailablePosition) {

    Position positionOfJob_j = Position(instance, n); // the position where add the j^th jobs on the current solution

    // first, we start from the right, we loop over mu_j - n_j jobs :
    unsigned int indexJob = n - n_j;
    unsigned int lastIndex = n - mu_j;

    while (indexJob < lastIndex) {
        // the position of the job at indexJob
        Position positionOfindexJob = Position(instance, indexJob + 1);

        // defined the index to loop over all machines in the block
        auto rangeIndexMachineInBlock = positionOfindexJob.range(instance);

        for (unsigned int indexMachineInBlock = rangeIndexMachineInBlock.first;
             indexMachineInBlock < rangeIndexMachineInBlock.second; ++indexMachineInBlock) {
            // the speed of the machine
            double speed = indexMachineInBlock < instance->getNbOfHighSpeedMachines() ? instance->getHighSpeed()
                                                                                      : instance->getLowSpeed();
            // check if the machine is available
            if (A[indexMachineInBlock]) {
                double jobCj = double(C[indexMachineInBlock]) / speed;
                // update Cj
                C[indexMachineInBlock] -= pj;
                // pass to the next job
                indexJob++;

                // if the completion time is negative then stop
                if (jobCj < 0) break;
                auto pairCjPositionOnBlockAndMachine = std::make_pair(jobCj, std::make_pair(indexMachineInBlock, positionOfindexJob.getBlockIndex(
                        indexMachineInBlock, instance) - 1));
                // move the created pair
                listCjAndAvailablePosition.insert(std::move(pairCjPositionOnBlockAndMachine));
                //set the machine unavailable
                A[indexMachineInBlock] = false;

            }
        }

        // the position of the next job at indexJob + 1
        Position positionOfNextJobFromindexJob = Position(instance, indexJob + 2);
        updateAvailableMachine(positionOfindexJob, positionOfNextJobFromindexJob, A);

    }

    // Loop on the number of machines in the last block and update A et C according the combination of chosen jobs
    auto rangeIndexLastBlock = positionOfJob_j.range(instance);

    // Then, we go to the last block's jobs
    indexJob = n - mu_j;
    lastIndex = n;
    for (unsigned int indexMachineInBlock = 0 + rangeIndexLastBlock.first; indexMachineInBlock <
                                                                           rangeIndexLastBlock.first +
                                                                           combinationAffectationLastBlock.size(); ++indexMachineInBlock) {
        // the position of the job at indexJob
        Position positionOfindexJob = Position(instance, indexJob + 1);
        // the speed of the machine
        double speed = indexMachineInBlock < instance->getNbOfHighSpeedMachines() ? instance->getHighSpeed()
                                                                                  : instance->getLowSpeed();
        // if we have not scheduled all jobs
        if (indexJob <= lastIndex) {

            // /!\ set the index of combinationAffectationLastBlock between 0 and its size,
            // here indexMachineInBlock can be in range [0, nb High Speed Machine] or [nb High Speed Machine, nb Low Speed Machine]

            // if the index of machine is in combination, i.e. the machine is chosen
            if (combinationAffectationLastBlock[indexMachineInBlock - rangeIndexLastBlock.first]) {
                double jobCj;
                // if the machine is available
                if (A[indexMachineInBlock]) {
                    jobCj = double(C[indexMachineInBlock]) / speed;
                    // update Cj
                    C[indexMachineInBlock] -= pj;
                    // if the completion time is negative then stop
                    if (jobCj < 0) break;
                    auto pairCjPositionOnBlockAndMachine = std::make_pair(jobCj, std::make_pair(indexMachineInBlock, positionOfindexJob.getBlockIndex(
                            indexMachineInBlock, instance) - 1));
                    // move the created pair
                    listCjAndAvailablePosition.insert(std::move(pairCjPositionOnBlockAndMachine));
                    //set the machine unavailable
                    A[indexMachineInBlock] = false;
                    // pass to the next job
                    ++indexJob;
                }
                    // else there is a job on the machine, pass to next machine
                else break;
            }
        }
    }
    // the position of the next job at n + 1
    Position positionOfNextJobFromindexJob = Position(instance, n + 1);
    updateAvailableMachine(positionOfJob_j, positionOfNextJobFromindexJob, A);

}

inline unsigned int DynamicPrograming::computeNumberOfJobsInLastBlock(unsigned int n, unsigned int n_j) {
    // Determine the number of job to put on the last block, called mu_j
    Position positionOfJob_J = Position(instance, n);
    Position positionOfJob_j_minus_n_j = Position(instance, n - n_j);
    unsigned int mu_j;
    // if job j and job n_j are in same machine and same block then mu_j equals n_j
    if (positionOfJob_J == positionOfJob_j_minus_n_j) mu_j = n_j;
    else {
        // Find the first index k such as k < j and f(k) = f(j) and f(k-1) != f(j)
        unsigned int k = n;
        while (Position(instance, k - 1) == positionOfJob_J) {
            k -= 1;
        }
        mu_j = n - k + 1;
    }
    return mu_j;
}

inline void DynamicPrograming::solveSubProblem(std::vector<bool> &A, std::vector<int> &C, std::vector<unsigned int> &indexIdenticalJobs, unsigned int n
                                               , std::vector<DynamicPrograming::State_recursion> &listStateRecursion) {

    /****************************/
    /*      Declaration         */
    /****************************/

    int pj = int(instance->getListJobs()[indexIdenticalJobs[0]].getPi());

    // Determine the number of job to put on the last block, called mu_j
    // with the j^th job that have be added, i.e. n selected job and n_j+1 identical jobs
    unsigned int mu_j = DynamicPrograming::computeNumberOfJobsInLastBlock(n, indexIdenticalJobs.size());

    // declare the binary tree of completion times that is not assigned with the corresponding position, i.e. key is the completion time and the value is (index of machine, index of block)
    std::multimap<double, std::pair<unsigned int, unsigned int>> listCjAndAvailablePosition;

    // declare the copy of the list of jobs. We store it on multimap for O(log n) operations. A key is a pair (due date, number of job)
    std::multimap<unsigned int, unsigned, std::greater<>> copyMapJobs;

    /*****************************/
    /*      Combinations         */
    /*****************************/

    // loop over each combination for mu_j assignments
    // the position where add the j^th jobs, i.e. with n selected jobs
    Position positionOfJob_j = Position(instance, n);
    auto rangeIndexLastBlock = positionOfJob_j.range(instance);
    unsigned int nbOfMachinesLastBlock = rangeIndexLastBlock.second - rangeIndexLastBlock.first;

    // reserve the size of listStateRecursion, by n choose k possibilities
    listStateRecursion.reserve(nChoosek(nbOfMachinesLastBlock, mu_j));

    std::vector<bool> mu_jAssigment = std::vector<bool>(nbOfMachinesLastBlock - mu_j);
    mu_jAssigment.resize(nbOfMachinesLastBlock, true);

    do {
        // the weighted for tardy jobs
        double minSumWj = 0.0;

        listCjAndAvailablePosition.clear();

        // declare the queue for max weighted early job
        auto queueEarlyJob = std::priority_queue<std::pair<double, unsigned int>, std::vector<std::pair<double, unsigned int>>>();
        // A copy of vector A and C
        std::vector<bool> copyA(A);
        std::vector<int> copyC(C);

        // Compute the list of completion times and the available corresponding positions
        // with the j^th job that have be added, i.e. n selected job and n_j+1 identical jobs
        computeListCompletionTimeAndAvailablePosition(pj, indexIdenticalJobs.size(), n, mu_j, copyA, copyC, mu_jAssigment, listCjAndAvailablePosition);

        // if we have not the same completion time as the number of identical job then break
        if (listCjAndAvailablePosition.size() != indexIdenticalJobs.size()) continue;

        // copy the list of jobs in map in order to remove the added jobs
        for (unsigned int j: indexIdenticalJobs) copyMapJobs.emplace(instance->getListJobs()[j].getDi(), j);

        /********************************************/
        /*      Compute weighted tardy Jobs         */
        /********************************************/

        // use iterator
        auto itCompletionTime = listCjAndAvailablePosition.rbegin();
        auto itJobs = copyMapJobs.begin();

        // compute the set of all early jobs for each completion times, we start from the end because we remove element
        while (itCompletionTime != listCjAndAvailablePosition.rend() && itJobs != copyMapJobs.end()) {
            // If the job can be schedule on time,i.e. dj > Ci, we add it
            if (double(itJobs->first) >= itCompletionTime->first) {
                queueEarlyJob.emplace(instance->getListJobs()[itJobs->second].getWi(), itJobs->second);
                ++itJobs;
            } else if (!queueEarlyJob.empty()) {
                // else if we have some jobs that early (means there are in the queue), so we affect the one with max weight
                affectEarlyJobsToPosition(queueEarlyJob, copyMapJobs);
                itCompletionTime = std::reverse_iterator(
                        listCjAndAvailablePosition.erase(std::next(itCompletionTime).base()));
            } else {
                // we pass to the next Cj
                ++itCompletionTime;
            }
        }

        // We clear out the queue if it not empty
        while (!queueEarlyJob.empty() && itCompletionTime != listCjAndAvailablePosition.rend()) {
            affectEarlyJobsToPosition(queueEarlyJob, copyMapJobs);
            itCompletionTime = std::reverse_iterator(
                    listCjAndAvailablePosition.erase(std::next(itCompletionTime).base()));
        }

        // now, we schedule not assigned jobs,i.e. the late ones, so iterate through the map of jobs
        itCompletionTime = listCjAndAvailablePosition.rbegin();
        itJobs = copyMapJobs.begin();
        unsigned int nbRemainderJobsToAffect = copyMapJobs.size();
        for (unsigned int i = 0; i < nbRemainderJobsToAffect; ++i) {
            minSumWj += instance->getListJobs()[itJobs->second].getWi();
            itCompletionTime = std::reverse_iterator(
                    listCjAndAvailablePosition.erase(std::next(itCompletionTime).base()));
            itJobs = copyMapJobs.erase(itJobs);
        }

        /******************************/
        /*      Schedule Jobs         */
        /******************************/


        // add the state recursion
        State_recursion newStateRecursion(std::move(copyA), std::move(copyC), minSumWj);
        listStateRecursion.push_back(std::move(newStateRecursion));

    } while (std::next_permutation(mu_jAssigment.begin(), mu_jAssigment.end()));

}

/**************************************/
/*     DYNAMIC PROGRAMING METHODS     */
/**************************************/



inline double
DynamicPrograming::solveDynamicPrograming(std::vector<bool> A, std::vector<int> C, unsigned int j, std::vector<unsigned int> indexIdenticalJobs, unsigned int n) {

    /**************************/
    /*      MEMOIZATION       */
    /**************************/

    std::vector<int> key;
    key.reserve(A.size() + C.size() + 2);

    // sort on both type machines
    for (char type = 0; type < 2; type++) {
        size_t nbElement = type == 0 ? instance->getNbOfHighSpeedMachines() : instance->getNbOfLowSpeedMachines();
        unsigned int startIndex = type == 0 ? 0 : instance->getNbOfHighSpeedMachines();
        // permutation to apply on A after sorted C
        std::vector<std::size_t> permutation(nbElement);
        std::iota(permutation.begin(), permutation.end(), 0);
        std::sort(permutation.begin(), permutation.end(), [&](std::size_t i, std::size_t j) { return C[i + startIndex] > C[j + startIndex]; });

        // add A element for type machine
        for (auto const &perm: permutation) {
            key.push_back(A[startIndex + perm]);
        }
        // add C element for type machine
        for (auto const &perm: permutation) {
            key.push_back(C[startIndex + perm]);
        }

    }

    // add identical jobs on key
    for (const auto identicalJob: indexIdenticalJobs) {
        key.push_back(identicalJob);
    }

    // append j and n_j
    key.push_back(int(j));
    key.push_back(int(n));

    double minValue = std::numeric_limits<double>::infinity();

    // if C is not positive vector
    if (std::any_of(C.cbegin(), C.cend(), [](const int &i) { return i < 0; })) {
        memoization[key] = std::numeric_limits<double>::infinity();
        return memoization[key];
    } else
        //if we have select more jobs than expected, and we have not identical jobs to schedule
    if (n == instance->getNbToSelectJob() && indexIdenticalJobs.empty() &&
        std::all_of(C.cbegin(), C.cend(), [](const int &i) { return i == 0; })) {
        memoization[key] = 0.0;
        return memoization[key];
    } else
        // if we have explored all jobs
    if (j >= instance->getNbJobs()) {
        memoization[key] = std::numeric_limits<double>::infinity();
        return memoization[key];
    }


    /************************/
    /*      VARIABLES       */
    /************************/


    // position of job j
    Position pos_j = Position(instance, n + 1);
    // position of job j
    Position pos_after_j = Position(instance, n + 2);

    //speeds
    auto v0 = instance->getLowSpeed();
    auto vMax = instance->getHighSpeed();

    /************************/
    /*      RECURSION       */
    /************************/

    auto itMemo = memoization.find(key);
    // if we have already solved this state
    if (itMemo != memoization.end()) {
        minValue = itMemo->second;
    } else {

        /*****************************/
        /*      Identical Jobs       */
        /*****************************/


        // if j and j+1 are identical jobs
        if (j < (instance->getNbJobs() - 1) &&
            int(instance->getListJobs()[j].getPi()) == int(instance->getListJobs()[j + 1].getPi())) {
            indexIdenticalJobs.push_back(j);
            // add the next job
            minValue = std::min(minValue, solveDynamicPrograming(A, C, j + 1, indexIdenticalJobs, n + 1));
            // remove the identical job
            indexIdenticalJobs.pop_back();
        }
            // else if we have already identical jobs, that means j is identical as the jobs before
        else if (!indexIdenticalJobs.empty()) {

            double minValue1 = minValue;
            double minValue2 = minValue;

            // the list of state of the recursion generate by solve the sub-problem
            std::vector<DynamicPrograming::State_recursion> listStateRecursion;

            // the list of identicalJobs
            indexIdenticalJobs.push_back(j);

            // a new list of identical jobs
            std::vector<unsigned int> newIdenticalJobs;
            newIdenticalJobs.reserve(instance->getNbToSelectJob());

            // solve the sub-problem with the list of identical jobs (so n+1 selected jobs)
            solveSubProblem(A, C, indexIdenticalJobs, n + 1, listStateRecursion);
            for (auto &state: listStateRecursion) {
                minValue1 = std::min(minValue1, solveDynamicPrograming(state.A, state.C, j + 1, newIdenticalJobs, n + 1) +
                                                state.minSumWj);
            }

            //we don't take j
            indexIdenticalJobs.pop_back();
            // clear list of recursion state
            listStateRecursion.clear();

            // solve the sub-problem with the list of identical jobs (so n selected jobs, because we don't take j)
            solveSubProblem(A, C, indexIdenticalJobs, n, listStateRecursion);
            for (auto &state: listStateRecursion) {
                minValue2 = std::min(minValue2, solveDynamicPrograming(state.A, state.C, j + 1, newIdenticalJobs, n) +
                                                state.minSumWj);
            }
            listStateRecursion.clear();

            // reduce n by the number of identical jobs
            n -= indexIdenticalJobs.size();
            indexIdenticalJobs.clear();
            minValue = std::min(minValue1, minValue2);

        } else {

            /****************************/
            /*      Schedule Jobs       */
            /****************************/

            // add the job j to the schedule, we need to take the block structure properties in consideration

            // if we are on low speed machine
            if (pos_j.machines[1] == -1 && pos_j.machines[0] == 1) {
                //make available all high speed machines
                std::replace(A.begin(), A.begin() + instance->getNbOfHighSpeedMachines(), false, true);

                // try to affect the job
                for (unsigned int i = instance->getNbOfHighSpeedMachines(); i < instance->getNbMachines(); ++i) {
                    // if the machine is available
                    if (A[i]) {
                        // make the machine unavailable and change the completion time
                        A[i] = false;
                        int completionTimeBeforeAddJob = C[i];
                        C[i] -= int(instance->getListJobs()[j].getPi());

                        // go next job
                        minValue = std::min(minValue, solveDynamicPrograming(A, C, j + 1, indexIdenticalJobs, n + 1) +
                                                      weightedTobeLate(completionTimeBeforeAddJob, v0, instance->getListJobs()[j]));

                        // undo the modification
                        A[i] = true;
                        C[i] += int(instance->getListJobs()[j].getPi());
                    }
                }
            }
                // if we are on high speed machine and the next job go on low speed machine
            else if (pos_j.machines[1] == -1 && pos_j.machines[0] == 0 && pos_after_j.machines[1] == -1 &&
                     pos_after_j.machines[0] == 1) {
                //make available all low speed machines
                std::replace(A.begin() + instance->getNbOfHighSpeedMachines(), A.end(), false, true);

                // try to affect the job
                for (unsigned int i = 0; i < instance->getNbOfHighSpeedMachines(); ++i) {
                    // if the machine is available
                    if (A[i]) {
                        // make the machine unavailable and change the completion time
                        A[i] = false;
                        int completionTimeBeforeAddJob = C[i];
                        C[i] -= int(instance->getListJobs()[j].getPi());

                        // go next job
                        minValue = std::min(minValue, solveDynamicPrograming(A, C, j + 1, indexIdenticalJobs, n + 1) +
                                                      weightedTobeLate(completionTimeBeforeAddJob, vMax, instance->getListJobs()[j]));

                        // undo the modification
                        A[i] = true;
                        C[i] += int(instance->getListJobs()[j].getPi());
                    }
                }

            }
                // if we are on high speed machine and next job go on new block of high speed machine
            else if (pos_j.blocks != pos_after_j.blocks && pos_j.machines == pos_after_j.machines &&
                     pos_j.machines[1] == -1 && pos_j.machines[0] == 1) {
                // use a copy of A
                std::vector<bool> newA(A);
                //make available all high speed machines
                std::replace(newA.begin(), newA.begin() + instance->getNbOfHighSpeedMachines(), false, true);
                //make unavailable all low speed machine
                std::replace(newA.begin() + instance->getNbOfHighSpeedMachines(), newA.end(), true, false);

                // try to affect the job
                for (unsigned int i = 0; i < instance->getNbOfHighSpeedMachines(); ++i) {
                    // if the machine is available
                    if (A[i]) {
                        // make the machine unavailable and change the completion time
                        A[i] = false;
                        int completionTimeBeforeAddJob = C[i];
                        C[i] -= int(instance->getListJobs()[j].getPi());
                        // go next job
                        minValue = std::min(minValue, solveDynamicPrograming(newA, C, j + 1, indexIdenticalJobs, n + 1) +
                                                      weightedTobeLate(completionTimeBeforeAddJob, vMax, instance->getListJobs()[j]));

                        // undo the modification
                        A[i] = true;
                        C[i] += int(instance->getListJobs()[j].getPi());
                    }
                }
            }
                // if we are on high speed machine and next job go on high speed and low speed machine
            else if (pos_j.machines[1] == -1 && pos_j.machines[0] == 0 && pos_after_j.machines[1] != -1) {
                // use a copy of A
                std::vector<bool> newA(A);
                //make available all machines
                std::replace(newA.begin(), newA.end(), false, true);

                // try to affect the job
                for (unsigned int i = 0; i < instance->getNbOfHighSpeedMachines(); ++i) {
                    // if the machine is available
                    if (A[i]) {
                        // make the machine unavailable and change the completion time
                        A[i] = false;
                        int completionTimeBeforeAddJob = C[i];
                        C[i] -= int(instance->getListJobs()[j].getPi());

                        // go next job
                        minValue = std::min(minValue, solveDynamicPrograming(newA, C, j + 1, indexIdenticalJobs, n + 1) +
                                                      weightedTobeLate(completionTimeBeforeAddJob, vMax, instance->getListJobs()[j]));

                        // undo the modification
                        A[i] = true;
                        C[i] += int(instance->getListJobs()[j].getPi());
                    }
                }
            }
                // if we are on high and low speed machines and the next job go on high speed machines
            else if (pos_j.machines[1] != -1 && pos_after_j.machines[1] == -1 && pos_after_j.machines[0] == 0) {
                // use a copy of A
                std::vector<bool> newA(A);
                //make available all high speed machines
                std::replace(newA.begin(), newA.begin() + instance->getNbOfHighSpeedMachines(), false, true);

                // try to affect the job
                for (unsigned int i = 0; i < instance->getNbMachines(); ++i) {
                    // if the machine is available
                    if (A[i]) {
                        // make the machine unavailable and change the completion time
                        A[i] = false;
                        int completionTimeBeforeAddJob = C[i];
                        C[i] -= int(instance->getListJobs()[j].getPi());

                        // go next job
                        double speed = i < instance->getNbOfHighSpeedMachines() ? vMax : v0;
                        minValue = std::min(minValue, solveDynamicPrograming(newA, C, j + 1, indexIdenticalJobs, n + 1) +
                                                      weightedTobeLate(completionTimeBeforeAddJob, speed, instance->getListJobs()[j]));

                        // undo the modification
                        A[i] = true;
                        C[i] += int(instance->getListJobs()[j].getPi());
                    }
                }
            }
                // else, we go one on recursion from the range of the position of job j
            else {
                // try to affect the job
                auto rangeJobJ = pos_j.range(instance);
                for (unsigned int i = rangeJobJ.first; i < rangeJobJ.second; ++i) {
                    // if the machine is available
                    if (A[i]) {
                        // make the machine unavailable and change the completion time
                        A[i] = false;
                        int completionTimeBeforeAddJob = C[i];
                        C[i] -= int(instance->getListJobs()[j].getPi());
                        // go next job
                        double speed = i < instance->getNbOfHighSpeedMachines() ? vMax : v0;
                        minValue = std::min(minValue, solveDynamicPrograming(A, C, j + 1, indexIdenticalJobs, n + 1) +
                                                      weightedTobeLate(completionTimeBeforeAddJob, speed, instance->getListJobs()[j]));

                        // undo the modification
                        A[i] = true;
                        C[i] += int(instance->getListJobs()[j].getPi());
                    }
                }
            }
        }
    }


    // add the next job, identical jobs does not change because we have not added job
    minValue = std::min(minValue, solveDynamicPrograming(A, C, j + 1, indexIdenticalJobs, n));

    memoization[key] = minValue;
    return minValue;
}


#endif //BILEVEL_SCHEDULING_DYNAMICPROGRAMING_H
