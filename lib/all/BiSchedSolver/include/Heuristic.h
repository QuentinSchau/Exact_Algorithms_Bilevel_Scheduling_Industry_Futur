// Copyright (C) 2024
// Laboratoire d'Informatique Fondamentale et Appliquée de Tours, Tours, France
//
// DIGEP, Politecnico di Torino, Corso Duca degli Abruzzi 24, Torino, Italy
// This file is part of bilevel-scheduling.
//
// bilevel-scheduling is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// bilevel-scheduling is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with bilevel-scheduling. If not, see <https://www.gnu.org/licenses/>.

//
// Created by schau on 1/22/25.
//

#ifndef BILEVEL_SCHEDULING_HEURISTIC_H
#define BILEVEL_SCHEDULING_HEURISTIC_H


#include "ISolver.h"
#include "HungarianAlgorithm.h"

class Heuristic : public ISolver {
private:
    bool use_predictor = false;
    bool trainWeight = false;
    double alpha0 = 1.0;
    double alpha1 = 1.0;
    double alpha2 = 1.0;
    double alpha3 = 1.0;

public:

    /************************/
    /*      CONSTRUCTOR     */
    /************************/

    explicit Heuristic();

    explicit Heuristic(Instance *instance);

    explicit Heuristic(Instance *instance, nlohmann::json &object);


    /***********************/
    /*      DESTRUCTOR     */
    /***********************/

    ~Heuristic() override;

    void solve() override;

    /**********************/
    /*      Heuristic     */
    /**********************/

    /** Local Search **/

    void localSearch();


    std::pair<unsigned int, unsigned int> exploreNeighbourhood_1_OPT(std::vector<Job> &initialListJobs, std::vector<bool> &alreadySelectedJobs, Solution *bestSolInNeighbourhood);

    /**
     * Method that upgrade the solution given in parameter with heuristic.
     * @param sol The solution to improve.
     * @param listJobsAvailable The list of jobs from which we try to fill the block that we have freed.
     */
    void upgradeSolutionWithHeuristic(Solution &sol, const std::vector<Job> &listJobsAvailable);

    /**
     * Method that release jobs in block.
     * @param blockStruct The block structure from where we release jobs.
     * @param indexBlock The index of block to release.
     * @param listJobsAvailable The list of jobs from which we try to fill the block that we have freed.
     */
    void freeAndAssignmentBlock(Solution::BlockStructure &blockStruct, unsigned int indexBlock, const std::vector<Job> *listOfJobsAvailable);

    /**
     * Method that computes the cost of the schedule at each block and for each machines.
     * @param costAtEachBlock The matrix of cost for each machine and each block. costAtEachBlock[i][k] is the cost of
     * the machine 'i' when the job ended at the block 'k'. The last element of 'costAtEachBlock[i]' is the cost of the
     * machine schedule.
     */
    static void computeAllWeightsAtEachBlock(std::vector<std::vector<double>> &costAtEachBlock, Solution::BlockStructure &blockStruct);


    /*********************/
    /*      Features     */
    /*********************/

    /**
     * Method that computes the average load on machine at position k, where the load refers to the average processing
     * times of tasks on machines at position k divided by sum of machine's speed.
     * @param solution The a solution of the problem.
     * @param positionInMachine The position in the machine.
     * @return The average load at a given position.
     */
    double computeAvgLoad(Solution &sol, unsigned int positionInMachine);

    /**
     * Method that computes the average due date in a given position in all machines.
     * @param solution The a solution of the problem.
     * @param positionInMachine The position in the machine.
     * @return The average due date in all machine at a given position.
     */
    double computeAvgDueDate(Solution &sol, unsigned int positionInMachine);

    /**
     * Method that computes the average weight in a given position in all machines.
     * @param solution The a solution of the problem.
     * @param positionInMachine The position in the machine.
     * @return The average weight in all machine at a given position.
     */
    double computeAvgWeight(Solution &sol, unsigned int positionInMachine);

    /**
     * Method that computes a feature from a solution. This method loops through each available position in machine schedule,
     * i.e., for a machine i we have M_i = k_1, k_2, ..., k_b where k is the position.
     * Some machines have more or less positions depending on the block structure.
     *
     * For a given position k, this method computes the sum of all average loads on positions before k, denoted as load_k,
     * which estimates the average starting time on position k. This method also computes the average of due date and weight
     * in all machines at the given position k.
     *
     * Then, the method returns the sum of average weights for all positions k where the average load is greater than the average due
     * date.
     *
     * @param sol The solution used to compute the feature.
     * @return The sum of average weights for all positions k where the average load was greater than the average due date
     */
    double feature1(Solution &sol);

    /**
     * Method that computes a feature from a solution. This method loops through each available position in machine schedule,
     * i.e., for a machine i we have M_i = k_1, k_2, ..., k_b where k is the position.
     * Some machines have more or less positions depending on the block structure.
     *
     * For a given position k, this method computes the sum of all average loads on positions before k, denoted as load_k.
     * This loads is an estimation of the average starting time on position k.
     *
     * By using this estimation, this method computes for all job j on position k, the difference D_{k,j} = (load_k + p_j / V_{max}) - d_j.
     * If for all jobs j in position k, the average of D_{k,j} is positive, then all differences will be used to compute the features;
     * otherwise we pass to the next position.
     *
     * Thus, at the end, this method returns the average of differences D_{k,j} where some D_{k,j} could be negative individually,
     * but the average on a position k is positive.

     * @param sol The solution used to compute the feature.
     * @return The average of difference between an estimation of the starting time of a job and its due date.
     */
    double feature2(Solution &sol);

    /**
     * Method that computes a feature from a solution. This method loops through each available position in machine schedule,
     * i.e., for a machine i we have M_i = k_1, k_2, ..., k_b where k is the position.
     * Some machines have more or less positions depending on the block structure.
     *
     * For a given position k, this method computes the sum of all average loads on positions before k, denoted as load_k.
     * This load is an estimation of the average starting time on position k.
     *
     * By using this estimation, this method computes for all job j on position k, the difference D_{k,j} = (load_k + p_j / V_{max}) - d_j.
     *
     * This method returns the sum of weights of all jobs with positive difference D_{k,j}
     *
     * @param sol The solution used to compute the feature.
     * @return The weighted number of jobs with positive difference D_{k,j}
     */
    double feature3(Solution &sol);


    /********************/
    /*      SETTER      */
    /********************/

    void setUsePredictor(bool usePredictor) {
        use_predictor = usePredictor;
    }

    void setWeight(std::vector<double> newWeights) {
        if (!newWeights.empty()) {
            alpha0 = newWeights[0];
        }
        if (newWeights.size() >= 2) {
            alpha1 = newWeights[1];
        }
        if (newWeights.size() >= 3) {
            alpha2 = newWeights[1];
        }
        if (newWeights.size() >= 4) {
            alpha3 = newWeights[2];
        }

    };

    void setTrainWeight(bool trainWeight) {
        Heuristic::trainWeight = trainWeight;
    }
    /********************/
    /*      GETTER      */
    /********************/

    bool isUsePredictor() const {
        return use_predictor;
    }

    /**
     * Method that saves the result of the instance in a file
     * @param fileOutputName The name of the file
     * @param outputFile The stream of the file
     */
    void printOutput(std::string &fileOutputName, std::ofstream &outputFile);
};

inline double
Heuristic::computeAvgLoad(Solution &sol, const unsigned int positionInMachine) {
    double sumC_j = 0.0;
    double sumSpeed = 0.0;
    // get the sum of completion of all jobs in position
    for (unsigned int indexMachine = 0; indexMachine < instance->getNbMachines(); indexMachine++) {
        auto &machine = sol[indexMachine];
        if (positionInMachine < machine.size()) {
            sumSpeed += machine.getSpeed();
            auto &job = machine[positionInMachine];
            sumC_j += job.getPi();
        }
    }
    return sumC_j / sumSpeed;
}

inline double Heuristic::feature1(Solution &sol) {

    double sumCharge = 0.0;
    double estimationWj = 0.0;
    // loop over each block
    for (unsigned int positionInMachine = 0; positionInMachine < instance->getE().size(); positionInMachine++) {
        sumCharge += computeAvgLoad(sol, positionInMachine);
        double avgDueDate = computeAvgDueDate(sol, positionInMachine);
        double avgWeight = computeAvgWeight(sol, positionInMachine);
        if (std::isless(avgDueDate, sumCharge)) estimationWj += avgWeight;
    }
    return estimationWj;
}

inline double Heuristic::computeAvgDueDate(Solution &sol, unsigned int positionInMachine) {
    double sumD_j = 0.0;
    double nbMachines = 0.0;
    // get the sum of due date of all jobs in position
    for (unsigned int indexMachine = 0; indexMachine < instance->getNbMachines(); indexMachine++) {
        auto &machine = sol[indexMachine];
        if (positionInMachine < machine.size()) {
            nbMachines += 1.0;
            auto &job = machine[positionInMachine];
            sumD_j += job.getDi();
        }
    }
    return sumD_j / nbMachines;
}

inline double Heuristic::computeAvgWeight(Solution &sol, unsigned int positionInMachine) {
    double sumW_j = 0.0;
    double nbMachines = 0.0;
    // get the sum of completion of all jobs in the block
    for (unsigned int indexMachine = 0; indexMachine < instance->getNbMachines(); indexMachine++) {
        auto &machine = sol[indexMachine];
        if (positionInMachine < machine.size()) {
            nbMachines += 1.0;
            auto &job = machine[positionInMachine];
            sumW_j += job.getWi();
        }
    }
    return sumW_j / nbMachines;
}

inline double Heuristic::feature2(Solution &sol) {
    auto blockStructure = sol.toBlockStruct(instance);
    double highSpeed = instance->getHighSpeed();
    double sumDiff_pj_dj = 0.0;
    double nbJobs = 0.0;
    double sumCharge = 0.0;
    // loop over each block
    for (unsigned int positionInMachine = 0; positionInMachine < instance->getE().size(); positionInMachine++) {
        sumCharge += (positionInMachine == 0) ? 0.0 : computeAvgLoad(sol, positionInMachine - 1);
        double sum_diff_in_position = 0.0;
        double nbJobsInBlock = 0.0;
        for (unsigned int indexMachine = 0; indexMachine < instance->getNbMachines(); indexMachine++) {
            auto &machine = sol[indexMachine];
            if (positionInMachine < machine.size()) {
                nbJobsInBlock += 1.0;
                auto &job = machine[positionInMachine];
                double diff_pj_dj = (sumCharge + job.getPi() / highSpeed - job.getDi());
                sum_diff_in_position += diff_pj_dj;
            }
        }
        if (sum_diff_in_position / nbJobsInBlock >= -EPSILON) {
            sumDiff_pj_dj += sum_diff_in_position;
            nbJobs += nbJobsInBlock;
        }
    }
    // if the nbJobs is 0.0, that's means we don't have any positive sum of diff, so let be equal to - max pj/Vmax
    if (nbJobs <= EPSILON)
        return -instance->getMaxPj() / highSpeed;
    else
        return sumDiff_pj_dj / nbJobs;
}

inline double Heuristic::feature3(Solution &sol) {
    auto blockStructure = sol.toBlockStruct(instance);
    double highSpeed = instance->getHighSpeed();
    double sum_weight_late_jobs = 0.0;
    double sumCharge = 0.0;
    // loop over each block
    for (unsigned int positionInMachine = 0; positionInMachine < instance->getE().size(); positionInMachine++) {
        sumCharge += (positionInMachine == 0) ? 0.0 : computeAvgLoad(sol, positionInMachine - 1);
        for (unsigned int indexMachine = 0; indexMachine < instance->getNbMachines(); indexMachine++) {
            auto &machine = sol[indexMachine];
            if (positionInMachine < machine.size()) {
                auto &job = machine[positionInMachine];
                double diff_pj_dj = (sumCharge + job.getPi() / highSpeed - job.getDi());
                if (diff_pj_dj >= -EPSILON) {
                    sum_weight_late_jobs += job.getWi();
                }
            }
        }
    }
    return sum_weight_late_jobs;
}

inline void Heuristic::upgradeSolutionWithHeuristic(Solution &sol, const std::vector<Job> &listJobsAvailable) {
    // declare the binary tree of completion times that is not assigned with the corresponding position,
    // i.e. key is the completion time and the value is (index of machine, index of block). This is sorted in increasing order
    typedef std::multimap<double, std::pair<unsigned int, unsigned int>> TreeCj;

    double objValBefore = sol.getSumWjUj(); // keep the value of solution before apply the heuristic
    auto blockStruct = sol.toBlockStruct(instance);
    // try to free and re-assign jobs from the RIGHT to LEFT
    unsigned int indexBlock = 0;
    // try now from the LEFT to RIGHT
    while (indexBlock < instance->getE().size()) {
        freeAndAssignmentBlock(blockStruct, indexBlock, &listJobsAvailable);
        ++indexBlock;
    }

    // first we determine all identical jobs in the schedule and we minimize the weight optimally
    // use a vector of (index in list grouped jobs, tree with completion time) we will solve optimally all element
    // of this vector
    std::vector<std::pair<unsigned int, TreeCj>> listIdenticalJobAndTheirTree;
    indexBlock = instance->getE().size();
    while (indexBlock) {
        --indexBlock;
        for (auto [indexMachine, indexBlockInStruct]: instance->getE()[indexBlock]) {
            auto job = blockStruct[indexMachine][indexBlockInStruct].first;
            // check if a job is scheduled
            if (job != nullptr) {
                // check if there exist identical jobs
                for (unsigned int indexInListOfGroupedJobs = 0;
                     indexInListOfGroupedJobs < instance->getListGrpJobs().size(); ++indexInListOfGroupedJobs) {
                    auto groupJobs = instance->getListGrpJobs()[indexInListOfGroupedJobs];
                    auto predIdenticalGroup = [job](const Job &jobInGroup) { return jobInGroup.getPi() == job->getPi(); };
                    // if the job in the schedule is identical job
                    if (std::find_if(groupJobs.begin(), groupJobs.end(), predIdenticalGroup) != groupJobs.end() && groupJobs.size() > 1) {
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
                        break;
                    }
                }
            }
        }
    }
    for (unsigned int indexLoopTreeCj = 0;
         indexLoopTreeCj < listIdenticalJobAndTheirTree.size(); ++indexLoopTreeCj) {
        auto &[indexListOfIdenticalJobs, listCjAndAvailablePosition] = listIdenticalJobAndTheirTree[indexLoopTreeCj];
        // use a copy of the list of grouped jobs, because it's modifying by the method
        std::vector<Job> listOfIdenticalJobs(instance->getListGrpJobs()[indexListOfIdenticalJobs]);
        solveProblemWithIdenticalJobs(nullptr, &blockStruct, listCjAndAvailablePosition, listOfIdenticalJobs);
    }
    sol.fromBlockStruct(blockStruct);
    if (verbose >= 3) std::cout << "UB before upgrade: " << objValBefore << " after: " << sol.getSumWjUj() << std::endl;
}

inline void Heuristic::freeAndAssignmentBlock(Solution::BlockStructure &blockStruct, unsigned int indexBlock, const std::vector<Job> *listOfJobsAvailable) {
    std::vector<std::vector<double>> costAtEachBlock;
    computeAllWeightsAtEachBlock(costAtEachBlock, blockStruct);
    auto &E = instance->getE();
    // get the max pj from the block before
    double maxPj = -1.0, minPj = std::numeric_limits<double>::infinity();
    if (indexBlock > 0) {
        for (auto &location: E[indexBlock - 1]) {
            auto [job, cj] = blockStruct[location.first][location.second];
            if (job != nullptr) {
                maxPj = std::max(maxPj, job->getPi());
            }
        }
    }
    // get the min pj from the block after
    if (indexBlock < E.size() - 1) {
        for (auto &location: E[indexBlock + 1]) {
            auto [job, cj] = blockStruct[location.first][location.second];
            if (job != nullptr) {
                minPj = std::min(minPj, job->getPi());
            }
        }
    }

    // compute the list of jobs that can be schedule on this block
    std::vector<unsigned int> listAvailableIndexJob;
    listAvailableIndexJob.reserve(instance->getNbJobs());

    if (listOfJobsAvailable) {
        // add jobs that can be scheduled on the block
        for (auto &job: *listOfJobsAvailable) {
            if (job.getPi() > maxPj
                && job.getPi() < minPj
                && std::find(listAvailableIndexJob.begin(), listAvailableIndexJob.end(), job.getIndex()) ==
                   listAvailableIndexJob.end()) {
                listAvailableIndexJob.emplace_back(job.getIndex());
            }
        }
    }
    // add the job from the block
    for (auto &location: E[indexBlock]) {
        auto [job, cj] = blockStruct[location.first][location.second];
        if (job != nullptr) {
            // check if we have not added the job
            auto itFindJob = std::find_if(listAvailableIndexJob.begin(), listAvailableIndexJob.end(), [job](const unsigned int indexJob) { return indexJob == job->getIndex(); });
            if (itFindJob == listAvailableIndexJob.end())
                listAvailableIndexJob.push_back(job->getIndex());
        }
    }

    //compute the matrix of cost
    std::vector<std::vector<double>> costMatrix;
    for (auto &[indexMachine, indexBlockInStruct]: E[indexBlock]) {
        double speed = (indexMachine < instance->getNbOfHighSpeedMachines()) ? instance->getHighSpeed()
                                                                             : instance->getLowSpeed();
        std::vector<double> costOfMachine; // the cost of the machine starting from indexBlock

        // loop over each jobs
        for (unsigned int indexJob: listAvailableIndexJob) {
            // get the completion time of the location before
            double completionTime = (indexBlockInStruct == 0) ? 0.0 : blockStruct[indexMachine][indexBlockInStruct - 1].second;
            // adjust the completion time with the current job
            auto &job = instance->getListJobs()[indexJob];
            completionTime += job.getPi() / speed;

            // initiate the cost of schedule with the sum of costs of others machines
            double costSchedule = 0.0;
            for (unsigned int indexOtherMachine = 0; indexOtherMachine < blockStruct.size(); ++indexOtherMachine) {
                // add the cost of the machine schedule of the block before
                if (indexOtherMachine == indexMachine) {
                    costSchedule += (indexBlock > 0) ? costAtEachBlock[indexMachine][indexBlock - 1] : 0.0;
                }// else add the total cost of the other machine schedule
                else
                    costSchedule += costAtEachBlock[indexOtherMachine].back();
            }

            // add the cost of the job if is late
            double costJob = (job.getDi() < completionTime) ? job.getWi() : 0.0;
            costSchedule += costJob;
            // compute the cost of the machine schedule because we change the completion time of other next blocks
            // So loop over the next block
            unsigned int indexNextBlock = indexBlock + 1;
            auto predFindIndexMachine = [&indexMachine](std::pair<unsigned int, unsigned int> locationInBlock) { return locationInBlock.first == indexMachine; };
            while (indexNextBlock < E.size()) {
                // if a job is schedule on the machine in the next block
                if (std::find_if(E[indexNextBlock].begin(), E[indexNextBlock].end(), predFindIndexMachine) != E[indexNextBlock].end()) {
                    auto &[_, indexNextBlockInStruct] = E[indexNextBlock][indexMachine];
                    auto nextJob = blockStruct[indexMachine][indexNextBlockInStruct].first;
                    completionTime += nextJob->getPi() / speed;
                    costJob = (nextJob->getDi() < completionTime) ? nextJob->getWi() : 0.0;
                    costSchedule += costJob;
                }
                ++indexNextBlock;
            }

            costOfMachine.emplace_back(costSchedule);
        }
        costMatrix.push_back(std::move(costOfMachine));
    }

    std::vector<int> assignment;
    auto H = HungarianAlgorithm();
    H.Solve(costMatrix, assignment);

    // if we are in the first block, then we have to select a subset of the assignment
    if (indexBlock == 0) {
        // the number of jobs that must be scheduled on the block
        unsigned int numJobsToScheduleOnBlock = instance->getNbJobsToScheduleOnFirstBlock();
        // create a list of pair (cost, indexMachine) given by the assignment
        std::vector<std::pair<double, unsigned int>> listPairCostIndexMachine;
        for (unsigned int indexLoopAssignment = 0; indexLoopAssignment < assignment.size(); ++indexLoopAssignment) {
            if (assignment[indexLoopAssignment] >= 0)
                listPairCostIndexMachine.emplace_back(costMatrix[indexLoopAssignment][assignment[indexLoopAssignment]], indexLoopAssignment);
        }
        // clear the block
        for (auto [indexMachineToClear, _]: E[0]) {
            blockStruct[indexMachineToClear][0].first = nullptr;
            blockStruct[indexMachineToClear][0].second = 0.0;
        }

        // get the numJobsToScheduleOnBlock smallest value
        std::nth_element(listPairCostIndexMachine.begin(), listPairCostIndexMachine.begin() + numJobsToScheduleOnBlock, listPairCostIndexMachine.end());
        for (unsigned int indexLoopListPair = 0; indexLoopListPair < numJobsToScheduleOnBlock; ++indexLoopListPair) {
            unsigned int indexMachine = listPairCostIndexMachine[indexLoopListPair].second;
            auto [_, indexBlockInStruct] = E[indexBlock][indexMachine];
            unsigned int indexAssignJob = listAvailableIndexJob[assignment[indexMachine]];
            // change the jobs of the structure
            blockStruct[indexMachine][indexBlockInStruct].first = &instance->getListJobs()[indexAssignJob];
            // update the completion time of the structure and for the next block
            double completionTime = (indexBlockInStruct == 0) ? 0.0 : blockStruct[indexMachine][indexBlockInStruct -
                                                                                                1].second;
            double speed = (indexMachine < instance->getNbOfHighSpeedMachines()) ? instance->getHighSpeed()
                                                                                 : instance->getLowSpeed();
            unsigned int indexNextBlock = indexBlock;
            auto predFindIndexMachine = [&indexMachine](std::pair<unsigned int, unsigned int> locationInBlock) { return locationInBlock.first == indexMachine; };
            while (indexNextBlock < E.size()) {
                // if a job is schedule on the machine in the next block
                if (std::find_if(E[indexNextBlock].begin(), E[indexNextBlock].end(), predFindIndexMachine) != E[indexNextBlock].end()) {
                    auto &[_, indexNextBlockInStruct] = E[indexNextBlock][indexMachine];
                    auto nextJob = blockStruct[indexMachine][indexNextBlockInStruct].first;
                    completionTime += nextJob->getPi() / speed;
                    blockStruct[indexMachine][indexNextBlockInStruct].second = completionTime;
                }
                ++indexNextBlock;
            }
        }
    }// apply the new assigment
    else {
        for (unsigned int indexMachine = 0; indexMachine < E[indexBlock].size(); indexMachine++) {
            auto [_, indexBlockInStruct] = E[indexBlock][indexMachine];
            unsigned int indexAssignJob = listAvailableIndexJob[assignment[indexMachine]];
            // change the jobs of the structure
            blockStruct[indexMachine][indexBlockInStruct].first = &instance->getListJobs()[indexAssignJob];
            // update the completion time of the structure and for the next block
            double completionTime = (indexBlockInStruct == 0) ? 0.0 : blockStruct[indexMachine][indexBlockInStruct -
                                                                                                1].second;
            double speed = (indexMachine < instance->getNbOfHighSpeedMachines()) ? instance->getHighSpeed()
                                                                                 : instance->getLowSpeed();
            unsigned int indexNextBlock = indexBlock;
            while (indexNextBlock < E.size()) {
                // if a job is schedule on the machine in the next block
                if (indexMachine < E[indexNextBlock].size()) {
                    auto &[_, indexNextBlockInStruct] = E[indexNextBlock][indexMachine];
                    auto nextJob = blockStruct[indexMachine][indexNextBlockInStruct].first;
                    completionTime += nextJob->getPi() / speed;
                    blockStruct[indexMachine][indexNextBlockInStruct].second = completionTime;
                }
                ++indexNextBlock;
            }
        }
    }
}

inline void Heuristic::computeAllWeightsAtEachBlock(std::vector<std::vector<double>> &costAtEachBlock, Solution::BlockStructure &blockStruct) {
    costAtEachBlock.clear();
    for (auto &machine: blockStruct) {
        costAtEachBlock.emplace_back();
        double costMachine = 0.0;
        for (auto &[job, CompletionTime]: machine) {
            if (job != nullptr) {
                costMachine += std::isless(job->getDi(), CompletionTime) ? job->getWi() : 0.0;
            }
            costAtEachBlock.back().emplace_back(costMachine);
        }
    }
}


#endif //BILEVEL_SCHEDULING_HEURISTIC_H
