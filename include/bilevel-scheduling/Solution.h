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

#ifndef BILEVEL_SCHEDULING_SOLUTION_H
#define BILEVEL_SCHEDULING_SOLUTION_H

#include <vector>
#include <ostream>
#include "Machine.h"
#include "Instance.h"
#include <queue>
#include "BiSchException.h"

class Solution {

private:

    std::vector<Machine> listHighSpeedMachines; // list of high speed machines
    std::vector<Machine> listLowSpeedMachines; // list of low speed machines
    double sum_wj_Uj; // the sum of the number weighted tardy jobs of the solution
    double sum_Cj; // the sum of the completion time of the solution
    size_t nbScheduledJobs;


public:

    // partial schedule [i][k] where i is the index of machine and k the index of the block. At position [i,k] we have a job and its completion time
    typedef std::vector<std::vector<std::pair<const Job *, double>>> BlockStructure;

    void static printB(const Solution::BlockStructure &blockStructure) {
        for (unsigned int indexMachine = 0; indexMachine < blockStructure.size(); ++indexMachine) {
            auto &machine = blockStructure[indexMachine];
            std::cout << "M" << indexMachine << ": " ;
            for (auto &job: machine) {
                if (job.first != nullptr) {
                    std::cout << "|" << job.first->getIndex() << " p:" << job.first->getPi() << " d:" << job.first->getDi() << "|";
                } else
                    std::cout << "\t X";
            }
            std::cout << std::endl;
        }
        std::cout << "----------------" << std::endl;
    }

    /************************/
    /*      CONSTRUCTOR     */
    /************************/

    /**
     * Default Constructor
     */
    Solution() : sum_wj_Uj(std::numeric_limits<double>::infinity()), sum_Cj(std::numeric_limits<double>::infinity()), nbScheduledJobs(0) {};


    /**
     * Constructor method for a Solution object. It initializes the list of machines by creating empty machines.
     * @param nbHighSpeedMachines The number of machines with high speed
     * @param nbLowSpeedMachines The number of machines with low speed
     */
    explicit Solution(Instance *instance) : sum_wj_Uj(instance->getSumWj()), sum_Cj(instance->getSumPj()), nbScheduledJobs(0) {
        listHighSpeedMachines = std::vector<Machine>(instance->getNbOfHighSpeedMachines(), Machine(instance->getHighSpeed()));
        listLowSpeedMachines = std::vector<Machine>(instance->getNbOfLowSpeedMachines(), Machine(instance->getLowSpeed()));
    }

    /**
     * Constructor method for a Solution object. It initializes the list of machines and number of scheduled jobs based
     * on an Instance object and BlockStructure.
     * @param instance The Instance object used to initialize the machines
     * @param blockStruct A BlockStructure object that defines the block structure property
     */
    explicit Solution(Instance *instance, const BlockStructure &blockStruct) {
        listHighSpeedMachines = std::vector<Machine>(instance->getNbOfHighSpeedMachines(), Machine(instance->getHighSpeed()));
        listLowSpeedMachines = std::vector<Machine>(instance->getNbOfLowSpeedMachines(), Machine(instance->getLowSpeed()));
        nbScheduledJobs = 0;
        fromBlockStruct(blockStruct);
    };

    /********************/
    /*      METHODS     */
    /********************/


    /**
     * This method exactly solves, for a given list of jobs, the criterion sum Cj by Brucker algorithm. The solution
     * is evaluate when it is returned.
     * @param listJobs The list of jobs to schedule. This list MUST BE SORTED ACCORDING LPT RULE
     * @param instance The instance from where we need to check if the solution is feasible
     */
    static Solution solveSumCjCriteria(const std::vector<Job> &listJobs, Instance *instance) {

        Solution sol = Solution(instance);
        // create a min heap of pair (weighted,index of machine)
        std::vector<std::pair<double, unsigned int>> weighted_indexes;
        weighted_indexes.reserve(instance->getNbJobs());
        for (unsigned int j = 0; j < instance->getNbMachines(); ++j) {
            double weight = j < instance->getNbOfHighSpeedMachines() ? (1 / instance->getHighSpeed()) : (1 /
                                                                                                         instance->getLowSpeed());
            weighted_indexes.emplace_back(weight, j);
        }

        // Use a priority queue, giving O(log n) performance for inserts and removals
        // (+ complexity of the containers, that here is in worst case amortized constant for vector push_back)
        // Using lambda to compare elements. If we have x=(a,b) and y=(c,d) then compare x>y we compare a>c or if a=c then look up after b<d
        // the idea is selected the smallest weight and the largest index
        auto cmp = [](std::pair<double, unsigned int> left, std::pair<double, unsigned int> right) {
            return (left.first == right.first) ? left.second < right.second : left.first > right.first;
        };
        std::priority_queue<std::pair<double, unsigned int>, std::vector<std::pair<double, unsigned int>>, decltype(cmp)> queue(
                weighted_indexes.begin(), weighted_indexes.end(), cmp);

        //loop over jobs to assign them
        for (const Job &job: listJobs) {
            // find the minimum weight with the largest index
            auto weighted_index = queue.top();
            queue.pop();
            sol.add_job(weighted_index.second, job);
            double weight =
                    weighted_index.second < instance->getNbOfHighSpeedMachines() ? (1 / instance->getHighSpeed())
                                                                                 : (1 / instance->getLowSpeed());
            double newWeight = weighted_index.first + weight;
            queue.emplace(newWeight, weighted_index.second);
        }
        sol.reverse();
        sol.evaluate();
        return sol;
    }


    /**
     * Method that reverses all jobs on all machines
     */
    void reverse() {
        for (Machine &machine: listHighSpeedMachines) {
            machine.reverse();
        }
        for (Machine &machine: listLowSpeedMachines) {
            machine.reverse();
        }
    }

    /**
     * Method that returns if the solution is empty
     * @return solution is empty
     */
    bool empty() {
        bool isEmpty = true;
        for (Machine &machine: listHighSpeedMachines) {
            isEmpty = isEmpty && machine.empty();
        }
        for (Machine &machine: listLowSpeedMachines) {
            isEmpty = isEmpty && machine.empty();
        }
        return isEmpty;
    }

    /**
     * Method that resets the current solution
     */
    void reset() {
        // reset all machines
        for (Machine &machine: listHighSpeedMachines) {
            machine.reset();
        }
        for (Machine &machine: listLowSpeedMachines) {
            machine.reset();
        }
        sum_wj_Uj = std::numeric_limits<double>::infinity();
        sum_Cj = std::numeric_limits<double>::infinity();
        nbScheduledJobs = 0;
    }

    /**
     * Method that evaluates the solution
     */
    void evaluate() {
        sum_wj_Uj = 0.0;
        sum_Cj = 0.0;
        // compute the objective function for the high speed machines
        for (auto &machine: listHighSpeedMachines) {
            machine.evaluate();
            sum_Cj += machine.getSumCj();
            sum_wj_Uj += machine.getSumWjUj();
        }

        // compute the objective function for the low speed machines
        for (auto &machine: listLowSpeedMachines) {
            machine.evaluate();
            sum_Cj += machine.getSumCj();
            sum_wj_Uj += machine.getSumWjUj();
        }
    }

    /**
     * Method that adds job to the corresponding idMachine at the position 'position'. This method update the attribute
     * nbScheduledJobs
     * @param idMachine The id of the machine where add the job
     * @param position The position in the machine
     * @param affectedJob The job to affect to the machine
     */
    void add_job(unsigned int idMachine, unsigned int position, const Job &affectedJob) {
        operator[](idMachine).add_job(position, affectedJob);
        ++nbScheduledJobs;
    }

    /**
     * Method that adds a job to the corresponding machine at its end.
     * @param idMachine The ID of the machine where the job will be added
     * @param affectedJob The job to be scheduled on the machine
     */
    void add_job(unsigned int idMachine, const Job &affectedJob) {
        operator[](idMachine).add_job(affectedJob);
        ++nbScheduledJobs;
    }

    /**
     * Method that retrieves the sorted list of jobs from the solution.The list is ordered based on
     * the Longest Processing Time (LPT) rule.
     * @return The list of jobs in the order of their longest processing time.
     */
    [[nodiscard]] std::vector<Job> extractListOfJobs() const {
        // Get the list of jobs
        std::vector<Job> listOfScheduledJobs;
        listOfScheduledJobs.reserve(nbScheduledJobs);
        // from the high speed
        for (const Machine &machine: listHighSpeedMachines) {
            for (const Job &job: machine.getAffectedJob()) {
                listOfScheduledJobs.push_back(job);
            }
        }
        // from the low speed
        for (const Machine &machine: listLowSpeedMachines) {
            for (const Job &job: machine.getAffectedJob()) {
                listOfScheduledJobs.push_back(job);
            }
        }
        // sort the jobs according to the LPT rule
        std::sort(listOfScheduledJobs.begin(), listOfScheduledJobs.end(), std::greater<>());
        return listOfScheduledJobs;
    }

    /**
     * Method that check if the solution is feasible for the leader, i.e. the right number of jobs have been selected
     * and the scheduling is optimal for the sum of Cj. The solution MUST HAVE BEEN already evaluate !
     * @param instance The problem instance from which the solution originates
     * @return The solution is feasible, i.e. optimal for the follower objective : Sum Cj
     */
    bool feasible(Instance *instance) {
        bool isFeasible = true;
        // check if the number of scheduled jobs is equal to the number of selected jobs
        if (nbScheduledJobs != instance->getNbToSelectJob()) return false;
        //check if all machines have enough and not more jobs for respect block structure
        // loop over high-speed machines
        for (auto &machine: listHighSpeedMachines) {
            if (machine.size() < instance->getMinNbJobsOnHS() || machine.size() > instance->getMaxNbJobsOnHS()) {
                isFeasible = false;
                return isFeasible;
            }
        }


        // loop over low-speed machines
        for (auto &machine: listLowSpeedMachines) {
            if (machine.size() < instance->getMinNbJobsOnLS() || machine.size() > instance->getMaxNbJobsOnLS()) {
                isFeasible = false;
                return isFeasible;
            }
        }

        BlockStructure blockStructure = toBlockStruct(instance);
        auto &E = instance->getE();
        unsigned int indexBlock = 0;
        while (indexBlock < E.size()) {
            // the number of jobs that must be scheduled on the block k
            unsigned int numJobsToScheduleOnBlock = (indexBlock == 0) ? instance->getNbJobsToScheduleOnFirstBlock() : E[indexBlock].size();
            double maxPjInBlock = -1.0, minPjNextBlock = std::numeric_limits<double>::infinity();
            // get the maximum of the pj in the block and reduced numJobsToScheduleOnBlock if there is a job scheduled on machine
            for (auto &location: E[indexBlock]) {
                auto [job, cj] = blockStructure[location.first][location.second];
                if (job != nullptr) {
                    maxPjInBlock = std::max(maxPjInBlock, job->getPi());
                    --numJobsToScheduleOnBlock;
                }
            }
            // if we have not scheduled all jobs in the block
            if (numJobsToScheduleOnBlock != 0) {
                isFeasible = false;
                break;
            }
            //pass to the next block
            ++indexBlock;
            if (indexBlock < E.size()) {
                // compute the min of the next block
                for (auto &location: E[indexBlock]) {
                    auto [job, cj] = blockStructure[location.first][location.second];
                    if (job != nullptr) {
                        minPjNextBlock = std::min(minPjNextBlock, job->getPi());
                    }
                }
            }
            // check if the max of pj from the block is smaller than the min of pj from the next block
            if (maxPjInBlock > minPjNextBlock) {
                isFeasible = false;
                break;
            }

        }
        return isFeasible;
    }

    bool feasibleOld(Instance *instance) const {
        // check if the number of scheduled jobs is equal to the number of selected jobs
        if (nbScheduledJobs != instance->getNbToSelectJob()) return false;
        auto listOfScheduledJobs = extractListOfJobs();
        Solution optForSumCj = solveSumCjCriteria(listOfScheduledJobs, instance);
        // check the objectives values
        if (sum_Cj != optForSumCj.sum_Cj) return false;
        else return true;
    }

    /********************/
    /*      GETTER      */
    /********************/

    [[nodiscard]] const std::vector<Machine> &getListHighSpeedMachines() const { return listHighSpeedMachines; }

    [[nodiscard]] const std::vector<Machine> &getListLowSpeedMachines() const { return listLowSpeedMachines; }

    [[nodiscard]] double getSumWjUj() const { return sum_wj_Uj; }

    [[nodiscard]] double getSumCj() const { return sum_Cj; }

    [[nodiscard]] size_t getNbScheduledJobs() const { return nbScheduledJobs; }

    /********************/
    /*      SETTER      */
    /********************/

    void setListHighSpeedMachines(
            const std::vector<Machine> &listHighSpeedMachines) { Solution::listHighSpeedMachines = listHighSpeedMachines; }

    void setListLowSpeedMachines(
            const std::vector<Machine> &listLowSpeedMachines) { Solution::listLowSpeedMachines = listLowSpeedMachines; }

    void setSumWjUj(double sumWjUj) { Solution::sum_wj_Uj = sumWjUj; }

    void setSumCj(double sumCj) { sum_Cj = sumCj; }


    /************************/
    /*      OPERATORS       */
    /************************/

    /**
     * Returns a reference to the element at specified location pos.
     * If pos is greater than the number of high speed machine then the reference is from the list of low speed machine
     * Else the reference is from the list of high speed machine. No bounds checking is performed.
     */
    Machine &operator[](size_t pos) {
        return pos < listHighSpeedMachines.size() ? listHighSpeedMachines[pos] : listLowSpeedMachines[pos -
                                                                                                      listHighSpeedMachines.size()];
    }

    void fromBlockStruct(const BlockStructure &blockStruct) {
        reset();
        for (unsigned int indexMachine = 0; indexMachine < blockStruct.size(); ++indexMachine) {
            for (auto pairJobCj: blockStruct[indexMachine]) {
                if (pairJobCj.first != nullptr) {
                    add_job(indexMachine, *pairJobCj.first);
                }
            }
        }
        evaluate();
    }

    /**
     * Method that constructs a schedule with block of job pairs (*Job, double) based on the solution's machines.
     * @param instance The instance to solve
     */
    BlockStructure toBlockStruct(Instance *instance) {
        if (listHighSpeedMachines.empty() && listLowSpeedMachines.empty())
            throw BiSchException("None machines have jobs");
        // schedule [i][k] where i is the index of machine and k the index of the block. At position [i,k] we have a job and the completion time
        BlockStructure blockStruc;
        // initialize completion time and block struct
        for (unsigned int i = 0; i < instance->getNbMachines(); ++i) {
            blockStruc.emplace_back(instance->getMaxNbJobsOnHS());
        }
        // loop over high-speed machines
        for (unsigned int indexLoopMachine = 0; indexLoopMachine < listHighSpeedMachines.size(); ++indexLoopMachine) {
            auto &machine = listHighSpeedMachines[indexLoopMachine];
            unsigned int indexBlock = (machine.size() < instance->getMaxNbJobsOnHS()) ? 1 : 0;
            double completionTime = 0.0;
            for (unsigned int indexLoopJobs = 0; indexLoopJobs < machine.size(); indexLoopJobs++) {
                auto &job = machine[indexLoopJobs];
                completionTime += job.getPi() / instance->getHighSpeed();
                assert(indexLoopMachine < blockStruc.size());
                assert(indexBlock < blockStruc[indexLoopMachine].size());
                blockStruc[indexLoopMachine][indexBlock] = std::make_pair(&instance->getListJobs()[job.getIndex()], completionTime);
                ++indexBlock;
            }
        }

        // loop over low-speed machines
        for (unsigned int indexLoopMachine = 0; indexLoopMachine < listLowSpeedMachines.size(); ++indexLoopMachine) {
            auto &machine = listLowSpeedMachines[indexLoopMachine];
            unsigned int indexBlock = (machine.size() < instance->getMaxNbJobsOnLS()) ? 1 : 0;
            double completionTime = 0.0;
            for (unsigned int indexLoopJobs = 0; indexLoopJobs < machine.size(); indexLoopJobs++) {
                auto &job = machine[indexLoopJobs];
                completionTime += job.getPi() / instance->getLowSpeed();
                blockStruc[indexLoopMachine + instance->getNbOfHighSpeedMachines()][indexBlock] = std::make_pair(&instance->getListJobs()[job.getIndex()], completionTime);
                ++indexBlock;
            }
        }
        return blockStruc;
    }
};

inline std::ostream &operator<<(std::ostream &os, const Solution &solution) {
    os << "Solution : " << std::endl
       << "high speed machines : [" << std::endl;
    for (const Machine &machine: solution.getListHighSpeedMachines()) {
        os << machine << std::endl;
    }
    os << "]" << std::endl
       << "low speed machines : [" << std::endl;
    for (const Machine &machine: solution.getListLowSpeedMachines()) {
        os << machine << std::endl;
    }
    os << "]" << std::endl << "Sum Cj : " << solution.getSumCj() << std::endl << "Sum wjUj : " << solution.getSumWjUj()
       << std::endl;
    return os;

}

#endif //BILEVEL_SCHEDULING_SOLUTION_H
