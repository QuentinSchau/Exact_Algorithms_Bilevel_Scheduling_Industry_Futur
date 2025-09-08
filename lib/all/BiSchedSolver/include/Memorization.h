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

//
// Created by schau on 11/12/24.
//

#ifndef BILEVEL_SCHEDULING_MEMORIZATION_H
#define BILEVEL_SCHEDULING_MEMORIZATION_H

#include <vector>
#include <unordered_map>
#include "Instance.h"
#include "Node.h"

class Memorization {

public:
    /**
     * \enum cleaning_t
     * \brief enum to know to chose the strategy of cleaning the memory
     */
    typedef enum {
        NOT_USED = 0,
    } cleaning_t;

    /**
     * Information for a node that allow to apply memorization in order to prune the node
     */
    struct Infos_t {
        #if defined DEBUG_BaB
        unsigned int nodeId; // Node ID for debugging purposes
        #endif
        std::bitset<MAX_NUMBER_JOB> setRemainingJobs; // Remaining jobs on this node
        std::vector<std::bitset<MAX_NUMBER_JOB>> encodingSelectedJobOnMachines; // Encoding of selected jobs on machines
        unsigned int nbJobToScheduleOnHSMachines; // Number of jobs to schedule on high-speed machines
        unsigned int nbJobToScheduleOnLSMachines; // Number of jobs to schedule on low-speed machines
        // List of completion times for each machine, sorted in descending order (high speed first)
        std::vector<double> listCompletionTime;
        double sum_wj_Uj; // Partial objective function value
        bool haveCut = false; // Whether the node has been cut

        inline Infos_t(std::bitset<MAX_NUMBER_JOB> &setRemainingJobs, const std::vector<std::bitset<MAX_NUMBER_JOB>> &encodingSelectedJobOnMachines, unsigned int nbJobToScheduleOnHSMachines
                       , unsigned int nbJobToScheduleOnLSMachines, const std::vector<double> &listCompletionTime, double sumWjUj)
                : setRemainingJobs(setRemainingJobs), encodingSelectedJobOnMachines(encodingSelectedJobOnMachines), nbJobToScheduleOnHSMachines(nbJobToScheduleOnHSMachines)
                  , nbJobToScheduleOnLSMachines(nbJobToScheduleOnLSMachines), listCompletionTime(listCompletionTime), sum_wj_Uj(sumWjUj) {}

        /**
         * Method that compares two lists of completion times: 'listCompletionTime' and 'otherListCompletionTime'.
          It checks if 'listCompletionTime' dominates 'otherListCompletionTime', which means all completion times in
          'listCompletionTime' are smaller than or equal to those in 'otherListCompletionTime'.
          @param otherListCompletionTime The list of completion times to compare with
          @return True if 'listCompletionTime' has all completion times smaller than ones in 'otherListCompletionTime', otherwise false.
         */
        inline bool compareListCompletionTime(const std::vector<double> &otherListCompletionTime) {
            // Check if lists have the same size
            if (otherListCompletionTime.size() != listCompletionTime.size()) return false;
            // Compare elements of the two lists
            for (unsigned int indexLoopCj = 0; indexLoopCj < listCompletionTime.size(); ++indexLoopCj) {
                // If an element in 'listCompletionTime' is greater than the corresponding element in 'otherListCompletionTime', return false
                if (listCompletionTime[indexLoopCj] > otherListCompletionTime[indexLoopCj]) return false;
            }
            // If no larger elements were found, return true (i.e., 'listCompletionTime' dominates 'otherListCompletionTime')
            return true;
        }
    };

    /**
     * The database where we store all node's information, each row `database_t[j]` contains the nodes information for the case where there are `j` remaining jobs
     */
    typedef std::vector<std::vector<Infos_t>> database_t;


private:

    database_t database;

    // The maximum capacity of the database.
    const unsigned int capacity = 1000000;

    // The strategy used to clean the database.
    const cleaning_t cleaning;

    const Instance *instance;

    unsigned int size_db = 0;
    unsigned int _rem;
    unsigned int _cut;
    unsigned int nbCleaning = 0;

public:

    Memorization() : cleaning(cleaning_t::NOT_USED), instance(nullptr), _rem(0), _cut(0) {}

    Memorization(unsigned int capacity, cleaning_t cleaning, const Instance *instance);

    /*******************/
    /*      GETTER     */
    /*******************/

    unsigned int getRem() const {
        return _rem;
    }

    unsigned int getCut() const {
        return _cut;
    }

    unsigned int getNbCleaning() const {
        return nbCleaning;
    }

    /**
     * Method that check if 'infoOfNode' dominate 'infoOfNodeToCheck'.
     * @param infoOfNode The node that we want to check if it dominate the other.
     * @param infoOfNodeToCheck The node that can be dominated
     * @return True if 'infoOfNode' dominate 'infoOfNodeToCheck' false otherwise
     */
    static bool isDominated(Infos_t &infoOfNode, const Infos_t &infoOfNodeToCheck) {
        if ((infoOfNodeToCheck.setRemainingJobs & infoOfNode.setRemainingJobs) == infoOfNodeToCheck.setRemainingJobs
            && isSmallerOrEqual(infoOfNode.sum_wj_Uj, infoOfNodeToCheck.sum_wj_Uj)
            && infoOfNode.nbJobToScheduleOnHSMachines <= infoOfNodeToCheck.nbJobToScheduleOnHSMachines
            && infoOfNode.nbJobToScheduleOnLSMachines <= infoOfNodeToCheck.nbJobToScheduleOnLSMachines
            && infoOfNode.compareListCompletionTime(infoOfNodeToCheck.listCompletionTime)) {
            // if we have same number of scheduled jobs, we need to check that it's not same partial solution
            if ((infoOfNode.nbJobToScheduleOnLSMachines + infoOfNode.nbJobToScheduleOnHSMachines) ==
                (infoOfNodeToCheck.nbJobToScheduleOnLSMachines + infoOfNodeToCheck.nbJobToScheduleOnHSMachines)) {
                // Use std::equal for fast pairwise comparison, if encoding are equals, then the node is not dominated
                return !std::equal(infoOfNode.encodingSelectedJobOnMachines.begin(), infoOfNode.encodingSelectedJobOnMachines.end(), infoOfNodeToCheck.encodingSelectedJobOnMachines.begin());
            }
            return true;
        }
        return false;

    }

    void cleanDatabase() {
        ++nbCleaning;
        switch (cleaning) {
            case NOT_USED:auto pred_remove_if_not_used = [](Infos_t &infosOfNode) { return !infosOfNode.haveCut; };
                for (auto itlistOfInfosOfNode = database.begin(); itlistOfInfosOfNode != database.end();) {
                    auto itRemovedElement = std::remove_if(itlistOfInfosOfNode->begin(), itlistOfInfosOfNode->end(), pred_remove_if_not_used);
                    unsigned int nbRemovedElements = (itlistOfInfosOfNode->end() - itRemovedElement);
                    itlistOfInfosOfNode->erase(itRemovedElement, itlistOfInfosOfNode->end());
                    if (nbRemovedElements > 0) {
                        size_db -= nbRemovedElements;
                        _rem += nbRemovedElements;
                    }
                    // pass next element
                    ++itlistOfInfosOfNode;
                }
                break;
        }
    }

    /**
     * Method that check if the node is dominated by an already explored node. This method update the database by deleting
     * all dominated nodes.
     * @param node The node to check if it is dominated.
     * @param addToDataBase Boolean to know if we need to add the node to the database. If it's false, that mean we only check if the node is dominated.
     * @return True if the node is dominated, in such case, we don't need to branch on it, otherwise false.
     */
    bool checkIfNodeIsDominated(Node &node, bool addToDataBase) {

        std::vector<double> listOfCompletionTime;
        node.getListOfCompletion(listOfCompletionTime);
        std::bitset<MAX_NUMBER_JOB> encodingAlreadyDecidedJobs;
        for (auto &encodingMachine: node.getEncodingSelectedJobOnMachine()) {
            encodingAlreadyDecidedJobs |= encodingMachine;
        }
        encodingAlreadyDecidedJobs |= node.getEncodingRemoveDecision();
        //get the set of remaining job
        std::bitset<MAX_NUMBER_JOB> remainingSet = (std::bitset<MAX_NUMBER_JOB>{}.set() ^ encodingAlreadyDecidedJobs);
        auto [nbJobToSelectOnHSMachines, nbJobToSelectOnLSMachines] = node.getNumberJobsToScheduleOnMachines();
        Infos_t infoOfNode = Infos_t(remainingSet, node.getEncodingSelectedJobOnMachine(), nbJobToSelectOnHSMachines, nbJobToSelectOnLSMachines, listOfCompletionTime, node.getPartialSumWjUj());

        #if defined DEBUG_BaB
        infoOfNode.nodeId = node.id;
        #endif
        // if we want to remove node that are dominated by the current one
        if (addToDataBase) {
            // loop over all subset of remaining jobs, and check if the node dominate other job
            for (unsigned int nbRemainingJob = instance->getNbToSelectJob();
                 nbRemainingJob < infoOfNode.setRemainingJobs.count(); ++nbRemainingJob) {
                if (nbRemainingJob >= database.size()) throw BiSchException("Error memo");
                auto &listNode = database[nbRemainingJob];
                // if we have node where the remaining set have the cardinality equals to 'nbRemainingJob'
                if (!listNode.empty()) {

                    // remove all node that are dominated by the current node
                    auto pred_remove_dominated_nodes = [&](const Infos_t &infoOfNodeToCheck) {
                        return isDominated(infoOfNode, infoOfNodeToCheck);
                    };
                    auto itRemoveElement = std::remove_if(listNode.begin(), listNode.end(), pred_remove_dominated_nodes);
                    // reduce the size of database by the number of removed elements
                    unsigned int nbRemovedElements = (listNode.end() - itRemoveElement);
                    if (nbRemovedElements > 0) {
                        size_db -= nbRemovedElements;
                        // number of node that are cut
                        _cut += nbRemovedElements;
                        infoOfNode.haveCut = true;
                    }
                    listNode.erase(itRemoveElement, listNode.end());
                }
            }
        }
        // loop over set of remaining job that equals or greater than the one given by the current node to check if the node is not dominated
        for (unsigned int nbRemainingJob = infoOfNode.setRemainingJobs.count(); nbRemainingJob <= MAX_NUMBER_JOB; nbRemainingJob++) {
            if (nbRemainingJob >= database.size()) throw BiSchException("Error memo");
            auto &listNode = database[nbRemainingJob];
            // if we have nodes with 'nbRemainingJob' jobs in the remaining set of jobs, i.e, the list of node is not empty
            if (!listNode.empty()) {
                // check if there exist a node that dominate the current node, i.e. node where nbRemainingJob = nbRemainingJob
                auto pred_check_if_node_is_dominated = [&](Infos_t &infoOfNodeToCheck) {
                    return isDominated(infoOfNodeToCheck, infoOfNode);
                };
                auto itFindDominantNode = std::find_if(listNode.begin(), listNode.end(), pred_check_if_node_is_dominated);
                // if there is no dominant node and nbRemainingJob is the cardinality of the set of remaining jobs of the current node, then add the infos_t only if is needed
                if (addToDataBase && itFindDominantNode == listNode.end()) {
                    listNode.push_back(infoOfNode);
                    ++size_db;
                    return false;
                } else if (itFindDominantNode != listNode.end()) {
                    //the node is dominated
                    #if defined DEBUG_BaB && defined DEBUG_DOT
                    std::ofstream dot(DEBUG_DOT, std::ios::app);
                        dot << node.id << DOT_MEMO_DOMINATED << ";" << std::endl;
                        node.stateDebug.back().append(",DomBy: ").append(std::to_string(itFindDominantNode->nodeId));
                    #endif
                    // set the node that have cut
                    itFindDominantNode->haveCut = true;
                    ++_cut;
                    return true;
                }
            }
                // if we are looking for the case where 'nbRemainingJob' is equal to the cardinality of the remaining job of the current node then add it is it's needed
            else if (addToDataBase && nbRemainingJob == remainingSet.count()) {
                database[remainingSet.count()].emplace_back(infoOfNode);
                ++size_db;
            }
            // check if we need to clean the database
            if (size_db > capacity) {
                cleanDatabase();
                return false;
            }
        }
        return false;
    }
};


#endif //BILEVEL_SCHEDULING_MEMORIZATION_H
