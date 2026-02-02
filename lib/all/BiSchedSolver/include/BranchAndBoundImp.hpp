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
// Created by schau on 11/25/24.
//

#ifndef BILEVEL_SCHEDULING_BRANCHANDBOUNDIMP_HPP
#define BILEVEL_SCHEDULING_BRANCHANDBOUNDIMP_HPP

inline void BranchAndBound::createChildrenNodeWithIdenticalJobs(Node &pNode) {
    ++nbSubProcessBaB;
    // case where we take no jobs:
    Node cNode = Node(pNode);
    #ifdef DEBUG_BaB
    std::string name = std::string("p=").append(
            std::to_string(cNode.getCurrentPj())).append(
            ",|s|=").append(std::to_string(cNode.getSelectedJobCount())).append(",L=").append(std::to_string(cNode.getIndexBlock())).append(",X");
    cNode.stateDebug.emplace_back(name);
    #endif
    unsigned int numberIdenticalJobsInGroups = cNode.removeGroupOfIdenticalJobs();
    #if defined DEBUG_BaB && defined DEBUG_DOT
    std::ofstream dot(DEBUG_DOT, std::ios::app);
    #endif
    addNode(cNode);
    #if defined DEBUG_BaB && defined DEBUG_DOT
    // if cNode have not the same id as its parents, means we add it to the list
    if (cNode.id != pNode.id)
        dot << pNode.id << " -- " << cNode.id << std::endl << cNode.id
            << "[label=\"" << "(id:" << cNode.id << "," << cNode.stateDebug.back().append(")").c_str() << "\"];"
            << std::endl;
    #endif

    // take at least one jobs
    for (unsigned int numberOfJobToSchedule = 1; numberOfJobToSchedule <= numberIdenticalJobsInGroups; ++numberOfJobToSchedule) {
        isWithinTimeLimit();
        // check if we respect the follower problem, i.e., we must select n jobs
        if (numberOfJobToSchedule + pNode.getSelectedJobCount() > instance->getNbToSelectJob()) break;
        // n_L the number of selected machine in the last block
        auto [indexLastBlock, n_L] = computeNumberJobsOnLastBlock(pNode, numberOfJobToSchedule);
        std::set<std::vector<unsigned int>> setOfAssignmentOnFirstBlock;
        // if the first block is the block 0, we need to fill it, that means we need to use different assignment for this block because it's the only one where there
        // the number of available machines a_0 > eta_0 the number of machines we must select. We have exactly a_0 + eta_0 - instance->getE()[0] available location
        if (pNode.getIndexBlock() == 0 and indexLastBlock > 0) {
            unsigned int numberAvailablePosition = pNode.getAvailableLocations() + pNode.getNumJobsToScheduleOnBlock() - instance->getE()[0].size();
            auto assignmentFirstBlock = std::vector<unsigned int>({});
            computeSetCombinationWithoutSymmetry(numberAvailablePosition, 0, pNode, assignmentFirstBlock, setOfAssignmentOnFirstBlock);
        } else setOfAssignmentOnFirstBlock.emplace();

        for (auto &assignmentFirstBlock: setOfAssignmentOnFirstBlock) {
            isWithinTimeLimit();
            std::set<std::vector<unsigned int>> setOfAssignmentOnLastBlock;
            computeSetCombinationWithoutSymmetry(n_L, indexLastBlock, pNode, assignmentFirstBlock, setOfAssignmentOnLastBlock);
            for (auto &assignmentLastBlock: setOfAssignmentOnLastBlock) {
                isWithinTimeLimit();
                createNodeWithAssignment(pNode, indexLastBlock, assignmentLastBlock, assignmentFirstBlock);
            }

        }
    }

}

inline std::pair<unsigned int, unsigned int> BranchAndBound::computeNumberJobsOnLastBlock(Node &pNode, unsigned int numberJobsToSchedule) {
    unsigned int startIndexBlock = pNode.getIndexBlock();
    std::vector<unsigned int> availableMachines;
    availableMachines.reserve(numberJobsToSchedule);
    unsigned int numJobsToScheduleOnBlock = pNode.getNumJobsToScheduleOnBlock(); // number of machine that must be selected in the block k
    unsigned int availableLocations = pNode.getAvailableLocations(); // number of available machine in block k
    unsigned int n_L = 0; // the number of job in last block
    while (numberJobsToSchedule > 0) {
        unsigned int numberAvailablePosition = availableLocations + numJobsToScheduleOnBlock - instance->getE()[startIndexBlock].size();
        if (numberJobsToSchedule > numberAvailablePosition) {
            numberJobsToSchedule -= numberAvailablePosition;
            ++startIndexBlock;
            numJobsToScheduleOnBlock = instance->getE()[startIndexBlock].size();
            availableLocations = instance->getE()[startIndexBlock].size();
        } else {
            n_L = numberJobsToSchedule;
            numberJobsToSchedule = 0;
        }
    }
    return {startIndexBlock, n_L};
}

inline void BranchAndBound::computeSetCombinationWithoutSymmetry(unsigned int nbSelect, unsigned int indexLastBlock, Node &node, const std::vector<unsigned int> &assigmentOnStartingBlock
                                                                 , std::set<std::vector<unsigned int>> &setOfAssignmentOnLastBlock) {
    isWithinTimeLimit();
    // create vector for the assignment on the first block
    std::vector<unsigned int> assignmentOnFirstBlock;
    assignmentOnFirstBlock.reserve(instance->getE()[node.getIndexBlock()].size());
    // get the processing time of identical jobs
    double pj = node.getCurrentPj();
    //if not assigment was defined for the first block, then construct it with available machine
    if (assigmentOnStartingBlock.empty()) {
        for (auto &location: instance->getE()[node.getIndexBlock()]) {
            // if not jobs where schedule, that means we have null value for the pointer
            if (!node.getBlockStruc()[location.first][location.second].first)
                assignmentOnFirstBlock.push_back(location.first);
        }
    } else
        assignmentOnFirstBlock.insert(assignmentOnFirstBlock.end(), assigmentOnStartingBlock.begin(), assigmentOnStartingBlock.end());


    // list of completion time of available machine use to compute the set of combinations
    std::vector<unsigned int> availableCompletionTime;
    // We store the list of elegant pair with completion time and the speed (it's bijection N^2 -> N)
    availableCompletionTime.reserve(instance->getE()[indexLastBlock].size());
    // use a map, where the key is elegant pair with completion time and the speed and the value is a vector of
    // index machine. The idea is to keep with this map the machines with same completion time on same speed.
    std::unordered_map<unsigned int, std::vector<unsigned int>> indexMachinesGroupByCompletionTime;
    for (auto &location: instance->getE()[indexLastBlock]) {
        isWithinTimeLimit();
        // location = (index Machine, index Block)
        // in (Job, Ci) in block structure /!\ Ci in the blockStructure are already divided by the speed
        double speed = (location.first < instance->getNbOfHighSpeedMachines()) ? instance->getHighSpeed() : instance->getLowSpeed();
        auto locationAssigned = node.getBlockStruc()[location.first][location.second];
        // if we have null value for the pointer that means no job was schedule, in other words the location is available
        if (!locationAssigned.first) {
            // compute the completion time of the machine.
            unsigned int Ci;
            // Check if we have several blocks.
            if (indexLastBlock > 0) {
                // compute the number of block that was filled,
                unsigned int numberOfFilledBlock = indexLastBlock - node.getIndexBlock() + 1;
                // if the machine was not present in the first block, then we have one less block to fulfill
                if (std::find(assignmentOnFirstBlock.begin(), assignmentOnFirstBlock.end(), location.first) == assignmentOnFirstBlock.end())
                    --numberOfFilledBlock;
                // We need to get the last block index where the machine is present and where there is a scheduled job that is
                unsigned int indexBlockWithSameMachine = indexLastBlock - 1;
                auto predFindIndexMachine = [&location](std::pair<unsigned int, unsigned int> locationInBlock) {
                    return locationInBlock.first == location.first;
                };
                bool foundLastBlockIndex = false;
                // create a location with index last block equal to 0
                auto locationLastScheduledJob = std::pair(location.first, 0);
                while (!foundLastBlockIndex && indexBlockWithSameMachine > 0) {
                    // looking for the machine's index in the last block
                    auto itLocationInBlockWithSameMachine = std::find_if(instance->getE()[indexBlockWithSameMachine].begin(), instance->getE()[indexBlockWithSameMachine].end(), predFindIndexMachine);
                    // if we have machine's index in the block, and we have a job that is already scheduled
                    if (itLocationInBlockWithSameMachine != instance->getE()[indexBlockWithSameMachine].end() &&
                        node.getBlockStruc()[itLocationInBlockWithSameMachine->first][itLocationInBlockWithSameMachine->second].first) {
                        locationLastScheduledJob = *itLocationInBlockWithSameMachine;
                        foundLastBlockIndex = true;
                    } else {
                        // go to the block before
                        --indexBlockWithSameMachine;
                    }
                }

                Ci = static_cast<unsigned int>(node.getBlockStruc()[locationLastScheduledJob.first][locationLastScheduledJob.second].second * speed +
                                               numberOfFilledBlock * pj); // get the completion time
            } else Ci = static_cast<unsigned int>(pj);
            // the key is the bijection N^2 -> N of (Ci, speed)
            unsigned int keyCi = (location.first < instance->getNbOfHighSpeedMachines()) ? elegantPair(Ci, 0u) : elegantPair(Ci, 1u);
            auto itIndexMachineGroupByCi = indexMachinesGroupByCompletionTime.find(keyCi);
            if (itIndexMachineGroupByCi == indexMachinesGroupByCompletionTime.end())
                indexMachinesGroupByCompletionTime.insert({keyCi, {location.first}});
            else itIndexMachineGroupByCi->second.emplace_back(location.first);
            availableCompletionTime.push_back(keyCi);
        }
    }

    //sort completion time in order to don't have duplicate combination
    std::sort(availableCompletionTime.begin(), availableCompletionTime.end());
    // compute the set of combinations with nbSelect elements
    std::set<std::vector<unsigned int>> setOfCombinationFromCompletionTime;
    findCombinations(availableCompletionTime, nbSelect, setOfCombinationFromCompletionTime);
    // Each combination contains elegant pair with completion time and the speed, we just want a combination of machine's index
    for (auto &combination: setOfCombinationFromCompletionTime) {
        isWithinTimeLimit();
        // we use 'indexListMachineWithCi' that correspond of the index in the list of machine's index
        unsigned int indexListMachineWithCi = 0;
        std::vector<unsigned int> newCombinationWithMachineIndex(combination);
        for (unsigned int indexInCombination = 0; indexInCombination < combination.size(); ++indexInCombination) {
            // if it's not the first element and the completion time is same that the previous element, we increase the index in the list of machine's index
            if (indexInCombination > 0 && combination[indexInCombination - 1] == combination[indexInCombination])
                ++indexListMachineWithCi;
            else
                // else we reset the index to the first element, i.e. 0
                indexListMachineWithCi = 0;
            // we get the real machine's index, by finding the completion time and get the index machine stored in the position 'indexListMachineWithCi'
            unsigned int indexMachine = indexMachinesGroupByCompletionTime.find(combination[indexInCombination])->second[indexListMachineWithCi];
            // change the completion time in combination by the machine's index
            newCombinationWithMachineIndex[indexInCombination] = indexMachine;
        }
        // add the new combination
        setOfAssignmentOnLastBlock.insert(newCombinationWithMachineIndex);
    }
}

inline void
BranchAndBound::createNodeWithAssignment(Node &pNode, unsigned int indexLastBlock, const std::vector<unsigned int> &assignmentLastBlock, const std::vector<unsigned int> &assignmentFirstBlock) {
    Node cNode = Node(pNode);

    // get the processing time of identical jobs
    double pj = cNode.getCurrentPj();

    // declare the binary tree of completion times that is not assigned with the corresponding position,
    // i.e. key is the completion time and the value is (index of machine, index of block). This is sorted in increasing order
    std::multimap<double, std::pair<unsigned int, unsigned int>> listCjAndAvailablePosition;

    if (!assignmentFirstBlock.empty()) {
        // compute completion of the first block and from assignmentFirstBlock std::vector
        for (unsigned int indexMachine: assignmentFirstBlock) {
            double machineSpeed = indexMachine < instance->getNbOfHighSpeedMachines() ? instance->getHighSpeed()
                                                                                      : instance->getLowSpeed();
            // use the real position in the machine with instance->getE()[index block][index machine]
            double newCj = (cNode.getIndexBlock() > 0) ?
                           cNode.getBlockStruc()[indexMachine][instance->getE()[cNode.getIndexBlock() -
                                                                                1][indexMachine].second].second +
                           pj / machineSpeed : pj / machineSpeed;
            // update the completion time in the partial scheduling
            unsigned int positionInMachine = instance->getE()[cNode.getIndexBlock()][indexMachine].second;
            cNode.setCompletionTimeOfBlockStructure(indexMachine, positionInMachine, newCj);
            cNode.updateNbJobToSchedule(indexMachine);
            listCjAndAvailablePosition.insert({newCj, {indexMachine, positionInMachine}});
        }
        cNode.updateChangeBlock(instance->getE()[cNode.getIndexBlock() + 1].size());
    }

    while (cNode.getIndexBlock() < indexLastBlock) {
        for (auto &location: instance->getE()[cNode.getIndexBlock()]) {
            isWithinTimeLimit();
            if (!cNode.getBlockStruc()[location.first][location.second].first) {
                double machineSpeed = location.first < instance->getNbOfHighSpeedMachines() ? instance->getHighSpeed()
                                                                                            : instance->getLowSpeed();
                double newCj = (location.second > 0) ?
                               cNode.getBlockStruc()[location.first][location.second - 1].second + pj / machineSpeed :
                               pj / machineSpeed;
                // update the completion time in the partial scheduling
                cNode.setCompletionTimeOfBlockStructure(location.first, location.second, newCj);
                cNode.updateNbJobToSchedule(location.first);
                listCjAndAvailablePosition.insert({newCj, location});
            }
        }
        cNode.updateChangeBlock(instance->getE()[cNode.getIndexBlock() + 1].size());
    }

    // compute completion of the last block and from assignmentLastBlock
    for (unsigned int indexMachine: assignmentLastBlock) {
        isWithinTimeLimit();
        auto predFindIndexMachine = [&indexMachine](std::pair<unsigned int, unsigned int> location) {
            return location.first == indexMachine;
        };
        double machineSpeed = indexMachine < instance->getNbOfHighSpeedMachines() ? instance->getHighSpeed()
                                                                                  : instance->getLowSpeed();
        // use the real position in the machine with instance->getE()[index block][index location] = (position in schedule, index machine), first find the last block (from the
        // left) where we have the machine
        double newCj;
        if (cNode.getIndexBlock() > 0) {
            // go on the left in order to find a block where the indexMachine appear for the first time
            unsigned indexBlockWithSameMachine = cNode.getIndexBlock() - 1;

            while (std::find_if(instance->getE()[indexBlockWithSameMachine].begin(), instance->getE()[indexBlockWithSameMachine].end(), predFindIndexMachine) ==
                   instance->getE()[indexBlockWithSameMachine].end()
                   && indexBlockWithSameMachine > 0)
                --indexBlockWithSameMachine;
            // if we indexBlockWithSameMachine == 0 and the machine 'indexMachine' is not in the block then the job appear for
            // the first time on this machine
            if (std::find_if(instance->getE()[indexBlockWithSameMachine].begin(), instance->getE()[indexBlockWithSameMachine].end(), predFindIndexMachine) ==
                instance->getE()[indexBlockWithSameMachine].end() && indexBlockWithSameMachine == 0)
                newCj = pj / machineSpeed;
                //else we can get the completion time of the given machine
            else {
                auto location = std::find_if(instance->getE()[indexBlockWithSameMachine].begin(), instance->getE()[indexBlockWithSameMachine].end(), predFindIndexMachine);
                newCj = cNode.getBlockStruc()[location->first][location->second].second + pj / machineSpeed;
            }
        } else newCj = pj / machineSpeed;
        unsigned int positionInMachine = std::find_if(instance->getE()[cNode.getIndexBlock()].begin(), instance->getE()[cNode.getIndexBlock()].end(), predFindIndexMachine)->second;
        listCjAndAvailablePosition.insert({newCj, {indexMachine, positionInMachine}});
        cNode.updateNbJobToSchedule(indexMachine);
    }

    /********************************************/
    /*      Compute weighted tardy Jobs         */
    /********************************************/

    // use a copy of the list of grouped jobs, because this list is sorted according SPT-EDD
    std::vector<Job> listOfIdenticalJobs(instance->getListGrpJobs()[cNode.getIndexGroup()]);
    solveProblemWithFixedCompletionTime(&cNode, nullptr, listCjAndAvailablePosition, listOfIdenticalJobs);

    // remove other jobs in the group
    for (auto &remainJob: listOfIdenticalJobs) {
        cNode.removeOneJob(remainJob);
    }

    cNode.passNextGroup();

    // remove the number of assignment in the last block to the number of available machines.
    cNode.setAvailableLocations(cNode.getAvailableLocations() - assignmentLastBlock.size());
    changeBlock(cNode);
    #ifdef DEBUG_BaB
    std::string name;
    name.append("p=").append(std::to_string(static_cast<unsigned int>(pj))).append(",");
    for (auto assignment: assignmentFirstBlock) {
        name.append(std::to_string(assignment));
    }
    name.append("i€|");
    for (auto assignment: assignmentLastBlock) {
        name.append(std::to_string(assignment));
    }
    name.append(",|s|=").append(std::to_string(cNode.getSelectedJobCount()));
    name.append(",L=").append(std::to_string(indexLastBlock));
    cNode.stateDebug.emplace_back(name);
    #endif
    #if defined DEBUG_BaB && defined DEBUG_DOT
    std::ofstream dot(DEBUG_DOT, std::ios::app);
    #endif
    addNode(cNode);
    #if defined DEBUG_BaB && defined DEBUG_DOT
    // if cNode have not the same id as its parents, means we add it to the list
    if (cNode.id != pNode.id)
        dot << pNode.id << " -- " << cNode.id << std::endl << cNode.id
            << "[label=\"" << "(id:" << cNode.id << "," << cNode.stateDebug.back().append(")").c_str() << "\"];"
            << std::endl;
    #endif


}

inline void BranchAndBound::addNode(Node &node) {
    #if defined DEBUG_BaB
    node.id = nbNodeLoc+1;
        #if defined DEBUG_DOT
        std::ofstream dot(DEBUG_DOT, std::ios::app);
        #endif
    #endif
    if (memorizationActivate) {
        bool isDominated = memorization.checkIfNodeIsDominated(node, false);
        if (isDominated) {
            ++nbCut;
            #if defined DEBUG_BaB
            ++nbNodeLoc;
            #endif
            return;
        }
    }
    // if not enough job or if we are not in the last group then cut
    if (instance->getNbJobs() - node.getRemovedJobsCount() < instance->getNbToSelectJob()) {
        ++nbCut;
        #if defined DEBUG_BaB
        ++nbNodeLoc;
            #if defined DEBUG_DOT
            dot << node.id << DOT_CUT << ";";
            dot << node.id << "[label=\""<<Node::debug_dot_node(node)<<",LB: NE_JOBS\"];" << std::endl;
            #endif
        #endif
        return;
    }
    ++nbNodeLoc; // incr the nb node explored
    // if we have to select and schedule jobs on the last block, we can solve it by assignment algorithm instead of branching
    if (node.getSelectedJobCount() == instance->getNbToSelectJob() - instance->getE().back().size()) {
        Solution solFromNode(instance);
        auto blockStruc = node.getBlockStruc();
        solFromNode.fromBlockStruct(blockStruc);
        std::vector<Job> listAvailableJob;
        listAvailableJob.reserve(instance->getNbJobs() - (node.getSelectedJobCount() + node.getRemovedJobsCount()));
        for (auto job: instance->getListJobs()) {
            if (not node.isScheduled(job.getIndex()) && not node.isRemoved(job.getIndex())) {
                listAvailableJob.push_back(job);
            }
        }

        //here the method freeAndAssignmentBlock, is optimal for the last block since we solve an assignment problem with all available jobs on the last block
        heuristicSolver.freeAndAssignmentBlock(blockStruc, node.getIndexBlock(), &listAvailableJob);
        solFromNode.fromBlockStruct(blockStruc);
        if (solFromNode.feasible(instance)) {
            if (solFromNode.getSumWjUj() < globalUB) {

                #if defined DEBUG_BaB && defined DEBUG_DOT
                dot << node.id << DOT_OPT << ";";
                node.stateDebug.emplace_back(std::string("LB:").append(std::to_string(solFromNode.getSumWjUj())));
                #endif
                *solution = solFromNode;
                globalUB = solFromNode.getSumWjUj();
            }
            #if defined DEBUG_BaB && defined DEBUG_DOT
            else {
                dot << node.id << DOT_OPT << ";";
                node.stateDebug.emplace_back(std::string("SOL:").append(std::to_string(solFromNode.getSumWjUj())));
            }
            #endif
            return;
        } else throw BiSchException("Select n job, unfeasible where as branching scheme make feasible");
    } else // if we have selected the right nb of job
    if (node.getSelectedJobCount() == instance->getNbToSelectJob()) {
        Solution solFromNode(instance);
        solFromNode.fromBlockStruct(node.getBlockStruc());
        if (solFromNode.feasible(instance)) {
            if (solFromNode.getSumWjUj() < globalUB) {
                #if defined DEBUG_BaB && defined DEBUG_DOT
                dot << node.id << DOT_OPT << ";";
                node.stateDebug.emplace_back(std::string("LB:").append(std::to_string(solFromNode.getSumWjUj())));
                #endif
                *solution = solFromNode;
                globalUB = solFromNode.getSumWjUj();
            }
            #if defined DEBUG_BaB && defined DEBUG_DOT
            else {
                dot << node.id << DOT_OPT << ";";
                node.stateDebug.emplace_back(std::string("SOL:").append(std::to_string(solFromNode.getSumWjUj())));
            }
            #endif
            return;
        } else throw BiSchException("Select n job, unfeasible where as branching scheme make feasible");
    }
    double lowerBound = 0.0;
    if (lBStrategy == LB_MIP) {
        lowerBound = lbFromMIP.computeLB(node);
        if (lbFromMIP.getStatus() == IloAlgorithm::Optimal && lowerBound < globalUB) {
            globalUB = lowerBound;
            lbFromMIP.computeSolution();
            *solution = *lbFromMIP.getSolution();
        }
    } else {
        columnGeneration.solve(node);
        lowerBound = columnGeneration.getSumWjUj();
    }

    if (node.getPartialSumWjUj() - lowerBound >= EPSILON) {
        throw BiSchException(std::string("LB is smaller than the partial objective function value for instance : ").append(instance->getInstancePath()).c_str());
    }

    // add the node depending on the strategy
    switch (walkStrategy) {
        #if defined DEBUG_BaB && defined DEBUG_DOT
        // copy the node to get same stateDebug after the call of the method
        case DEPTH_FIRST:
            stackActiveNode.push_back(NodeWithLB{lowerBound,node});
            break;
        case BREADTH_FIRST:
            queueActiveNode.push(NodeWithLB{lowerBound,node});
            break;
        case BEST_FIRST:
            heap.push(NodeWithLB{lowerBound, node});
            break;
        #else
        case DEPTH_FIRST:stackActiveNode.push_back(NodeWithLB{lowerBound, std::move(node)});
            break;
        case BREADTH_FIRST:queueActiveNode.push(NodeWithLB{lowerBound, std::move(node)});
            break;
        case BEST_FIRST:heap.push(NodeWithLB{lowerBound, std::move(node)});
            break;
            #endif
    }
}

inline void BranchAndBound::branchingLocation(NodeWithLB &nodeWithLb) {

    #if defined DEBUG_BaB && defined DEBUG_DOT
    std::ofstream dot(DEBUG_DOT, std::ios::app);
    #endif
    //index use to sort added node according LB
    unsigned int indexBeforeAddingChildNodes = 0;
    if (walkStrategy == DEPTH_FIRST) {
        indexBeforeAddingChildNodes = stackActiveNode.size();
    }

    double nodeLB = nodeWithLb.lowerBound;
    Node &pNode = nodeWithLb.node;

    if (verbose >= 2) {
        unsigned int sizeActive = (walkStrategy == DEPTH_FIRST) ?
                                  stackActiveNode.size()
                                                                : (walkStrategy == BREADTH_FIRST) ?
                                                                  queueActiveNode.size()
                                                                                                  : heap.size();
        std::cout << "Active size:" << sizeActive << " UB: " << globalUB << " LB: " << nodeLB << " nbCleaningSetCol: " << columnGeneration.getNbCleaningSetCol() << " setCol size: "
                  << columnGeneration.getXs().getSize() << " nbNode: " << nbNodeLoc << std::endl;
    }

    /**
     * Check if the node is dominated by memorization
     */

    if (memorizationActivate) {
        bool isDominated = memorization.checkIfNodeIsDominated(pNode, true);
        if (isDominated) {
            ++nbCut;
            #if defined DEBUG_BaB && defined DEBUG_DOT
            dot << pNode.id << DOT_MEMO_DOMINATED << ";" << std::endl <<
                pNode.id << "[label=\"" << "(id:" << pNode.id << "," << pNode.stateDebug.back().append(")").c_str()
                << "\"];" << std::endl;
            #endif
            return;
        }
    }


    // if the LB is equal to the UB then stop
    if (std::fabs(nodeLB - globalUB) <= EPSILON) {
        #if defined DEBUG_BaB && defined DEBUG_DOT
        dot << pNode.id << "[style=filled,fillcolor=blue];" << std::endl;
        dot << pNode.id << "[label=\"" << Node::debug_dot_node(pNode) << ",LB:" << nodeLB << "=UB\"];" << std::endl;
        #endif
        return;
    }
        // if the LB is smaller than the UB
    else if (nodeLB < globalUB) {

        // get the index of the group, it's the max between the indexGroup of both kind of machine.
        unsigned int indexGroup = pNode.getIndexGroup();
        auto &listGroupedJobs = instance->getListGrpJobs();
        isWithinTimeLimit();
        if (listGroupedJobs[indexGroup].size() > 1) {
            createChildrenNodeWithIdenticalJobs(pNode);
        } else {
            // use a boolean to know if the current job is on time whatever its position.
            bool jobIsEarly = true;
            // fixe a location
            for (auto itLocation = instance->getE()[pNode.getIndexBlock()].begin(); itLocation != instance->getE()[pNode.getIndexBlock()].end(); ++itLocation) {
                isWithinTimeLimit();
                // create child node
                Node cNode = Node(pNode);

                // if the location is available
                if (!cNode.getBlockStruc()[itLocation->first][itLocation->second].first) {

                    // the machine speed
                    char machineSpeed = itLocation->first < instance->getNbOfHighSpeedMachines() ? 0 : 1;
                    double speed = (machineSpeed == 0) ? instance->getHighSpeed() : instance->getLowSpeed();

                    // we add the jobs to the fixed location
                    const Job &scheduledJob = listGroupedJobs[indexGroup].back();

                    double C_j = (itLocation->second > 0) ?
                                 cNode.getBlockStruc()[itLocation->first][itLocation->second - 1].second
                                 + scheduledJob.getPi() / speed : scheduledJob.getPi() / speed;
                    // the job is late ?
                    if (jobIsEarly)
                        jobIsEarly = isSmallerOrEqual(C_j, scheduledJob.getDi());
                    //add the pointer from the list of job (of instance)
                    cNode.scheduleOneJob(itLocation->first, itLocation->second, &instance->getListJobs()[scheduledJob.getIndex()], C_j);
                    cNode.updateNbJobToSchedule(itLocation->first);

                    // update constant that count the number of available location
                    cNode.setAvailableLocations(cNode.getAvailableLocations() - 1);
                    changeBlock(cNode);
                    // add the node to the list of active node
                    #if defined DEBUG_BaB
                    std::string name = std::string("J").append(std::to_string(scheduledJob.getIndex())).append(
                            ",i").append(
                            std::to_string(itLocation->first)).append(",k=").append(
                            std::to_string(itLocation->second));
                    cNode.stateDebug.emplace_back(name);
                    #endif
                    addNode(cNode);
                    #if defined DEBUG_BaB && defined DEBUG_DOT
                    // if cNode have not the same id as its parents, means we add it to the list
                    if (cNode.id != pNode.id)
                        dot << pNode.id << " -- " << cNode.id << std::endl << cNode.id
                            << "[label=\"" << "(id:" << cNode.id << "," << cNode.stateDebug.back().append(")").c_str()
                            << "\"];" << std::endl;
                    #endif
                    // if we are in the first block, then  don't try for the same job to schedule it on identical
                    // machine.
                    if (pNode.getIndexBlock() == 0) {
                        // So if we are on high speed machine, then go to low speed machine
                        if (itLocation->first < instance->getNbOfHighSpeedMachines()) {
                            // it's the nb of high speed machine - the index of current machine - 1 (-1 because the loop
                            //increase of one)
                            unsigned int nbJump = instance->getNbOfHighSpeedMachines() - itLocation->first - 1;
                            std::advance(itLocation, nbJump);
                        }//else we are on low speed machine then we can break
                        else {
                            break;
                        }
                    }
                }
            }
            // if we are in first block and there is no late jobs and the current job is on time, then don't create the node where we remove a job
            if (jobIsEarly && pNode.getIndexBlock() == 0 && isSmallerOrEqual(pNode.getPartialSumWjUj(), 0.0)) {
                return;
            } else {
                // don't select job
                Node cNode = Node(pNode);
                #if defined DEBUG_BaB
                unsigned int indexOfRemovedJob =
                #endif
                cNode.removeCurrentJob();
                #if defined DEBUG_BaB
                std::string name = std::to_string(indexOfRemovedJob).append("X");
                cNode.stateDebug.emplace_back(name);
                #endif
                addNode(cNode);
                #if defined DEBUG_BaB && defined DEBUG_DOT
                // if cNode have not the same id as its parents, means we add it to the list
                if (cNode.id != pNode.id)
                    dot << pNode.id << " -- " << cNode.id << std::endl << cNode.id << "[label=\""
                        << "(id:" << cNode.id << "," << cNode.stateDebug.back().append(")").c_str() << ","
                        << indexOfRemovedJob << "X\"];" << std::endl;
                #endif
            }
        }
    } else if (nodeLB > globalUB) {
        ++nbCut;
        #if defined DEBUG_BaB && defined DEBUG_DOT
        dot << pNode.id << DOT_CUT << ";";
        dot << pNode.id << "[label=\"" << Node::debug_dot_node(pNode) << ",LB:" << nodeLB << "\"];" << std::endl;
        #endif
        return;
    }
    if (walkStrategy == DEPTH_FIRST) {
        std::sort(stackActiveNode.begin() + indexBeforeAddingChildNodes, stackActiveNode.end());
    }
}

inline void BranchAndBound::changeBlock(Node &node) {
    // if we have totally filled the block, then we prepare the next (if is not last, i.e., indexBlock < |E|-1).
    if (node.getAvailableLocations() + node.getNumJobsToScheduleOnBlock() ==
        instance->getE()[node.getIndexBlock()].size()
        &&
        node.getIndexBlock() < instance->getE().size() - 1) {
        //update other jobs with release date that is the minimum completion time of the block
        double minCj = std::numeric_limits<double>::infinity();
        // for all machine, compute the completion time of each ones and update the release date of other jobs
        for (unsigned int i = 0; i < instance->getNbMachines(); ++i) {
            //if there is no jobs, then go block before until we find a job
            unsigned int indexBlock = node.getIndexBlock();
            while (!node.getBlockStruc()[i][indexBlock].first && indexBlock > 0)
                --indexBlock;
            double C_ik = node.getBlockStruc()[i][indexBlock].second;
            if (C_ik < minCj) {
                minCj = C_ik;
            }
        }
        node.updateChangeBlock(instance->getE()[node.getIndexBlock() + 1].size());
    }
}


#endif //BILEVEL_SCHEDULING_BRANCHANDBOUNDIMP_HPP
