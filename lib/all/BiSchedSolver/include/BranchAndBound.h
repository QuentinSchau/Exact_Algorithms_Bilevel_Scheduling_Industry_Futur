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
// Created by schau on 3/26/24.
//

#ifndef BILEVEL_SCHEDULING_BRANCHANDBOUND_H
#define BILEVEL_SCHEDULING_BRANCHANDBOUND_H

#include "ISolver.h"
#include "MIP.h"
#include "Heuristic.h"
#include "Node.h"
#include "ColumnGeneration.h"
#include "LowerBoundMIP.h"
#include "Memorization.h"


class BranchAndBound : public ISolver{

public:

    struct NodeWithLB
    {
        double lowerBound{-std::numeric_limits<double>::infinity()};
        Node node;

        inline friend bool operator<(NodeWithLB const& lhs, NodeWithLB const& rhs)
        {
            return lhs.lowerBound > rhs.lowerBound;
        }

        inline friend std::ostream& operator<<(std::ostream& os, NodeWithLB const& e)
        {
            Solution solOfNode;
            solOfNode.fromBlockStruct(e.node.getBlockStruc());
            return os << '{' << e.lowerBound << ", '" << solOfNode  << "'}";
        }
    };

    enum Strategy { DEPTH_FIRST , BREADTH_FIRST, BEST_FIRST };
    enum Scheme {LOCATION}; // if we want try other branching scheme
    enum LowerBound {CG,LB_MIP};

    /************************/
    /*      CONSTRUCTOR     */
    /************************/

    explicit BranchAndBound();

    explicit BranchAndBound(Instance *instance);

    explicit BranchAndBound(Instance *instance, nlohmann::json &object);


    /***********************/
    /*      DESTRUCTOR     */
    /***********************/

    ~BranchAndBound() override;

    /*******************/
    /*      GETTER     */
    /*******************/

    bool isOPT() const {return isOptimal;}

    [[nodiscard]] std::string getStrategy(){
        std::string strat;
        switch (walkStrategy) {
            case DEPTH_FIRST:
                strat = "depth-first";
                break;
            case BREADTH_FIRST:
                strat = "breadth-first";
                break;
            case BEST_FIRST:
                strat = "best-first";
                break;
        }
        return strat;
    }

    [[nodiscard]] std::string getLowerBoundStrategy(){
        std::string lB_strat;
        switch (lBStrategy) {
            case LB_MIP:
                lB_strat = "LB_from_MIP";
                break;
            case CG:
                lB_strat = "Columns_Generation";
                break;
        }
        return lB_strat;
    }
    /*******************/
    /*      SETTER     */
    /*******************/



    void setStrategy(std::string strat){
        if (strat == "breadth-first") walkStrategy = BREADTH_FIRST;
        else if (strat == "depth-first") walkStrategy = DEPTH_FIRST;
        else if (strat == "best-first") walkStrategy = BEST_FIRST;
        else throw BiSchException("Strategy is not known for the branch-and-bound algorithm, read \"README\" file for more details on which strategy to use.");
    }

    void setMemorizationActivate(bool memorizationActivate) {
        BranchAndBound::memorizationActivate = memorizationActivate;
    };

    void setLBStrategy(std::string strat){
        if (strat == "LB_MIP") lBStrategy = LB_MIP;
        else if (strat == "CG") lBStrategy = CG;
    };

    void setBranchingScheme(std::string scheme) {
        if (scheme == "location") branchingScheme = LOCATION;
    };

    /********************/
    /*      METHODS     */
    /********************/


    /***************************************/
    /*      Sub-Problem Identical Jobs     */
    /***************************************/

    /**
     * Method that creates children nodes.
     * @param pNode The parent node.
     * @return True if a node have been added False otherwise
     */
    void createChildrenNodeWithIdenticalJobs(Node &pNode);

    /**
     * Method that computes the number of jobs on last block if we need to schedule 'numberJobsToSchedule' of identical jobs.
     * @param pNode The parent node
     * @param numberJobsToSchedule The number of jobs that we must scheduled
     * @return (listAvailableMachines,k',n_L) the list of available machines from where we can assign n_L jobs on last last block k'
     */
    std::pair<unsigned int,unsigned int> computeNumberJobsOnLastBlock(Node &pNode, unsigned int numberJobsToSchedule);

    /**
     * Computes a set combination without symmetry.
     *
     * @param nbSelect The number of elements to select from the set.
     * @param indexLastBlock The index of the last block.
     * @param node A reference to a Node object.
     * @param assigmentOnStartingBlock A vector of unsigned integers representing the assignment on the starting block.
     * @param setOfAssignmentOnLastBlock A set of combination, represented as vectors of machine indexes, representing the set of assignments on the last block.
     */
    void computeSetCombinationWithoutSymmetry(unsigned int nbSelect,unsigned int indexLastBlock, Node &node,const std::vector<unsigned int> &assigmentOnStartingBlock,
                                              std::set<std::vector<unsigned int>> &setOfAssignmentOnLastBlock);

    /**
    * Method that creates a child node from a given parent node `pNode` with identical jobs using an assignment.
    * @param pNode The parent node from which we create child nodes.
    * @param indexLastBlock The index of the last block where scheduled jobs are located.
    * @param assignmentLastBlock A vector of machine indices where jobs must be scheduled on the last block.
    * @param assignmentFirstBlock A vector of machine indices where jobs must be scheduled on the first block. This vector is not empty when the first block is the leftmost block (block 0).
    *
    */
    void createNodeWithAssignment(Node &pNode, unsigned int indexLastBlock, const std::vector<unsigned int> &assignmentLastBlock, const std::vector<unsigned int> &assignmentFirstBlock);

    /********************/
    /*      General     */
    /********************/

    /**
     * Method that adds the node to the right structure depending on the walk strategy.
     * @param node The node to add.
     */
    void addNode(Node &node);

    /**
     * Method that generates new nodes by branching from a given node on the all types of machine. The branching scheme
     * use the location. The new nodes are then added to the active list.
     * @param nodeWithLb The parent node from where all children nodes are created with its lower bound.
     */
    void branchingLocation(NodeWithLB &nodeWithLb);

    /**
     * Checks if a block has been filled. If so, it updates all relevant constants.
     * @param node The current node in the tree of branch and bound.
     */
    void changeBlock(Node &node);

    /**
     * Method that compute an upper bound thanks to the MIP.
     * @return if the solution is optimal thanks to the MIP
     */
    bool computeUpperBound();


    /**
     * Method that initialize the tree branch and bound algorithm
     */
    void initialize();

    /**
     * Method that solves the model by using a branch and bound algorithm
     */
    void solve() override;

    /**
     * Method that saves the result of the instance in a file
     * @param fileOutputName The name of the file
     * @param outputFile The stream of the file
     */
    void printOutput(std::string &fileOutputName,std::ofstream &outputFile);

private:

    unsigned long nbNodeLoc; // number of nodes that was created by location scheme
    unsigned long nbCut; // number of nodes that was cut
    unsigned long nbSubProcessBaB=0;
    bool isOptimal = true;
    bool memorizationActivate = true;
    // The stack to store active node in deep-first strategy
    std::vector<NodeWithLB> stackActiveNode;

    // The queue to store active node in breadth-first strategy
    std::queue<NodeWithLB> queueActiveNode;

    // The heap to store active node in best-first strategy
    std::priority_queue<BranchAndBound::NodeWithLB,std::vector<BranchAndBound::NodeWithLB>> heap;

    // column generation model to solve the relaxation of a node
    ColumnGeneration columnGeneration;
    Heuristic heuristicSolver;
    LowerBoundMIP lbFromMIP;
    Memorization memorization;

    // best know solution
    double globalUB;

    Strategy walkStrategy;
    Scheme branchingScheme;
    LowerBound lBStrategy;

};

#ifndef BILEVEL_SCHEDULING_BRANCHANDBOUND_IMP_H
// DO NOT INCLUDE IN OTHER FILE ! IT'S IMPLEMENTATION OF BranchAndBound.h METHODS INLINED
#include "BranchAndBoundImp.hpp"
#endif


#endif //BILEVEL_SCHEDULING_BRANCHANDBOUND_H
