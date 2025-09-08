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
// Created by schau on 12/4/24.
//

#ifndef BILEVEL_SCHEDULING_LOWERBOUNDMIP_H
#define BILEVEL_SCHEDULING_LOWERBOUNDMIP_H

#include "ISolver.h"
#include "ilconcert/ilomodel.h"
#include "ilcplex/ilocplex.h"
#include "Node.h"

class LowerBoundMIP : public ISolver {
public:
    /************************/
    /*      CONSTRUCTOR     */
    /************************/

    explicit LowerBoundMIP();

    explicit LowerBoundMIP(Instance *instance);

    explicit LowerBoundMIP(Instance *instance, nlohmann::json &object);

    /***********************/
    /*      DESTRUCTOR     */
    /***********************/

    ~LowerBoundMIP() override;

    /********************/
    /*      METHODS     */
    /********************/

    /**
     * Method that computes the lower bound (LB) of a given node in the branch-and-bound algorithm.
     * @param node The node whose LB is to be computed
     * @return The computed LB value
     */
    double computeLB(Node &node);

    /**
     * Method that solves the Bilevel problem.
     */
    void solve() override;

    /**
     * Method that solves the model using a given node, that holds lists of job indices and jobs.
     * @param node The node that contains all lists of indexes and jobs.
     */
    void solve(Node &node);

    /**
     * Method that save the result of the instance in a file
     * @param fileOutputName The name of the file
     * @param outputFile The stream of the file
     */
    void printOutput(std::string &fileOutputName, std::ofstream &outputFile);

    /**
     * Method that updates the MIP model with new information from a given node in the branch-and-bound algorithm.
     * @param node The node whose information is used to update the model
     */
    void updateModelFromNode(Node &node);

    /**
     * Method that initializes the basic MIP formulation
     */
    void initializeMIP();

    /**
     * Method that solves the model that is active in memory.
     */
    void solveMIP();

    /**
     * Method that parametrizes the solver.
     */
    void parametrizeMIP();

    /**
     * Method that computes the solution for the current model and returns it.
     * This method is typically called after the model has been solved..
     */
    void computeSolution();

    /**
     * Method that saves the current MIP model to a file in LP format.
     *
     * @param modelPath The path where the model will be saved
     */
    void saveModel(const std::string &modelPath) {
        cplex.exportModel(modelPath.c_str());
    };

    /********************/
    /*      GETTER      */
    /********************/

    /**
     * The model has been solved
     * @return The value of the objective function value
     */
    double getObjt() { return cplex.getBestObjValue(); }

    IloAlgorithm::Status getStatus() { return cplex.getStatus(); }

    /**
     * The model has been solved
     * @return The value of the number of nodes
     */
    long getNodes() { return cplex.getNnodes(); }

    /********************/
    /*      SETTER      */
    /********************/

    void setDebug(bool isDebugged) { LowerBoundMIP::debug = isDebugged; }

    void setParameterFromJSON(nlohmann::json &object);

private:
    IloEnv env;
    IloModel model;
    IloNumVar obj;
    IloCplex cplex;
    std::vector<std::vector<IloNumVarArray>> x, t;
    std::vector<IloNumVarArray> C;
    bool debug;
    bool haveBeenInitiated = false;

    struct indexVarJob {
        unsigned int i = 0;
        unsigned int j = 0;
        unsigned int k = 0;

        inline indexVarJob(unsigned int i, unsigned int j, unsigned int k) : i(i), j(j), k(k) {}
    };

    struct indexVarCompletionTime {
        unsigned int i = 0;
        unsigned int k = 0;

        inline indexVarCompletionTime(unsigned int i, unsigned int k) : i(i), k(k) {}
    };
    // keep the list of var that was converted to bool instead of float. We update this list at each call of computeLB

    std::vector<indexVarJob> listVarJobModified;
    std::vector<indexVarCompletionTime> listVarCompletionTimeModified;
};

void inline LowerBoundMIP::updateModelFromNode(Node &node) {

    // reset all bounds
    for (auto &indexVarJob: listVarJobModified) {
        t[indexVarJob.i][indexVarJob.j][indexVarJob.k].setBounds(0.0, 1.0);
        x[indexVarJob.i][indexVarJob.j][indexVarJob.k].setBounds(0.0, 1.0);
    }

    for (auto &indexVarCompletionTime: listVarCompletionTimeModified) {
        C[indexVarCompletionTime.i][indexVarCompletionTime.k].setBounds(0.0, instance->getSumPj());
    }

    listVarJobModified.clear();
    listVarCompletionTimeModified.clear();

    // set upper bound and lower bound of the variable according the removed jobs from node

    for (size_t indexLoopEncodingRemoved = 0; indexLoopEncodingRemoved < instance->getNbJobs() - 1; ++indexLoopEncodingRemoved) {
        auto &removedJob = instance->getListJobs()[indexLoopEncodingRemoved];
        for (unsigned int i = 0; i < t.size(); ++i) {
            for (unsigned int k = 0; k < t[i][removedJob.getIndex()].getSize(); ++k) {
                t[i][removedJob.getIndex()][k].setBounds(0.0, 0.0);
                x[i][removedJob.getIndex()][k].setBounds(0.0, 0.0);
                listVarJobModified.emplace_back(i, removedJob.getIndex(), k);
            }
        }
    }
    // set upper bound and lower bound of the variable according the partial scheduling of the given node
    for (unsigned int indexLoopMachine = 0; indexLoopMachine < instance->getNbMachines(); ++indexLoopMachine) {
        auto &machineSchedule = node.getBlockStruc()[indexLoopMachine];
        //loop over block in machine
        for (unsigned int indexBlock = 0; indexBlock < machineSchedule.size(); ++indexBlock) {
            if (machineSchedule[indexBlock].first != nullptr) {
                double Cik = machineSchedule[indexBlock].second;
                const Job *currentJob = machineSchedule[indexBlock].first;
                // set UB of the corresponding x variable to 1
                x[indexLoopMachine][currentJob->getIndex()][indexBlock].setBounds(1.0, 1.0);
                // set the value of the correspond t variable, depending on if the job is late or note
                if (currentJob->getDi() < Cik)
                    t[indexLoopMachine][currentJob->getIndex()][indexBlock].setBounds(1.0, 1.0);
                else
                    t[indexLoopMachine][currentJob->getIndex()][indexBlock].setBounds(0.0, 0.0);
                listVarJobModified.emplace_back(indexLoopMachine, currentJob->getIndex(), indexBlock);
                //set the value of the current completion time
                C[indexLoopMachine][indexBlock].setBounds(Cik, Cik);
                listVarCompletionTimeModified.emplace_back(indexLoopMachine, indexBlock);
            } else {
                // if we are not in the first block and there is no job on the machine then break because we have not more jobs
                if (indexBlock > 0) break;
            }
        }
    }
}

double inline LowerBoundMIP::computeLB(Node &node) {
    try {
        if (!haveBeenInitiated) {
            initializeMIP();
            parametrizeMIP();
        }
        // set bounds of variables depending on the node
        updateModelFromNode(node);
        solveMIP();
        return getObjt();
    } catch (const IloException &e) {
        env.error() << e << std::endl;
        throw; // if you like
    }
}

void inline LowerBoundMIP::solveMIP() {
    try {
        cplex.solve();
    } catch (const IloException &e) {
        env.error() << e;
        throw; // if you like
    }
}

#endif //BILEVEL_SCHEDULING_LOWERBOUNDMIP_H
