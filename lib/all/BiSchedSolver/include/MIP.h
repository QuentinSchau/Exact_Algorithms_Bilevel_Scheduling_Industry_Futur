//  Copyright (C) 2024
//  Laboratoire d'Informatique Fondamentale et Appliquée de Tours, Tours, France
//  DIGEP, Politecnico di Torino, Corso Duca degli Abruzzi 24, Torino, Italy
//
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

#ifndef BILEVEL_SCHEDULING_MIP_H
#define BILEVEL_SCHEDULING_MIP_H


#include <utility>
#include <set>
#include "Math.h"
#include "ISolver.h"
#include "ilconcert/ilomodel.h"
#include "ilcplex/ilocplex.h"
#include "BranchAndBound.h"


class MIP : public ISolver {
public:

    /************************/
    /*      CONSTRUCTOR     */
    /************************/

    explicit MIP();

    explicit MIP(Instance *instance);

    explicit MIP(Instance *instance, nlohmann::json &object);

    /***********************/
    /*      DESTRUCTOR     */
    /***********************/

    ~MIP() override;

    /********************/
    /*      METHODS     */
    /********************/


    /**
     * Method that solves the Bilevel problem.
     */
    void solve() override;

    /**
     * Method that save the result of the instance in a file
     * @param fileOutputName The name of the file
     * @param outputFile The stream of the file
     */
    void printOutput(std::string &fileOutputName, std::ofstream &outputFile);


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

    /********************/
    /*      GETTER      */
    /********************/

    /**
     * The model has been solved
     * @return The value of the objective function value
     */
    double getObjt() { return cplex.getObjValue(); }

    /**
     * Method that retrieves the status of the optimization algorithm.
     * @return The current status of the optimization algorithm, which can be one of several possible values (e.g. IloAlgorithm::Optimal, IloAlgorithm::Infeasible, etc.)
     */
    IloAlgorithm::Status getStatus() { return cplex.getStatus(); }

    /**
     * The model has been solved
     * @return The value of the number of nodes
     */
    long getNodes() { return cplex.getNnodes(); }

    /********************/
    /*      SETTER      */
    /********************/

    void setDebug(bool isDebugged) { MIP::debug = isDebugged; }

    void setRelaxation(bool isRelaxed) { MIP::relaxation = isRelaxed; }

private:
    IloEnv env;
    IloModel model;
    IloNumVar obj;
    IloCplex cplex;
    std::vector<std::vector<IloNumVarArray>> x, t;
    std::vector<IloNumVarArray> C;
    bool debug;
    bool relaxation;
};

#endif //BILEVEL_SCHEDULING_MIP_H
