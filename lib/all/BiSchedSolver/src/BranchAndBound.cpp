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

#include "BranchAndBound.h"

BranchAndBound::BranchAndBound() = default;

BranchAndBound::BranchAndBound(Instance *instance) :
        ISolver(instance), nbNodeLoc(0), nbCut(0), columnGeneration(instance), lbFromMIP(instance), memorization(1000000, Memorization::cleaning_t::NOT_USED, instance) {}

BranchAndBound::BranchAndBound(Instance *instance, nlohmann::json &object) :
        ISolver(instance), nbNodeLoc(0), nbCut(0), columnGeneration(instance),lbFromMIP(instance), memorization(1000000, Memorization::cleaning_t::NOT_USED, instance), globalUB(
        instance->getSumWj()) {
    if (object.contains("name")) {
        if (object["name"] == "BaB") {
            // first set the column generation if we have parameters
            bool CG_haveParameters = false;
            bool LB_MIP_haveParameters = false;
            if (object.contains("LB_parameters")) {
                if (object["LB_parameters"].contains("name")) {
                    if (object["LB_parameters"]["name"].is_string()) {
                        if (object["LB_parameters"]["name"] == "CG") {
                            columnGeneration.setParameters(object["LB_parameters"]);
                            CG_haveParameters = true;
                            setLBStrategy("CG");
                        } else if (object["LB_parameters"]["name"] == "LB_MIP") {
                            lbFromMIP.setParameterFromJSON(object["LB_parameters"]);
                            LB_MIP_haveParameters = true;
                            setLBStrategy("LB_MIP");
                        } else throw std::invalid_argument(R"(The JSON object for setting lower bound ,i.e. 'LB_parameters', have 'name' that not corresponding to current implemented lower bound.)");
                    } else throw std::invalid_argument(R"(The JSON object for setting lower bound ,i.e. 'LB_parameters', have 'name' attribute that is not a string)");
                } else throw std::invalid_argument(R"(The JSON object for setting lower bound ,i.e. 'LB_parameters', have not a name, check so documentation to set the lower bound)");
            } else setLBStrategy("CG");

            // check if we have verbose mode
            if (object.contains("verbose")) {
                if (object["verbose"].is_number_unsigned()) {
                    setVerbose(object["verbose"].template get<char>());
                    if (!CG_haveParameters) columnGeneration.setVerbose(object["verbose"].template get<char>());
                    if (!LB_MIP_haveParameters) lbFromMIP.setVerbose(object["verbose"].template get<char>());
                    if (!CG_haveParameters && !LB_MIP_haveParameters && verbose >= 1)
                        std::cout << "Using default LB, i.e., Column Generation" << std::endl;
                } else
                    throw std::invalid_argument(
                            R"(The attribute "verbose" of JSON object must be an "unsigned integer" value for the constructor of BaB solver)");
            }

            if (object.contains("timeLimits")) {
                if (object["timeLimits"].is_number_unsigned()) {
                    setTimeLimit(object["timeLimits"]);
                    columnGeneration.setTimeLimit(object["timeLimits"]);
                } else
                    throw std::invalid_argument(
                            R"(The attribute "timeLimits" of JSON object must be an "unsigned int" value for the constructor of BaB solver)");
            }

            if (object.contains("strategy")) {
                if (object["strategy"].is_string()) setStrategy(object["strategy"]);
                else
                    throw std::invalid_argument(
                            R"(The attribute "strategy" of JSON object must be an "string" value for the constructor of BaB solver)");
            } else walkStrategy = DEPTH_FIRST;

            if (object.contains("scheme")) {
                if (object["scheme"].is_string()) setBranchingScheme(object["scheme"]);
                else
                    throw std::invalid_argument(
                            R"(The attribute "scheme" of JSON object must be an "string" value for the constructor of BaB solver)");
            } else branchingScheme = LOCATION;

            if (object.contains("memorization")) {
                if (object["memorization"].is_boolean()) setMemorizationActivate(object["memorization"]);
                else
                    throw std::invalid_argument(
                            R"(The attribute "memorization" of JSON object must be an "boolean" value for the constructor of BaB solver)");
            };


        } else throw std::invalid_argument(R"(The JSON object have not the right name to instance a BaB solver)");
    } else throw std::invalid_argument(R"(The JSON object have not a name to instance a BaB solver)");
}

BranchAndBound::~BranchAndBound() = default;

bool BranchAndBound::computeUpperBound() {
    Heuristic heuristic(instance);
    MIP optSolMIP = MIP(instance);
    optSolMIP.setTimeLimit(20);
    optSolMIP.solve();
    auto modelStatus = optSolMIP.getStatus();
    if (modelStatus != IloAlgorithm::Feasible && modelStatus != IloAlgorithm::Optimal) {
        if (verbose >= 1) std::cout << "Failed to find a solution with MIP" << std::endl;
        std::vector<Job> listJobs;
        for (unsigned int i = 0; i < instance->getNbToSelectJob(); ++i) {
            listJobs.push_back(instance->getListJobs()[instance->getNbJobs() - i - 1]);
        }
        auto sol = Solution::solveSumCjCriteria(listJobs, instance);
        *solution = sol;
    } else {
        *solution = *optSolMIP.getSolution();
        globalUB = optSolMIP.getSolution()->getSumWjUj();
        if (verbose >= 1)
            std::cout << "MIP found a solution UB=" << globalUB << " with status: " << optSolMIP.getStatus() << std::endl;
    }
    if (modelStatus != IloAlgorithm::Optimal) {
        bool continueUpgrading = true;
        do {
            globalUB = solution->getSumWjUj();
            Solution newSol = *solution;
            heuristic.upgradeSolutionWithHeuristic(newSol, instance->getListJobs());
            if (isSmaller(newSol.getSumWjUj(), globalUB)) *solution = newSol;
            else continueUpgrading = false;
        } while (continueUpgrading);
    }

    //resort instance according SPT
    instance->sort_by_SPT();
    return modelStatus == IloAlgorithm::Optimal;
}

void BranchAndBound::initialize() {
    const auto start = std::chrono::steady_clock::now();
    if (verbose >= 2) std::cout << "Computing the upper bound" << std::endl;
    isOptimal = computeUpperBound();
    // if the solution is optimal than stop the branching scheme
    if (isOptimal) return;
        //else set isOptimal to true, because we will set to false if we use too much time than expected
    else isOptimal = true;
    if (verbose >= 1) std::cout << "Upper bound value : " << globalUB << std::endl;
    // creating the root node
    if (verbose >= 2) std::cout << "Creating the root node" << std::endl;
    Node rootNode = Node(instance);
    columnGeneration.generateStartingColumns(rootNode);
    columnGeneration.clearConstraintOfModel();
    columnGeneration.initializeModel(rootNode);
    columnGeneration.updateValueOfLmax(rootNode, 0);
    columnGeneration.updateValueOfLmax(rootNode, 1);

    #ifdef DEBUG_BaB
    std::string name = "rootNode";
    rootNode.stateDebug.emplace_back(name);
    #endif
    const auto endInit{std::chrono::steady_clock::now()};
    time_elapsed = endInit - start;
    addNode(rootNode);
}

void BranchAndBound::solve() {

    #ifdef DEBUG_DOT
    std::ofstream dot(DEBUG_DOT, std::ios::app);
    dot << "graph search {" << std::endl;
    #endif
    try {
        const auto start = std::chrono::steady_clock::now();
        initialize();
        const auto endInit{std::chrono::steady_clock::now()};
        time_elapsed = std::chrono::duration<double>{endInit - start};
        bool haveActiveNode = (walkStrategy == DEPTH_FIRST) ? !stackActiveNode.empty() : (walkStrategy == BREADTH_FIRST) ? !queueActiveNode.empty() : !heap.empty();
        while (haveActiveNode) {
            NodeWithLB nodeWithLb;
            switch (walkStrategy) {
                case DEPTH_FIRST:nodeWithLb = stackActiveNode.back();
                    stackActiveNode.pop_back();
                    break;
                case BREADTH_FIRST:nodeWithLb = queueActiveNode.front();
                    queueActiveNode.pop();
                    break;
                case BEST_FIRST:nodeWithLb = heap.top();
                    heap.pop();
                    break;
            }
            //if we are in Best_First, and UB - floor(LB) < 1, then we can stop because the smallest value reachable from the UB is UB - 1 < LB
            if (walkStrategy == BEST_FIRST && isSmaller(globalUB - std::floor(nodeWithLb.lowerBound), 1.0))
                break;

            switch (branchingScheme) {
                case LOCATION:branchingLocation(nodeWithLb);
                    break;
            }
            const auto end{std::chrono::steady_clock::now()};
            time_elapsed = std::chrono::duration<double>{end - start};
            haveActiveNode = (walkStrategy == DEPTH_FIRST) ? !stackActiveNode.empty() : (walkStrategy == BREADTH_FIRST) ? !queueActiveNode.empty() : !heap.empty();
            if (time_elapsed.count() > time_limits.count()) {
                throw BiSchTimeOutException();
            }
        }
    } catch (const BiSchTimeOutException &e) {
        isOptimal = false;
    }
    #ifdef DEBUG_DOT
    dot << "}" << std::endl;
    #endif

    solution->evaluate();
    if (verbose >= 1)
        std::cout << "Branch and Bound is over after " << time_elapsed.count() << " seconds" << std::endl
                  << "The objective value is : " << solution->getSumWjUj() << std::endl;
    if (verbose >= 2)
        std::cout << "Nb nodes Loc: " << nbNodeLoc << std::endl
                  << "Nb cuts in memo : " << memorization.getCut() << std::endl
                  << "Nb cleaning of memo : " << memorization.getNbCleaning() << std::endl
                  << "Nb remove in memo : " << memorization.getRem() << std::endl
                  << "Nb cuts : " << nbCut << std::endl;
}

void BranchAndBound::printOutput(std::string &fileOutputName, std::ofstream &outputFile) {
    bool fileExists = std::filesystem::exists(fileOutputName);
    auto filePath = std::filesystem::path(fileOutputName);
    std::filesystem::create_directories(filePath.lexically_normal().parent_path());
    outputFile.open(fileOutputName, std::ios::out | std::ios::app | std::ios::ate);
    // print header
    if (!fileExists)
        outputFile << "InstanceName"
                   << "\t" << "InstancePath"
                   << "\t" << "N"
                   << "\t" << "n"
                   << "\t" << "m_Max"
                   << "\t" << "m_0"
                   << "\t" << "V_max"
                   << "\t" << "V_0"
                   << "\t" << "Method"
                   << "\t" << "Time"
                   << "\t" << "LimitTime"
                   << "\t" << "NBMinColum"
                   << "\t" << "Strategy"
                   << "\t" << "LowerBound"
                   << "\t" << "GenerateCol"
                   << "\t" << "maxNbCallHeuristics"
                   << "\t" << "Memorization"
                   << "\t" << "NBNodes"
                   << "\t" << "NBCut"
                   << "\t" << "NBCleaning"
                   << "\t" << "NBCutInMemo"
                   << "\t" << "NBRemInMemo"
                   << "\t" << "NBComputeCost"
                   << "\t" << "NBCallsHeu"
                   << "\t" << "NBCallDP"
                   << "\t" << "NBCallSubProcessCG"
                   << "\t" << "NBCallSubProcessBaB"
                   << "\t" << "isOptimal"
                   << "\t" << "Objective" << std::endl;

    solution->evaluate();
    // write value
    outputFile << instance->getInstanceName()
               << "\t" << instance->getInstancePath().string()
               << "\t" << instance->getNbJobs()
               << "\t" << instance->getNbToSelectJob()
               << "\t" << instance->getNbOfHighSpeedMachines()
               << "\t" << instance->getNbOfLowSpeedMachines()
               << "\t" << instance->getHighSpeed()
               << "\t" << instance->getLowSpeed()
               << "\t" << "BaB"
               << "\t" << time_elapsed.count()
               << "\t" << time_limits.count()
               << "\t" << columnGeneration.getNbMinStateDp()
               << "\t" << getStrategy()
               << "\t" << getLowerBoundStrategy();
    if (lBStrategy == LowerBound::CG) {
        outputFile << "\t" << static_cast<unsigned int>(columnGeneration.getGenerateColumn()) << "\t" << columnGeneration.getMaxNbCallHeuristic();
    } else {
        outputFile << "\t null \t null \t";
    }
    outputFile << "\t" << memorizationActivate
               << "\t" << nbNodeLoc
               << "\t" << nbCut
               << "\t" << memorization.getNbCleaning()
               << "\t" << memorization.getCut()
               << "\t" << memorization.getRem()
               << "\t" << columnGeneration.getNbCallComputeCost()
               << "\t" << columnGeneration.getNbCallsHeu()
               << "\t" << columnGeneration.getNbCallsDp()
               << "\t" << columnGeneration.getNbCallSubProcessCG()
               << "\t" << nbSubProcessBaB
               << "\t" << isOptimal
               << "\t" << solution->getSumWjUj() << std::endl;
    outputFile.close();
}
