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

#include "MIP.h"
#include "ilcplex/ilocplex.h"

ILOSTLBEGIN

MIP::MIP() : model(env), obj(env, 0.0, IloInfinity, ILOFLOAT), cplex(model), debug(false), relaxation(false) {}

MIP::MIP(Instance *instance) : ISolver(instance), model(env), obj(env, 0.0, instance->getSumWj(), ILOFLOAT), cplex(model), debug(false), relaxation(false) {}

MIP::MIP(Instance *instance, nlohmann::json &object) : ISolver(instance, object), model(env), obj(env, 0.0, instance->getSumWj(), ILOFLOAT), cplex(model), debug(false), relaxation(false) {
    if (object.contains("name")) {
        if (object["name"] == "MIP" || (object["name"] == "LP" && object["relaxation"].is_boolean() && object["relaxation"])) {
            // check if we have verbose mode
            if (object.contains("verbose")) {
                if (object["verbose"].is_number_unsigned()) setVerbose(object["verbose"].template get<char>());
                else throw std::invalid_argument(R"(The attribute "verbose" of JSON object must be an "unsigned integer" value for the constructor of MIP solver)");
            }
            if (object.contains("debug")) {
                if (object["debug"].is_boolean()) setDebug(object["debug"]);
                else throw std::invalid_argument(R"(The attribute "debug" of JSON object must be an "boolean" value for the constructor of MIP solver)");
            }
            if (object.contains("relaxation")) {
                if (object["relaxation"].is_boolean()) setRelaxation(object["relaxation"]);
                else throw std::invalid_argument(R"(The attribute "relaxation" of JSON object must be an "boolean" value for the constructor of MIP solver)");
            }
            if (object.contains("timeLimits")) {
                if (object["timeLimits"].is_number_unsigned()) setTimeLimit(object["timeLimits"]);
                else throw std::invalid_argument(R"(The attribute "timeLimits" of JSON object must be an "unsigned int" value for the constructor of MIP solver)");
            }
        } else
            throw std::invalid_argument(
                    R"(The JSON object have not the right name to instance a MIP solver, or if you use the LP solver, you must provide the property relaxation to true. Cf the documentation to how use LP.)");
    } else throw std::invalid_argument(R"(The JSON object have not a name to instance a MIP solver)");
}

MIP::~MIP() { env.end(); }


void MIP::solve() {

    // start time to measure performance
    const auto start = std::chrono::steady_clock::now();

    if (verbose >= 2) std::cout << "Initializing MIP" << std::endl;
    initializeMIP();
    parametrizeMIP();
    if (verbose >= 2) std::cout << "Solving MIP" << std::endl;
    solveMIP();

    const auto end{std::chrono::steady_clock::now()};

    // stop time to measure performance
    time_elapsed = std::chrono::duration<double>{end - start};
}

void MIP::initializeMIP() {

    /*******************/
    /*      DATAS      */
    /*******************/

    unsigned int n = instance->getNbToSelectJob();
    unsigned int N = instance->getNbJobs();
    unsigned int m = instance->getNbMachines();
    unsigned int maxPj = static_cast<unsigned int>(instance->getMaxPj());

    // the set of all locations to be optimal for sum Cj
    auto &E = instance->getE();

    /***********************/
    /*      VARIABLES      */
    /***********************/

    // We create x_{i,j,k} variables
    // first the x_i variables
    for (unsigned int i = 0; i < m; ++i) {
        x.emplace_back();
    }

    // loop over each batches
    for (auto const &batch: E) {
        // loop over each position
        for (auto const &position: batch) {
            // loop over all jobs
            for (unsigned int j = 0; j < N; ++j) {
                if (x[position.first].size() < N) x[position.first].emplace_back(env);
                if (relaxation)
                    x[position.first][j].add(IloNumVar(env, 0.0, 1.0, ILOFLOAT));
                else
                    x[position.first][j].add(IloNumVar(env, 0.0, 1.0, ILOBOOL));
            }
        }
    }

    // We create C_{i,k} variables

    // first the C_i variables
    for (unsigned int i = 0; i < m; ++i) {
        C.emplace_back(env);
    }

    for (auto const &batch: E) {
        // loop over each position
        for (auto const &position: batch) {
            C[position.first].add(IloNumVar(env, 0.0, IloInfinity, ILOFLOAT));
        }
    }

    // We create t_{i,j,k} variables

    // first the t_i variables
    for (unsigned int i = 0; i < m; ++i) {
        t.emplace_back();
    }

    // loop over each batches
    for (auto const &batch: E) {
        // loop over each position
        for (auto const &position: batch) {
            // loop over all jobs
            for (unsigned int j = 0; j < N; ++j) {
                if (j <= t[position.first].size()) t[position.first].emplace_back(env);
                if (relaxation) t[position.first][j].add(IloNumVar(env, 0.0, 1.0, ILOFLOAT));
                else t[position.first][j].add(IloNumVar(env, 0.0, 1.0, ILOBOOL));
            }
        }
    }

    /***********************/
    /*      Objective      */
    /***********************/

    // The last variable Obj contains the definition of the objective function value: it has been created as an attribute

    // Step 1: Declaring the objective function
    model.add(IloMinimize(env, obj));
    IloExpr ExprObj(env);
    // loop over each batches
    for (auto const &batch: E) {
        // loop over each position
        for (auto const &position: batch) {
            // loop over all jobs
            for (unsigned int j = 0; j < N; ++j)
                ExprObj += (t[position.first][j][position.second] * int(instance->getListJobs()[j].getWi()));
        }
    }

    ExprObj -= obj;
    model.add(ExprObj == 0);
    ExprObj.end();

    /*************************/
    /*      Constraints      */
    /*************************/

    // Step 2: Declaring the constraints

    // Constraints 1 : One job for each location
    for (unsigned int j = 0; j < N; ++j) {
        IloExpr ExprC1(env);
        // loop over each batches
        for (auto const &batch: E) {
            // loop over each position
            for (auto const &position: batch)
                ExprC1 += x[position.first][j][position.second];
        }
        ExprC1 -= 1;
        model.add(ExprC1 <= 0);
        ExprC1.end();
    }

    // Constraints 2 : Each location can have only one job
    for (auto const &batch: E) {
        // loop over each position
        for (auto const &position: batch) {
            IloExpr ExprC2(env);
            // loop over all jobs
            for (unsigned int j = 0; j < N; ++j)
                ExprC2 += x[position.first][j][position.second];
            ExprC2 -= 1;
            model.add(ExprC2 <= 0);
            ExprC2.end();
        }
    }

    // Constraints 3 : Each location, except maybe the ones on first batch, are affected
    IloExpr ExprC3(env);
    for (unsigned int j = 0; j < N; ++j) {
        for (auto itE = E.begin() + 1; itE != E.end(); ++itE) {
            // loop over each position
            for (auto const &position: *itE)
                ExprC3 += x[position.first][j][position.second];

        }
    }
    ExprC3 -= int(instance->getMaxNbLocation() - E[0].size());
    model.add(ExprC3 == 0);
    ExprC3.end();

    // Constraints 4 : Affect job on location in first batch. there are exactly |B1| + n - |E| jobs inside B1
    IloExpr ExprC4(env);
    for (unsigned int j = 0; j < N; ++j) {
        // loop over each position on first batch
        for (auto const &position: E[0])
            ExprC4 += x[position.first][j][position.second];
    }
    ExprC4 -= int(E[0].size() + n - instance->getMaxNbLocation());
    model.add(ExprC4 == 0);
    ExprC4.end();

    // Constraints 5 : The processing times must be non-decreasing between 2 consecutive batches.
    for (unsigned int h = 0; h < E.size() - 1; ++h) {
        // loop over each position on the batch h
        for (auto const &positionBatchH: E[h]) {
            // loop over each position on the batch h+1
            for (auto const &positionNextBathH: E[h + 1]) {
                IloExpr ExprC5(env);
                for (unsigned int j = 0; j < N; ++j) {
                    ExprC5 += (x[positionBatchH.first][j][positionBatchH.second] * int(instance->getListJobs()[j].getPi()));
                    ExprC5 -= (x[positionNextBathH.first][j][positionNextBathH.second] * int(instance->getListJobs()[j].getPi()));
                }
                model.add(ExprC5 <= 0);
                ExprC5.end();
            }
        }
    }

    // Constraints 6 : Compute completion times
    // loop over each batches
    for (auto const &batch: E) {
        // loop over each position
        for (auto const &position: batch) {
            IloExpr ExprC6(env);
            // If were are not at first batch
            if (position.second > 0) {
                ExprC6 = C[position.first][position.second - 1];
            }
            double speed =
                    position.first < instance->getNbOfHighSpeedMachines() ? instance->getHighSpeed() : instance->getLowSpeed();
            for (unsigned int j = 0; j < N; ++j) {
                auto pj = instance->getListJobs()[j].getPi();
                ExprC6 += (x[position.first][j][position.second] * pj / speed);
            }
            ExprC6 -= C[position.first][position.second];
            model.add(ExprC6 == 0);
            ExprC6.end();
        }
    }

    // Constraints 7 : One job can be tardy
    for (auto const &batch: E) {
        // loop over each position
        for (auto const &position: batch) {
            // loop over all jobs
            for (unsigned int j = 0; j < N; ++j) {
                IloExpr ExprC7(env);
                ExprC7 += (t[position.first][j][position.second] - x[position.first][j][position.second]);
                model.add(ExprC7 <= 0);
                ExprC7.end();
            }
        }
    }

    // Constraints 8 : verify completion time <= due date
    for (auto const &batch: E) {
        // loop over each position
        for (auto const &position: batch) {
            IloExpr ExprC8(env);
            double speed =
                    position.first < instance->getNbOfHighSpeedMachines() ? instance->getHighSpeed() : instance->getLowSpeed();
            double HV = double(position.second + 1) * double(maxPj) / speed;
            // loop over all jobs
            for (unsigned int j = 0; j < N; ++j)
                ExprC8 -= t[position.first][j][position.second];
            ExprC8 *= HV;
            for (unsigned int j = 0; j < N; ++j)
                ExprC8 -= (x[position.first][j][position.second] * int(instance->getListJobs()[j].getDi()));
            ExprC8 += C[position.first][position.second];
            model.add(ExprC8 <= 0);
            ExprC8.end();
        }
    }
}

void MIP::solveMIP() {
    // start time to measure performance
    const auto start = std::chrono::steady_clock::now();
    try {
        if (!cplex.solve()) env.error() << "Failed to optimise MIP" << std::endl;
        else {
            if (verbose >= 2) std::cout << "Status : " << cplex.getStatus() << std::endl;
            std::vector<std::tuple<unsigned int, unsigned int, unsigned int>> listOfMachineJobPositionTuple;

            Solution sol(instance);
            for (unsigned int i = 0; i < x.size(); ++i) {
                for (unsigned int j = 0; j < x[i].size(); ++j) {
                    for (unsigned int k = 0; k < x[i][j].getSize(); ++k) {
                        if (debug)
                            std::cout << "x[" << i << "," << j << "," << k << "]=" << cplex.getValue(x[i][j][k])
                                      << std::endl;
                        // check if x_i,j,k equals 1.0 with 10^-6 precision
                        if ((1.0 - std::abs(cplex.getValue(x[i][j][k]))) < EPSILON)
                            listOfMachineJobPositionTuple.emplace_back(i, j, k);
                    }
                }
            }

            // sort the jobs according to machine and the position
            std::sort(listOfMachineJobPositionTuple.begin(), listOfMachineJobPositionTuple.end(), [](auto &lhs, auto &rhs) {
                return std::get<0>(lhs) == std::get<0>(rhs) ? std::get<2>(lhs) < std::get<2>(rhs) : std::get<0>(lhs) < std::get<0>(rhs);
            });

            for (auto const &[idMachine, idJobs, idPosition]: listOfMachineJobPositionTuple) {
                sol.add_job(idMachine, instance->getListJobs()[idJobs]);
            }
            sol.evaluate();
            if (int(sol.getSumWjUj()) != int(getObjt())) {
                std::cout << "ERROR with objective value = " << getObjt() << " and the one of solution from x_{i,j,k} = " << sol.getSumWjUj() << std::endl;
                sol.setSumWjUj(getObjt());
            }
            *solution = std::move(sol);
        }
        // stop time to measure performance
        const auto end{std::chrono::steady_clock::now()};
        time_elapsed = std::chrono::duration<double>{end - start};
        std::string methodName = (relaxation) ? "LP" : "MIP";
        if (verbose >= 1)
            std::cout << methodName << " is over after " << time_elapsed.count() << " seconds" << std::endl << "The objective value is : " << solution->getSumWjUj() << std::endl;
    } catch (const IloException &e) {
        env.error() << e;
        throw; // if you like
    }
}

void MIP::parametrizeMIP() {

    cplex.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, 0);
    cplex.setParam(IloCplex::Param::Threads, 1);
    cplex.setParam(IloCplex::Param::TimeLimit, time_limits.count());
    cplex.setParam(IloCplex::Param::RandomSeed, 1);
    if (!debug) cplex.setParam(IloCplex::Param::MIP::Display, 0);
    if (verbose <= 1) cplex.setOut(env.getNullStream());

}

void MIP::printOutput(std::string &fileOutputName, std::ofstream &outputFile) {
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
                   << "\t" << "Objective"
                   << "\t" << "isOptimal"
                   << "\t" << "NBNodes" << std::endl;

    // write value
    std::string method = (relaxation) ? "LP" : "MIP";
    outputFile << instance->getInstanceName()
               << "\t" << instance->getInstancePath().string()
               << "\t" << instance->getNbJobs()
               << "\t" << instance->getNbToSelectJob()
               << "\t" << instance->getNbOfHighSpeedMachines()
               << "\t" << instance->getNbOfLowSpeedMachines()
               << "\t" << instance->getHighSpeed()
               << "\t" << instance->getLowSpeed()
               << "\t" << method
               << "\t" << time_elapsed.count()
               << "\t" << time_limits.count()
               << "\t" << solution->getSumWjUj()
               << "\t" << (cplex.getStatus() == IloAlgorithm::Optimal ? 1 : 0)
               << "\t" << getNodes() << std::endl;
    outputFile.close();
}

