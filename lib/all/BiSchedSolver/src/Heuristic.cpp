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
// Created by schau on 1/22/25.
//

#include "Heuristic.h"

Heuristic::Heuristic() = default;

Heuristic::Heuristic(Instance *instance) :
        ISolver(instance) {}

Heuristic::Heuristic(Instance *instance, nlohmann::json &object) :
        ISolver(instance) {
    if (object.contains("name")) {
        if (object["name"] == "Heuristic") {
            if (object.contains("timeLimits")) {
                if (object["timeLimits"].is_number_unsigned()) setTimeLimit(object["timeLimits"]);
                else throw std::invalid_argument(R"(The attribute "timeLimits" of JSON object must be an "unsigned int" value for the constructor of Heuristic solver)");
            }
            if (object.contains("verbose")) {
                if (object["verbose"].is_number_unsigned()) setVerbose(object["verbose"].template get<char>());
                else throw std::invalid_argument(R"(The attribute "verbose" of JSON object must be an "unsigned integer" value for the constructor of Heuristic solver)");
            }
            if (object.contains("predictor")) {
                if (object["predictor"].is_boolean()) setUsePredictor(object["predictor"].template get<bool>());
                else throw std::invalid_argument(R"(The attribute "predictor" of JSON object must be an "boolean" value for the constructor of Heuristic solver)");
            }
            if (object.contains("trainWeight")) {
                if (object["trainWeight"].is_boolean()) setTrainWeight(object["trainWeight"].template get<bool>());
                else throw std::invalid_argument(R"(The attribute "trainWeight" of JSON object must be an "boolean" value for the constructor of Heuristic solver)");
            }
            if (object.contains("weights")) {
                if (object.contains("predictor")) {
                    if (object["weights"].is_array()) setWeight(object["weights"].template get<std::vector<double>>());
                    else
                        throw std::invalid_argument(
                                R"(The attribute "weights" of JSON object must be an "array" of float value for the constructor of Heuristic solver)");
                } else throw std::invalid_argument(R"(You define weights attribute in JSON object for the heuristic, but you do not define attribute "predictor": true)");
            }
        } else throw std::invalid_argument(R"(The JSON object have not the right name to instance a Heuristic solver)");
    } else throw std::invalid_argument(R"(The JSON object have not a name to instance a Heuristic solver)");
}

Heuristic::~Heuristic() = default;

void Heuristic::localSearch() {
    // start time to measure performance
    const auto start = std::chrono::steady_clock::now();

    // First step, found an initial solution
    // sort jobs according to (p_j - d_j) / w_j rules
    std::vector<Job> listOfJobs = instance->getListJobs();
    auto sortingRule = [](Job &jobLeft, Job &jobRight) { return std::isless((jobLeft.getPi() - jobLeft.getDi()) / jobLeft.getWi(), (jobRight.getPi() - jobRight.getDi()) / jobRight.getWi()); };
    std::sort(listOfJobs.begin(), listOfJobs.end(), sortingRule);
    //select the first n jobs from this list and use a vector of bool to know which jobs that have been already selected
    std::vector<bool> alreadySelectedJobs(instance->getNbJobs(), false);
    std::vector<Job> jobsSelected;
    jobsSelected.reserve(instance->getNbToSelectJob());
    for (unsigned int loopJobSelection = 0; loopJobSelection < instance->getNbToSelectJob(); ++loopJobSelection) {
        jobsSelected.push_back(listOfJobs[loopJobSelection]);
        alreadySelectedJobs[jobsSelected.back().getIndex()].flip();
    }
    std::sort(jobsSelected.begin(), jobsSelected.end(), std::greater<>());
    Solution bestKnowSolution = Solution::solveSumCjCriteria(jobsSelected, instance);
    upgradeSolutionWithHeuristic(bestKnowSolution, jobsSelected);

    const auto end{std::chrono::steady_clock::now()};
    // stop time to measure performance
    time_elapsed = std::chrono::duration<double>{end - start};
    // make local search, by exploring the neighbourhood
    unsigned int maxIter = 20;
    unsigned int nbIter = 0;
    if (verbose >= 2)
        std::cout << "Begin Local search " << std::endl;
    while (nbIter < maxIter && time_elapsed < time_limits) {
        if (verbose >= 2)
            std::cout << "Iteration: " << nbIter << " best know value: " << bestKnowSolution.getSumWjUj() << " time elapsed: " << time_elapsed.count() << " seconds" << std::endl;
        // find the best solution using a heuristic of the neighbourhood
        Solution bestNeighbourhood = Solution(instance);
        // if we don't use predictor, then for each neighbourhood we compute the solution with the heuristic
        std::pair<unsigned int, unsigned int> bestSwap = exploreNeighbourhood_1_OPT(jobsSelected, alreadySelectedJobs, &bestNeighbourhood);
        const auto end_time_exploration{std::chrono::steady_clock::now()};
        // stop time to measure performance
        time_elapsed = std::chrono::duration<double>{end_time_exploration - start};
        // if we use prediction, so we can create solution from swap and update it with heuristic
        if (use_predictor) {
            // swap jobs and compute a solution
            std::vector<Job> newListOfJobs = jobsSelected;
            newListOfJobs[bestSwap.first] = instance->getListJobs()[bestSwap.second];
            std::sort(newListOfJobs.begin(), newListOfJobs.end(), std::greater<>());
            bestNeighbourhood = Solution::solveSumCjCriteria(newListOfJobs, instance);
            upgradeSolutionWithHeuristic(bestNeighbourhood, newListOfJobs);


        }
        // check if the bestNeighbourhood solution is better than the best known solution
        if (std::isless(bestNeighbourhood.getSumWjUj(), bestKnowSolution.getSumWjUj())) {
            bestKnowSolution = bestNeighbourhood;
            // update already selected jobs
            alreadySelectedJobs[jobsSelected[bestSwap.first].getIndex()].flip();
            alreadySelectedJobs[bestSwap.second].flip();
            // swap the jobs from the neighbourhood
            jobsSelected[bestSwap.first] = instance->getListJobs()[bestSwap.second];
            std::sort(jobsSelected.begin(), jobsSelected.end(), std::greater<>()); // sort according LPT

            nbIter++;
        } else break; // we don't need to loop again since we don't have change the list of jobs
    }
    *solution = bestKnowSolution;
}

std::pair<unsigned int, unsigned int> Heuristic::exploreNeighbourhood_1_OPT(std::vector<Job> &initialListJobs, std::vector<bool> &alreadySelectedJobs, Solution *bestSolInNeighbourhood) {
    std::pair<unsigned int, unsigned int> bestSwap{initialListJobs.back().getIndex(), initialListJobs.back().getIndex()}; // initial best swap as swap same element => do nothing

    // Create a best neighbourhood solution if we don't use prediction
    Solution bestNeighbourhoodSolution;
    bestNeighbourhoodSolution.setSumWjUj(instance->getSumWj());

    // Create the best value of predictor if we use prediction
    double bestPrediction = std::numeric_limits<double>::infinity();
    const auto start_time_exploration{std::chrono::steady_clock::now()}; // measure time
    // We use the 1-OPT operator, which consists of swapping a non-selected job with a selected job.
    for (unsigned int loopSelectedJobs = 0; loopSelectedJobs < instance->getNbToSelectJob(); ++loopSelectedJobs) {
        for (unsigned int loopNotSelectedJobs = 0;
             loopNotSelectedJobs < instance->getNbJobs(); ++loopNotSelectedJobs) {
            const auto end_time_exploration{std::chrono::steady_clock::now()};
            // stop time to measure performance
            auto time_elapsed_in_exploration = std::chrono::duration<double>{end_time_exploration - start_time_exploration};
            if (time_elapsed + time_elapsed_in_exploration > time_limits) {
                // if we have bestSolInNeighbourhood != null
                if (bestSolInNeighbourhood) *bestSolInNeighbourhood = bestNeighbourhoodSolution;
                return bestSwap;
            }
            // first check if the jobs is not already selected
            if (!alreadySelectedJobs[loopNotSelectedJobs]) {
                // try the set with swap jobs
                std::vector<Job> newListOfJobs = initialListJobs;
                newListOfJobs[loopSelectedJobs] = instance->getListJobs()[loopNotSelectedJobs];
                std::sort(newListOfJobs.begin(), newListOfJobs.end(), std::greater<>());
                Solution newSol = Solution::solveSumCjCriteria(newListOfJobs, instance);

                // if we don't use prediction then try to get better solution with heuristic
                if (!use_predictor) {
                    upgradeSolutionWithHeuristic(newSol, newListOfJobs);
                    if (std::isless(newSol.getSumWjUj(), bestNeighbourhoodSolution.getSumWjUj())) {
                        bestNeighbourhoodSolution = newSol;
                        bestSwap.first = loopSelectedJobs;
                        bestSwap.second = loopNotSelectedJobs;
                    }
                    if (bestNeighbourhoodSolution.getSumWjUj() < EPSILON) {
                        // if we have bestSolInNeighbourhood != null
                        if (bestSolInNeighbourhood) *bestSolInNeighbourhood = bestNeighbourhoodSolution;
                        return bestSwap;
                    }
                } else {
                    // make prediction
                    double prediction = alpha0 * newSol.getSumWjUj()
                                        + alpha1 * feature1(newSol)
                                        + alpha2 * feature2(newSol)
                                        + alpha3 * feature3(newSol);
                    if (std::isless(prediction, bestPrediction)) {
                        bestPrediction = prediction;
                        bestSwap.first = loopSelectedJobs;
                        bestSwap.second = loopNotSelectedJobs;
                    }
                }
            }
        }
    }
    // if we have bestSolInNeighbourhood != null
    if (bestSolInNeighbourhood) *bestSolInNeighbourhood = bestNeighbourhoodSolution;
    return bestSwap;
}


void Heuristic::solve() {

    // start time to measure performance
    const auto start = std::chrono::steady_clock::now();
    if (trainWeight) {
        auto listOfJobs = instance->getListJobs();
        std::sort(listOfJobs.begin(), listOfJobs.end(), std::greater<>());
        Solution newSol = Solution::solveSumCjCriteria(listOfJobs, instance);
        // make prediction
        double prediction = alpha0 * newSol.getSumWjUj()
                            + alpha1 * feature1(newSol)
                            + alpha2 * feature2(newSol)
                            + alpha3 * feature3(newSol);
        solution->setSumWjUj(prediction);
    } else {
        localSearch();
    }
    const auto end{std::chrono::steady_clock::now()};

    // stop time to measure performance
    time_elapsed = std::chrono::duration<double>{end - start};
    if (verbose >= 1)
        std::cout << "Heuristic is over after " << time_elapsed.count() << " seconds" << std::endl << "The objective value is : " << solution->getSumWjUj() << std::endl;
}


void Heuristic::printOutput(std::string &fileOutputName, std::ofstream &outputFile) {
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
                   << "\t" << "Objective" << std::endl;

    if (not trainWeight) solution->evaluate();
    // write value
    outputFile << instance->getInstanceName()
               << "\t" << instance->getInstancePath().string()
               << "\t" << instance->getNbJobs()
               << "\t" << instance->getNbToSelectJob()
               << "\t" << instance->getNbOfHighSpeedMachines()
               << "\t" << instance->getNbOfLowSpeedMachines()
               << "\t" << instance->getHighSpeed()
               << "\t" << instance->getLowSpeed()
               << "\t" << "Heuristic"
               << "\t" << time_elapsed.count()
               << "\t" << time_limits.count()
               << "\t" << solution->getSumWjUj() << std::endl;
    outputFile.close();
}


