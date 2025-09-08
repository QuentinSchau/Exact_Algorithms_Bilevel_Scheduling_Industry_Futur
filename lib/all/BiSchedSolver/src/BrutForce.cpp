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
// Created by schau on 4/11/24.
//
#include "BrutForce.h"

BrutForce::BrutForce(Instance *instance, nlohmann::json &object) : ISolver(instance, object) {
    if (object.contains("name")) {
        if (object["name"] == "BF") {
            // check if we have verbose mode
            if (object.contains("verbose")) {
                if (object["verbose"].is_number_unsigned()) setVerbose(object["verbose"].template get<char>());
                else throw std::invalid_argument(R"(The attribute "verbose" of JSON object must be an "unsigned integer" value for the constructor of Brut Force solver)");
            }
            if (object.contains("timeLimits")) {
                if (object["timeLimits"].is_number_unsigned()) setTimeLimit(object["timeLimits"]);
                else throw std::invalid_argument(R"(The attribute "timeLimits" of JSON object must be an "unsigned int" value for the constructor of Brut Force solver)");
            }
        } else throw std::invalid_argument(R"(The JSON object have not the right name to instance a Brut Force solver)");
    } else throw std::invalid_argument(R"(The JSON object have not a name to instance a Brut Force solver)");
}

void BrutForce::solve() {

    auto start = std::chrono::steady_clock::now();
    if (verbose >= 2) std::cout << "Compute of the weight matrix ..." << std::endl;
    //It's the list that contains all batches
    std::vector<std::vector<std::pair<unsigned int, unsigned int>>> listOfLocationForSumCj;

    instance->computeListAvailableLocation(listOfLocationForSumCj);

    if (verbose >= 2) std::cout << "Compute permutations of possible positions for the jobs ..." << std::endl;
    std::vector<std::vector<std::pair<unsigned int, unsigned int>>> listOfPermutedLocationForSelectedJob;

    generateConcatenatedPermutations(listOfPermutedLocationForSelectedJob, listOfLocationForSumCj);

    // solve for all subVector
    solution->setSumWjUj(std::numeric_limits<double>::infinity());
    std::vector<unsigned int> subVector;
    if (verbose >= 2) std::cout << "Iterate over all feasible subset of jobs ..." << std::endl;

    //solve problem
    solveForAllSubVector(subVector, instance->getListJobs(), listOfPermutedLocationForSelectedJob, instance->getNbToSelectJob(), -1, 0);

    // stop time to measure performance
    const auto end{std::chrono::steady_clock::now()};
    time_elapsed = std::chrono::duration<double>{end - start};

    if (verbose >= 1) std::cout << "Brut Force is over after " << time_elapsed.count() << " seconds" << std::endl << "The objective value is : " << solution->getSumWjUj() << std::endl;
}

void BrutForce::printOutput(std::string &fileOutputName, std::ofstream &outputFile) {

    bool fileExists = std::filesystem::exists(fileOutputName);
    outputFile.open(fileOutputName, std::ios::out | std::ios::app | std::ios::ate);
    // print header, no metrics to add
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
    // write value
    outputFile << instance->getInstanceName()
               << "\t" << instance->getInstancePath().string()
               << "\t" << instance->getNbJobs()
               << "\t" << instance->getNbToSelectJob()
               << "\t" << instance->getNbOfHighSpeedMachines()
               << "\t" << instance->getNbOfLowSpeedMachines()
               << "\t" << instance->getHighSpeed()
               << "\t" << instance->getLowSpeed()
               << "\t" << "BF"
               << "\t" << time_elapsed.count()
               << "\t" << time_limits.count()
               << "\t" << solution->getSumWjUj() << std::endl;
    outputFile.close();
}


