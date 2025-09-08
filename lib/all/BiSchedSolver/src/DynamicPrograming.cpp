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
#include "DynamicPrograming.h"


DynamicPrograming::DynamicPrograming(Instance *instance) : ISolver(instance), absMemoizationSize(computeMemorySize(500)), nbCleaning(0) {}


DynamicPrograming::DynamicPrograming(Instance *instance, nlohmann::json &object) : ISolver(instance, object), absMemoizationSize(
        computeMemorySize(500)), nbCleaning(0) {
    if (object.contains("name")) {
        if (object["name"] == "DP") {
            // check if we have verbose mode
            if (object.contains("verbose")) {
                if (object["verbose"].is_number_unsigned()) setVerbose(object["verbose"].template get<char>());
                else
                    throw std::invalid_argument(
                            R"(The attribute "verbose" of JSON object must be an "unsigned integer" value for the constructor of Dynamic Programing solver)");
            }
            if (object.contains("timeLimits")) {
                if (object["timeLimits"].is_number_unsigned()) setTimeLimit(object["timeLimits"]);
                else
                    throw std::invalid_argument(
                            R"(The attribute "timeLimits" of JSON object must be an "unsigned int" value for the constructor of Dynamic Programing solver)");
            }
            if (object.contains("memorySize")) {
                if (object["memorySize"].is_number_unsigned())
                    absMemoizationSize = computeMemorySize(object["memorySize"]);
                else
                    throw std::invalid_argument(
                            R"(The attribute "memorySize" of JSON object must be an "unsigned int" value for the constructor of Dynamic Programing solver)");
            }
        } else
            throw std::invalid_argument(
                    R"(The JSON object have not the right name to instance a Dynamic Programing solver)");
    } else
        throw std::invalid_argument(R"(The JSON object have not a name to instance a Dynamic Programing solver)");
}


void DynamicPrograming::solve() {

    // start time to measure performance
    start = std::chrono::steady_clock::now();
    if (verbose >= 2) std::cout << "Compute the minimum and maximum for the vector of Sum of Cj" << std::endl;

    std::vector<std::size_t> indexOfJobs(instance->getNbToSelectJob());

    // compute the min SumCj for the n smallest jobs
    std::iota(indexOfJobs.begin(), indexOfJobs.end(), instance->getNbJobs() - instance->getNbToSelectJob());
    auto minMaxCjFromSmallestJob = DynamicPrograming::computeMinMaxSumCompletionTimes(indexOfJobs);

    // compute the min SumCj for the n largest jobs
    std::iota(indexOfJobs.begin(), indexOfJobs.end(), 0);
    auto minMaxCjFromLargestJob = DynamicPrograming::computeMinMaxSumCompletionTimes(indexOfJobs);

    std::array<unsigned int, 4> minMaxCj = {minMaxCjFromSmallestJob[0], minMaxCjFromLargestJob[1],
                                            minMaxCjFromSmallestJob[2], minMaxCjFromLargestJob[3]};
    std::vector<int> C(instance->getNbMachines());

    if (verbose >= 2) {
        std::cout << "Try all possibilities over each vector of Sum of Cj " << std::endl;
        std::cout << "min-max-Cj for High speed = [ " << minMaxCj[0] << " " << minMaxCj[1] << " ]"
                  << std::endl;
        std::cout << "min-max-Cj for low speed = [ " << minMaxCj[2] << " " << minMaxCj[3] << " ]"
                  << std::endl;
    }

    memoization.reserve(absMemoizationSize);
    // try all possibility for the vector of Sum Cj
    solution->setSumWjUj(tryAllSumCj(C, 0, std::numeric_limits<double>::infinity(), minMaxCj));
    memoization.clear();

    // stop time to measure performance
    const auto end{std::chrono::steady_clock::now()};
    time_elapsed = std::chrono::duration<double>{end - start};
    if (verbose >= 1)
        std::cout << "Dynamic Programing is over after " << time_elapsed.count() << " seconds" << std::endl
                  << "The objective value is : " << solution->getSumWjUj() << std::endl;

}


void DynamicPrograming::printOutput(std::string &fileOutputName, std::ofstream &outputFile) {

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
                   << "\t" << "NbCleaning"
                   << "\t" << "MaxNbNode" << std::endl;


    // write value
    outputFile << instance->getInstanceName()
               << "\t" << instance->getInstancePath().string()
               << "\t" << instance->getNbJobs()
               << "\t" << instance->getNbToSelectJob()
               << "\t" << instance->getNbOfHighSpeedMachines()
               << "\t" << instance->getNbOfLowSpeedMachines()
               << "\t" << instance->getHighSpeed()
               << "\t" << instance->getLowSpeed()
               << "\t" << "DP"
               << "\t" << time_elapsed.count()
               << "\t" << time_limits.count()
               << "\t" << solution->getSumWjUj()
               << "\t" << std::to_string(nbCleaning)
               << "\t" << absMemoizationSize << std::endl;

    outputFile.close();

}

size_t DynamicPrograming::computeMemorySize(size_t nbMegabytesAllowedByUser) {
    // estimate the size of unordered map in order to use the memory allow by user
    // the size of a key (maximum size)
    size_t sizeKey = sizeof(int) * (instance->getNbMachines() * 2 + instance->getNbToSelectJob() + 3);
    size_t estimationSizeMap = static_cast<size_t>(
            ((sizeKey + sizeof(void *)) + // data list
             (sizeof(void *) + sizeof(size_t)) / memoization.max_load_factor()) // bucket index
            * 1.5);// estimated allocation overheads

    return (nbMegabytesAllowedByUser * 1000000) / estimationSizeMap;;
}