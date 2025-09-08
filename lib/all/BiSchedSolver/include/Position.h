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

#ifndef BILEVEL_SCHEDULING_POSITION_H
#define BILEVEL_SCHEDULING_POSITION_H

#include "Math.h"
#include <vector>
#include "Instance.h"
#include "BiSchException.h"

/**
 * Position of a job on machine
 */
struct Position {
    // set -1 for undefined machine and block
    std::array<char, 2> machines{-1, -1}; // machines where the job can be scheduled
    std::array<int, 2> blocks{-1, -1}; // blocks where the job can be scheduled

    /**
     * Method that returns the position of the block given the machine index and the corresponding instance.
     * This method throw an error, if the index of machine is out of range.
     * @param indexOfMachine The index of the machine.
     * @param instance The instance use to determine the kind of machine from the index.
     * @return the block position
     */
    inline int getBlockIndex(unsigned int indexOfMachine, const Instance *instance) {
        if (instance == nullptr) throw std::invalid_argument("There is no instance to compute indexes");
        if (machines[0] == -1 || blocks[0] == -1) throw BiSchException("Position is empty");
        // if we are on both kind of machines
        if (machines[1] != -1) {
            return indexOfMachine < instance->getNbOfHighSpeedMachines() ? blocks[0] : blocks[1];
        }
        // we are on one high speed machines
        if (indexOfMachine < instance->getNbOfHighSpeedMachines()) {
            return machines[0] == 0 ? blocks[0] : throw BiSchException(
                    "Out of range : the index of machine is for high speed machine, whereas position is on low speed machine");
        } else if (instance->getNbOfHighSpeedMachines() <= indexOfMachine &&
                   indexOfMachine < instance->getNbMachines()) {
            return machines[0] == 1 ? blocks[0] : throw BiSchException(
                    "Out of range : the index of machine is for low speed machine, whereas position is on high speed machine");
        } else {
            throw BiSchException("Out of range : the index of machine is out of range");
        }

    }

    /**
     * Method that returns a pair of the start and the end index of machines given current position.
     * @param instance The instance we use
     * @return (start index, end index). If the position is on high speed machine then, start = 0
     * and end = nb of high speed machine, else if the position is on low speed machine,
     * then start = nb of high speed machine and end = nb of machine. Else, if position if on both kind machines,
     * then start = 0 and end = nb of machines.
     */
    inline std::pair<size_t, size_t> range(const Instance *instance) {
        if (instance == nullptr) throw std::invalid_argument("There is no instance to compute indexes");
        if (machines[0] == -1 || blocks[0] == -1) throw BiSchException("Position is empty");
        size_t startIndex = machines[1] != -1 ?
                            0 :
                            machines[0] == 0 ?
                            0 : instance->getNbOfHighSpeedMachines();
        size_t endIndex = machines[1] != -1 ?
                          instance->getNbMachines() :
                          machines[0] == 0 ?
                          instance->getNbOfHighSpeedMachines() : instance->getNbMachines();
        return std::make_pair(startIndex, endIndex);
    }

    /**
     * Constructor of a position. It is based on the implementation of a function f(n).
     * The entry 'n' is in range 1 to N the number of jobs. This method computes where the n^th job must be scheduled.
     * The machine can be 0 for high speed or 1 for low speed. If job can be scheduled on 2 kind of machine, then we will get
     * two elements in machines and blocks. The blocks vector contains the current position of the job on the machines.
     * If 'n'=0 then the function add 0 on machines and 1 on blocks.
     */
    inline Position(const Instance *instance, unsigned int indexJob) {
        if (indexJob == 0) {
            machines[0] = 0;
            blocks[0] = 1;
        }

        // compute euclidean division
        auto q = static_cast<unsigned int>(instance->getHighSpeed() / instance->getLowSpeed());
        auto r = static_cast<unsigned int>(instance->getHighSpeed()) % static_cast<unsigned int>(instance->getLowSpeed());

        auto alpha = int(indexJob / static_cast<unsigned int>(q * instance->getNbOfHighSpeedMachines() + instance->getNbOfLowSpeedMachines()));
        auto delta = indexJob % static_cast<unsigned int>((q * instance->getNbOfHighSpeedMachines() + instance->getNbOfLowSpeedMachines()));

        auto gamma = delta / instance->getNbOfHighSpeedMachines();
        auto epsilon = delta % instance->getNbOfHighSpeedMachines();

        if (gamma == 0 && epsilon == 0) {
            if (r > 0) {
                machines[0] = 1;
                blocks[0] = alpha;
            }
            if (r == 0) {
                machines[0] = 0;
                machines[1] = 1;
                blocks[0] = alpha * q;
                blocks[1] = alpha;
            }
        } else if (gamma >= q) {
            if (r > 0) {
                if ((gamma - q) * instance->getNbOfHighSpeedMachines() + epsilon >= 1) {
                    machines[0] = 1;
                    blocks[0] = alpha + 1;
                } else {
                    machines[0] = 0;
                    blocks[0] = q * (alpha + 1);
                }
            } else if (r == 0) {
                machines[0] = 0;
                machines[1] = 1;
                blocks[0] = (alpha + 1) * q;
                blocks[1] = alpha + 1;
            }
        } else if (gamma == (q - 1)) {
            if (r > 0) {
                machines[0] = 0;
                blocks[0] = (alpha + 1) * q - 1 +
                            static_cast<int>(std::ceil(double(epsilon) / double(instance->getNbOfHighSpeedMachines())));
            }
            if (r == 0) {
                if (epsilon == 0) {
                    machines[0] = 0;
                    blocks[0] = (alpha + 1) * q - 1;
                } else {
                    machines[0] = 0;
                    machines[1] = 1;
                    blocks[0] = (alpha + 1) * q;
                    blocks[1] = alpha + 1;
                }
            }
        } else if (gamma < (q - 1)) {
            machines[0] = 0;
            blocks[0] = alpha * q + gamma +
                        static_cast<int>(std::ceil(double(epsilon) / double(instance->getNbOfHighSpeedMachines())));
        }
            //if we are here, there is an error
        else throw std::domain_error("A position must have been created");
    }

    inline bool operator==(const Position &rhs) const {
        return machines == rhs.machines &&
               blocks == rhs.blocks;
    }

    inline bool operator!=(const Position &rhs) const {
        return !(rhs == *this);
    }


};

#endif //BILEVEL_SCHEDULING_POSITION_H
