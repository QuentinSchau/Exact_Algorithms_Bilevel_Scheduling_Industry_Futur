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
// Created by schau on 11/12/24.
//

#include "Memorization.h"

Memorization::Memorization(unsigned int capacity, Memorization::cleaning_t cleaning, const Instance *instance) :
        capacity(capacity), cleaning(cleaning), instance(instance), _rem(0), _cut(0) {
    std::vector<Infos_t> listOfJob;
    listOfJob.reserve(capacity / MAX_NUMBER_JOB);
    database.resize(MAX_NUMBER_JOB + 1, listOfJob);
}

