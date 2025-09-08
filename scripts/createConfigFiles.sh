#!/bin/bash

#
# Copyright (C) 2024
# Laboratoire d'Informatique Fondamentale et Appliquée de Tours, Tours, France
#
# DIGEP, Politecnico di Torino, Corso Duca degli Abruzzi 24, Torino, Italy
# This file is part of bilevel-scheduling.
#
# bilevel-scheduling is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published
# by the Free Software Foundation, either version 3 of the License,
# or (at your option) any later version.
#
# bilevel-scheduling is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with bilevel-scheduling. If not, see <https://www.gnu.org/licenses/>.
#

# Define the range of values for the variables
seed=0
nbInstGen=10
methods=("MIP" "BaB_CG")
n_values=(40 50 60 70 80)
frac_n=(0.25 0.5 0.75)
gen_col=(0 1) # methods for generating column
maxNbCallHeuristic=(1 3 5)
strategies=("depth-first" "breadth-first" "best-first") # Strategies to use for specific methods
nbMinStateDP=(1 3 5)
memo_value=(0 1)
timeLimit=300

# Get the absolute directory of the script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Generation of instances
# Loop over each method to generate json file configuration
for n in "${n_values[@]}"; do
  for method in "${methods[@]}"; do
    # Check if method is "BaB_CG" or "BaB_MIP"
    if [[ "$method" == "BaB_CG" || "$method" == "BaB_MIP" ]]; then
      for memo in "${memo_value[@]}"; do
        # Loop over each strategy
        for strategy in "${strategies[@]}"; do
          if [[ "$method" == "BaB_CG" ]]; then
            for genMethod in "${gen_col[@]}"; do
              for dp in "${nbMinStateDP[@]}"; do
                if [[ "$genMethod" != "0" ]]; then
                  for maxCallHeur in "${maxNbCallHeuristic[@]}"; do
                    echo "Running with method=$method, n=$n, heuristic=$genMethod, maxNbCall=$maxCallHeur, nb backward of dp=$dp and strategy=$strategy"
                    # First instance generation and execution
                    python3 "$(realpath $SCRIPT_DIR/generateInstances.py)" \
                       --memorization $memo\
                      --frac-of-n "${frac_n[@]}" --timeLimit $timeLimit\
                      --m 1,1 --nb-instance-to-generate $nbInstGen \
                      --n-list $n --output-result /instances/N$n/2M_ --methods $method \
                      --strategy $strategy \
                      --use-heuristic-gen-col ${genMethod} --maxNbCallHeuristic ${maxCallHeur} --nbMinStateDP $dp\
                      --config-file-name-solve config_solve/config_solve_N${n}_M2_${method}_${strategy}_memo_${memo}_genCol_${genMethod}_maxNbCalls_${maxCallHeur}_nbBackW_${dp}_ \
                      --path-save-instance /instances/N$n/instances/

                    # Second instance generation and execution
                    python3 "$(realpath $SCRIPT_DIR/generateInstances.py)" \
                       --memorization $memo\
                      --frac-of-n "${frac_n[@]}" --timeLimit $timeLimit \
                      --m 2,2 --nb-instance-to-generate $nbInstGen \
                      --n-list $n --output-result /instances/N$n/4M_ --methods $method \
                      --strategy $strategy \
                      --use-heuristic-gen-col ${genMethod} --maxNbCallHeuristic ${maxCallHeur} --nbMinStateDP $dp\
                      --config-file-name-solve config_solve/config_solve_N${n}_M4_${method}_${strategy}_memo_${memo}_genCol_${genMethod}_maxNbCalls_${maxCallHeur}_nbBackW_${dp}_ \
                      --path-save-instance /instances/N$n/instances/
                  done
                else
                  echo "Running with method=$method, n=$n, heuristic=$genMethod, nb backward of dp=$dp and strategy=$strategy"
                  # First instance generation and execution
                  python3 "$(realpath $SCRIPT_DIR/generateInstances.py)" \
                     --memorization $memo\
                    --frac-of-n "${frac_n[@]}" --timeLimit $timeLimit\
                    --m 1,1 --nb-instance-to-generate $nbInstGen \
                    --n-list $n --output-result /instances/N$n/2M_ --methods $method \
                    --strategy $strategy \
                    --use-heuristic-gen-col ${genMethod} --nbMinStateDP $dp \
                    --config-file-name-solve config_solve/config_solve_N${n}_M2_${method}_${strategy}_memo_${memo}_genCol_${genMethod}_nbBackW_${dp}_ \
                    --path-save-instance /instances/N$n/instances/

                  # Second instance generation and execution
                  python3 "$(realpath $SCRIPT_DIR/generateInstances.py)" \
                     --memorization $memo\
                    --frac-of-n "${frac_n[@]}" --timeLimit $timeLimit \
                    --m 2,2 --nb-instance-to-generate $nbInstGen \
                    --n-list $n --output-result /instances/N$n/4M_ --methods $method \
                    --strategy $strategy \
                    --use-heuristic-gen-col ${genMethod} --nbMinStateDP $dp \
                    --config-file-name-solve config_solve/config_solve_N${n}_M4_${method}_${strategy}_memo_${memo}_genCol_${genMethod}_nbBackW_${dp}_ \
                    --path-save-instance /instances/N$n/instances/
                fi
              done
            done
          else
            echo "Running with method=$method, n=$n, and strategy=$strategy"
            # First instance generation and execution
            python3 "$(realpath $SCRIPT_DIR/generateInstances.py)" \
              --memorization $memo\
              --frac-of-n "${frac_n[@]}" --timeLimit $timeLimit\
              --m 1,1 --nb-instance-to-generate $nbInstGen \
              --n-list $n --output-result /instances/N$n/2M_ --methods $method \
              --strategy $strategy \
              --config-file-name-solve config_solve/config_solve_N${n}_M2_${method}_${strategy}_memo_${memo}_ \
              --path-save-instance /instances/N$n/instances/

            # Second instance generation and execution
            python3 "$(realpath $SCRIPT_DIR/generateInstances.py)" \
              --memorization $memo\
              --frac-of-n "${frac_n[@]}" --timeLimit $timeLimit \
              --m 2,2 --nb-instance-to-generate $nbInstGen \
              --n-list $n --output-result /instances/N$n/4M_ --methods $method \
              --strategy $strategy \
              --config-file-name-solve config_solve/config_solve_N${n}_M4_${method}_${strategy}_memo_${memo}_ \
              --path-save-instance /instances/N$n/instances/
          fi
          echo "Completed for method=$method, n=$n, and strategy=$strategy"
        done
      done
    else
      # If method is "MIP", no strategy is used
      echo "Running with method=$method and n=$n"

      # First instance generation and execution
      python3 "$(realpath $SCRIPT_DIR/generateInstances.py)" \
        --frac-of-n "${frac_n[@]}" --timeLimit $timeLimit \
        --m 1,1 --nb-instance-to-generate $nbInstGen \
        --n-list $n --output-result /instances/N$n/2M_ --methods $method \
        --strategy depth-first \
        --config-file-name-solve config_solve/config_solve_N${n}_M2_${method}_ \
        --path-save-instance /instances/N$n/instances/

      # Second instance generation and execution
      python3 "$(realpath $SCRIPT_DIR/generateInstances.py)" \
        --frac-of-n "${frac_n[@]}" --timeLimit $timeLimit \
        --m 2,2 --nb-instance-to-generate $nbInstGen \
        --n-list $n --output-result /instances/N$n/4M_ --methods $method \
        --strategy depth-first \
        --config-file-name-solve config_solve/config_solve_N${n}_M4_${method}_ \
        --path-save-instance /instances/N$n/instances/

      echo "Completed for method=$method and n=$n"
    fi
  done
done
