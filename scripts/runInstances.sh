#!/bin/bash

# Number of parallel jobs (set to the number of CPUs or desired concurrency)

# Get the absolute directory of the script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [ "$#" -eq 0 ]
then
	echo "Not enough argument." >&2
	echo "Usage: $0 MODE={local,slurm}" >&2
	exit 1
fi

# Check if "slurm" parameter is provided
if [ "$1" == "slurm" ]; then
    # Path for the tasks file
    tasks_file="$HOME/tasks"
    outDir="$HOME/outputs/"
    mkdir -p $outDir
    # Initialize counter
    counter=1
    
    # Create or overwrite the tasks file
    > "$tasks_file"
    
    # Find all JSON files in the directory
    json_files=$(find "$(realpath $SCRIPT_DIR/../instances/configs/config_solve/)" -type f -name "*.json")
    total_files=$(echo "$json_files" | wc -l)  # Count the total number of files

    # Notify the user
    echo "Found $total_files JSON files to process."
    echo "Creating the tasks file at: $tasks_file"

    counter=1
    bar_length=40

    progress_bar() {
        local progress=$1
        local total=$2
        local filled=$((progress * bar_length / total))
        local empty=$((bar_length - filled))
        printf "\r[%s%s] %d/%d" \
            "$(printf '#%.0s' $(seq 1 $filled))" \
            "$(printf ' %.0s' $(seq 1 $empty))" \
            "$progress" "$total"
    }

    echo "$json_files" | while read -r file; do
        # Show progress
        # Update progress bar
        progress_bar "$counter" "$total_files"

        # Generate the Slurm command and append it to the tasks file
        cmd="$(realpath $SCRIPT_DIR/../bin/Release/bilevel-scheduling) $file > $(realpath $outDir${counter}.txt)"
        echo "$cmd" >> "$tasks_file"
        ((counter++))
    done

    echo "\n Tasks file creation complete. Processed $total_files files."

elif [ "$1" == "local" ]; then
    num_cpus=12
    # Default behavior: Use the provided find method
    find "$(realpath "$SCRIPT_DIR/../instances/configs/config_solve/")" -type f -name "*.json" | \
    parallel -j "$num_cpus" --bar --joblog execution.log \
    "$(realpath "$SCRIPT_DIR/../bin/Release/bilevel-scheduling") {} >> bilevel_scheduling.log 2>&1"

else
    echo "The mode you provided is not defined."
    echo "Accepted values are \"local\" to run on your device. You can adjust the number of CPUs used to run in parallel by modifying the script."
    echo "The other mode is \"slurm\", which generates a file \"$HOME/tasks\" containing all commands. By using the slurm file \"$(realpath $SCRIPT_DIR/../campagne.job)\", it will run in parallel on a Slurm server."
    exit 1
fi

