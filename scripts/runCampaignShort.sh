#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <input_file> "
    exit 1
fi

INPUT_FILE=$1
OUTPUT_PREFIX=${INPUT_FILE}_split
LINES_PER_FILE=1248

# Check if the input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: File '$INPUT_FILE' not found!"
    exit 1
fi

# Split the file into chunks of 1200 lines each
# The --numeric-suffixes ensures numeric ordering of the output files
# The -l flag specifies the number of lines per file
split -l $LINES_PER_FILE --suffix-length=4 --numeric-suffixes "$INPUT_FILE" "$OUTPUT_PREFIX"
mkdir splits
mv ${OUTPUT_PREFIX}* ./splits/
# Provide feedback to the user
echo "File '$INPUT_FILE' has been split into files with prefix '$OUTPUT_PREFIX'."

# Loop through the generated files and run the sbatch command
for FILE in ./splits/${OUTPUT_PREFIX}*; do
    # Extract the numeric suffix from the file name
    SUFFIX=$(echo "$FILE" | grep -o -E '[0-9]+')

    # Construct the task and log file names
    LOG_FILE="log$SUFFIX.txt"

    # Run the sbatch command
    sbatch --mail-user=mail.to.use@email.com --qos=short -p short campagne.job "$FILE" "$LOG_FILE"

    # Provide feedback for each command
    echo "Submitted sbatch for $TASKS_FILE with log file $LOG_FILE."
done
