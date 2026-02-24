#!/bin/bash

# Get the directory of the current script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Set the parent directory of the script directory
PARENT_DIR="$(dirname "$SCRIPT_DIR")"

# Set the PATH for sratoolkit using an environment variable or default path
SRATOOLKIT_PATH="${SRATOOLKIT_PATH:-$PARENT_DIR/softwares/sratoolkit.3.2.1-ubuntu64/bin}"
if [ ! -d "$SRATOOLKIT_PATH" ]; then
  echo "Error: sratoolkit path $SRATOOLKIT_PATH does not exist."
  exit 1
fi
export PATH="$PATH:$SRATOOLKIT_PATH"

# Configure sratoolkit
vdb-config -i
if [ $? -ne 0 ]; then
  echo "Error: vdb-config failed."
  exit 1
fi

# Check for the correct number of arguments
if [ $# -lt 1 ] || [ $# -gt 2 ]; then
  echo "Usage: $0 <input_file> [number_of_threads]"
  exit 1
fi

# Input file containing the list of SRR accessions
INPUT_FILE="$1"

# Set the number of threads (default to 20 if not specified)
THREADS="${2:-20}"

# Check if the input file exists
if [ ! -f "$INPUT_FILE" ]; then
  echo "Error: Input file $INPUT_FILE not found."
  exit 1
fi

# Function to process each line in the input file
process_line() {
  local accession=$1
  local threads=$2
  echo "Processing $accession with $threads threads"
  
  # Run prefetch
  prefetch -v "$accession"
  if [ $? -ne 0 ]; then
    echo "Error: prefetch failed for $accession."
    return 1
  fi

  # Run fasterq-dump with the specified number of threads
  fasterq-dump "$accession" -e "$threads"
  if [ $? -ne 0 ]; then
    echo "Error: fasterq-dump failed for $accession."
    return 1
  fi

  return 0
}

# Read the input file line by line and process each accession
while read -r line; do 
  process_line "$line" "$THREADS"
done < "$INPUT_FILE"

echo "Processing completed."

#move files to the folder 
mv *.fastq ../../data/fastq_files/

#remove .sra files 
rm -rf "$accession".sra
