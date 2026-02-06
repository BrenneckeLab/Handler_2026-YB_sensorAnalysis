#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem=4g
#SBATCH -e "%x.e.%j.txt"
#SBATCH -o "%x.o.%j.txt" 
#SBATCH --time=30:00  
#SBATCH --qos=short

# Print the hostname of the compute node
hostname

# Exit if an unset variable is used
set -u

###################################################################################################
###################################################################################################
# Extract variables from the first argument (comma-separated key=value pairs)
VARI=$(echo $1 | sed 's/,/\t/g;s/"//g')
eval $VARI
VARI=$1
# Print each variable on a new line for logging/debugging
echo $VARI | sed 's/,/\n/g'

# Source utility functions from the tools script in SCRIPTdir
source ${SCRIPTdir}tools

# Determine number of threads (double the SLURM_CPUS_PER_TASK)
CORES=$(( $SLURM_CPUS_PER_TASK * 2 ))
# Calculate available memory for the job (subtract 5GB from allocated memory)
MEM=$(scontrol show job $SLURM_JOBID | grep TRES | awk '{ split($NF, X, /,|=|G/);{print X[5]-5}}' | head -n 1)

###################################################################################################
# Setup phase

# Set the raw data directory
rawOPENdir="$OPENdir"
echo $OPENdir
# Create required directories for plots and raw data
mkdir -p ${OPENdir}/plots/
mkdir ${OPENdir}/raw/

# Optionally remove all files in the plots directory (currently commented out)
#@ rm -rf ${OPENdir}/plots/*

# Initiate stats file with header
echo NAME NORM tj sh ratio > ${OPENdir}plots/stats.tmp

#---------------------------------------------------------------------------------------------------------
# Map sRNA data to the sensor sequences

COMMAND=${SCRIPTdir}map-sRNA.sh

if [[ $COMPUTING == C ]]
then
  # If running on a cluster, submit an array job for each library and wait for completion
  sbatch --array=1-$nLIBs --wait $COMMAND ${VARI}
else 
  if [[ $localALL == Y ]]; then
    # If running locally and processing all libraries, loop through each library
    for i in $(seq 1 $nLIBs); do
      currVARI="${VARI},SLURM_ARRAY_TASK_ID=${i}"
      $COMMAND ${currVARI}
    done
  else
    # Otherwise, run the command for the first library only
    VARI="${VARI},SLURM_ARRAY_TASK_ID=1"
    $COMMAND ${VARI}
  fi
fi

#---------------------------------------------------------------------------------------------------------
# Concatenate all sensor background files, keeping only the first header

# Remove any existing ALL.bg file
rm -rf ${OPENdir}raw/ALL.bg
# Loop through each library name listed in library-fasta-pairing.txt
for libNAME in $(cat ${OPENdir}library-fasta-pairing.txt | tr ' ' '\t' | cut -f 1 | tr '\n' '\t'  ); do
  echo $libNAME
  if [[ ! -s ${OPENdir}raw/ALL.bg  ]]; then
    # For the first library, copy the full files (including header)
    cat ${TMPdir}TMP_map-data/${libNAME}/${libNAME}_sense.norm.bg > ${OPENdir}raw/ALL.bg
    cat ${TMPdir}TMP_map-data/${libNAME}/${libNAME}_sense.norm.UTRonly.bg > ${OPENdir}raw/ALL.UTRonly.bg
    echo libNAME miRNA_NORM piRNA_NORM shRNA_NORM QPCRnorm UTRstart UTRend UTRcount | tr ' ' '\t' >${OPENdir}raw/ALL.stats.txt
    cat ${TMPdir}TMP_map-data/${libNAME}/stats.tmp >> ${OPENdir}raw/ALL.stats.txt
    cat ${TMPdir}TMP_map-data/${libNAME}/U-content.txt > ${TMPdir}U-content.txt
  else
    # For subsequent libraries, append data without the header (skip first line)
    tail -n +2 ${TMPdir}TMP_map-data/${libNAME}/${libNAME}_sense.norm.bg >> ${OPENdir}raw/ALL.bg
    tail -n +2 ${TMPdir}TMP_map-data/${libNAME}/${libNAME}_sense.norm.UTRonly.bg >> ${OPENdir}raw/ALL.UTRonly.bg
    cat ${TMPdir}TMP_map-data/${libNAME}/stats.tmp >> ${OPENdir}raw/ALL.stats.txt
    tail -n +2 ${TMPdir}TMP_map-data/${libNAME}/U-content.txt >> ${TMPdir}U-content.txt
  fi
done

# Remove duplicate lines from U-content.txt and save to final output
cat ${TMPdir}U-content.txt | uniq > ${OPENdir}raw/U-content.txt

# Exit the script
exit
