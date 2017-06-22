#!/bin/bash

while getopts ":s:p:" opt; do
  case $opt in
    s) step="$OPTARG"
    ;;
    p) params="$OPTARG"
    ;;
  esac
done



## Help display of function

if [ "$1" == "-h" ]; then
    
    cat `dirname $0`/README.txt

  echo "Usage: `basename $0` [-h] [-s n] [-p paramfile.txt]

Pimp My Rads: a program to perform assembly of RAD loci 
          where:
            -s (step)  perform the PmR pipeline according to step number 
               (all: perform all step)
            -p (params)  set the params file

            * Documentation:  https://"
  exit 0
fi




## Processing PimpMyRad using STACKS program

source $params
source `dirname $0`/tools.sh

set -e

if (( $VERBOSE==1 )); then
      set -x
fi

echo ""  
echo ""
echo "Processing PimpMyRads program"
echo ""

   if (( $step == 1 )); then
      echo ""
      echo "############### Step 1 : Adapters trimming ###############"
      echo ""
      source `dirname $0`/src/step_1_filterReads.sh
      exit 0

   elif (( $step == 2 )); then
      echo ""
      echo "############### Step 2 : Data demultiplexing and reads filtering ###############"
      echo ""
      source `dirname $0`/src/step_2_demultiplexing.sh
      exit 0

   elif (( $step == 3 )); then
      echo ""
      echo "############### Step 3 : Samples files preparation before loci reconstruction ###############"
      echo ""
      source `dirname $0`/src/step_3_Prepare_Samples.sh
      exit 0

   elif (( $step == 4 )); then
      echo ""
      echo "############### Step 4 : Range of M threshold for loci reconstruction ###############"
      echo ""
      source `dirname $0`/src/step_4_assembling_loci_gamme.sh
      exit 0

   elif (( $step == 5 )); then
      echo ""
      echo "############### Step 5 : Loci reconstruction for each individuals using a M value of " ${M_CHOOSEN} " ###############"
      echo "" 
      source `dirname $0`/src/step_5_AssemblyStacks.sh
      exit 0
   
   elif (( $step == 6 )); then
      echo ""
      echo "############### Step 6 : Build the catalog of loci ###############"
      echo ""
      source `dirname $0`/src/step_6_build_catalog.sh
      exit 0

   elif (( $step == 7 )); then
      echo ""
      echo "############### Step 7 : Match of individual loci to the catalog ###############"
      echo ""
      source `dirname $0`/src/step_7_match_catalog.sh
      exit 0

   elif (( $step == 8 )); then
      echo ""
      echo "############### Step 8 : Genetic dataset export ###############"
      echo ""
      source `dirname $0`/src/step_8_dataset_exportation.sh
      exit 0

   elif (( $step == all )); then
      
      echo ""
      echo "############### Step 1 : Adapters trimming ###############"
      echo "" 
      source `dirname $0`/src/step_1_filterReads.sh
      echo ""
      echo "############### Step 2 : Data demultiplexing and reads filtering ###############"
      echo ""      
      source `dirname $0`/src/step_2_demultiplexing.sh
      echo ""      
      echo "############### Step 3 : Samples files preparation before loci reconstruction ###############"
      echo ""      
      source `dirname $0`/src/step_3_Prepare_Samples.sh
      echo ""
      echo "############### Step 4 : Range of M threshold for loci reconstruction ###############"
      echo ""
      source `dirname $0`/src/step_4_assembling_loci_gamme.sh
      echo ""
      echo "############### Step 5 : Loci reconstruction for each individuals using a M value of " ${M_CHOOSEN} " ###############"
      echo ""
      source `dirname $0`/src/step_5_AssemblyStacks.sh
      echo ""      
      echo "############### Step 6 : Build the catalog of loci ###############"
      echo ""      
      source `dirname $0`/src/step_6_build_catalog.sh
      echo ""      
      echo "############### Step 7 : Match of individual loci to the catalog ###############"
      echo ""      
      source `dirname $0`/src/step_7_match_catalog.sh
      echo ""      
      echo "############### Step 8 : Genetic dataset export ###############"
      echo ""      
      source `dirname $0`/src/step_8_dataset_exportation.sh
      exit 0

   fi

exit 0

