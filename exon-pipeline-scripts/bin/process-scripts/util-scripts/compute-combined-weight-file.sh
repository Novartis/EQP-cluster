#!/bin/sh

## Copyright 2015 Novartis Institutes for BioMedical Research
## Inc.Licensed under the Apache License, Version 2.0 (the "License"); you
## may not use this file except in compliance with the License. You may
## obtain a copy of the License at
##
## http://www.apache.org/licenses/LICENSE-2.0
##
## Unless required by applicable law or agreed to in writing,
## software distributed under the License is distributed on an "AS IS"
## BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
## implied. See the License for the specific language governing
## permissions and limitations under the License. 

#$ -cwd
#$ -j y
#$ -S /bin/sh
#$ -V

PROG_NAME=`basename $0`
PROG_DIR=`dirname $0`
VERSION=2.0

## Ensure that pipes report the exit status of the first failed job
set -o pipefail

## echo $JAVA, echo "$TOOLS ... ", echo rm


################################################################################
##
##  Read options
##
################################################################################

PRINT_HELP="FALSE"
MAX_MEMORY="12G"
while [ "$1" = "-h" -o "$1" = "-M" ]
do

  if [ "$1" = "-h" ]
  then
    shift
    PRINT_HELP="TRUE"
  fi

  if [ "$1" = "-M" ]
  then
    shift
    MAX_MEMORY="$1"
    shift
  fi

done


################################################################################
##
##  Print help
##
################################################################################

if [ "$PRINT_HELP" = "TRUE" -o "$2" = "" ]
then
  echo "Usage: $PROG_NAME <options> <project dir> <combined SAM file>"
  echo
  echo "where <options> is"
  echo "   [-h]"
  echo
  echo " project dir: directory where the relevant files are stored"
  echo " combined SAM file: the file containing the combined alignments of the reads"
  exit
fi


################################################################################
##
##  Read arguments
##
################################################################################

PROJECT_DIR=$1
shift

SETUP_FILE=$PROJECT_DIR/bin/setup.sh
if [ -f $SETUP_FILE ]
then
  source $SETUP_FILE
else
  echo "File $SETUP_FILE not found ... exiting."
  exit 1  
fi


GENE_MODEL_PREFIX=$FILE_BASE


COMBINED_SAM_PATH=$1
COMBINED_SAM_FILE=`basename $COMBINED_SAM_PATH`
COMBINED_SAM_FILE_BASE=`echo $COMBINED_SAM_FILE | sed -e 's/.gz$//' | sed -e 's/.[bs]am$//'` 
COMBINED_SAM_FILE_EXT=`echo $COMBINED_SAM_FILE | sed -e "s/$COMBINED_SAM_FILE_BASE[.]//"` 
shift


################################################################################
##
##  Set and check directories
##
################################################################################

SAM_DIR=`dirname $COMBINED_SAM_PATH`
OLD_WD=`pwd`
cd $SAM_DIR
SAM_DIR=`pwd`
cd $OLD_WD

if [ ! -d $PROJECT_DIR/exon-pipeline-scripts ]
then
  echo "Directory $PROJECT_DIR/exon-pipeline-scripts does not exist ... exiting"
  exit 1
fi

export JAVA_DIR=$PROJECT_DIR/exon-pipeline-scripts/java
if [ ! -d $JAVA_DIR ]
then
  echo "Directory $JAVA_DIR does not exist ... exiting"
  exit 1
fi

export BIN_DIR=$PROJECT_DIR/exon-pipeline-scripts/bin/process-scripts
if [ ! -d $BIN_DIR ]
then
  echo "Directory $BIN_DIR does not exist ... exiting"
  exit 1
fi

export UTIL_BIN_DIR=$BIN_DIR/util-scripts
if [ ! -d $UTIL_BIN_DIR ]
then
  echo "Directory $UTIL_BIN_DIR does not exist ... exiting"
  exit 1
fi

HALF_MAX_MEMORY=`echo $MAX_MEMORY | sed -e 's;G$; / 2;' | bc | sed -e 's/$/G/'`

JAVA_CLASS_DIR="$JAVA_DIR/classes:$JAVA_DIR"
JAVA="java -oss8M -ss8M -ms$HALF_MAX_MEMORY -mx$HALF_MAX_MEMORY -cp ${JAVA_CLASS_DIR}:${CLASSPATH}"
echo "JAVA=$JAVA"



################################################################################
##
##  Create result directories
##
################################################################################

WEIGHT_DIR=`echo $SAM_DIR | sed -e "s;/sam-files;/weight-files;"`
if [ ! -d $WEIGHT_DIR ]
then
  echo "Making directory $WEIGHT_DIR"
  mkdir $WEIGHT_DIR
fi

date

################################################################################
##
##  Compute weights for the combined alignments.
##
################################################################################

COMBINED_SAM_EDIT_DISTANCE_FILE="$COMBINED_SAM_FILE_BASE-sam.wgt"
echo "Creating combined SAM weights file $COMBINED_SAM_EDIT_DISTANCE_FILE"

CMD="ComputeReadWeightsSam -s - -o $WEIGHT_DIR/$COMBINED_SAM_EDIT_DISTANCE_FILE"
echo $CMD
zcat $COMBINED_SAM_PATH | $JAVA $CMD
if [ $? -ne 0 ]
then
  echo "Computation of the combined read weight file failed ... exiting."
  exit 1
fi

date
