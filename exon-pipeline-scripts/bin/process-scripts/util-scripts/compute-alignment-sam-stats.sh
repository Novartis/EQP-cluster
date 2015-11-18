#!/bin/bash

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
#$ -S /bin/bash
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
HEADER_OPTION=""
SAMPLE=""
CHUNK=""
while [ "$1" = "-h" -o "$1" = "-s" -o "$1" = "-c" -o "$1" = "-M" -o "$1" = "-H" ]
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

  if [ "$1" = "-H" ]
  then
    shift
    HEADER_OPTION="-H"
  fi
  
  if [ "$1" = "-s" ]
  then
    shift
    SAMPLE="$1"
    shift
  fi
  
  if [ "$1" = "-c" ]
  then
    shift
    CHUNK="$1"
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
  echo "   [-h] [-H]"
  echo
  echo " project dir: directory where the relevant files are stored"
  echo " combined SAM file: the file containing the combined alignments of the reads"
  echo "-H: print header"
  echo "-s STRING: sample id"
  echo "-c STRING: chunk id"
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
  . $SETUP_FILE
else
  echo "File $SETUP_FILE not found ... exiting."
  exit 1  
fi


GENE_MODEL_PREFIX=$FILE_BASE

INPUT_SAM_PATH=$1
INPUT_SAM_FILE=`basename $INPUT_SAM_PATH`
INPUT_SAM_FILE_BASE=`echo $INPUT_SAM_FILE | sed -e 's/.gz$//' | sed -e 's/.[bs]am$//'` 
INPUT_SAM_FILE_EXT=`echo $INPUT_SAM_FILE | sed -e "s/$INPUT_SAM_FILE_BASE[.]//"` 
shift


################################################################################
##
##  Set and check directories
##
################################################################################

SAM_DIR=`dirname $INPUT_SAM_PATH`
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

PROJECT_EXON_DIR=$PROJECT_DIR/exon-pipeline-files
if [ ! -d $PROJECT_EXON_DIR ]
then
  echo "Directory $PROJECT_EXON_DIR does not exist ... exiting"
  exit 1
fi

PROJECT_GENOME_DIR=$PROJECT_EXON_DIR/genome-files
if [ ! -d $PROJECT_GENOME_DIR ]
then
  echo "Directory $PROJECT_GENOME_DIR does not exist ... exiting"
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

if [ "$SAMPLE" = "" ]
then
  SAMPLE=`echo $INPUT_SAM_FILE_BASE | sed -e 's/-C[0-9][0-9][0-9]-.*//'`
fi

if [ "$CHUNK" = "" ]
then
  CHUNK=`echo $INPUT_SAM_FILE_BASE | sed -e 's/.*-\(C[0-9][0-9][0-9]\)-.*/\1/'`
fi

INPUT_SAM_STATS_FILE="$INPUT_SAM_FILE_BASE-sam.stats"
echo "Creating SAM stats file $INPUT_SAM_STATS_FILE"

CMD="ComputeSamAlignmentStatistics -s - -c $PROJECT_GENOME_DIR/genome.fa.fai $HEADER_OPTION -S $SAMPLE -C $CHUNK -o $WEIGHT_DIR/$INPUT_SAM_STATS_FILE"
echo $CMD
zcat $INPUT_SAM_PATH | $JAVA $CMD
if [ $? -ne 0 ]
then
  echo "Computation of the SAM alignment statistics failed ... exiting."
  exit 1
fi

date
