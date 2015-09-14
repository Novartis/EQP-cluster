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

PAIRED_END="TRUE"
PRINT_HELP="FALSE"
MAX_MEMORY="25G"
USE_WEIGHT_FILE="FALSE"
WEIGHT_FILE=""
QUANTIFY_OPTION=""
QUIET_OPTION=""
while [ "$1" = "-h" -o "$1" = "--help" -o "$1" = "-M" -o "$1" = "-Q" -o "$1" = "-q" -o "$1" = "-W" -o "$1" = "-w" ]
do

  if [ "$1" = "-h" -o "$1" = "--help" ]
  then
    shift
    PRINT_HELP="TRUE"
  fi

  if [ "$1" = "-w" ]
  then
    shift
    WEIGHT_FILE="$1"
    shift
    USE_WEIGHT_FILE="TRUE"
  fi
  
  if [ "$1" = "-W" ]
  then
    shift
    USE_WEIGHT_FILE="TRUE"
  fi
  
  if [ "$1" = "-Q" ]
  then
    shift
    QUANTIFY_OPTION="-Q"
  fi
  
  if [ "$1" = "-q" ]
  then
    shift
    QUIET_OPTION="-q"
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
  echo "Usage: $PROG_NAME <options> <project dir> <genome SAM file>"
  echo
  echo "where <options> is"
  echo "   [-W|-w <weight file>] [-Q]"
  echo
  echo " project dir: directory where the relevant files are stored"
  echo " genome SAM file: the file containing the spliced alignments of the reads"
  echo "    against the genome."
  echo " -W: Use weight file"
  echo " -w STRING: Use weight file STRING"
  echo " -Q: Quantify intron spanning reads"
  exit
fi


################################################################################
##
##  Read arguments
##
################################################################################

PROJECT_DIR=$1
shift

#set ORGANISM, ORGANISM_SHORT, READ_LENGTH, RADIUS, FILE_BASE, PAIRED_END, and STRAND_SPECIFIC
PAIRED_END="TRUE"
STRAND_SPECIFIC="FALSE"
SETUP_FILE=$PROJECT_DIR/bin/setup.sh
if [ -f $SETUP_FILE ]
then
  source $SETUP_FILE
else
  echo "File $SETUP_FILE not found ... exiting."
  exit 1  
fi
GENE_MODEL_PREFIX=$FILE_BASE


########################### Genome file ######################################
GENOME_SAM_PATH=$1
if [ ! -f $GENOME_SAM_PATH ]
then
  echo "Genome SAM file $GENOME_SAM_PATH not found ... exiting"
  exit 1
fi

GENOME_SAM_FILE=`basename $GENOME_SAM_PATH`
GENOME_SAM_FILE_BASE=`echo $GENOME_SAM_FILE | sed -e 's/.gz$//' | sed -e 's/.[bs]am$//'` 
GENOME_SAM_FILE_EXT=`echo $GENOME_SAM_FILE | sed -e "s/$GENOME_SAM_FILE_BASE[.]//"` 
shift

if [ "$1" != "" ]
then
  echo "Unused arguments: $*."
fi


################################################################################
##
##  Check naming convention
##
################################################################################

if [ "$PAIRED_END" = "TRUE" ]
then
  SAM_FILE_BASE=`echo $GENOME_SAM_FILE_BASE | sed -e "s/-pe//" | sed -e "s/-mixed//"`
  if [ "$SAM_FILE_BASE" = "$GENOME_SAM_FILE_BASE" ]
  then
    echo "$GENOME_SAM_FILE_BASE does not conform to naming convention (-pe|-mixed) ... exiting."
    exit 1
  fi
else
  SAM_FILE_BASE=`echo $GENOME_SAM_FILE_BASE | sed -e "s/-sr//"`
  if [ "$SAM_FILE_BASE" = "$GENOME_SAM_FILE_BASE" ]
  then
    echo "$GENOME_SAM_FILE_BASE does not conform to naming convention (-sr) ... exiting."
    exit 1
  fi
fi


################################################################################
##
##  Set and check directories
##
################################################################################

SAM_DIR=`dirname $GENOME_SAM_PATH`
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

export TOOLS_DIR=$PROJECT_DIR/exon-pipeline-scripts/tools
if [ ! -d $TOOLS_DIR ]
then
  echo "Directory $TOOLS_DIR does not exist ... exiting"
  exit 1
fi
PATH="$TOOLS_DIR:$PATH"

PROJECT_EXON_DIR=$PROJECT_DIR/exon-pipeline-files
if [ ! -d $PROJECT_EXON_DIR ]
then
  echo "Directory $PROJECT_EXON_DIR does not exist ... exiting"
  exit 1
fi

PROJECT_MAP_DIR=$PROJECT_EXON_DIR/map-files
if [ ! -d $PROJECT_MAP_DIR ]
then
  echo "Directory $PROJECT_MAP_DIR does not exist ... exiting"
  exit 1
fi

PROJECT_GTF_DIR=$PROJECT_EXON_DIR/gtf-files
if [ ! -d $PROJECT_GTF_DIR ]
then
  echo "Directory $PROJECT_GTF_DIR does not exist ... exiting"
  exit 1
fi

GTF_FILE="$PROJECT_GTF_DIR/$GENE_MODEL_PREFIX.gtf"
if [ ! -f $GTF_FILE ]
then
  echo "File $GTF_FILE not found ... exiting"
  exit 1
fi

JAVA_CLASS_DIR="$JAVA_DIR/classes:$JAVA_DIR"
JAVA="java -oss8M -ss8M -ms$MAX_MEMORY -mx$MAX_MEMORY -cp ${JAVA_CLASS_DIR}:${CLASSPATH}"
echo "JAVA=$JAVA"


################################################################################
##
##  Create result directories
##
################################################################################

BED_DIR=`echo $SAM_DIR | sed -e "s;/sam-files;/bed-files;"`
if [ ! -d $BED_DIR ]
then
  echo "Making directory $BED_DIR"
  mkdir $BED_DIR
fi
OUTPUT_FILE="$BED_DIR/$GENOME_SAM_FILE_BASE-junction.map"

COUNT_DIR=`echo $SAM_DIR | sed -e "s;/sam-files;/count-files;"`
if [ ! -d $COUNT_DIR ]
then
  echo "Making directory $COUNT_DIR"
  mkdir $COUNT_DIR
fi

WEIGHT_DIR=`echo $SAM_DIR | sed -e "s;/sam-files;/weight-files;"`
if [ "$USE_WEIGHT_FILE" = "TRUE" -a ! -d $WEIGHT_DIR ]
then
  echo "Directory $WEIGHT_DIR not found ... exiting"
  exit 1
fi

if [ "$WEIGHT_FILE" = "" ]
then
  WEIGHT_FILE="$WEIGHT_DIR/$GENOME_SAM_FILE_BASE.wgt"
fi

if [ "$USE_WEIGHT_FILE" = "TRUE" ]
then
  if [ ! -f $WEIGHT_FILE ]
  then
    echo "File $WEIGHT_FILE not found ... exiting"
    exit 1
  else
    WEIGHT_OPTION="-w $WEIGHT_FILE"
  fi
else
  WEIGHT_OPTION=""
fi

date


################################################################################
##
##  Create combined SAM file name
##
################################################################################

JAVA_CMD="ExtractSplicedExonExonIds $QUANTIFY_OPTION $QUIET_OPTION $WEIGHT_OPTION -g $GTF_FILE -o $OUTPUT_FILE -s $GENOME_SAM_PATH"
echo "Java call: $JAVA_CMD"
$JAVA $JAVA_CMD
if [ $? -ne 0 ]
then
  echo "ERROR: Problem with $JAVA $JAVA_CMD ... exiting."
  if [ -f $OUTPUT_FILE ]
  then
    rm $OUTPUT_FILE
  fi
  exit 1
fi
date
echo "Compute intron counts successfully completed."
