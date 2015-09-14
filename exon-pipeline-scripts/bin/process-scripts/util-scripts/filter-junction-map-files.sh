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
#$ -v PATH
#$ -pe smp 4 -R y

PROG_NAME=`basename $0`
PROG_DIR=`dirname $0`
VERSION=2.0

set -o pipefail


################################################################################
##
##  Read options
##
################################################################################

CREATE_JUNCTION_JUNCTION_MAP_FILE="FALSE"
OUTPUT_PREFIX=""
FILTERED_JUNCTION_JUNCTION_MAP_FILE=""
while [ "$1" = "-j" -o "$1" = "-o" -o "$1" = "-T" -o "$1" = "--help" -o "$1" = "-h" ]
do
  if [ "$1" = "-h" -o "$1" = "--help" ]
  then
    shift
    HELP="TRUE"
  fi
  
  if [ "$1" = "-o" ]
  then
    shift
    OUTPUT_PREFIX=$1
    shift
  fi
  
  if [ "$1" = "-j" ]
  then
    shift
    FILTERED_JUNCTION_JUNCTION_MAP_FILE="$1"
    shift
  fi

  if [ "$1" = "-T" ]
  then
    shift
    TEMP_DIR="$1"
    shift
  fi

done


################################################################################
##
##  Check number of arguments and show help
##
################################################################################

if [ "$HELP" = "TRUE" -o "$3" = "" ]
then
  echo "Usage: $PROG_NAME [-j <filtered junction junction map file>]"
  echo "  [-o <output prefix>] <project dir> <combined SAM reference id file>"
  echo "  <filtered exon junction map file>"
  echo
  echo "project dir: directory where the relevant files are stored"
  echo "combined SAM file: file with the reference ids of the combined alignments"
  echo "filtered exon junction map file: a file with the exon junction mappings that"
  echo "  contains all possible exon junction mappings that can occur given the"
  echo "  junction SAM file. These junctions are extended to cover junctions that arise"
  echo "  from the fact that the junction Fasta entries are sometimes composed of more"
  echo "  than two exons, the transcript exon junctions, and junction between"
  echo "  identical exons that belong to different genes - the latter need to be"
  echo "  included because again the junction Fasta entries may consist of more than"
  echo "  two exons."
  echo 
  echo " -j STRING: compute the filtered junction junction map file STRING"
  echo " -o STRING: A directory STRING is created in the current working directory"
  echo "     and all intermediate files are stored in a directory structure under"
  echo "     this directory; furthermore, the count files are generated in the current"
  echo "     working directory with the prefix STRING."
  
  exit
fi


################################################################################
##
##  Read arguments PROJECT_DIR, COMBINED_SAM_FILE, FILTERED_EXON_JUNCTION_MAP_FILE
##
################################################################################

PROJECT_DIR=$1
shift

## Read ORGANISM, ORGANISM_SHORT, READ_LENGTH, FILE_BASE, GENE_MODEL, and RADIUS
SETUP_FILE=$PROJECT_DIR/bin/setup.sh
if [ -f $SETUP_FILE ]
then
  source $SETUP_FILE
else
  echo "File $SETUP_FILE not found ... exiting."
  exit 1  
fi
GENE_MODEL_PREFIX=$FILE_BASE

COMBINED_SAM_REFERENCE_ID_FILE="$1"
if [ ! -f "$COMBINED_SAM_REFERENCE_ID_FILE" ]
then
  echo "Combined reference id file $COMBINED_SAM_REFERENCE_ID_FILE not found ... exiting"
  exit 1
fi
COMBINED_SAM_FILE_BASE=`basename $COMBINED_SAM_REFERENCE_ID_FILE | sed -e 's/-reference.ids$//'`
shift


FILTERED_EXON_JUNCTION_MAP_FILE="$1"
shift


################################################################################
##
##  Set and check directories
##
################################################################################

export BIN_DIR=$PROJECT_DIR/exon-pipeline-scripts/bin/process-scripts
if [ ! -d $BIN_DIR ]
then
  echo "Directory $BIN_DIR does not exist ... exiting"
  exit 1
fi

export UTIL_BIN_DIR=$PROJECT_DIR/exon-pipeline-scripts/bin/util-lib
if [ ! -d $UTIL_BIN_DIR ]
then
  echo "Directory $UTIL_BIN_DIR does not exist ... exiting"
  exit 1
fi

PROJECT_MAP_DIR=$PROJECT_DIR/exon-pipeline-files/map-files
if [ ! -d $PROJECT_MAP_DIR ]
then
  echo "Directory $PROJECT_MAP_DIR does not exist ... exiting"
  exit 1
fi


if [ "$OUTPUT_PREFIX" = "" ]
then
  SAM_DIR=`dirname $COMBINED_SAM_REFERENCE_ID_FILE`
  OLD_WD=`pwd`
  cd $SAM_DIR
  SAM_DIR=`pwd`
  cd $OLD_WD
else
  echo "Using directory $OUTPUT_PREFIX for intermediate results"
  if [ ! -d $OUTPUT_PREFIX ]
  then
    echo "Making directory $OUTPUT_PREFIX"
    mkdir $OUTPUT_PREFIX
  fi
  SAM_DIR="$OUTPUT_PREFIX/sam-files"
  if [ ! -d $SAM_DIR ]
  then
    echo "Making directory $SAM_DIR"
    mkdir $SAM_DIR
  fi
fi

REMOVE_TEMP_DIR="FALSE"
if [ "$TEMP_DIR" = "" ]
then
  ## TEMP_DIR is used for sort
  TMP_NAME=`mktemp | sed -e 's;.*/;;'`
  export TEMP_DIR=$SAM_DIR/.$TMP_NAME
  if [ ! -d $TEMP_DIR ]
  then
    mkdir $TEMP_DIR
  fi
  REMOVE_TEMP_DIR="TRUE"
fi


################################################################################
##
##  Set and check files
##
################################################################################

EXON_JUNCTION_MAP_FILE=$PROJECT_MAP_DIR/${GENE_MODEL_PREFIX}_exon_junction.map.gz
if [ ! -f $EXON_JUNCTION_MAP_FILE ]
then
  echo "File $EXON_JUNCTION_MAP_FILE not found ... exiting"
  exit 1
fi

EXON_JUNCTION_TRANSCRIPT_MAP_FILE=$PROJECT_MAP_DIR/${GENE_MODEL_PREFIX}_exon_junction_transcript.map
if [ ! -f $EXON_JUNCTION_TRANSCRIPT_MAP_FILE ]
then
  echo "File $EXON_JUNCTION_TRANSCRIPT_MAP_FILE not found ... exiting"
  exit 1
fi


JUNCTION_JUNCTION_MAP_FILE=$PROJECT_MAP_DIR/${GENE_MODEL_PREFIX}_${RADIUS}_junction_junction_coordinate.map.gz
if [ ! -f $JUNCTION_JUNCTION_MAP_FILE ]
then
  echo "$JUNCTION_JUNCTION_MAP_FILE does not exist ... exiting"
  exit 1
fi


if [ "$CREATE_JUNCTION_JUNCTION_MAP_FILE" = "TRUE" ]
then
  JUNCTION_JUNCTION_TRANSCRIPT_MAP_FILE=$PROJECT_MAP_DIR/${GENE_MODEL_PREFIX}_${RADIUS}_junction_junction_transcript.map.gz
  if [ ! -f $JUNCTION_JUNCTION_TRANSCRIPT_MAP_FILE ]
  then
    echo "File  $JUNCTION_JUNCTION_TRANSCRIPT_MAP_FILE not found ... exiting."
    exit 1
  fi
fi

EQUIVALENT_JUNCTIONS_FILE=$PROJECT_MAP_DIR/${GENE_MODEL_PREFIX}_equivalent_junctions.map
if [ ! -f $EQUIVALENT_JUNCTIONS_FILE ]
then
  echo "$EQUIVALENT_JUNCTIONS_FILE does not exist ... exiting"
  exit 1
fi


################################################################################
##
##  Filter the files
##
################################################################################

echo "Creating extended reference ids"
date
$BIN_DIR/util-scripts/createExtendedJunctionIds.py -j $COMBINED_SAM_REFERENCE_ID_FILE -o $SAM_DIR/$COMBINED_SAM_FILE_BASE-extended-reference.ids
if [ $? -ne 0 ]
then
  echo "Problem with createExtendedJunctionIds.py -j $COMBINED_SAM_REFERENCE_ID_FILE -o $COMBINED_SAM_FILE_BASE-extended-reference.ids"
  echo "... exiting"
  exit 1
fi


echo "Creating external junction ids"
## filterFile.py options: -I -> index of key field, -S -> separator string to the split off a part of an entry
$UTIL_BIN_DIR/filterFile.py -I 0 -S ":" -f $SAM_DIR/$COMBINED_SAM_FILE_BASE-extended-reference.ids -i $JUNCTION_JUNCTION_MAP_FILE | \
     tee $SAM_DIR/$COMBINED_SAM_FILE_BASE-junction-junction.map | cut -f 2 | sort -u -S 4G -T $TEMP_DIR > $SAM_DIR/$COMBINED_SAM_FILE_BASE-external-junction.ids
if [ $? -ne 0 ]
then
  echo "Problem with creating $COMBINED_SAM_FILE_BASE-external-junction.ids ... exiting"
  echo "$UTIL_BIN_DIR/filterFile.py -I 0 -S \":\" -f $SAM_DIR/$COMBINED_SAM_FILE_BASE-extended-reference.ids -i $JUNCTION_JUNCTION_MAP_FILE | \
     tee $SAM_DIR/$COMBINED_SAM_FILE_BASE-junction-junction.map | cut -f 2 | sort -u -S 4G -T $TEMP_DIR > $SAM_DIR/$COMBINED_SAM_FILE_BASE-external-junction.ids"
  exit 1
fi


echo "Adding equivalent junction ids"
$UTIL_BIN_DIR/filterFile.py -f $SAM_DIR/$COMBINED_SAM_FILE_BASE-external-junction.ids -i $EQUIVALENT_JUNCTIONS_FILE > $SAM_DIR/$COMBINED_SAM_FILE_BASE-equivalent-junction.ids
EQUIVALENT_JUNCTION_NUM=`cat $SAM_DIR/$COMBINED_SAM_FILE_BASE-equivalent-junction.ids | wc -l`
if [ $EQUIVALENT_JUNCTION_NUM -gt 0 ]
then
  cut -f 2 $SAM_DIR/$COMBINED_SAM_FILE_BASE-equivalent-junction.ids | cat $SAM_DIR/$COMBINED_SAM_FILE_BASE-external-junction.ids - | \
    sort -u -S 4G -T $TEMP_DIR > $SAM_DIR/$COMBINED_SAM_FILE_BASE-external-junction-all.ids
  if [ $? -ne 0 ]
  then
    echo "Problem with creating $COMBINED_SAM_FILE_BASE-external-junction-all.ids ... exiting"
    exit 1
  fi
  rm $SAM_DIR/$COMBINED_SAM_FILE_BASE-external-junction.ids
else
  mv $SAM_DIR/$COMBINED_SAM_FILE_BASE-external-junction.ids $SAM_DIR/$COMBINED_SAM_FILE_BASE-external-junction-all.ids
fi
rm $SAM_DIR/$COMBINED_SAM_FILE_BASE-equivalent-junction.ids


echo "Combining the mapping junction ids and the exon junction transcript entries"
$UTIL_BIN_DIR/filterFile.py -f $SAM_DIR/$COMBINED_SAM_FILE_BASE-external-junction-all.ids -i $EXON_JUNCTION_MAP_FILE -I 1 | cat $EXON_JUNCTION_TRANSCRIPT_MAP_FILE - | \
     sort -u -S 4G -T $TEMP_DIR | sort -S 4G -k 2 -T $TEMP_DIR > $FILTERED_EXON_JUNCTION_MAP_FILE
if [ $? -ne 0 ]
then
  echo "Problem with creating $FILTERED_EXON_JUNCTION_MAP_FILE ... exiting"
  exit 1
fi


if [ "$FILTERED_JUNCTION_JUNCTION_MAP_FILE" != "" ]
then
  zcat $JUNCTION_JUNCTION_TRANSCRIPT_MAP_FILE | cat - $SAM_DIR/$COMBINED_SAM_FILE_BASE-junction-junction.map | sort -u -S 4G -T $TEMP_DIR > $FILTERED_JUNCTION_JUNCTION_MAP_FILE
  if [ $? -ne 0 ]
  then
    echo "Problem with creating $FILTERED_JUNCTION_JUNCTION_MAP_FILE ... exiting"
    exit 1
  fi
fi

#rm -f $SAM_DIR/$COMBINED_SAM_FILE_BASE-extended-reference.ids $SAM_DIR/$COMBINED_SAM_FILE_BASE-external-junction-all.ids $SAM_DIR/$COMBINED_SAM_FILE_BASE-junction-junction.map

echo "Extended junction ids successfully created."
date


################################################################################
##
##  Cleaning up
##
################################################################################

if [ "$REMOVE_TEMP_DIR" = "TRUE" -a -d $TEMP_DIR ]
then
  rm -r $TEMP_DIR
fi

exit 0
