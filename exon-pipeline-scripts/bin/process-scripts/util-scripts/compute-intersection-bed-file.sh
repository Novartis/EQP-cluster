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
VERSION=1.0

## Ensure that pipes report the exit status of the first failed job
set -o pipefail

## echo $JAVA, echo "$TOOLS ... ", echo rm


################################################################################
##
##  waitPid
##
################################################################################

waitPid ()
{
  PID="$1"
  shift
  
  OUTPUT_FILE="$1"
  shift
    
  CMD="$1"
  shift

  MESSAGE="$1"
  shift
  if [ "$MESSAGE" = "" ]
  then
    MESSAGE="$CMD"
  fi

  LOG_FILE=""
  if [ "$1" != "" ]
  then
    LOG_FILE="$1"
    shift
  fi

  if [ "$PID" != "" ]
  then
    CONFIRMED_PID=`pgrep -P $$ | grep "^$PID"`
    if [ "$CONFIRMED_PID" = "" ]
    then
      echo "No waiting required: process $PID is no longer a child process of $$"
    fi
    PID="$CONFIRMED_PID"
  fi

  PID_EXIT_STATUS=0
  if [ "$PID" != "" ]
  then
    echo "$MESSAGE (pid $PID)"
    wait $PID
    if [ $? -ne 0 ]
    then
      PID_EXIT_STATUS=1
    fi
  fi
  
  if [ "$LOG_FILE" != "" ]
  then
    if [ -f "$LOG_FILE" ]
    then
      cat $LOG_FILE
      rm  $LOG_FILE
    else
      echo "WARNING: File $LOG_FILE not found."
    fi
  fi

  if [ "$PID_EXIT_STATUS" != "0" ]
  then
    echo "ERROR: Problem with $CMD ... exiting."
    if [ -f $OUTPUT_FILE ]
    then
      rm $OUTPUT_FILE
    fi
    EXIT_STATUS=1
  else
    NUM_LINES=0
    if [ -f "$OUTPUT_FILE" ]
    then
      NUM_LINES=`cat $OUTPUT_FILE | head -100 | wc -l`
    fi
    if [ "$NUM_LINES" = "0" ]
    then
      echo "ERROR: Problem with $CMD."
      echo "File $OUTPUT_FILE does not exist or is empty ... removing and exiting."
      if [ -f "$OUTPUT_FILE" ]
      then
        rm $OUTPUT_FILE
      fi
      PID_EXIT_STATUS=1
      EXIT_STATUS=1
    fi
  fi

  if [ "$PID_EXIT_STATUS" = "0" ]
  then
    echo "Done"
  fi
}

################################################################################
##
## Check version
##
################################################################################

checkVersion ()
{
  TOOL_NAME="$1"
  VERSION1="$2"
  VERSION2="$3"

  VERSION1_FORMAT=`echo $VERSION1 | sed -e 's/^\([0-9][0-9]*\)[.]\([0-9][0-9]*\)[.]\([0-9][0-9]*\).*/okay/'`
  if [ "$VERSION1_FORMAT" != "okay" ]
  then
    echo "Version format: $VERSION1 for $TOOL_NAME not recognized."
    exit 1
  fi
  VERSION2_FORMAT=`echo $VERSION2 | sed -e 's/^\([0-9][0-9]*\)[.]\([0-9][0-9]*\)[.]\([0-9][0-9]*\).*/okay/'`
  if [ "$VERSION2_FORMAT" != "okay" ]
  then
    echo "Version format: $VERSION2 for $TOOL_NAME not recognized."
    exit 1
  fi

  VERSION1_FIRST=`echo $VERSION1 | sed -e 's/^\([0-9]*\).*/\1/'`
  VERSION2_FIRST=`echo $VERSION2 | sed -e 's/^\([0-9]*\).*/\1/'`

  ## echo "VERSION1_FIRST: $VERSION1_FIRST, VERSION2_FIRST: $VERSION2_FIRST"

  if [ "$VERSION1_FIRST" = "" -o "$VERSION1_FIRST" -lt "$VERSION2_FIRST" ]
  then
    echo "Please use $TOOL_NAME version >= $VERSION2. Found version: $VERSION1 ... exiting"
    exit 1
  fi
  
  if [ "$VERSION1_FIRST" -gt "$VERSION2_FIRST" ]
  then
    return
  fi
  
  VERSION1_SECOND=`echo $VERSION1 | sed -e 's/^\([0-9]*\)[.]\([0-9]*\).*/\2/'`
  VERSION2_SECOND=`echo $VERSION2 | sed -e 's/^\([0-9]*\)[.]\([0-9]*\).*/\2/'`

  ## echo "VERSION1_SECOND: $VERSION1_SECOND, VERSION2_SECOND: $VERSION2_SECOND"
  
  if [ "$VERSION1_SECOND" = "" -o "$VERSION1_SECOND" -lt "$VERSION2_SECOND" ]
  then
    echo "Please use $TOOL_NAME version >= $VERSION2. Found version: $VERSION1 ... exiting"
    exit 1
  fi

  if [ "$VERSION1_SECOND" -gt "$VERSION2_SECOND" ]
  then
    return
  fi

  VERSION1_THIRD=`echo $VERSION1 | sed -e 's/^\([0-9]*\)[.]\([0-9]*\)[.]\([0-9]*\).*/\3/'`
  VERSION2_THIRD=`echo $VERSION2 | sed -e 's/^\([0-9]*\)[.]\([0-9]*\)[.]\([0-9]*\).*/\3/'`

  if [ "$VERSION1_THIRD" = "" -o "$VERSION1_THIRD" -lt "$VERSION2_THIRD" ]
  then
    echo "Please use $TOOL_NAME version >= $VERSION2. Found version: $VERSION1 ... exiting"
    exit 1
  fi
}

################################################################################
##
## Check version
##
################################################################################

checkTool ()
{
  TOOL_DIR="$1"
  TOOL_NAME="$2"
  TOOL_VERSION="$3"
  TOOL_EXE="$TOOL_DIR/$TOOL_NAME"
  if [ ! -f $TOOL_EXE ]
  then
    TOOL_EXE=`which $TOOL_NAME 2>&1`
    TOOL_EXE_COMMAND_NOT_FOUND=`echo $TOOL_EXE | sed -e 's/.*Command not found.$/Command not found./' | sed -e "s/^which: no $TOOL_NAME in .*/Command not found./"`
    if [ "$TOOL_EXE_COMMAND_NOT_FOUND" = "Command not found." ]
    then
      echo "$TOOL_NAME not found in PATH. Please make sure that $TOOL_NAME (>=$TOOL_VERSION) is"
      echo "available ... exiting"
      exit 1
    else
      if [ "$TOOL_NAME" = "bedtools" ]
      then
        ACTUAL_TOOL_VERSION=`$TOOL_EXE --version | sed -e "s/$TOOL_NAME v*//"`
      else
        echo "Cannot check tool $TOOL_NAME"
	return
      fi
      ## echo "ACTUAL_TOOL_VERSION: $ACTUAL_TOOL_VERSION"
      checkVersion $TOOL_NAME $ACTUAL_TOOL_VERSION $TOOL_VERSION
    fi
  fi
}


################################################################################
##
##  Read options
##
################################################################################

RECOMPUTE="FALSE"
USE_GTF=FALSE
EDIT_DISTANCE_OPTION=
STRAND_SPECIFIC_OPTION="FALSE"
STRAND_SPECIFIC_DIRECTION_OPTION=""
OUTPUT_PREFIX=""
OUTPUT_PREFIX_OPTION=""
COMPUTE_COUNTS="FALSE"
WAIT="FALSE"
while [ "$1" = "-r" -o "$1" = "-L" -o "$1" = "-gtf" -o "$1" = "-d" -o "$1" = "-G" -o "$1" = "-u" -o \
        "$1" = "-s" -o "$1" = "-o" -o "$1" = "-c" -o "$1" = "-w" ]
do

  if [ "$1" = "-r" ]
  then
    shift
    RECOMPUTE="TRUE"
  fi
  
  if [ "$1" = "-gtf" ]
  then
    shift
    USE_GTF=TRUE
  fi

  if [ "$1" = "-d" ]
  then
    shift
    EDIT_DISTANCE_OPTION="-e $1"
    shift
  fi

  if [ "$1" = "-c" ]
  then
    shift
    COMPUTE_COUNTS="TRUE"
  fi

  if [ "$1" = "-s" ]
  then
    shift
    STRAND_SPECIFIC_OPTION="TRUE"
    STRAND_SPECIFIC_DIRECTION_OPTION="$1"
    shift
  fi

  if [ "$1" = "-o" ]
  then
    shift
    OUTPUT_PREFIX=$1
    shift
    OUTPUT_PREFIX_OPTION="-o $OUTPUT_PREFIX"
  fi
  
  if [ "$1" = "-w" ]
  then
    shift
    WAIT="TRUE"
  fi

done


################################################################################
##
##  Read arguments
##
################################################################################

if [ "$4" = "" ]
then
  echo "Usage: $PROG_NAME <options> <project dir> <transcript SAM file>"
  echo "        <junction SAM file> <genome SAM-file>"
  echo
  echo "where <options> is"
  echo "   [-r] [-gtf] [-d <edit distance>] [-s <direction>]"
  echo "   [-o <output-prefix>]"
  echo
  echo " project dir: directory where the relevant files are stored"
  echo " transcript SAM file: the file containing the alignments of the"
  echo "     reads against the transcripts."
  echo " junction SAM file: the file containing the alignments of the"
  echo "     reads against the junctions."
  echo " genome SAM file: the file containing the alignments of the reads against"
  echo "    the genome."
  echo " -r: force the recomputation of SAM and BED files even if they exist"
  echo " -gtf: use GTF transcript file instead of the transcript file downloaded"
  echo "       from $GENE_MODEL"
  echo " -s STRING: direction of how to process the reads as strand-specific: forward"
  echo "     or backward"
  echo " -d INT: edit distance threshold in conversion of SAM to BED file (only SAM"
  echo "     records with an edit distance of at most INT are converted)"
  echo " -o STRING: A directory STRING is created in the current working directory"
  echo "     and all intermediate files are stored in a directory structure under"
  echo "     this directory; furthermore, the count files are generated in the current"
  echo "     working directory with the prefix STRING."
  echo " -w: Wait until the compression of the combined SAM file is finished"
  exit
fi

PROJECT_DIR=$1
shift

## Read ORGANISM, ORGANISM_SHORT, READ_LENGTH, FILE_BASE, GENE_MODEL, and RADIUS
SETUP_FILE=$PROJECT_DIR/bin/setup.sh
if [ -f $SETUP_FILE ]
then
  . $SETUP_FILE
else
  echo "File $SETUP_FILE not found ... exiting."
  exit 1  
fi

GENE_MODEL_PREFIX=$FILE_BASE

BED_STRAND_SPECIFIC_OPTION=""
CONVERT_SAM_BED_STRAND_SPECIFIC_OPTION=""
if [ "$STRAND_SPECIFIC" = "TRUE" -o "$STRAND_SPECIFIC_OPTION" = "TRUE" ]
then
  if [ "$STRAND_SPECIFIC_OPTION" = "TRUE" ]
  then
    STRAND_SPECIFIC="TRUE"
    STRAND_SPECIFIC_DIRECTION="$STRAND_SPECIFIC_DIRECTION_OPTION"
  fi
  if [ "$STRAND_SPECIFIC_DIRECTION" != "forward" -a "$STRAND_SPECIFIC_DIRECTION" != "backward" ]
  then
    echo "Unknown direction for strand-specific processing: $STRAND_SPECIFIC_DIRECTION"
    echo "Please set the direction to forward or backward (option -s) ... exiting."
    exit 1
  fi
  CONVERT_SAM_BED_STRAND_SPECIFIC_OPTION="-S $STRAND_SPECIFIC_DIRECTION"

  ## Note that in the conversion to the BED file the read alignments are changed to the
  ## correct orientation; hence, we can use the intersectBed option for the correct
  ## orientation (-s). The option for the reverse orientation (-S) is not necessary.
  BED_STRAND_SPECIFIC_OPTION="-s"
fi


######################### Transcript file ######################################
TRANSCRIPT_SAM_PATH=$1
shift
if [ $RECOMPUTE = "TRUE" -a ! -f "$TRANSCRIPT_SAM_PATH" ]
then
  echo "Transcript SAM file $TRANSCRIPT_SAM_PATH not found ... exiting"
  exit 1
fi

TRANSCRIPT_SAM_FILE=`basename $TRANSCRIPT_SAM_PATH`
TRANSCRIPT_SAM_FILE_BASE=`echo $TRANSCRIPT_SAM_FILE | sed -e 's/.gz$//' | sed -e 's/.[bs]am$//'` 
TRANSCRIPT_SAM_FILE_EXT=`echo $TRANSCRIPT_SAM_FILE | sed -e "s/$TRANSCRIPT_SAM_FILE_BASE[.]//"` 


########################### Junction file ######################################
JUNCTION_SAM_PATH=$1
shift
if [ $RECOMPUTE = "TRUE" -a ! -f "$JUNCTION_SAM_PATH" -a "$JUNCTION_FILE_MISSING" != "TRUE" ]
then
  echo "Junction SAM file $JUNCTION_SAM_PATH not found ... exiting"
  exit 1
fi

if [ "$JUNCTION_FILE_MISSING" != "TRUE" ]
then
  JUNCTION_SAM_FILE=`basename $JUNCTION_SAM_PATH`
  JUNCTION_SAM_FILE_BASE=`echo $JUNCTION_SAM_FILE | sed -e 's/.gz$//' | sed -e 's/.[bs]am$//'` 
  JUNCTION_SAM_FILE_EXT=`echo $JUNCTION_SAM_FILE | sed -e "s/$JUNCTION_SAM_FILE_BASE[.]//"` 
fi


########################### Genome file ######################################
GENOME_SAM_PATH=$1
shift
if [ $RECOMPUTE = "TRUE" -a ! -f "$GENOME_SAM_PATH" -a "$GENOME_FILE_MISSING" != "TRUE" ]
then
  echo "Genome SAM file $GENOME_SAM_PATH not found ... exiting"
  exit 1
fi

if [ "$GENOME_FILE_MISSING" != "TRUE" ]
then
  GENOME_SAM_FILE=`basename $GENOME_SAM_PATH`
  GENOME_SAM_FILE_BASE=`echo $GENOME_SAM_FILE | sed -e 's/.gz$//' | sed -e 's/.[bs]am$//'` 
  GENOME_SAM_FILE_EXT=`echo $GENOME_SAM_FILE | sed -e "s/$GENOME_SAM_FILE_BASE[.]//"` 
fi

if [ "$1" != "" ]
then
  echo "Unused arguments: $*."
fi


################################################################################
##
##  Set and check directories
##
################################################################################

if [ "$OUTPUT_PREFIX" = "" ]
then
  SAM_DIR=`dirname $TRANSCRIPT_SAM_PATH`
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

export UTIL_LIB_DIR=$PROJECT_DIR/exon-pipeline-scripts/bin/util-lib
if [ ! -d $UTIL_LIB_DIR ]
then
  echo "Directory $UTIL_LIB_DIR does not exist ... exiting"
  exit 1
fi


export TOOLS_DIR=$PROJECT_DIR/exon-pipeline-scripts/tools
if [ -d $TOOLS_DIR ]
then
  PATH="${TOOLS_DIR}:$PATH"
fi
checkTool $TOOLS_DIR bedtools 2.24.0
BEDTOOLS_EXE=$TOOL_EXE

PROJECT_EXON_DIR=$PROJECT_DIR/exon-pipeline-files
if [ ! -d $PROJECT_EXON_DIR ]
then
  echo "Directory $PROJECT_EXON_DIR does not exist ... exiting"
  exit 1
fi

PROJECT_BED_DIR=$PROJECT_EXON_DIR/bed-files
if [ ! -d $PROJECT_BED_DIR ]
then
  echo "Directory $PROJECT_BED_DIR does not exist ... exiting"
  exit 1
fi

PROJECT_MAP_DIR=$PROJECT_EXON_DIR/map-files
if [ ! -d $PROJECT_MAP_DIR ]
then
  echo "Directory $PROJECT_MAP_DIR does not exist ... exiting"
  exit 1
fi

PROJECT_GENOME_DIR=$PROJECT_EXON_DIR/genome-files
if [ ! -d $PROJECT_GENOME_DIR ]
then
  echo "Directory $PROJECT_GENOME_DIR does not exist ... exiting"
  exit 1
fi

JAVA_CLASS_DIR="$JAVA_DIR/classes:$JAVA_DIR"
JAVA="java -oss8M -ss8M -ms3G -mx3G -cp ${JAVA_CLASS_DIR}:${CLASSPATH}"
echo "Using $JAVA"


################################################################################
##
##  Set BED files and map files
##
################################################################################

CHROMOSOME_ID_FILE=$PROJECT_GENOME_DIR/genome.fa.fai
if [ ! -f $CHROMOSOME_ID_FILE ]
then
  echo "File $CHROMOSOME_ID_FILE not found ... exiting"
  exit 1
fi

if [ "$USE_GTF" = "FALSE" ]
then
  COMBINED_EXON_BED_FILE=$PROJECT_BED_DIR/${GENE_MODEL_PREFIX}_${RADIUS}_combined_exons.bed
else
  COMBINED_EXON_BED_FILE=$PROJECT_BED_DIR/${GENE_MODEL_PREFIX}_${RADIUS}_combined_exons_gtf.bed
fi
GENOME_EXON_BED_FILE=$PROJECT_BED_DIR/${GENE_MODEL_PREFIX}_genome_exons.bed

if [ ! -f $COMBINED_EXON_BED_FILE ]
then
  echo "$COMBINED_EXON_BED_FILE not found ... exiting"
  exit 1
fi

TRANSCRIPT_BED_FILE=$PROJECT_BED_DIR/${GENE_MODEL_PREFIX}_transcript_exons.bed
if [ ! -f $TRANSCRIPT_BED_FILE ]
then
  echo "$TRANSCRIPT_BED_FILE not found ... exiting"
  exit 1
fi


################################################################################
##
##  Create directories for intermediate and final results
##
################################################################################

BED_DIR=`echo $SAM_DIR | sed -e "s/sam-files/bed-files/"`
if [ ! -d $BED_DIR ]
then
  echo "Making directory $BED_DIR"
  mkdir -p $BED_DIR
fi

WEIGHT_DIR=`echo $SAM_DIR | sed -e "s/sam-files/weight-files/"`
if [ ! -d $WEIGHT_DIR ]
then
  echo "Making directory $WEIGHT_DIR"
  mkdir -p $WEIGHT_DIR
fi


################################################################################
##
##  Start the computation of SAM files, weight files and counts files
##
################################################################################

TRANSCRIPT_SAM_EDIT_DISTANCE_FILE="$TRANSCRIPT_SAM_FILE_BASE-sam.wgt"
## Test if GENOME_SAM_WEIGHT_FILE is already set by the calling script
if [ "$GENOME_SAM_WEIGHT_FILE" != "" ]
then
  GENOME_SAM_EDIT_DISTANCE_FILE="$GENOME_SAM_WEIGHT_FILE"
else
  GENOME_SAM_EDIT_DISTANCE_FILE="$GENOME_SAM_FILE_BASE-sam.wgt"
  GENOME_SAM_WEIGHT_FILE="$GENOME_SAM_EDIT_DISTANCE_FILE"
fi
JUNCTION_SAM_EDIT_DISTANCE_FILE="$JUNCTION_SAM_FILE_BASE-sam.wgt"

if [ "$RECOMPUTE" = "TRUE" ]
then
  rm -f $WEIGHT_DIR/$TRANSCRIPT_SAM_EDIT_DISTANCE_FILE
  rm -f $WEIGHT_DIR/$GENOME_SAM_EDIT_DISTANCE_FILE
  rm -f $WEIGHT_DIR/$JUNCTION_SAM_EDIT_DISTANCE_FILE
fi


################################################################################
##
##  Compute weights for the transcript and genome alignments. Note that the
##  chromosome id file ensures that we only consider alignments against
##  chromosomes even if the alignment file contains other alignments, e.g.
##  against transcripts or junctions. This file is later needed for the
##  computation of the read weights.
##
################################################################################

EXIT_STATUS="0"

echo "Start of weight computation"
if [ ! -r $WEIGHT_DIR/$TRANSCRIPT_SAM_EDIT_DISTANCE_FILE ]
then
  echo "Creating transcript SAM weights $TRANSCRIPT_SAM_EDIT_DISTANCE_FILE" > $WEIGHT_DIR/transcript.log
  CMD="ComputeReadWeightsSam -s - -o $WEIGHT_DIR/$TRANSCRIPT_SAM_EDIT_DISTANCE_FILE"
  echo $CMD > $WEIGHT_DIR/transcript.log
  zcat $TRANSCRIPT_SAM_PATH | $JAVA $CMD >> $WEIGHT_DIR/transcript.log 2>&1 &
  TRANSCRIPT_WEIGHT_PID=$!
fi

if [ ! -r $WEIGHT_DIR/$JUNCTION_SAM_EDIT_DISTANCE_FILE -a "$JUNCTION_FILE_MISSING" != "TRUE" ]
then
  echo "Creating junction SAM weights $JUNCTION_SAM_EDIT_DISTANCE_FILE" > $WEIGHT_DIR/junction.log
  CMD="ComputeReadWeightsSam -s - -o $WEIGHT_DIR/$JUNCTION_SAM_EDIT_DISTANCE_FILE"
  echo $CMD >> $WEIGHT_DIR/junction.log
  zcat $JUNCTION_SAM_PATH | $JAVA $CMD >> $WEIGHT_DIR/junction.log 2>&1 &
  JUNCTION_WEIGHT_PID=$!
fi

if [ ! -r $WEIGHT_DIR/$GENOME_SAM_WEIGHT_FILE -a "$GENOME_FILE_MISSING" != "TRUE" ]
then
  echo "Creating genome SAM weight file $GENOME_SAM_WEIGHT_FILE" > $WEIGHT_DIR/genome.log
  CMD="ComputeReadWeightsSam -c $CHROMOSOME_ID_FILE -s - -o $WEIGHT_DIR/$GENOME_SAM_WEIGHT_FILE"
  echo $CMD >> $WEIGHT_DIR/genome.log
  zcat $GENOME_SAM_PATH | $JAVA $CMD >> $WEIGHT_DIR/genome.log 2>&1  &
  GENOME_WEIGHT_PID=$!
fi


################################################################################
##
##  Wait for the weight processes to finish and get their exit status
##
################################################################################

## Get exit status of transcript weight file computation
if [ "$TRANSCRIPT_WEIGHT_PID" != "" ]
then
  waitPid $TRANSCRIPT_WEIGHT_PID $WEIGHT_DIR/$TRANSCRIPT_SAM_EDIT_DISTANCE_FILE "Compute transcript SAM edit distance" "" $WEIGHT_DIR/transcript.log
fi

## Junction weight file computation
if [ "$JUNCTION_WEIGHT_PID" != "" ]
then
  waitPid $JUNCTION_WEIGHT_PID $WEIGHT_DIR/$JUNCTION_SAM_EDIT_DISTANCE_FILE "Compute junction SAM edit distance" "" $WEIGHT_DIR/junction.log
fi

## Genome weight file computation
if [ "$GENOME_WEIGHT_PID" != "" ]
then
  waitPid $GENOME_WEIGHT_PID $WEIGHT_DIR/$GENOME_SAM_EDIT_DISTANCE_FILE "Compute genome SAM edit distance" "" $WEIGHT_DIR/genome.log
fi

if [ "$EXIT_STATUS" != "0" ]
then
  echo "Computation of read weight files failed."
  exit 1
  date
fi


################################################################################
##
##  Combine transcript, junction, and genome SAM files
##
##  Step 1: set the name of the combined file
##
################################################################################

if [ "$COMBINED_SAM_FILE_BASE" = "" ]
then
  COMBINED_SAM_FILE_BASE=`echo $TRANSCRIPT_SAM_FILE_BASE | sed -e "s/-transcript-/-combined-/"`
  if [ "$COMBINED_SAM_FILE_BASE" = "$TRANSCRIPT_SAM_FILE_BASE" ]
  then
    COMBINED_SAM_FILE_BASE="$TRANSCRIPT_SAM_FILE_BASE-combined"
  fi
fi

COMBINED_SAM_FILE="$COMBINED_SAM_FILE_BASE.sam.gz"
if [ "$RECOMPUTE" = "TRUE" ]
then
  rm -f $SAM_DIR/$COMBINED_SAM_FILE
fi


################################################################################
##
##  Step 2: Combine transcript and junction weight files
##
################################################################################

TRANSCRIPT_JUNCTION_SAM_FILE_BASE=`echo $COMBINED_SAM_FILE_BASE | sed -e "s/-combined-/-trans-junction-/"`
TRANSCRIPT_JUNCTION_SAM_FILE="$TRANSCRIPT_JUNCTION_SAM_FILE_BASE.$TRANSCRIPT_SAM_FILE_EXT"
TRANSCRIPT_JUNCTION_SAM_EDIT_DISTANCE_FILE="$TRANSCRIPT_JUNCTION_SAM_FILE_BASE-sam.wgt"

if [ "$RECOMPUTE" = "TRUE" ]
then
  rm -f $WEIGHT_DIR/$TRANSCRIPT_JUNCTION_SAM_EDIT_DISTANCE_FILE
fi

if [ ! -r $WEIGHT_DIR/$TRANSCRIPT_JUNCTION_SAM_EDIT_DISTANCE_FILE -a "$JUNCTION_FILE_MISSING" != "TRUE"  ]
then
  echo "Creating transcript+junction weight file $TRANSCRIPT_JUNCTION_SAM_EDIT_DISTANCE_FILE"
  date
  COMBINE_READ_WEIGHT_CMD="CombineReadWeightFiles -e -1 $WEIGHT_DIR/$TRANSCRIPT_SAM_EDIT_DISTANCE_FILE -2 $WEIGHT_DIR/$JUNCTION_SAM_EDIT_DISTANCE_FILE \
       -o $WEIGHT_DIR/$TRANSCRIPT_JUNCTION_SAM_EDIT_DISTANCE_FILE"
  echo "$COMBINE_READ_WEIGHT_CMD"
  $JAVA $COMBINE_READ_WEIGHT_CMD
fi


################################################################################
##
##  Creating the commands to combine the transcript, junction, and genome SAM files
##
################################################################################

ZCAT_TRANSCRIPT_SAM_FILE_CMD="zcat $TRANSCRIPT_SAM_PATH"

if [ "$JUNCTION_FILE_MISSING" = "TRUE" ]
then
  JUNCTION_COMBINE_SAM_FILE_CMD="cat"
  JUNCTION_COMBINE_SAM_FILE_CMD_FULL="cat"
else
  JUNCTION_COMBINE_SAM_FILE_CMD="CombineSamFiles -1 - -2 $JUNCTION_SAM_PATH \
    -w $WEIGHT_DIR/$TRANSCRIPT_SAM_EDIT_DISTANCE_FILE -W $WEIGHT_DIR/$JUNCTION_SAM_EDIT_DISTANCE_FILE -o -"
  JUNCTION_COMBINE_SAM_FILE_CMD_FULL="$JAVA $JUNCTION_COMBINE_SAM_FILE_CMD"
fi
      
if [ "$GENOME_FILE_MISSING" = "TRUE" ]
then
  GENOME_COMBINE_SAM_FILE_CMD="cat"
  GENOME_COMBINE_SAM_FILE_CMD_FULL="cat"
else
  GENOME_COMBINE_SAM_FILE_CMD="CombineSamFiles -1 - -2 $GENOME_SAM_PATH  \
    -w $WEIGHT_DIR/$TRANSCRIPT_JUNCTION_SAM_EDIT_DISTANCE_FILE -W $WEIGHT_DIR/$GENOME_SAM_EDIT_DISTANCE_FILE -o -"
  GENOME_COMBINE_SAM_FILE_CMD_FULL="$JAVA $GENOME_COMBINE_SAM_FILE_CMD"
fi


################################################################################
##
##  Combine transcript, junction, and genome SAM file
##
################################################################################

if [ "$RECOMPUTE" = "TRUE" ]
then
  rm -f $SAM_DIR/$COMBINED_SAM_FILE
fi

if [ ! -r $SAM_DIR/$COMBINED_SAM_FILE ]
then
  echo "Combining the transcript, junction, and genome SAM files"
  date

  echo "Command pipeline:"
  echo "  $ZCAT_TRANSCRIPT_SAM_FILE_CMD |"
  echo "  $JUNCTION_COMBINE_SAM_FILE_CMD |"
  echo "  $GENOME_COMBINE_SAM_FILE_CMD | gzip > $SAM_DIR/$COMBINED_SAM_FILE"

  ## Excute the commands
  $ZCAT_TRANSCRIPT_SAM_FILE_CMD | $JUNCTION_COMBINE_SAM_FILE_CMD_FULL | $GENOME_COMBINE_SAM_FILE_CMD_FULL | gzip > $SAM_DIR/$COMBINED_SAM_FILE
  if [ $? -ne 0 ]
  then
    echo "ERROR: Problem with generation of combined SAM file"
    OUTPUT_FILE=$SAM_DIR/$COMBINED_SAM_FILE
    if [ -f $OUTPUT_FILE ]
    then
      rm $OUTPUT_FILE
    fi
  fi

fi


################################################################################
##
##  Create the preliminary combined reference id file: We only need to include
##  the reference ids of the combined SAM file in the exon BED file.
##
################################################################################

if [ "$RECOMPUTE" = "TRUE" ]
then
  rm -f $SAM_DIR/$COMBINED_SAM_FILE_BASE-reference.ids
fi

COMPUTE_JUNCTION_REFERENCE_ID_FILE="FALSE"
if [ ! -r $SAM_DIR/$COMBINED_SAM_FILE_BASE-reference.ids ]
then
  echo "Extracting reference ids of file $SAM_DIR/$COMBINED_SAM_FILE_BASE.sam.gz"
  zcat $SAM_DIR/$COMBINED_SAM_FILE | cut -f 3 | sort -u -S 4G -T $TEMP_DIR > $SAM_DIR/$COMBINED_SAM_FILE_BASE-reference.ids
  if [ $? -ne 0 ]
  then
    echo "ERROR: Problem with zcat $SAM_DIR/$COMBINED_SAM_FILE | cut -f 3 | sort -u -S 4G > $SAM_DIR/$COMBINED_SAM_FILE_BASE-reference.ids ... exiting."
    exit 1
  fi

  COMPUTE_JUNCTION_REFERENCE_ID_FILE="TRUE"
fi


################################################################################
##
##  Filter the exon BED file
##
################################################################################

if [ "$INTERSECTION_BED_FILE" = "" ]
then
  INTERSECTION_BED_FILE="$COMBINED_SAM_FILE_BASE-intersection.bed.gz"
  echo "Using $INTERSECTION_BED_FILE as intersection BED file"
fi

if [ "$RECOMPUTE" = "TRUE" ]
then
  rm -f $BED_DIR/$INTERSECTION_BED_FILE
fi

COMBINED_EXON_FILTERED_BED_FILE="$BED_DIR/${GENE_MODEL_PREFIX}_${RADIUS}_combined_exons_filtered.bed"
if [ ! -r $BED_DIR/$INTERSECTION_BED_FILE ]
then
  echo "Filtering $COMBINED_EXON_BED_FILE"
  ## filterFile.py options: -I -> index of key field, -s -> separator string for key field (needed for transcripts ids with a ".<cover version>" extension)
  $UTIL_LIB_DIR/filterFile.py -I 0 -s "." -f $SAM_DIR/$COMBINED_SAM_FILE_BASE-reference.ids -i $COMBINED_EXON_BED_FILE -o $COMBINED_EXON_FILTERED_BED_FILE
  if [ $? -ne 0 ]
  then
    echo "ERROR: Problem with $UTIL_LIB_DIR/filterFile.py -I 0 -s \".\" -f $SAM_DIR/$COMBINED_SAM_FILE_BASE-reference.ids -i $COMBINED_EXON_BED_FILE \
      -o $COMBINED_EXON_FILTERED_BED_FILE ... exiting."
    exit 1
  fi
  date
fi


################################################################################
##
##  Create the final combined reference id file: Now we also consider all
##  mapped junction ids. The reason is that there may be exons that are not
##  neighbouring exons in the GTF file but become neighbouring exons via the
##  exon to transcript mapping of EQP (of course, this is a rare case). By
##  including all mapped junction ids it should be possible to catch the majority
##  of these cases.
##
################################################################################

if [ "$COMPUTE_JUNCTION_REFERENCE_ID_FILE" = "TRUE" -a "$JUNCTION_FILE_MISSING" != "TRUE"  ]
then
  echo "Extracting reference ids of file $JUNCTION_SAM_PATH"
  zcat $JUNCTION_SAM_PATH | cut -f 3 | sort -u -S 4G -T $TEMP_DIR > $SAM_DIR/$JUNCTION_SAM_FILE_BASE-reference.ids
  if [ $? -ne 0 ]
  then
    echo "ERROR: Problem with zcat $JUNCTION_SAM_PATH | cut -f 3 | sort -u -S 4G -T $TEMP_DIR > $SAM_DIR/$JUNCTION_SAM_FILE_BASE-reference.ids ... exiting."
    exit 1
  fi
  date

  echo "Merging $COMBINED_SAM_FILE_BASE-reference.ids and $JUNCTION_SAM_FILE_BASE-reference.ids"
  mv $SAM_DIR/$COMBINED_SAM_FILE_BASE-reference.ids $SAM_DIR/$COMBINED_SAM_FILE_BASE-reference.ids.bak
  cat $SAM_DIR/$COMBINED_SAM_FILE_BASE-reference.ids.bak $SAM_DIR/$JUNCTION_SAM_FILE_BASE-reference.ids | sort -u -S 4G -T $TEMP_DIR > $SAM_DIR/$COMBINED_SAM_FILE_BASE-reference.ids

  rm $SAM_DIR/$JUNCTION_SAM_FILE_BASE-reference.ids $SAM_DIR/$COMBINED_SAM_FILE_BASE-reference.ids.bak
  date
fi


################################################################################
##
##  Creating the command to convert the combined SAM file to a BED file
##
################################################################################

ZCAT_COMBINED_SAM_FILE_CMD="zcat $SAM_DIR/$COMBINED_SAM_FILE"
CONVERT_SAM_BED_CMD="ConvertSamBed $EDIT_DISTANCE_OPTION $CONVERT_SAM_BED_STRAND_SPECIFIC_OPTION -T $TRANSCRIPT_BED_FILE -s - -o -"


################################################################################
##
##  Convert the combined SAM file to a BED file and intersect with filtered
##  exon BED file
##
################################################################################

COMBINED_BED_FILE="$COMBINED_SAM_FILE_BASE.bed"
if [ ! -r $BED_DIR/$INTERSECTION_BED_FILE ]
then

  echo "Converting to BED format file and intersecting with the exon BED file"
  date

  echo "Convert SAM/BAM file $SAM_FILE_BASE to a BED file and intersect with $EXON_BED_FILE_BASE"
  echo "Command pipeline:"
  echo "  $ZCAT_COMBINED_SAM_FILE_CMD |"
  echo "  $CONVERT_SAM_BED_CMD |"
  echo "  tee $BED_DIR/$COMBINED_BED_FILE |"
  echo "  $BEDTOOLS_EXE intersect -wo $BED_STRAND_SPECIFIC_OPTION -a stdin -b $COMBINED_EXON_FILTERED_BED_FILE |"
  echo "  awk '($1 == $7) {print}' | gzip > $BED_DIR/$INTERSECTION_BED_FILE"

  
  $ZCAT_COMBINED_SAM_FILE_CMD | $JAVA $CONVERT_SAM_BED_CMD | tee $BED_DIR/$COMBINED_BED_FILE | \
  $BEDTOOLS_EXE intersect -wo $BED_STRAND_SPECIFIC_OPTION -a stdin -b $COMBINED_EXON_FILTERED_BED_FILE | \
  awk '($1 == $7) {print}' | gzip > $BED_DIR/$INTERSECTION_BED_FILE
  if [ $? -ne 0 ]
  then
    echo "ERROR: Problem with generation of intersect BED file"
    OUTPUT_FILE="$BED_DIR/$COMBINED_BED_FILE"
    if [ -f $OUTPUT_FILE ]
    then
      rm $OUTPUT_FILE
    fi
    
    OUTPUT_FILE="$BED_DIR/$INTERSECTION_BED_FILE"
    if [ -f $OUTPUT_FILE ]
    then
      rm $OUTPUT_FILE
    fi
    exit 1
  fi
  echo "Intersection BED file computed"
  date

  if [[ -z "$OUTPUT_PREFIX" ]]
  then
    gzip -f $BED_DIR/$COMBINED_BED_FILE &
    ## We do not need to wait for the completion of the following gzips since
    ## compute-intersection-bed-file.sh is only an intermediate script
    if [ "$WAIT" = "TRUE" ]
    then
      echo "Waiting for compression of file $BED_DIR/$COMBINED_BED_FILE to finish."
      wait
      date
    fi
  fi
elif [ "$COMPUTE_COUNTS" = "TRUE" ]
then
  echo "Computing the number of mapped reads."
  NUMBER_MAPPED_READS=`zcat $BED_DIR/$COMBINED_BED_FILE.gz | cut -f 4 | tr "-" "\t" | cut -f 1 | sort -u -S 4G -T $TEMP_DIR | wc -l`
  echo "NUMBER_MAPPED_READS=$NUMBER_MAPPED_READS"
  date
fi


################################################################################
##
##  Cleaning up
##
################################################################################

## Cleaning up
# rm -f $COMBINED_EXON_FILTERED_BED_FILE

if [ "$REMOVE_TEMP_DIR" = "TRUE" -a -d $TEMP_DIR ]
then
  rm -r $TEMP_DIR
fi

echo "Computation of BED intersection file complete."
date
