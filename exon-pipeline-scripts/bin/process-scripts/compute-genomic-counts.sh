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
PROG_DIR=`cd "$PROG_DIR" && pwd`
VERSION=1.0.0
PID="$$"

## Ensure that pipes report the exit status of the first failed job
set -o pipefail


################################################################################
##
## wait for sub process to finish
##
################################################################################

waitForPid ()
{
  WAIT_PID="$1"
  shift

  WAIT_PPID=`ps -p $WAIT_PID -o ppid | tail -1`
  if [ "$WAIT_PPID" = "$$" ]
  then
    wait $WAIT_PID
  fi
}

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
      echoVerbose "No waiting required: process $PID is no longer a child process of $$"
    fi
    PID="$CONFIRMED_PID"
  fi

  PID_EXIT_STATUS=0
  if [ "$PID" != "" ]
  then
    echoVerbose "$MESSAGE (pid $PID)"
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
      echoVerbose "WARNING: File $LOG_FILE not found."
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
    echoVerbose "Done"
  fi
}


################################################################################
##
## Output if VERBOSE is TRUE
##
################################################################################

echoVerbose ()
{
  if [ "$VERBOSE" = "TRUE" ]
  then
    echo $1
  fi
}

dateVerbose ()
{
  if [ "$VERBOSE" = "TRUE" ]
  then
    date
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
## Check tool
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
      elif [ "$TOOL_NAME" = "samtools" ]
      then
        ACTUAL_TOOL_VERSION=`$TOOL_EXE 2>&1 | fgrep "Version:" | sed -e "s/Version: \([0-9]*[.][0-9]*[.][0-9]*\).*/\1/"`
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

COMPUTE_EXON_COUNT="FALSE"
COMPUTE_GENE_COUNT="FALSE"
COMPUTE_JUNCTION_COUNT="FALSE"
COMPUTE_OLD_JUNCTION_COUNTS="FALSE"
EXON_OVERLAP=5
JUNCTION_OVERLAP_DEFAULT="8"
JUNCTION_OVERLAP_INPUT=""
READ_WEIGHT_THRESHOLD=0.01
COMPUTE_UNWEIGHTED="FALSE"
STRAND_SPECIFIC_OPTION=""
STRAND_SPECIFIC_DIRECTION_OPTION=""
OUTPUT_PREFIX=""
NON_ZERO_JUNCTION_COUNT_OPTION=""
UNAMBIGUOUS_OPTION=""
CONTAINMENT_OPTION=""
NONSPLICE_CONFORMING_OPTION=""
OVERLAP_MODE="FALSE"
READ_WEIGHT_FILE=""
SORT_BAM_FILE="TRUE"
RECOMPUTE="TRUE"
KEEP_TEMPORARY_FILES="FALSE"
QSUB="FALSE"
PRINT_HELP="FALSE"
while [ "$1" = "-h" -o "$1" = "-d" -o "$1" = "-a" -o "$1" = "--all" -o "$1" = "-g" -o "$1" = "-e" -o "$1" = "-j" -o "$1" = "-E" -o \
        "$1" = "-J" -o "$1" = "-W" -o "$1" = "-sr" -o "$1" = "--no-recompute" -o "$1" = "--unweighted" -o "$1" = "-unweighted" -o "$1" = "-s" -o \
	"$1" = "-N" -o "$1" = "-o" -o "$1" = "-nj" -o "$1" = "--non-zero-junction-counts-only" -o "$1" = "--unambig" -o "$1" = "-unambig" -o \
	"$1" = "-w" -o "$1" = "--nosort" -o "$1" = "-C" -o "$1" = "-O" -o "$1" = "--keep-tmp" -o "$1" = "--verbose" -o "$1" = "-v" -o \
	"$1" = "--qsub" -o "$1" = "-q" -o "$1" = "--old-junction-counts" ]
do

  if [ "$1" = "-h" ]
  then
    shift
    PRINT_HELP="TRUE"
  fi
  
  if [ "$1" = "-d" ]
  then
    shift
    PROJECT_DATA_DIR="$1"
    shift
  fi
  
  if [ "$1" = "-a" -o "$1" = "--all" ]
  then
    shift
    COMPUTE_GENE_COUNT="TRUE"
    COMPUTE_EXON_COUNT="TRUE"
    COMPUTE_JUNCTION_COUNT="TRUE"
  fi
  
  if [ "$1" = "-g" ]
  then
    shift
    COMPUTE_GENE_COUNT="TRUE"
  fi

  if [ "$1" = "-e" ]
  then
    shift
    COMPUTE_EXON_COUNT="TRUE"
  fi

  if [ "$1" = "-j" ]
  then
    shift
    COMPUTE_JUNCTION_COUNT="TRUE"
  fi

  if [ "$1" = "--old-junction-counts" ]
  then
    shift
    COMPUTE_OLD_JUNCTION_COUNTS="TRUE"
  fi

  if [ "$1" = "-E" ]
  then
    shift
    EXON_OVERLAP=$1
    shift
  fi
  
  if [ "$1" = "-J" ]
  then
    shift
    JUNCTION_OVERLAP_INPUT=$1
    shift
  fi

  if [ "$1" = "-W" ]
  then
    shift
    READ_WEIGHT_THRESHOLD=$1
    shift
  fi

  if [ "$1" = "--unweighted" -o "$1" = "-unweighted" ]
  then
    shift
    COMPUTE_UNWEIGHTED="TRUE"
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
  fi

  if [ "$1" = "-nj" -o "$1" = "--non-zero-junction-counts-only" ]
  then
    shift
    NON_ZERO_JUNCTION_COUNT_OPTION="-n"
  fi
  
  if [ "$1" = "--unambig" -o "$1" = "-unambig" ]
  then
    shift
    UNAMBIGUOUS_OPTION="-U"
  fi
  
  if [ "$1" = "-N" ]
  then
    shift
    NONSPLICE_CONFORMING_OPTION="-N"
  fi
  
  if [ "$1" = "-w" ]
  then
    shift
    READ_WEIGHT_FILE=$1
    shift
  fi

  if [ "$1" = "--nosort" ]
  then
    shift
    SORT_BAM_FILE="FALSE"
  fi

  if [ "$1" = "--no-recompute" ]
  then
    shift
    RECOMPUTE="FALSE"
  fi
  
  if [ "$1" = "-C" ]
  then
    shift
    CONTAINMENT_OPTION="-C"
  fi

  if [ "$1" = "-O" ]
  then
    shift
    OVERLAP_MODE="TRUE"
  fi

  if [ "$1" = "--keep-tmp" ]
  then
    shift
    KEEP_TEMPORARY_FILES="TRUE"
  fi

  if [ "$1" = "--verbose" -o "$1" = "-v" ]
  then
    shift
    VERBOSE="TRUE"
  fi
  
  if [ "$1" = "--qsub" -o "$1" = "-q" ]
  then
    shift
    QSUB="TRUE"
  fi
done

if [ "$RECOMPUTE" = "TRUE" ]
then
  BED_OPTION_LIST="$BED_OPTION_LIST -r"
fi

if [ "$COMPUTE_GENE_COUNT" = "FALSE" -a "$COMPUTE_EXON_COUNT" = "FALSE" -a "$COMPUTE_JUNCTION_COUNT" = "FALSE" ]
then
  COMPUTE_GENE_COUNT="TRUE"
  COMPUTE_EXON_COUNT="TRUE"
  COMPUTE_JUNCTION_COUNT="TRUE"
fi

PROJECT_DIR_NAME="project dir"
if [ "$PROG_NAME" = "eqp-quantify.sh" ]
then
  PROJECT_DIR_NAME="output dir"
fi

################################################################################
##
##  Print help
##
################################################################################

if [ "$PRINT_HELP" = "TRUE" -o "$1" = "" -o "$2" = "" ]
then
  echo "Usage: $PROG_NAME <options> <$PROJECT_DIR_NAME> <SAM/BAM file>"
  echo
  echo "where <options> is"
  echo "   [-d <setup dir>] [-g] [-e] [-j] [-E <exon overlap>] [-J <junction overlap>]"
  echo "   [-W <min read weight>] [-s <direction>] [-o <output prefix>]"
  echo "   [--nosort] [--unambig] [--unweighted] [-w <weight file>]"
  echo
  if [ "$PROG_NAME" != "eqp-quantify.sh" ]
  then
    echo "$PROJECT_DIR_NAME: directory where the relevant files are stored"
  else
    echo "$PROJECT_DIR_NAME: directory used for the count files"
  fi
  echo "BAM file: the file containing the alignments of the"
  echo "     reads against the genome with the aligner."
  echo " -d STRING: Use STRING as the directory that contains the auxilliary files"
  echo " -g: compute gene counts [computed by default, without -g,-e,-j] "
  echo " -e: compute exon counts [computed by default, without -g,-e,-j]"
  echo " -j: compute junction counts [computed by default, without -g,-e,-j]"
  echo " -E INT: Minimal overlap of a read with an exon [$EXON_OVERLAP]"
  echo " -J INT: Minimal overlap of a read with both exons on a junction [$JUNCTION_OVERLAP_DEFAULT]"
  echo " -W FLOAT: Minimal weight of a read; reads with a lower weight are"
  echo "           disregarded [$READ_WEIGHT_THRESHOLD]"
  echo " -s STRING: process reads as strand-specific in direction STRING (forward for"
  echo "     orientation fr or backward for orientation rf)"
  echo " -o STRING: A directory STRING is created in the current working directory"
  echo "     and all intermediate files are stored in a directory structure under"
  echo "     this directory; furthermore, the count files are generated in the current"
  echo "     working directory with the prefix STRING."
  echo " --nosort: the alignment file is already sorted by names; do not sort it."
  echo " --unambig: count only reads that can be assigned unambiguously to a single gene"
  echo "    or exon when creating the gene or exon counts"
  echo " --unweighted: do not use read weights in the generation of counts"
  echo " -w STRING: use STRING as weight file instead of computing it from the genomic"
  echo "    alignments"
  exit
fi


################################################################################
##
##  Read arguments
##
################################################################################

PROJECT_DIR=$1
shift

if [ ! -d $PROJECT_DIR ]
then
  mkdir -p $PROJECT_DIR
fi
PROJECT_DIR=`cd "$PROJECT_DIR" && pwd`

STRAND_SPECIFIC="FALSE"
GENE_MODEL_PREFIX=
SETUP_SCRIPT=$PROJECT_DIR/bin/setup.sh
if [[ -x $SETUP_SCRIPT ]]
then
  source $SETUP_SCRIPT
  GENE_MODEL_PREFIX=$FILE_BASE
elif [ -e $PROJECT_DIR/exon-pipeline-files ]
then
  echo "File $PROJECT_DIR/bin/setup.sh not found."
  echo "Please make sure that all files of the project instance are copied."
  exit 1
elif [ ! -e $PROJECT_DIR ]
then
  echo "Directory $PROJECT_DIR not found."
  exit 1
fi

SAM_FILE_PATH="$1"
shift

if [ ! -f $SAM_FILE_PATH ]
then
  echo "File $SAM_FILE_PATH not found ... exiting."
  exit 1
fi

## SAM_FILE_PATH=`readlink -f $SAM_FILE_PATH`

SAM_FILE=`basename $SAM_FILE_PATH`
SAM_FILE_BASE=`echo $SAM_FILE | sed -e 's/.gz$//' | sed -e 's/.[bs]am$//'` 
SAM_FILE_EXT=`echo $SAM_FILE | sed -e "s/$SAM_FILE_BASE[.]//"`


################################################################################
##
##  Set JUNCTION_OVERLAP, COUNT_STRAND_SPECIFIC_OPTION, and
##  CONVERT_SAM_BED_STRAND_SPECIFIC_OPTION
##
################################################################################

if [ "$JUNCTION_OVERLAP_INPUT" != "" ]
then
  JUNCTION_OVERLAP="$JUNCTION_OVERLAP_INPUT"
else
  if [ "$READ_LENGTH" != "" -a "$RADIUS" != "" ]
  then
    if [ $READ_LENGTH -gt $RADIUS ]
    then
      JUNCTION_OVERLAP=`echo $READ_LENGTH - $RADIUS | bc`
    else
      JUNCTION_OVERLAP="$JUNCTION_OVERLAP_DEFAULT"
    fi
  else
    JUNCTION_OVERLAP="$JUNCTION_OVERLAP_DEFAULT"
  fi
fi

if [ "$JUNCTION_OVERLAP" -lt 1 ]
then
  echo "Please specify a positive junction overlap (option -J). The radius $RADIUS is"
  echo "larger than the read length $READ_LENGTH ... exiting."
  exit 1
fi

COUNT_STRAND_SPECIFIC_OPTION=""
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
  COUNT_STRAND_SPECIFIC_OPTION="-s"

  ## Note that in the conversion to the BED file the read alignments are changed to the
  ## correct orientation; hence, we can use the intersectBed option for the correct
  ## orientation (-s). The option for the reverse orientation (-S) is not necessary.
  BED_STRAND_SPECIFIC_OPTION="-s"
fi


################################################################################
##
##  Set SAM_DIR which is used to specify the location of intermediate files
##  e.g. BED or weight files. If the input SAM file needs to be sorted, the
##  tempory SAM will be stored in the directory $SAM_DIR as well.
##
################################################################################

PROG_DIR_PARENT=`dirname $PROG_DIR`
if [ "$OUTPUT_PREFIX" = "" ]
then
  SAM_DIR=`dirname $SAM_FILE_PATH`

  ## Check if we are calling compute-genomic-counts.sh with the standard directory structure
  if [ "${SAM_DIR#$PROJECT_DIR/*/C[0-9][0-9][0-9]/sam-files}" != "" ]
  then
    OUTPUT_PREFIX="$PROJECT_DIR/$SAM_FILE_BASE"
    SAM_DIR="$OUTPUT_PREFIX/sam-files"
  fi
else
  SAM_DIR="$OUTPUT_PREFIX/sam-files"
fi

if [ "$OUTPUT_PREFIX" != "" ]
then
  echo "Using directory $OUTPUT_PREFIX for intermediate results"
  if [ -d "$OUTPUT_PREFIX" ]
  then
    rm -r $OUTPUT_PREFIX
  fi
  mkdir -p $OUTPUT_PREFIX
fi

if [ ! -d $SAM_DIR ]
then
  mkdir -p $SAM_DIR
fi


################################################################################
##
##  Set directory variables and check directories
##
################################################################################

BIN_DIR=$PROJECT_DIR/exon-pipeline-scripts
if [ ! -d $BIN_DIR -a ! -h $BIN_DIR ]
then
  if [ "$QSUB" = "FALSE" ]
  then
    BIN_DIR=$PROG_DIR
  else
    echo "Directory $BIN_DIR not found ... exiting."
    exit 1
  fi
fi

if [ ! -d $BIN_DIR ]
then
  echo "Directory $BIN_DIR not found ... exiting"
  exit 1
fi
echoVerbose "Using $BIN_DIR as scripts directory"


if [ "$PROJECT_DATA_DIR" = "" ]
then
  if [ -d $PROJECT_DIR/exon-pipeline-files ]
  then
    PROJECT_DATA_DIR=$PROJECT_DIR/exon-pipeline-files
  else
    PROJECT_DATA_DIR=$PROJECT_DIR
  fi
fi

if [ ! -d $PROJECT_DATA_DIR ]
then
  echo "Directory $PROJECT_DATA_DIR not found ... exiting"
  exit 1
fi
echoVerbose "Using directory $PROJECT_DATA_DIR as data directory"


export JAVA_DIR=$BIN_DIR/java
if [ ! -d $JAVA_DIR ]
then
  echo "Directory $JAVA_DIR not found ... exiting"
  exit 1
fi

export TOOLS_DIR=$BIN_DIR/tools
if [ -d $TOOLS_DIR ]
then
  PATH="${TOOLS_DIR}:$PATH"
fi
checkTool $TOOLS_DIR bedtools 2.24.0
BEDTOOLS_EXE=$TOOL_EXE
checkTool $TOOLS_DIR samtools 0.1.17
SAMTOOLS_EXE=$TOOL_EXE


PROJECT_BED_DIR=$PROJECT_DATA_DIR/bed-files
if [ ! -d $PROJECT_BED_DIR ]
then
  echo "Directory $PROJECT_BED_DIR not found ... exiting"
  exit 1
fi

PROJECT_GTF_DIR=$PROJECT_DATA_DIR/gtf-files
if [ ! -d $PROJECT_GTF_DIR ]
then
  echo "Directory $PROJECT_GTF_DIR not found ... exiting"
  exit 1
fi

if [ "$GENE_MODEL_PREFIX" = "" ]
then
  GENE_MODEL_PREFIX=`ls -1 $PROJECT_GTF_DIR | fgrep ".gtf" | tail -1 | sed -e 's/[.]gtf$//'`
fi

PROJECT_MAP_DIR=$PROJECT_DATA_DIR/map-files
if [ ! -d $PROJECT_MAP_DIR ]
then
  echo "Directory $PROJECT_MAP_DIR not found ... exiting"
  exit 1
fi

JAVA_CLASS_DIR="$JAVA_DIR/classes"
JAVA="java -oss8M -ss8M -ms6G -mx6G -cp ${JAVA_CLASS_DIR}"
echoVerbose "JAVA=$JAVA"


################################################################################
##
##  Set gene exon files and BED file
##
################################################################################

EXON_BED_FILE=$PROJECT_BED_DIR/${GENE_MODEL_PREFIX}_genome_exons.bed
if [ "$COMPUTE_GENE_COUNT" = "TRUE" -a "$CONTAINMENT_OPTION" = "-C" ]
then
  EXON_BED_FILE=$PROJECT_BED_DIR/${GENE_MODEL_PREFIX}_genome_exons_combined.bed
fi
if [ ! -f $EXON_BED_FILE ]
then
  echo "File $EXON_BED_FILE not found. Please create the file with the -setup option"
  echo "... exiting"
  exit 1
fi

if [ "$CONTAINMENT_OPTION" = "-C" ]
then
  EXON_GENE_MAP_FILE=$PROJECT_MAP_DIR/${GENE_MODEL_PREFIX}_exon_gene_combined.map
else
  EXON_GENE_MAP_FILE=$PROJECT_MAP_DIR/${GENE_MODEL_PREFIX}_exon_gene.map
fi

if [ "$COMPUTE_GENE_COUNT" = "TRUE" -a ! -f $EXON_GENE_MAP_FILE ]
then
  echo "File $EXON_GENE_MAP_FILE not found ... exiting"
  exit 1
fi

EXON_EXON_MAP_FILE=$PROJECT_MAP_DIR/${GENE_MODEL_PREFIX}_exon_exon.map
if [ "$COMPUTE_EXON_COUNT" = "TRUE" -a ! -f $EXON_EXON_MAP_FILE ]
then
  echo "File $EXON_EXON_MAP_FILE not found ... exiting"
  exit 1
fi

EXON_JUNCTION_MAP_FILE=$PROJECT_MAP_DIR/${GENE_MODEL_PREFIX}_exon_junction.map.gz
JUNCTION_EXON_MAP_FILE=$PROJECT_MAP_DIR/${GENE_MODEL_PREFIX}_junction_exon.map.gz
GTF_FILE=$PROJECT_GTF_DIR/${GENE_MODEL_PREFIX}.gtf
if [ "$COMPUTE_JUNCTION_COUNT" = "TRUE" ]
then
  if [ ! -f $EXON_JUNCTION_MAP_FILE ]
  then
    echo "File $EXON_JUNCTION_MAP_FILE not found ... exiting"
    exit 1
  fi
  if [ ! -f $GTF_FILE ]
  then
    echo "File $GTF_FILE not found ... exiting"
    exit 1
  fi
fi


################################################################################
##
##  Output meta information to log file
##
################################################################################

echoVerbose "QUANTIFIER=EQP $VERSION"
echoVerbose "READ_WEIGHT_THRESHOLD=$READ_WEIGHT_THRESHOLD"
echoVerbose "EXON_OVERLAP=$EXON_OVERLAP"
echoVerbose "JUNCTION_OVERLAP=$JUNCTION_OVERLAP"
echoVerbose "STRAND_SPECIFIC=$STRAND_SPECIFIC"
echoVerbose "STRAND_SPECIFIC_DIRECTION=$STRAND_SPECIFIC_DIRECTION"


################################################################################
##
##  Create directories for intermediate results
##
################################################################################

BED_DIR=`echo $SAM_DIR | sed -e 's;/sam-files$;/bed-files;'`
if [ ! -d $BED_DIR ]
then
  echoVerbose "Making directory $BED_DIR"
  mkdir $BED_DIR
  if [ $? -ne 0 ]
  then
    echo "Cannot create $BED_DIR in $PWD"
    exit 1
  fi
fi

WEIGHT_DIR=`echo $SAM_DIR | sed -e 's;/sam-files$;/weight-files;'`
if [ ! -d $WEIGHT_DIR ]
then
  echoVerbose "Making directory $WEIGHT_DIR"
  mkdir $WEIGHT_DIR
  if [ $? -ne 0 ]
  then
    echo "Cannot create $WEIGHT_DIR in $PWD"
    exit 1
  fi
fi

COUNT_DIR=`echo $SAM_DIR | sed -e 's;/sam-files$;/count-files;'`
if [ ! -d $COUNT_DIR ]
then
  echoVerbose "Making directory $COUNT_DIR"
  mkdir $COUNT_DIR
  if [ $? -ne 0 ]
  then
    echo "Cannot create $WEIGHT_DIR in $PWD"
    exit 1
  fi
fi

dateVerbose


################################################################################
##
##  Sort BAM files by read names
##
################################################################################

IS_BAM_FILE="FALSE"
BAM_FILE=`echo $SAM_FILE_PATH | sed -e 's/[.]gz$//' | sed -e 's/[.]sam$/.bam/'`
if [ "$BAM_FILE" = "$SAM_FILE_PATH" ]
then
  IS_BAM_FILE="TRUE"
fi

BAM_FILE_BASE=`basename $BAM_FILE`
if [ "$IS_BAM_FILE" = "FALSE" ]
then
  ## Create BAM files in $SAM_DIR (which may be in $OUTPUT_PREFIX)
  BAM_FILE=$SAM_DIR/$BAM_FILE_BASE

  if [ "$SAM_FILE_EXT" = "sam.gz" ]
  then
    CAT="zcat"
  elif [ "$SAM_FILE_EXT" = "sam" ]
  then
    CAT="cat"
  else
    echo "Unknown extension $SAM_FILE_EXT for SAM file $SAM_FILE_PATH ... exiting"
    exit 1
  fi

  echo "Converting file $SAM_FILE to a BAM file"
  $CAT $SAM_FILE_PATH | $SAMTOOLS_EXE view -Sb - > $BAM_FILE

  if [ $? -ne 0 ]
  then
    echo "ERROR: Problem with $CAT $SAM_FILE_PATH | $SAMTOOLS_EXE view -Sb - > $BAM_FILE ... exiting."
    exit 1
  fi
  echo "Conversion completed"

  TEMP_FILES="$BAM_FILE"
  IS_BAM_FILE="TRUE"
fi

if [ "$SORT_BAM_FILE" = "TRUE" ]
then
  
  SORTED_NAMES_BAM_FILE=`echo $SAM_DIR/$BAM_FILE_BASE | sed -e 's/[.]bam$/_sorted_names/'`
  if [ ! -f $SORTED_NAMES_BAM_FILE -o "$RECOMPUTE" = "TRUE" ]
  then
    echo "Sorting BAM file by name"
    $SAMTOOLS_EXE sort -n -m 5000000000 $BAM_FILE $SORTED_NAMES_BAM_FILE

    if [ $? -ne 0 ]
    then
      echo "ERROR: Problem with samtools sort -n -m 5000000000 $BAM_FILE $SORTED_NAMES_BAM_FILE ... exiting."
      exit 1
    fi

    echoVerbose "Done"
    dateVerbose
  fi

  SAM_FILE_PATH=$SORTED_NAMES_BAM_FILE.bam
  TEMP_FILES="$TEMP_FILES $SORTED_NAMES_BAM_FILE.bam"
else
  SAM_FILE_PATH=$BAM_FILE
fi

$SAMTOOLS_EXE view -H $SAM_FILE_PATH | grep "^@SQ" | cut -f 2 | sed -e "s/^SN://" | sed -e 's/^/#/' | sed -e 's/$/#/' > $SAM_DIR/chromsomes.txt
if [ $? -ne 0 ]
then
  echo "ERROR: Problem with $SAMTOOLS_EXE view -H $SAM_FILE_PATH ... exiting."
  exit 1
fi

NUM_CHROMOSOMES=`cat $SAM_DIR/chromsomes.txt | wc -l`
if [ "$NUM_CHROMOSOMES" == 0 ]
then
  echo "WARNING: The SAM input file does not contain a header - please make sure that"
  echo "the genome used for the generation of the SAM file and the GTF file used for"
  echo "the setup of EQP contain the same set of chromosome ids."
else
  NUM_CHROMOSOMES_IN_MISSED=`fgrep -v -f $PROJECT_GTF_DIR/${GENE_MODEL_PREFIX}-chromosomes.txt $SAM_DIR/chromsomes.txt | wc -l`
  if [ "$NUM_CHROMOSOMES_IN_MISSED" = "$NUM_CHROMOSOMES" ]
  then
    echo "ERROR: The chromosomes used for the alignment of the SAM file and the ones"
    echo "used in the setup of EQP have different identifiers. Please setup EQP with"
    echo "a GTF files that uses the same chromosome identifiers as used in the genome"
    echo "Fasta file used for alignment ... aborting."
    exit 1
  elif [ "$NUM_CHROMOSOMES_IN_MISSED" != "0" ]
  then
    echo "WARNING: The chromosomes used for the alignment of the SAM file and the ones"
    echo "used in the setup of EQP have different identifiers. Alignments against the"
    echo "following chromosomes will not be considered:"
    fgrep -v -f $PROJECT_GTF_DIR/${GENE_MODEL_PREFIX}-chromosomes.txt $SAM_DIR/chromsomes.txt | sed -e 's/^#//'  | sed -e 's/#$//'
  fi
fi

################################################################################
##
##  Compute weights from alignment file
##
################################################################################

if [[ -z "$READ_WEIGHT_FILE" ]]
then
  WEIGHT_FILE=$WEIGHT_DIR/${SAM_FILE_BASE}.wgt
else
  WEIGHT_FILE="$READ_WEIGHT_FILE"
fi

if [ "$RECOMPUTE" = "TRUE" -o ! -f $WEIGHT_FILE ]
then
  echo "Computing read weights" > $WEIGHT_DIR/${SAM_FILE_BASE}.log
  COMPUTE_READ_WEIGHT_JAVA_CMD="ComputeReadWeightsSam -o $WEIGHT_FILE"
  if [ "$IS_BAM_FILE" = "TRUE" ]
  then
    echoVerbose "Command line call: $SAMTOOLS_EXE view $SAM_FILE_PATH | $JAVA $COMPUTE_READ_WEIGHT_JAVA_CMD" >> $WEIGHT_DIR/${SAM_FILE_BASE}.log
    $SAMTOOLS_EXE view $SAM_FILE_PATH | $JAVA $COMPUTE_READ_WEIGHT_JAVA_CMD >> $WEIGHT_DIR/${SAM_FILE_BASE}.log 2>&1 &
    COMPUTE_READ_WEIGHT_PID=$!      
  else
    echoVerbose "Command line call: $JAVA $COMPUTE_READ_WEIGHT_JAVA_CMD -s $SAM_FILE_PATH" >> $WEIGHT_DIR/${SAM_FILE_BASE}.log
    $JAVA $COMPUTE_READ_WEIGHT_JAVA_CMD -s $SAM_FILE_PATH >> $WEIGHT_DIR/${SAM_FILE_BASE}.log 2>&1 &
    COMPUTE_READ_WEIGHT_PID=$!  
  fi
  dateVerbose
fi

WEIGHT_OPTION="-w $WEIGHT_FILE -W $READ_WEIGHT_THRESHOLD"
if [ "$COMPUTE_UNWEIGHTED" = "TRUE" ]
then
  WEIGHT_OPTION="-u"
fi


################################################################################
##
##  Intersect bed files
##
################################################################################

INTERSECTION_BED_FILE_GZIP="$BED_DIR/${SAM_FILE_BASE}-intersection.bed.gz"
EXON_BED_FILE_BASE=`basename $EXON_BED_FILE`

if [ "$RECOMPUTE" = "TRUE" -o ! -f $INTERSECTION_BED_FILE_GZIP ]
then
  echo "Converting SAM/BAM file:"
  echo $SAM_FILE_BASE
  echo "to a BED file and intersecting it with file:"
  echo $EXON_BED_FILE_BASE
  if [ "$IS_BAM_FILE" = "TRUE" ]
  then
    echoVerbose "Command:"
    echoVerbose "$SAMTOOLS_EXE view $SAM_FILE_PATH | $JAVA ConvertSamBed | \
       $BEDTOOLS_EXE intersect -wo $BED_STRAND_SPECIFIC_OPTION -a stdin -b $EXON_BED_FILE |  gzip > $INTERSECTION_BED_FILE_GZIP"

    $SAMTOOLS_EXE view $SAM_FILE_PATH | $JAVA ConvertSamBed $CONVERT_SAM_BED_STRAND_SPECIFIC_OPTION | \
       $BEDTOOLS_EXE intersect -wo $BED_STRAND_SPECIFIC_OPTION -a stdin -b $EXON_BED_FILE |  gzip > $INTERSECTION_BED_FILE_GZIP
  
    if [ $? -ne 0 ]
    then
      echo "ERROR: Problem with samtools view $SAM_FILE_PATH | $JAVA ConvertSamBed | \
       $BEDTOOLS_EXE intersect -wo $BED_STRAND_SPECIFIC_OPTION -a stdin -b $EXON_BED_FILE |  gzip > $INTERSECTION_BED_FILE_GZIP ... exiting."
      rm -f $INTERSECTION_BED_FILE_GZIP
      exit 1
    fi
  else
    echoVerbose "Command:"
    echoVerbose "ConvertSamBed -s $SAM_FILE_PATH | \ "
    echoVerbose "  $BEDTOOLS_EXE intersect -wo $BED_STRAND_SPECIFIC_OPTION -a stdin -b $EXON_BED_FILE | \ "
    echoVerbose "  gzip > $INTERSECTION_BED_FILE_GZIP"
    $JAVA ConvertSamBed $CONVERT_SAM_BED_STRAND_SPECIFIC_OPTION -s $SAM_FILE_PATH | \
       $BEDTOOLS_EXE intersect -wo $BED_STRAND_SPECIFIC_OPTION -a stdin -b $EXON_BED_FILE | gzip > $INTERSECTION_BED_FILE_GZIP
    if [ $? -ne 0 ]
    then
      echo "ERROR: Problem with ConvertSamBed -s $SAM_FILE_PATH | \
       $BEDTOOLS_EXE intersect -wo $BED_STRAND_SPECIFIC_OPTION -a stdin -b $EXON_BED_FILE | gzip > $INTERSECTION_BED_FILE_GZIP ... exiting."
      rm -f $INTERSECTION_BED_FILE_GZIP
      exit 1
    fi
  fi
fi
dateVerbose


################################################################################
##
##  Wait for the computation of read weights to complete
##
################################################################################

if [ "$COMPUTE_READ_WEIGHT_PID" != "" ]
then
  waitPid $COMPUTE_READ_WEIGHT_PID $WEIGHT_FILE "$COMPUTE_READ_WEIGHT_JAVA_CMD" "Waiting for completion of read weights" $WEIGHT_DIR/${SAM_FILE_BASE}.log
fi

if [ ! -f $WEIGHT_FILE ]
then
  echo "File $WEIGHT_FILE not found ... exiting."
  exit 1
else
  echoVerbose "Using file $WEIGHT_FILE as weight file."
fi


################################################################################
##
##  Compute the number of expressed reads
##
################################################################################

# zcat $INTERSECTION_BED_FILE_GZIP | cut -f 4 | tr "-" "\t" | cut -f 1 |  sort -u | wc -l > $COUNT_DIR/num-expressed-reads.txt &
# NUMBER_EXPRESSED_READS_PID=$!


################################################################################
##
##  Compute gene counts
##
################################################################################

if [ "$OUTPUT_PREFIX" = "" ]
then
  GENE_COUNT_FILE="$COUNT_DIR/${SAM_FILE_BASE}-gene.cnt"
else
  GENE_COUNT_FILE="$OUTPUT_PREFIX-gene.cnt"
fi


if [ "$COMPUTE_GENE_COUNT" = "TRUE" ]
then
  GENE_COUNT_OPTION="-g"
  if [ "$OVERLAP_MODE" = "TRUE" ]
  then
    GENE_COUNT_OPTION="-e -O 1"
  fi

  echo "Computing gene counts"
  echoVerbose "Read weight threshold: $READ_WEIGHT_THRESHOLD"
  echoVerbose "Weight file: $WEIGHT_FILE"
  GENE_COUNT_JAVA_CMD="ComputeCounts $GENE_COUNT_OPTION $WEIGHT_OPTION $COUNT_STRAND_SPECIFIC_OPTION $UNAMBIGUOUS_OPTION \
         $CONTAINMENT_OPTION $NONSPLICE_CONFORMING_OPTION -m $EXON_GENE_MAP_FILE -b - -o $GENE_COUNT_FILE"
  echoVerbose "Java cmd: $GENE_COUNT_JAVA_CMD"
  zcat $INTERSECTION_BED_FILE_GZIP | $JAVA $GENE_COUNT_JAVA_CMD &
  GENE_COUNT_PID=$!
  
fi


################################################################################
##
##  Compute exon counts
##
################################################################################

if [ "$OUTPUT_PREFIX" = "" ]
then
  EXON_COUNT_FILE="$COUNT_DIR/${SAM_FILE_BASE}-exon.cnt"
else
  EXON_COUNT_FILE="$OUTPUT_PREFIX-exon.cnt"
fi

if [ "$COMPUTE_EXON_COUNT" = "TRUE" ]
then

  echo "Computing exon counts with minimum overlap $EXON_OVERLAP" > $COUNT_DIR/genomic-exon-count.log
  EXON_COUNT_JAVA_CMD="ComputeCounts -e $WEIGHT_OPTION $COUNT_STRAND_SPECIFIC_OPTION $NONSPLICE_CONFORMING_OPTION $UNAMBIGUOUS_OPTION \
      -O $EXON_OVERLAP -m $EXON_EXON_MAP_FILE -b - -o $EXON_COUNT_FILE"
  echoVerbose "$EXON_COUNT_JAVA_CMD" >> $COUNT_DIR/genomic-exon-count.log
  zcat $INTERSECTION_BED_FILE_GZIP | $JAVA $EXON_COUNT_JAVA_CMD >> $COUNT_DIR/genomic-exon-count.log 2>&1 &
  EXON_COUNT_PID=$!
  
fi


################################################################################
##
##  Wait for the computation of the number of expressed reads to finish
##
################################################################################

if [ "$NUMBER_EXPRESSED_READS_PID" != "" ]
then
  waitForPid $NUMBER_EXPRESSED_READS_PID
  if [ $? -ne 0 ]
  then
    echo "WARNING: The number of expressed reads could not be computed."
    NUMBER_EXPRESSED_READS=-1
  else
    NUMBER_EXPRESSED_READS=`cat $COUNT_DIR/num-expressed-reads.txt | head -1`
  fi

  echoVerbose "NUMBER_EXPRESSED_READS=$NUMBER_EXPRESSED_READS"

  if [ -f $COUNT_DIR/num-expressed-reads.txt ]
  then
    rm $COUNT_DIR/num-expressed-reads.txt
  fi
  
  dateVerbose
fi


################################################################################
##
##  Wait for the computation of the gene counts to finish
##
################################################################################

EXIT_STATUS=0
if [ "$GENE_COUNT_PID" != "" ]
then
  waitPid $GENE_COUNT_PID $GENE_COUNT_FILE "$GENE_COUNT_JAVA_CMD" "Waiting for computation of gene counts to finish"
fi


################################################################################
##
##  Wait for the computation of the exon counts to finish
##
################################################################################

if [ "$EXON_COUNT_PID" != "" ]
then
  waitPid $EXON_COUNT_PID $EXON_COUNT_FILE "$EXON_COUNT_JAVA_CMD" "Waiting for computation of exon counts to finish" $COUNT_DIR/genomic-exon-count.log
fi


################################################################################
##
##  Abort if the computation of counts failed
##
################################################################################

if [ "$EXIT_STATUS" != 0 ]
then
  exit 1
fi


################################################################################
##
##  Compute junction counts using EXON_JUNCTION_MAP_FILE from the BED file
##
################################################################################

if [ "$OUTPUT_PREFIX" = "" ]
then
  JUNCTION_COUNT_FILE="$COUNT_DIR/${SAM_FILE_BASE}-junction.cnt"
else
  JUNCTION_COUNT_FILE="$OUTPUT_PREFIX-junction.cnt"
fi

if [ "$COMPUTE_JUNCTION_COUNT" = "TRUE" ]
then

  dateVerbose
  echoVerbose "Extracting the exon-exon pairs which border the gaps (introns) of spliced reads"
  FILTERED_EXON_EXON_FILE="$SAM_DIR/${SAM_FILE_BASE}-exon-exon.jnc"
  EXTRACT_EXON_EXON_CMD="ExtractSplicedExonExonIds -g $GTF_FILE -s - -o $FILTERED_EXON_EXON_FILE"

  if [ "$IS_BAM_FILE" = "TRUE" ]
  then
    echoVerbose "$SAMTOOLS_EXE view $SAM_FILE_PATH | $EXTRACT_EXON_EXON_CMD"
    $SAMTOOLS_EXE view $SAM_FILE_PATH | $JAVA $EXTRACT_EXON_EXON_CMD
    if [ $? -ne 0 ]
    then
      echo "Command $SAMTOOLS_EXE view $SAM_FILE_PATH | $EXTRACT_EXON_EXON_CMD failed ... exiting."
      exit 1
    fi
  else
    echoVerbose "zcat $SAM_FILE_PATH | $EXTRACT_EXON_EXON_CMD"
    zcat $SAM_FILE_PATH | $JAVA $EXTRACT_EXON_EXON_CMD
    if [ $? -ne 0 ]
    then
      echo "Command zcat $SAM_FILE_PATH | $EXTRACT_EXON_EXON_CMD failed ... exiting."
      exit 1
    fi
  fi

  dateVerbose
  echoVerbose "Filter the junction exon map file"
  FILTERED_JUNCTION_EXON_MAP_FILE="$SAM_DIR/${SAM_FILE_BASE}_junction_exon.map"
  $BIN_DIR/bin/util-lib/filterFile.py -f $FILTERED_EXON_EXON_FILE -F "0:1" -i $JUNCTION_EXON_MAP_FILE -I "1:2" -o $FILTERED_JUNCTION_EXON_MAP_FILE
  if [ $? -ne 0 ]
  then
    echo "Command $BIN_DIR/bin/util-lib/filterFile.py -f $FILTERED_EXON_EXON_FILE -F \"0:1\" -i $JUNCTION_EXON_MAP_FILE \
      -I \"1:2\" -o $FILTERED_JUNCTION_EXON_MAP_FILE failed ... exiting."
    exit 1
  fi

  dateVerbose
  echoVerbose "Filter the exon junction map file"
  FILTERED_EXON_JUNCTION_MAP_FILE="$SAM_DIR/${SAM_FILE_BASE}_exon_junction.map"
  $BIN_DIR/bin/util-lib/filterFile.py -f $FILTERED_JUNCTION_EXON_MAP_FILE -i $EXON_JUNCTION_MAP_FILE -I "1" -o $FILTERED_EXON_JUNCTION_MAP_FILE
  if [ $? -ne 0 ]
  then
    echo "Command $BIN_DIR/bin/util-lib/filterFile.py -f $FILTERED_JUNCTION_EXON_MAP_FILE -i $EXON_JUNCTION_MAP_FILE -I \"1\" \
       -o $FILTERED_EXON_JUNCTION_MAP_FILE failed ... exiting."
    exit 1
  fi

  dateVerbose
  ## Scale up the memory requirements of Java
  JAVA="java -oss8M -ss8M -ms8G -mx8G -cp ${JAVA_CLASS_DIR}"
  echoVerbose "JAVA=$JAVA"

  echo "Computing junction counts with minimum overlap $JUNCTION_OVERLAP"
  JUNCTION_COUNT_JAVA_CMD="ComputeCounts -j $NON_ZERO_JUNCTION_COUNT_OPTION $WEIGHT_OPTION -m $FILTERED_EXON_JUNCTION_MAP_FILE \
     -O $JUNCTION_OVERLAP $NONSPLICE_CONFORMING_OPTION -b - -o $JUNCTION_COUNT_FILE  $COUNT_STRAND_SPECIFIC_OPTION"
  
  echoVerbose "$JUNCTION_COUNT_JAVA_CMD"
  zcat $INTERSECTION_BED_FILE_GZIP | $JAVA $JUNCTION_COUNT_JAVA_CMD
  if [ $? -ne 0 ]
  then
    echo "Command $INTERSECTION_BED_FILE_GZIP | $JAVA $JUNCTION_COUNT_JAVA_CMD failed ... exiting."
    exit 1
  fi

  if [ "$COMPUTE_OLD_JUNCTION_COUNTS" = "TRUE" ]
  then
    ## Scale up the memory requirements of Java once more
    JAVA="java -oss8M -ss8M -ms16G -mx16G -cp ${JAVA_CLASS_DIR}"
    echoVerbose "JAVA=$JAVA"
    JUNCTION_COUNT_JAVA_CMD="ComputeCounts -j $NON_ZERO_JUNCTION_COUNT_OPTION $WEIGHT_OPTION -m $EXON_JUNCTION_MAP_FILE \
       -O $JUNCTION_OVERLAP $NONSPLICE_CONFORMING_OPTION -b - -o $JUNCTION_COUNT_FILE.old  $COUNT_STRAND_SPECIFIC_OPTION"
       
    echoVerbose "$JUNCTION_COUNT_JAVA_CMD"
    zcat $INTERSECTION_BED_FILE_GZIP | $JAVA $JUNCTION_COUNT_JAVA_CMD
    if [ $? -ne 0 ]
    then
      echo "Command $INTERSECTION_BED_FILE_GZIP | $JAVA $JUNCTION_COUNT_JAVA_CMD failed ... exiting."
      exit 1
    fi
  fi

fi


################################################################################
##
##  Clean-up and exit
##
################################################################################

if [ "$EXIT_STATUS" != "0" ]
then
  echo "Compute genomic counts failed."
  exit $EXIT_STATUS
fi

if [ "$KEEP_TEMPORARY_FILES" = "FALSE" ]
then
  for FILE in $TEMP_FILES
  do
    echoVerbose "Removing file $FILE"
    rm -r $FILE
  done
  
  
  if [ "$OUTPUT_PREFIX" != "" ]
  then
    echoVerbose "Removing intermediate results in directory $OUTPUT_PREFIX"
    rm -r $OUTPUT_PREFIX
  fi
fi

dateVerbose
echo "Compute genomic counts successfully completed."
exit $EXIT_STATUS
