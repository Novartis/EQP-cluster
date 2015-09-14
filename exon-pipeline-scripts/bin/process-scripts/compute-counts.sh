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

## echo $JAVA, echo rm


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
##  Read options
##
################################################################################

COMPUTE_EXON_COUNT="FALSE"
COMPUTE_GENE_COUNT="FALSE"
COMPUTE_JUNCTION_COUNT="FALSE"
EXON_OVERLAP=5
JUNCTION_OVERLAP_DEFAULT="8"
JUNCTION_OVERLAP_INPUT=""
READ_WEIGHT_THRESHOLD=0.01
PAIRED_END_OPTION=""
COMPUTE_UNWEIGHTED="FALSE"
STRAND_SPECIFIC_OPTION=""
STRAND_SPECIFIC_DIRECTION_OPTION=""
OUTPUT_PREFIX=""
BED_OPTION_LIST=""
NON_ZERO_JUNCTION_COUNT_OPTION=""
UNAMBIGUOUS_OPTION=""
CONTAINMENT_OPTION=""
NONSPLICE_CONFORMING_OPTION=""
OVERLAP_MODE="FALSE"
RECOMPUTE="FALSE"
PRINT_HELP="FALSE"
MAX_MEMORY="12G"
while [ "$1" = "-h" -o "$1" = "-a" -o "$1" = "-g" -o "$1" = "-e" -o "$1" = "-j" -o "$1" = "-js" -o "$1" = "-E" -o \
        "$1" = "-J" -o "$1" = "-W" -o "$1" = "-sr" -o "$1" = "--unweighted" -o "$1" = "-s" -o "$1" = "-N" -o \
	"$1" = "-o" -o "$1" = "-nj" -o "$1" = "--unambig" -o "$1" = "-r" -o "$1" = "-C" -o "$1" = "-O" -o "$1" = "-M" ]
do

  if [ "$1" = "-h" ]
  then
    shift
    PRINT_HELP="TRUE"
  fi
  
  if [ "$1" = "-a" ]
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

  if [ "$1" = "-js" ]
  then
    shift
    COMPUTE_JUNCTION_SAM_COUNT="TRUE"
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

  if [ "$1" = "-sr" ]
  then
    shift
    PAIRED_END_OPTION="FALSE"
  fi

  if [ "$1" = "--unweighted" ]
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
    BED_OPTION_LIST="$BED_OPTION_LIST -o $1 -w"
    OUTPUT_PREFIX_OPTION="-o $1"
    shift
  fi

  if [ "$1" = "-nj" ]
  then
    shift
    NON_ZERO_JUNCTION_COUNT_OPTION="-n"
  fi

  if [ "$1" = "--unambig" ]
  then
    shift
    UNAMBIGUOUS_OPTION="-U"
  fi
  
  if [ "$1" = "-r" ]
  then
    shift
    RECOMPUTE="TRUE"
    BED_OPTION_LIST="$BED_OPTION_LIST -r"
  fi

  if [ "$1" = "-C" ]
  then
    shift
    CONTAINMENT_OPTION="-C"
  fi
  
  if [ "$1" = "-N" ]
  then
    shift
    NONSPLICE_CONFORMING_OPTION="-N"
  fi

  if [ "$1" = "-O" ]
  then
    shift
    OVERLAP_MODE="TRUE"
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

if [ "$PRINT_HELP" = "TRUE" -o "$4" = "" ]
then
  echo "Usage: $PROG_NAME <options> <project dir> <transcript SAM file>"
  echo "        <junction SAM file> <genome SAM-file>"
  echo
  echo "where <options> is"
  echo "   [-a] [-g] [-e] [-j] [-E <min exon overlap>] [-J <min junction overlap]"
  echo "   [-W <min read weight>] [-sr] [-u] [-s <direction>] [-o <output prefix>]"
  echo "   [-nj] [-r] [-unambig]"
  echo
  echo " project dir: directory where the relevant files are stored"
  echo " transcript SAM file: the file containing the alignments of the"
  echo "     reads against the transcripts."
  echo " junction SAM file: the file containing the alignments of the"
  echo "     reads against the junctions."
  echo " genome SAM file: the file containing the alignments of the reads against"
  echo "    the genome."
  echo " -a: short hand for the combination of -g, -e, and -j"
  echo " -g: compute gene counts"
  echo " -e: compute exon counts"
  echo " -j: compute junction counts (the junction SAM file needs to be specified"
  echo "     for this option to have an effect)"
  echo " -E INT: Minimal overlap of a read with an exon [$EXON_OVERLAP]"
  echo " -J INT: Minimal overlap of a read with a junction [<read length> - <radius>]"
  echo " -W FLOAT: Minimal weight of a read; reads with a lower weight are"
  echo "           disregarded [$READ_WEIGHT_THRESHOLD]"
  echo " -sr: the input files are single read [the value from setup.sh]"
  echo " -u: do not use read weights and compute unweighted read counts only"
  echo "     (the read weight file is not computed and the minimal weight of a read"
  echo "     is set to -1)"
  echo " -s STRING: direction of how to process the reads as strand-specific: forward"
  echo "     or backward [the value of setup.sh]"
  echo " -o STRING: A directory STRING is created in the current working directory"
  echo "     and all intermediate files are stored in a directory structure under"
  echo "     this directory; furthermore, the count files are generated in the current"
  echo "     working directory with the prefix STRING."
  echo " -nj: output only non-zero junction counts"
  echo " -r: force the recomputation of SAM and BED files even if they exist"
  echo " -unambig: count only reads that can be assigned to a single gene"
  echo "    when creating the gene counts"
  exit
fi


################################################################################
##
##  Read arguments
##
################################################################################

PROJECT_DIR=$1
shift

BED_ARGUMENT_LIST="$PROJECT_DIR"

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

if [ "$PAIRED_END_OPTION" != "" ]
then
  PAIRED_END="$PAIRED_END_OPTION"
fi

COUNT_STRAND_SPECIFIC_OPTION=""
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
  COUNT_STRAND_SPECIFIC_OPTION="-s"
  BED_OPTION_LIST="$BED_OPTION_LIST -s $STRAND_SPECIFIC_DIRECTION"
fi


if [ "$JUNCTION_OVERLAP_INPUT" != "" ]
then
  JUNCTION_OVERLAP="$JUNCTION_OVERLAP_INPUT"
else
  if [ $READ_LENGTH -gt $RADIUS ]
  then
    JUNCTION_OVERLAP=`echo $READ_LENGTH - $RADIUS | bc`
  else
    JUNCTION_OVERLAP="$JUNCTION_OVERLAP_DEFAULT"
  fi
fi

if [ "$JUNCTION_OVERLAP" -lt 1 ]
then
  echo "Please specify a positive junction overlap (option -J). The radius $RADIUS is"
  echo "no less than the read length $READ_LENGTH ... exiting."
  exit 1
fi

GENE_MODEL_PREFIX=$FILE_BASE


TRANSCRIPT_SAM_PATH=$1
shift
if [ $RECOMPUTE = "TRUE" -a ! -f $TRANSCRIPT_SAM_PATH ]
then
  echo "Transcript SAM file $TRANSCRIPT_SAM_PATH not found ... exiting"
  exit 1
fi
BED_ARGUMENT_LIST="$BED_ARGUMENT_LIST $TRANSCRIPT_SAM_PATH"

######################### Transcript file ######################################
TRANSCRIPT_SAM_FILE=`basename $TRANSCRIPT_SAM_PATH`
TRANSCRIPT_SAM_FILE_BASE=`echo $TRANSCRIPT_SAM_FILE | sed -e 's/.gz$//' | sed -e 's/.[bs]am$//'` 
TRANSCRIPT_SAM_FILE_EXT=`echo $TRANSCRIPT_SAM_FILE | sed -e "s/$TRANSCRIPT_SAM_FILE_BASE[.]//"` 


########################### Junction file ######################################
JUNCTION_SAM_PATH=$1
shift
if [ $RECOMPUTE = "TRUE" -a ! -f "$JUNCTION_SAM_PATH" ]
then
  echo "WARNING: Junction SAM file not found ... no junction alignments will be used"
  export JUNCTION_FILE_MISSING="TRUE"
fi
BED_ARGUMENT_LIST="$BED_ARGUMENT_LIST $JUNCTION_SAM_PATH"

if [ "$JUNCTION_FILE_MISSING" != "TRUE" ]
then
  JUNCTION_SAM_FILE=`basename $JUNCTION_SAM_PATH`
  JUNCTION_SAM_FILE_BASE=`echo $JUNCTION_SAM_FILE | sed -e 's/.gz$//' | sed -e 's/.[bs]am$//'` 
  JUNCTION_SAM_FILE_EXT=`echo $JUNCTION_SAM_FILE | sed -e "s/$JUNCTION_SAM_FILE_BASE[.]//"`
fi


########################### Genome file ######################################
GENOME_SAM_PATH=$1
shift
if [ $RECOMPUTE = "TRUE" -a ! -f "$GENOME_SAM_PATH" ]
then
  echo "WARNING: Genome SAM file not found ... no genome alignments will be used"
  export GENOME_FILE_MISSING="TRUE"
fi
BED_ARGUMENT_LIST="$BED_ARGUMENT_LIST $GENOME_SAM_PATH"

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
##  Check aligner
##
################################################################################

if [ "$PAIRED_END" = "TRUE" ]
then
  SAM_FILE_BASE=`echo $TRANSCRIPT_SAM_FILE_BASE | sed -e "s/-pe//" | sed -e "s/-mixed//"`
  if [ "$SAM_FILE_BASE" = "$TRANSCRIPT_SAM_FILE_BASE" ]
  then
    echo "$TRANSCRIPT_SAM_FILE_BASE does not conform to naming convention (-pe|-mixed) ... exiting."
    exit 1
  fi
else
  SAM_FILE_BASE=`echo $TRANSCRIPT_SAM_FILE_BASE | sed -e "s/-sr//"`
  if [ "$SAM_FILE_BASE" = "$TRANSCRIPT_SAM_FILE_BASE" ]
  then
    echo "$TRANSCRIPT_SAM_FILE_BASE does not conform to naming convention (-sr) ... exiting."
    exit 1
  fi
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

HALF_MAX_MEMORY=`echo $MAX_MEMORY | sed -e 's;G$; / 2;' | bc | sed -e 's/$/G/'`

JAVA_CLASS_DIR="$JAVA_DIR/classes:$JAVA_DIR"
JAVA="java -oss8M -ss8M -ms$HALF_MAX_MEMORY -mx$HALF_MAX_MEMORY -cp ${JAVA_CLASS_DIR}:${CLASSPATH}"
echo "JAVA=$JAVA"


################################################################################
##
##  Set BED files and map files
##
################################################################################

if [ "$COMPUTE_GENE_COUNT" = "TRUE" ]
then
  EXON_GENE_MAP_FILE=$PROJECT_MAP_DIR/${GENE_MODEL_PREFIX}_exon_gene.map
  if [ $COMPUTE_GENE_COUNT = "TRUE" -a ! -f $EXON_GENE_MAP_FILE ]
  then
    echo "$EXON_GENE_MAP_FILE not found ... exiting"
    exit 1
  fi
fi

if [ "$COMPUTE_EXON_COUNT" = "TRUE" ]
then
  EXON_EXON_MAP_FILE=$PROJECT_MAP_DIR/${GENE_MODEL_PREFIX}_exon_exon.map
  if [ ! -f $EXON_EXON_MAP_FILE ]
  then
    echo "$EXON_EXON_MAP_FILE not found ... exiting"
    exit 1
  fi
fi

EXON_JUNCTION_MAP_FILE=$PROJECT_MAP_DIR/${GENE_MODEL_PREFIX}_exon_junction.map.gz
FILTERED_EXON_JUNCTION_MAP_FILE="$EXON_JUNCTION_MAP_FILE"
if [ "$COMPUTE_JUNCTION_COUNT" = "TRUE" ]
then
  ## This file depends on the junction alignments of the chunk and needs to be computed
  if [ "$NON_ZERO_JUNCTION_COUNT_OPTION" != "" ]
  then 
    FILTERED_EXON_JUNCTION_MAP_FILE="$SAM_DIR/${GENE_MODEL_PREFIX}_exon_junction.map"
  elif [ ! -f $EXON_JUNCTION_MAP_FILE ]
  then
    echo "$EXON_JUNCTION_MAP_FILE not found ... exiting"
    exit 1
  fi
fi


JUNCTION_JUNCTION_MAP_FILE=$PROJECT_MAP_DIR/${GENE_MODEL_PREFIX}_${RADIUS}_junction_junction_coordinate.map.gz
FILTERED_JUNCTION_JUNCTION_MAP_FILE="$JUNCTION_JUNCTION_MAP_FILE"
if [ "$COMPUTE_JUNCTION_SAM_COUNT" = "TRUE" ]
then
  ## This file also depends on the junction alignments of the chunk and needs to be computed
  if [ "$NON_ZERO_JUNCTION_COUNT_OPTION" != "" ]
  then 
    FILTERED_JUNCTION_JUNCTION_MAP_FILE="$SAM_DIR/${GENE_MODEL_PREFIX}_${RADIUS}_junction_junction.map"
  elif [ ! -f $JUNCTION_JUNCTION_MAP_FILE ]
  then
    echo "$JUNCTION_JUNCTION_MAP_FILE does not exist ... exiting"
    exit 1
  fi
fi



################################################################################
##
##  Output meta information to log file
##
################################################################################

echo "READ_WEIGHT_THRESHOLD=$READ_WEIGHT_THRESHOLD"
echo "EXON_OVERLAP=$EXON_OVERLAP"
echo "JUNCTION_OVERLAP=$JUNCTION_OVERLAP"
echo "STRAND_SPECIFIC=$STRAND_SPECIFIC"
echo "STRAND_SPECIFIC_DIRECTION=$STRAND_SPECIFIC_DIRECTION"
echo "PAIRED_END=$PAIRED_END"


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

WEIGHT_DIR=`echo $SAM_DIR | sed -e "s;/sam-files;/weight-files;"`
if [ ! -d $WEIGHT_DIR ]
then
  echo "Making directory $WEIGHT_DIR"
  mkdir $WEIGHT_DIR
fi

COUNT_DIR=`echo $SAM_DIR | sed -e "s;/sam-files;/count-files;"`
if [ ! -d $COUNT_DIR ]
then
  echo "Making directory $COUNT_DIR"
  mkdir $COUNT_DIR
fi

date



################################################################################
##
##  Create combined SAM file name
##
################################################################################

if [ "$GENOME_SAM_FILE_BASE" != "" ]
then
  export GENOME_SAM_WEIGHT_FILE="$GENOME_SAM_FILE_BASE-sam.wgt"
fi

export COMBINED_SAM_FILE_BASE=`echo $TRANSCRIPT_SAM_FILE_BASE | sed -e "s/-transcript/-combined/"`
if [ "$COMBINED_SAM_FILE_BASE" = "$TRANSCRIPT_SAM_FILE_BASE" ]
then
  COMBINED_SAM_FILE_BASE="$TRANSCRIPT_SAM_FILE_BASE-combined"
fi
COMBINED_SAM_FILE="$COMBINED_SAM_FILE_BASE.sam.gz"


################################################################################
##
##  Check for reference id files
##
################################################################################

##  Check for the existence of the combined reference id file
##  (it is created by compute-intersection-bed-file.sh)
if [ ! -r $SAM_DIR/$COMBINED_SAM_FILE_BASE-reference.ids -o ! -r $SAM_DIR/$COMBINED_SAM_FILE ]
then
  if [ "$RECOMPUTE" = "FALSE" ]
  then
    RECOMPUTE="TRUE"
    BED_OPTION_LIST="$BED_OPTION_LIST -r"
  fi
fi


################################################################################
##
##  Compute intersection BED file
##
################################################################################

export INTERSECTION_BED_FILE="$COMBINED_SAM_FILE_BASE-intersection.bed.gz"
## Call script to compute INTERSECTION_BED_FILE
if [ "$RECOMPUTE" = "TRUE" ]
then
  rm -f $BED_DIR/$INTERSECTION_BED_FILE
fi

if [ ! -f $BED_DIR/$INTERSECTION_BED_FILE ]
then
  echo "Call script $UTIL_BIN_DIR/compute-intersection-bed-file.sh"
  $UTIL_BIN_DIR/compute-intersection-bed-file.sh -c $BED_OPTION_LIST $BED_ARGUMENT_LIST
  if [ $? -ne 0 ]
  then
    echo "ERROR: Problem with compute-intersection-bed-file.sh -c $BED_OPTION_LIST $BED_ARGUMENT_LIST ... exiting."
    OUTPUT_FILE="$BED_DIR/$INTERSECTION_BED_FILE"
    if [ -f $OUTPUT_FILE ]
    then
      rm $OUTPUT_FILE
    fi	  
    exit 1
  fi
  echo "Done"
  date
fi


################################################################################
##
##  Filter the junction map files
##
################################################################################

if [ "$RECOMPUTE" = "TRUE" ]
then
  if [ "$FILTERED_EXON_JUNCTION_MAP_FILE" != "$EXON_JUNCTION_MAP_FILE" ]
  then
    rm -f $FILTERED_EXON_JUNCTION_MAP_FILE
  fi
  if [ "$FILTERED_JUNCTION_JUNCTION_MAP_FILE" != "$JUNCTION_JUNCTION_MAP_FILE" ]
  then
    rm -f $FILTERED_JUNCTION_JUNCTION_MAP_FILE
  fi
fi

RECOMPUTE_FILTER_FILES="FALSE"
if [ "$COMPUTE_JUNCTION_COUNT" = "TRUE" -a ! -f $FILTERED_EXON_JUNCTION_MAP_FILE ]
then
  RECOMPUTE_FILTER_FILES="TRUE"
fi

JUNCTION_JUNCTION_FILTER_OPTION=""
if [ "$COMPUTE_JUNCTION_SAM_COUNT" = "TRUE" -a ! -f $FILTERED_JUNCTION_JUNCTION_MAP_FILE ]
then
  RECOMPUTE_FILTER_FILES="TRUE"
  JUNCTION_JUNCTION_FILTER_OPTION="-j $FILTERED_JUNCTION_JUNCTION_MAP_FILE"
fi

if [ "$RECOMPUTE_FILTER_FILES" = "TRUE" ]
then
  echo "Filtering the junction map files with new script" > $SAM_DIR/filter-junction-map.log
  $UTIL_BIN_DIR/filter-junction-map-files.sh $JUNCTION_JUNCTION_FILTER_OPTION $OUTPUT_PREFIX_OPTION $PROJECT_DIR \
    $SAM_DIR/$COMBINED_SAM_FILE_BASE-reference.ids $FILTERED_EXON_JUNCTION_MAP_FILE >> $SAM_DIR/filter-junction-map.log 2>&1 &
  FILTER_JUNCTION_MAP_PID=$!
fi

################################################################################
##
##  Compute Bed file reads weights and combine them with read weights of
##  genomic SAM file
##
##  Note that though we may keep a read alignment from the transcript/junction
##  alignments the read may align in the genome and have many more alignments.
##
################################################################################

COMBINED_BED_WEIGHT_FILE="$COMBINED_SAM_FILE_BASE-bed.wgt"

FINAL_WEIGHT_FILE_BASE=`echo $COMBINED_SAM_FILE_BASE | sed -e "s/-combined.*/-final/"`
FINAL_WEIGHT_FILE="$FINAL_WEIGHT_FILE_BASE.wgt"

############################### Combine the weight files #######################

if [ "$GENOME_SAM_WEIGHT_FILE" != "" -a "$GENOME_FILE_MISSING" != "TRUE" ]
then
  if [ -f $WEIGHT_DIR/$GENOME_SAM_WEIGHT_FILE ]
  then
    if [ "$RECOMPUTE" = "TRUE" -o ! -f $WEIGHT_DIR/$FINAL_WEIGHT_FILE ]
    then

      if [ "$JUNCTION_FILE_MISSING" != "TRUE" ]
      then
        TRANSCRIPT_JUNCTION_SAM_FILE_BASE=`echo $COMBINED_SAM_FILE_BASE | sed -e "s/-combined-/-trans-junction-/"`
      else
        TRANSCRIPT_JUNCTION_SAM_FILE_BASE=`echo $COMBINED_SAM_FILE_BASE | sed -e "s/-combined-/-transcript-/"`
      fi
      TRANSCRIPT_JUNCTION_SAM_EDIT_DISTANCE_FILE="$TRANSCRIPT_JUNCTION_SAM_FILE_BASE-sam.wgt"

      if [ ! -f $WEIGHT_DIR/$TRANSCRIPT_JUNCTION_SAM_EDIT_DISTANCE_FILE ]
      then
        echo "File $WEIGHT_DIR/$TRANSCRIPT_JUNCTION_SAM_EDIT_DISTANCE_FILE not found ... exiting."
	exit 1
      fi

      ## Java command to compute BED read weights
      JAVA_BED_READ_WEIGHT_CMD="ComputeReadWeights -b - -o -"

      ## Java command to combine the genome and BED read weights
      JAVA_COMBINE_READ_WEIGHT_CMD="CombineReadWeightFiles -1 - -2 $WEIGHT_DIR/$GENOME_SAM_WEIGHT_FILE \
           -a $WEIGHT_DIR/$TRANSCRIPT_JUNCTION_SAM_EDIT_DISTANCE_FILE -o $WEIGHT_DIR/$FINAL_WEIGHT_FILE"

      echo "Combining BED read weights with $GENOME_SAM_WEIGHT_FILE"
      echo "using $TRANSCRIPT_JUNCTION_SAM_EDIT_DISTANCE_FILE writing to $FINAL_WEIGHT_FILE"
      echo "zcat $BED_DIR/$INTERSECTION_BED_FILE |"
      echo "$JAVA_BED_READ_WEIGHT_CMD |"
      echo "$JAVA_COMBINE_READ_WEIGHT_CMD"

      zcat $BED_DIR/$INTERSECTION_BED_FILE | $JAVA $JAVA_BED_READ_WEIGHT_CMD | $JAVA $JAVA_COMBINE_READ_WEIGHT_CMD
  
      if [ $? -ne 0 ]
      then
        echo "ERROR: Problem with $JAVA_BED_READ_WEIGHT_CMD | $JAVA_COMBINE_READ_WEIGHT_CMD ... exiting."
  	  OUTPUT_FILE="$WEIGHT_DIR/$FINAL_WEIGHT_FILE"
        if [ -f $OUTPUT_FILE ]
        then
          rm $OUTPUT_FILE
        fi
        exit 1
      fi
      echo "Done"
      date
    fi
  else
    echo "File $GENOME_SAM_WEIGHT_FILE not found ... exiting."
    exit 1
  fi
else
  if [ ! -f $FINAL_WEIGHT_FILE -o "$RECOMPUTE" = "TRUE" ]
  then
    echo "Computing BED read weights"
    JAVA_BED_READ_WEIGHT_CMD=" ComputeReadWeights -b - -o $FINAL_WEIGHT_FILE"
    echo $JAVA_BED_READ_WEIGHT_CMD
    zcat $BED_DIR/$INTERSECTION_BED_FILE | $JAVA $JAVA_BED_READ_WEIGHT_CMD
  
    if [ $? -ne 0 ]
    then
      echo "ERROR: Problem with $JAVA_BED_READ_WEIGHT_CMD ... exiting."
      OUTPUT_FILE="$FINAL_WEIGHT_FILE"
      if [ -f $OUTPUT_FILE ]
      then
        rm $OUTPUT_FILE
      fi
      exit 1
    fi
    echo "Done"
    date
  fi
fi


################################################################################
##
##  Set weight option
##
################################################################################

WEIGHT_OPTION="-w $WEIGHT_DIR/$FINAL_WEIGHT_FILE -W $READ_WEIGHT_THRESHOLD"
if [ "$COMPUTE_UNWEIGHTED" = "TRUE" ]
then
  WEIGHT_OPTION="-u"  
fi


################################################################################
##
##
##  Compute counts
##
##
################################################################################

## We always recompute the counts
RECOMPUTE="TRUE"


################################################################################
##
##  Compute gene counts
##
################################################################################

if [ "$OUTPUT_PREFIX" = "" ]
then
  GENE_COUNT_FILE="$COUNT_DIR/$COMBINED_SAM_FILE_BASE-gene.cnt"
else
  GENE_COUNT_FILE="$OUTPUT_PREFIX-gene.cnt"
fi

if [ "$COMPUTE_GENE_COUNT" = "TRUE" ]
then

  if [ "$RECOMPUTE" = "TRUE" -o ! -f $GENE_COUNT_FILE ]
  then
    GENE_COUNT_OPTION="-g"
    if [ "$OVERLAP_MODE" = "TRUE" ]
    then
      GENE_COUNT_OPTION="-e -O 1"
    fi
    
    echo "Computing gene counts"
    GENE_COUNT_JAVA_CMD="ComputeCounts $GENE_COUNT_OPTION $WEIGHT_OPTION $COUNT_STRAND_SPECIFIC_OPTION $UNAMBIGUOUS_OPTION \
        $CONTAINMENT_OPTION $NONSPLICE_CONFORMING_OPTION -m $EXON_GENE_MAP_FILE -b - -o -"
    echo "Java cmd: $GENE_COUNT_JAVA_CMD"
    zcat $BED_DIR/$INTERSECTION_BED_FILE | $JAVA $GENE_COUNT_JAVA_CMD > $GENE_COUNT_FILE &
    GENE_COUNT_PID=$!
  fi
fi


################################################################################
##
##  Compute exon counts
##
################################################################################

if [ "$OUTPUT_PREFIX" = "" ]
then
  EXON_COUNT_FILE="$COUNT_DIR/$COMBINED_SAM_FILE_BASE-exon.cnt"
else
  EXON_COUNT_FILE="$OUTPUT_PREFIX-exon.cnt"
fi

if [ "$COMPUTE_EXON_COUNT" = "TRUE" ]
then 
  if [ "$RECOMPUTE" = "TRUE" -o ! -f $EXON_COUNT_FILE ]
  then
    echo "Computing exon counts" > $COUNT_DIR/exon-count.log
    EXON_COUNT_JAVA_CMD="ComputeCounts -e $WEIGHT_OPTION $COUNT_STRAND_SPECIFIC_OPTION $NONSPLICE_CONFORMING_OPTION $UNAMBIGUOUS_OPTION \
       -O $EXON_OVERLAP -m $EXON_EXON_MAP_FILE -b - -o $EXON_COUNT_FILE"
    echo "$EXON_COUNT_JAVA_CMD" >> $COUNT_DIR/exon-count.log
    zcat $BED_DIR/$INTERSECTION_BED_FILE | $JAVA $EXON_COUNT_JAVA_CMD >> $COUNT_DIR/exon-count.log 2>&1 &
    EXON_COUNT_PID=$!
  fi
fi


################################################################################
##
##  Wait for the processes to finish and get their exit status
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
  waitPid $EXON_COUNT_PID $EXON_COUNT_FILE "$EXON_COUNT_JAVA_CMD" "Waiting for computation of exon counts to finish" $COUNT_DIR/exon-count.log
fi


################################################################################
##
##  Wait for the filtering of the junction map files to finish
##
################################################################################

if [ "$FILTER_JUNCTION_MAP_PID" != "" ]
then
  waitPid $FILTER_JUNCTION_MAP_PID $FILTERED_EXON_JUNCTION_MAP_FILE "compute the filtering of the junction map files" \
    "Waiting for the completion of the filtering of the junction map files" $SAM_DIR/filter-junction-map.log
fi

if [ "$EXIT_STATUS" != "0" ]
then
  echo "Compute counts failed."
  exit $EXIT_STATUS
fi


################################################################################
##
##  Compute junction counts
##
################################################################################

## Scale up the memory requirements of JAVA
JAVA="java -oss8M -ss8M -ms$MAX_MEMORY -mx$MAX_MEMORY -cp ${JAVA_CLASS_DIR}:${CLASSPATH}"
echo "JAVA=$JAVA"

if [ "$OUTPUT_PREFIX" = "" ]
then
  JUNCTION_COUNT_FILE="$COUNT_DIR/$COMBINED_SAM_FILE_BASE-junction.cnt"
else
  JUNCTION_COUNT_FILE="$OUTPUT_PREFIX-junction.cnt"
fi

if [ "$COMPUTE_JUNCTION_COUNT" = "TRUE" ]
then 
  if [ "$RECOMPUTE" = "TRUE" -o ! -f $JUNCTION_COUNT_FILE ]
  then
    echo "Computing junction counts"
    JUNCTION_COUNT_JAVA_CMD="ComputeCounts -j $NON_ZERO_JUNCTION_COUNT_OPTION $COUNT_STRAND_SPECIFIC_OPTION -m $FILTERED_EXON_JUNCTION_MAP_FILE \
       -O $JUNCTION_OVERLAP -b - $WEIGHT_OPTION -o $JUNCTION_COUNT_FILE"
    echo "$JUNCTION_COUNT_JAVA_CMD"
    zcat $BED_DIR/$INTERSECTION_BED_FILE | $JAVA $JUNCTION_COUNT_JAVA_CMD
  fi
fi


################################################################################
##
##  Compute junction SAM counts
##
################################################################################

if [ "$OUTPUT_PREFIX" = "" ]
then
  JUNCTION_SAM_COUNT_FILE="$COUNT_DIR/$COMBINED_SAM_FILE_BASE-junction-sam.cnt"
else
  JUNCTION_SAM_COUNT_FILE="$OUTPUT_PREFIX-junction-sam.cnt"
fi

if [ "$COMPUTE_JUNCTION_SAM_COUNT" = "TRUE" ]
then
  if [ "$JUNCTION_SAM_FILE" = "" ]
  then
    echo "Please specify the junction SAM file to compute the junction counts based on"
    echo "the junction alignments ... exiting"
    exit 1
  fi

  if [ "$JUNCTION_SAM_FILE" != "" -a \( "$RECOMPUTE" = "TRUE" -o ! -f $JUNCTION_SAM_COUNT_FILE \) ]
  then
    echo "Computing junction SAM counts"
    JUNCTION_SAM_COUNT_JAVA_CMD="ComputeGeneCountsSam $NON_ZERO_JUNCTION_COUNT_OPTION -m $FILTERED_JUNCTION_JUNCTION_MAP_FILE $WEIGHT_OPTION -s - \
       -O $JUNCTION_OVERLAP -o $JUNCTION_SAM_COUNT_FILE"
    echo "$JUNCTION_SAM_COUNT_JAVA_CMD"
    zcat $JUNCTION_SAM_PATH | $JAVA $JUNCTION_SAM_COUNT_JAVA_CMD
  fi
fi


################################################################################
##
##  Clean-up and exit
##
################################################################################


if [ "$CLEAN_UP" = "TRUE" ]
then
  if [ "$FILTERED_EXON_JUNCTION_MAP_FILE" != "$EXON_JUNCTION_MAP_FILE" ]
  then
    rm -f $FILTERED_EXON_JUNCTION_MAP_FILE
  fi
  if [ "$FILTERED_JUNCTION_JUNCTION_MAP_FILE" != "$JUNCTION_JUNCTION_MAP_FILE" ]
  then
    rm -f $FILTERED_JUNCTION_JUNCTION_MAP_FILE
  fi
  
  rm -r $BED_DIR
  rm -r $WEIGHT_DIR
  rm -r $SAM_DIR/*.ids $SAM_DIR/*.ids $SAM_DIR/$COMBINED_SAM_FILE
fi

if [ "$OUTPUT_PREFIX" != "" ]
then
  wait
  echo "Removing intermediate results in directory $OUTPUT_PREFIX"
  rm -r $OUTPUT_PREFIX
fi

date
if [ "$EXIT_STATUS" = "0" ]
then
  echo "Compute counts successfully completed."
else
  echo "Compute counts failed."
fi

exit $EXIT_STATUS

