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
#$ -pe smp 4 -R y

PROG_NAME=`basename $0`
PROG_DIR=`dirname $0`
VERSION=1.0

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

HELP="FALSE"
PICARD_MAX_RECORDS_IN_RAM="5000000"
SORT_ORDER="queryname"
INDEX_BAM_FILE="FALSE"
MERGE_BAM_FILES="FALSE"
while [ "$1" = "-M" -o "$1" = "-c" -o "$1" = "-m" -o "$1" = "-i" -o "$1" = "-h" -o "$1" = "--help" ]
do
  if [ "$1" = "-M" ]
  then
    shift
    PICARD_MAX_RECORDS_IN_RAM="$1"
    shift
  fi

  if [ "$1" = "-c" ]
  then
    shift
    SORT_ORDER="coordinate"
  fi

  if [ "$1" = "-i" ]
  then
    shift
    INDEX_BAM_FILE="TRUE"
  fi

  if [ "$1" = "-m" ]
  then
    shift
    MERGE_BAM_FILES="TRUE"
  fi  


  if [ "$1" = "-h" -o "$1" = "--help" ]
  then
    shift
    HELP="TRUE"
  fi

done

################################################################################
##
##  Check number of arguments and show help
##
################################################################################

if [ "$HELP" = "TRUE" -o "$1" = "" -o "$2" = "" -o "$3" = "" ]
then
  echo "Usage: $PROG_NAME <options> <project dir> <SAM file> <sorted SAM file> "
  echo
  echo "where <options> is"
  echo "   [-M <max records>] [-c]"
  echo
  echo "project dir: directory where the relevant files are stored"
  echo "SAM file: file to be sorted (use /dev/stdin for stdin)"
  echo "sorted SAM file: the sorted output SAM file"
  echo
  echo " -M INT: maximum number of records allowed in RAM [5000000]"
  echo " -c: sort by coordinate (default is sort by queryname)"
  exit
fi


################################################################################
##
##  Read arguments
##
################################################################################

PROJECT_DIR=$1
shift

if [ ! -d "$PROJECT_DIR" ]
then
  echo "Directory $PROJECT_DIR not found ... exiting"
  exit 1
fi

if [ ! -d $PROJECT_DIR/.temp ]
then
  mkdir $PROJECT_DIR/.temp
fi
export TEMP_DIR=$PROJECT_DIR/.temp

# Use /dev/stdout to write to stdout (BAM format only)
SORTED_SAM_FILE=$1
shift

# Use /dev/stdin to read from stdin
SAM_FILE=$1
shift

if [ "$SAM_FILE" != "/dev/stdin" -a ! -f "$SAM_FILE" ]
then
  echo "$SAM_FILE not found ... exiting."
  exit 1
fi


FIRST_GZIP=`echo "$SAM_FILE" | sed -e 's/.*gz$/gz/'`

SAM_FILES=$SAM_FILE
while [ "$SAM_FILE" != "" -a "$SAM_FILE" != "/dev/stdin" ]
do
  # Use /dev/stdin to read from stdin
  SAM_FILE=$1
  shift

  if [ "$SAM_FILE" != "" ]
  then
    if [ ! -f "$SAM_FILE" ]
    then
      echo "$SAM_FILE not found ... exiting."
      exit 1
    fi

    GZIP=`echo "$SAM_FILE" | sed -e 's/.*gz$/gz/'`

    if [ \( "$GZIP" = "gz" -a "$FIRST_GZIP" != "gz" \) -o \( "$GZIP" != "gz" -a "$FIRST_GZIP" = "gz" \) ]
    then
      echo "Mix of gzipped and umcompressed files as submitted input ... exiting."
      exit 1
    fi

    BAM=`echo "$SAM_FILE" | sed -e 's/.*bam$/bam/'`

    if [ "$MERGE_BAM_FILES" = "TRUE" -a "$BAM" != "bam" ]
    then
      echo "None BAM submitted for merging ... exiting."
      exit 1
    fi

    if [ "$MERGE_BAM_FILES" = "FALSE" -a "$BAM" = "bam" ]
    then
      echo "BAM submitted for sorting ... exiting."
      exit 1
    fi

    SAM_FILES="$SAM_FILES $SAM_FILE"
  fi

done


if [ "$SAM_FILES" = "/dev/stdin" ]
then
  STREAM_INPUT=""
elif [ "$FIRST_GZIP" = "gz" ]
then
  STREAM_INPUT="zcat $SAM_FILES |"
else
  STREAM_INPUT="cat $SAM_FILES |"
fi
  


################################################################################
##
##  Set java environment
##
################################################################################

if [ "$PICARD_JAR_DIR" = "" ]
then
  PICARD_JAR_DIR=$PROJECT_DIR/picard-scripts/jar-files
fi

if [ ! -d $PICARD_JAR_DIR ]
then
  echo "Jar directory $PICARD_JAR_DIR not found ... exiting"
  exit 1
fi

if [ "$TOOLS_DIR" = "" ]
then
  TOOLS_DIR=$PROJECT_DIR/exon-pipeline-scripts/tools
fi
checkTool $TOOLS_DIR samtools 0.1.17
SAMTOOLS_EXE=$TOOL_EXE


if [ ! -d $TOOLS_DIR ]
then
  echo "Directory $TOOLS_DIR not found ... exiting"
  exit 1
fi


JAVA_MEM_GB=`echo "22 * $PICARD_MAX_RECORDS_IN_RAM / 5000000" | bc`
JAVA="java -oss8M -ss8M -ms${JAVA_MEM_GB}G -mx${JAVA_MEM_GB}G"

echo $JAVA


################################################################################
##
##  Run SortSam
##
################################################################################

if [ "$MERGE_BAM_FILES" != "TRUE" ]
then
  CMD="$STREAM_INPUT \
      $JAVA -jar $PICARD_JAR_DIR/SortSam.jar \
      INPUT=/dev/stdin \
      OUTPUT=$SORTED_SAM_FILE \
      SORT_ORDER=$SORT_ORDER \
      MAX_RECORDS_IN_RAM=$PICARD_MAX_RECORDS_IN_RAM \
      VALIDATION_STRINGENCY=LENIENT \
      TMP_DIR=$TEMP_DIR \
      QUIET=TRUE"

  echo "Sorting SAM files:"
  for SAM_FILE in $SAM_FILES
  do
    echo $SAM_FILE
  done
  echo "into $SORTED_SAM_FILE"
  
  echo "$CMD"

  eval $CMD

  if [ $? -ne 0 ]
  then
    echo "Problem with Picard SortSam ... exiting"
    echo "Input: $SAM_FILE"
    echo "Output: $SORTED_SAM_FILE"
    echo "Sort order: $SORT_ORDER"
    exit 1
  fi

else
  echo "Merging BAM files:"
  for BAM_FILE in $SAM_FILES
  do
    echo "$BAM_FILE"
  done
  echo "into $SORTED_SAM_FILE"

  if [ -f $SORTED_SAM_FILE ]
  then
    rm $SORTED_SAM_FILE
  fi
  
  $SAMTOOLS_EXE merge $SORTED_SAM_FILE $SAM_FILES

  if [ $? -ne 0 ]
  then
    echo "Problem with samtools merge ... exiting"
    exit 1
  fi

fi

BAM=`echo "$SORTED_SAM_FILE" | sed -e 's/.*bam$/bam/'`

if [ "$INDEX_BAM_FILE" = "TRUE" ]
then
  if [ "$BAM" = "bam" ]
  then
    echo "Indexing BAM file $SORTED_SAM_FILE"
    $SAMTOOLS_EXE index $SORTED_SAM_FILE

    if [ $? -ne 0 ]
    then
      echo "Problem with samtools index ... exiting"
      exit 1
    fi
  else
    echo "Result file $SORTED_SAM_FILE is not a BAM file and cannot be indexed."
    exit 1
  fi

fi

echo "Done"
date
