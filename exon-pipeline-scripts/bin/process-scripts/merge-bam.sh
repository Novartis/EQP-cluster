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
KEEP_BAM_FILES="FALSE"
MAPPED_READS_OPTION=""
while [ "$1" = "-k" -o "$1" = "-m" -o "$1" = "--help" -o "$1" = "-h" ]
do
  if [ "$1" = "-h" -o "$1" = "--help" ]
  then
    shift
    HELP="TRUE"
  fi
  
  if [ "$1" = "-k" ]
  then
    shift
    KEEP_BAM_FILES="TRUE"
  fi

  if [ "$1" = "-m" ]
  then
    shift
    MAPPED_READS_OPTION="-F 4"
  fi
 
done

################################################################################
##
##  Check number of arguments and show help
##
################################################################################

if [ "$HELP" = "TRUE" -o "$1" = "" -o "$2" = "" -o "$3" = "" ]
then
  echo "Usage: $PROG_NAME [-k] [-m] <project dir> <merged BAM file> <SAM/BAM file 1>"
  echo "         [<SAM/BAM file 2> ...]"
  echo
  echo "-k: keep the sorted BAM files that are generated (default: remove them)"
  echo "-m: discard unmapped reads and keep only mapped reads"
  echo
  echo "project dir: directory where the relevant files are stored"
  echo "merged BAM file: result of merging input BAM files"
  echo "SAM/BAM file 1,2,...: the input SAM or BAM files"
  exit
fi


################################################################################
##
##  Read arguments
##
################################################################################

PROJECT_DIR=$1
shift

MERGED_BAM_FILE=$1
shift

OUTPUT_DIR=`dirname $MERGED_BAM_FILE`

INPUT_FILES=$*

echo "Merging and sorting file(s):"
NUM_INPUT_FILES=0
for INPUT_FILE in $INPUT_FILES
do
  echo $INPUT_FILE
  NUM_INPUT_FILES=`echo $NUM_INPUT_FILES + 1 | bc`
done
echo "into $MERGED_BAM_FILE"

if [ -f $MERGED_BAM_FILE ]
then
  rm $MERGED_BAM_FILE
fi


################################################################################
##
##  Set directories and files
##
################################################################################

TOOLS_DIR=$PROJECT_DIR/exon-pipeline-scripts/tools
checkTool $TOOLS_DIR samtools 0.1.17
SAMTOOLS_EXE=$TOOL_EXE

GENOME_DIR=$PROJECT_DIR/exon-pipeline-files/genome-files
GENOME_FASTA_INDEX_FILE=$GENOME_DIR/genome.fa.fai


################################################################################
##
##  Convert to BAM format and sort input files
##
################################################################################

# Change to the output directory to take care of temp files
cd $OUTPUT_DIR

SORTED_BAM_FILES=""
REMOVE_BAM_FILE_LIST=""
for INPUT_FILE_PATH in $INPUT_FILES
do
  INPUT_DIR=`dirname $INPUT_FILE_PATH`
  INPUT_FILE=`basename $INPUT_FILE_PATH`
  INPUT_FILE_BASE=`echo $INPUT_FILE | sed -e 's/.gz$//' | sed -e 's/.[bs]am$//'` 
  INPUT_FILE_EXT=`echo $INPUT_FILE | sed -e "s/$INPUT_FILE_BASE[.]//"`
  OUTPUT_DIR=`dirname $MERGED_BAM_FILE`
  SORTED_BAM_FILE_BASE=$OUTPUT_DIR/$INPUT_FILE_BASE-sorted
  if [ "$INPUT_FILE_EXT" = "bam" ]
  then
  
    echo "Sorting file $INPUT_FILE"
    echo "Pipeline calls:"
    echo "  $SAMTOOLS_EXE view -b $MAPPED_READS_OPTION -t $GENOME_FASTA_INDEX_FILE $INPUT_FILE_PATH | \ "
    echo "  $SAMTOOLS_EXE sort -m 6000000000 - $SORTED_BAM_FILE_BASE"
    
    $SAMTOOLS_EXE view -b $MAPPED_READS_OPTION -t $GENOME_FASTA_INDEX_FILE $INPUT_FILE_PATH | \
      $SAMTOOLS_EXE sort -m 6000000000 - $SORTED_BAM_FILE_BASE
    if [ $? -ne 0 ]
    then
      echo "Problem with samtools view/sort $INPUT_FILE_PATH ... exiting"
      exit 1
    fi
  else
    if [ "$INPUT_FILE_EXT" = "sam" ]
    then
      CAT="cat"
    elif [ "$INPUT_FILE_EXT" = "sam.gz" ]
    then
      CAT="zcat"
    else
      echo "Unknown extension: $INPUT_FILE_EXT for file $INPUT_FILE_PATH ... exiting"
      exit 1
    fi

    echo "Converting to BAM format and sorting file $INPUT_FILE"
    echo "Pipeline calls:"
    echo "  $CAT $INPUT_FILE_PATH | \ "
    echo "  $SAMTOOLS_EXE view -bS $MAPPED_READS_OPTION -t $GENOME_FASTA_INDEX_FILE - | \ "
    echo "  $SAMTOOLS_EXE sort -m 6000000000 - $SORTED_BAM_FILE_BASE"
    
    $CAT $INPUT_FILE_PATH | \
    $SAMTOOLS_EXE view -bS $MAPPED_READS_OPTION -t $GENOME_FASTA_INDEX_FILE - | \
    $SAMTOOLS_EXE sort -m 6000000000 - $SORTED_BAM_FILE_BASE
    if [ $? -ne 0 ]
    then
      echo "Problem with samtools view/sort $INPUT_FILE_PATH ... exiting"
      exit 1
    fi
    REMOVE_BAM_FILE_LIST="$REMOVE_BAM_FILE_LIST $SORTED_BAM_FILE_BASE.bam"
  fi
  SORTED_BAM_FILES="$SORTED_BAM_FILES $SORTED_BAM_FILE_BASE.bam"
done

if [ $NUM_INPUT_FILES -gt 1 ]
then
  $SAMTOOLS_EXE merge $MERGED_BAM_FILE $SORTED_BAM_FILES
  if [ $? -ne 0 ]
  then
    echo "Problem with samtools merge $MERGED_BAM_FILE $SORTED_BAM_FILES ... exiting"
    exit 1
  fi
elif [ -f $SORTED_BAM_FILES ]
then
  mv $SORTED_BAM_FILES $MERGED_BAM_FILE
else
  echo "File $SORTED_BAM_FILES does not exist ... exiting."
  exit 1
fi

echo "Indexing BAM file $MERGED_BAM_FILE"
$SAMTOOLS_EXE index $MERGED_BAM_FILE
if [ $? -ne 0 ]
then
  echo "Problem with samtools index $MERGED_BAM_FILE ... exiting"
  exit 1
fi

if [ "$KEEP_BAM_FILES" = "FALSE" -a $NUM_INPUT_FILES -gt 1 ]
then
  echo "Removing sorted BAM files:"
  echo $REMOVE_BAM_FILE_LIST
  rm -f $REMOVE_BAM_FILE_LIST
fi

echo "Merging of BAM files successfully completed."
date
