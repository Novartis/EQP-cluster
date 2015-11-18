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

RECOMPUTE="FALSE"
USE_GTF=FALSE
EDIT_DISTANCE_OPTION=
OUTPUT_PREFIX=""
QUIET_OPTION=""
USE_PIPED_COMPUTATION="TRUE"
BED_OPTION_LIST=""
CREATE_BAM_FILE="FALSE"
while [ "$1" = "-r" -o "$1" = "-gtf" -o "$1" = "-d" -o "$1" = "-o" -o $1 = "-p" -o $1 = "-q" -o "$1" = "-bam" ]
do

  if [ "$1" = "-r" ]
  then
    shift
    RECOMPUTE="TRUE"
    echo "Compute genomic alignments: RECOMPUTE=TRUE"
    BED_OPTION_LIST="$BED_OPTION_LIST -r"
  fi
  
  if [ "$1" = "-gtf" ]
  then
    shift
    USE_GTF=TRUE
    BED_OPTION_LIST="$BED_OPTION_LIST -gtf"
  fi

  if [ "$1" = "-d" ]
  then
    shift
    EDIT_DISTANCE_OPTION="-e $1"
    BED_OPTION_LIST="$BED_OPTION_LIST -d $1"
    shift
  fi
  
  if [ "$1" = "-o" ]
  then
    shift
    OUTPUT_PREFIX=$1
    BED_OPTION_LIST="$BED_OPTION_LIST -o $1 -w"
    shift
  fi

  if [ "$1" = "-p" ]
  then
    shift
    USE_PIPED_COMPUTATION="TRUE"
  fi

  if [ "$1" = "-q" ]
  then
    shift
    QUIET_OPTION="TRUE"
  fi

  if [ "$1" = "-bam" ]
  then
    shift
    CREATE_BAM_FILE="TRUE"
  fi

done


################################################################################
##
##  Read arguments
##
################################################################################

if [ "$5" = "" ]
then
  echo "Usage: $PROG_NAME <options> <project dir> <transcript SAM file>"
  echo "  <junction SAM file> <genome SAM-file> <EQP genome SAM file>"
  echo
  echo "where <options> is"
  echo "   [-r] [-gtf] [-d <edit distance>] [-s <direction>]"
  echo "   [-L <read length>] [-o <output-prefix>] [-p] [-q]"
  echo
  echo " project dir: directory where the relevant files are stored"
  echo " transcript SAM file: the file containing the alignments of the"
  echo "     reads against the transcripts."
  echo " junction SAM file: the file containing the alignments of the"
  echo "     reads against the junctions."
  echo " genome SAM file: the file containing the alignments of the reads against"
  echo "    the genome."
  echo " EQP genome SAM file: the file containing the genomic spliced alignments"
  echo "     of the reads as computed by EQP."
  echo " -r: force the recomputation of SAM and BED files even if they exist"
  echo " -gtf: use GTF transcript file instead of the transcript file downloaded"
  echo "       from $GENE_MODEL"
  echo " -s STRING: direction of how to process the reads as strand-specific: forward"
  echo "     or backward"
  echo " -L INT: set read length to INT; otherwise the value in"
  echo "    <project dir>/bin/setup.sh is taken."
  echo " -d INT: edit distance threshold in conversion of SAM to BED file (only SAM"
  echo "     records with an edit distance of at most INT are converted)"
  echo " -o STRING: A directory STRING is created in the current working directory"
  echo "     and all intermediate files are stored in a directory structure under"
  echo "     this directory; furthermore, the count files are generated in the current"
  echo "     working directory with the prefix STRING."
  echo " -p: use a piped computation (instead of a sequential one)"
  echo " -q: do not output messages to stderr"
  exit
fi

PROJECT_DIR=$1
shift

BED_ARGUMENT_LIST="$PROJECT_DIR"

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


TRANSCRIPT_SAM_PATH=$1
shift
if [ ! -f $TRANSCRIPT_SAM_PATH ]
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
if [ ! -f $JUNCTION_SAM_PATH ]
then
  echo "WARNING: Junction SAM file not found ... no junction alignments will be used"
  export JUNCTION_FILE_MISSING="TRUE"
fi
BED_ARGUMENT_LIST="$BED_ARGUMENT_LIST $JUNCTION_SAM_PATH"


########################### Genome file ######################################
GENOME_SAM_PATH=$1
shift
if [ ! -f $GENOME_SAM_PATH ]
then
  echo "WARNING: Genome SAM file not found ... no genome alignments will be used"
  export GENOME_FILE_MISSING="TRUE"
  GENOME_FILE_OPTION=""
else
  GENOME_FILE_OPTION="-g $GENOME_SAM_PATH"
fi
BED_ARGUMENT_LIST="$BED_ARGUMENT_LIST $GENOME_SAM_PATH"


########################### EQP genome SAM file ######################################
EQP_GENOME_SAM_PATH="$1"
shift

EQP_GENOME_SAM_FILE=`basename $EQP_GENOME_SAM_PATH`
EQP_GENOME_SAM_FILE_BASE=`echo $EQP_GENOME_SAM_FILE | sed -e 's/.gz$//' | sed -e 's/.[bs]am$//'` 
EQP_GENOME_SAM_FILE_EXT=`echo $EQP_GENOME_SAM_FILE | sed -e "s/$EQP_GENOME_SAM_FILE_BASE[.]//"` 
EQP_GENOME_SAM_DIR=`dirname $EQP_GENOME_SAM_PATH`

if [ ! -d $EQP_GENOME_SAM_DIR ]
then
  echo "Creating directory $EQP_GENOME_SAM_DIR"
  mkdir $EQP_GENOME_SAM_DIR
fi

if [ "$1" != "" ]
then
  echo "Unused arguments: $*."
fi


################################################################################
##
##  Set and check directories
##
## SAM_DIR and BED_DIR are used to store intermediate results whereas
## EQP_GENOME_SAM_DIR is used to store the end result. Note that usually these
## variables point to the same directory but if, for instance, OUTPUT_PREFIX
## is set, then they point to different directories.
##
################################################################################

if [ "$OUTPUT_PREFIX" = "" ]
then
  SAM_DIR=`dirname $TRANSCRIPT_SAM_PATH`
  OLD_WD=`pwd`
  cd $SAM_DIR
  SAM_DIR=`pwd`
  cd $OLD_WD
  EQP_SAM_FILE_BASE="$EQP_GENOME_SAM_DIR/$EQP_GENOME_SAM_FILE_BASE"
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
  EQP_SAM_FILE_BASE="$OUTPUT_PREFIX-eqp-genome"
fi

PROJECT_SCRIPTS_DIR="$PROJECT_DIR/exon-pipeline-scripts"
if [ ! -d $PROJECT_SCRIPTS_DIR ]
then
  echo "Directory $PROJECT_SCRIPTS_DIR does not exist ... exiting"
  exit 1
fi

export JAVA_DIR=$PROJECT_SCRIPTS_DIR/java
if [ ! -d $JAVA_DIR ]
then
  echo "Directory $JAVA_DIR does not exist ... exiting"
  exit 1
fi

export BIN_DIR=$PROJECT_SCRIPTS_DIR/bin/process-scripts
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

export TOOLS_DIR=$PROJECT_SCRIPTS_DIR/tools
if [ -d $TOOLS_DIR ]
then
  PATH="${TOOLS_DIR}:$PATH"
fi
checkTool $TOOLS_DIR bedtools 2.24.0
BEDTOOLS_EXE=$TOOL_EXE
checkTool $TOOLS_DIR samtools 0.1.17
SAMTOOLS_EXE=$TOOL_EXE

export UTIL_LIB_DIR=$PROJECT_SCRIPTS_DIR/bin/process-scripts/util-scripts
if [ ! -d $UTIL_LIB_DIR ]
then
  echo "Directory $UTIL_LIB_DIR not found ... exiting"
  exit 1
fi

PROJECT_GENOME_DIR=$PROJECT_DIR/exon-pipeline-files/genome-files
if [ ! -d $PROJECT_GENOME_DIR ]
then
  echo "Directory $PROJECT_GENOME_DIR does not exist ... exiting"
  exit 1
fi

################################################################################
##
##  Set genome files
##
################################################################################

CHROMOSOME_ID_FILE=$PROJECT_GENOME_DIR/genome.fa.fai
if [ ! -f $CHROMOSOME_ID_FILE ]
then
  echo "$CHROMOSOME_ID_FILE not found ... exiting"
  exit 1
fi

GENOME_SAM_HEADER_FILE=$PROJECT_GENOME_DIR/genome-header.sam
if [ ! -f $GENOME_SAM_HEADER_FILE ]
then
  echo "File $GENOME_SAM_HEADER_FILE not found ... exiting"
  exit 1
fi

GENOME_FILE=$PROJECT_GENOME_DIR/genome.fa
if [ ! -f $GENOME_FILE ]
then
  echo "File $GENOME_FILE not found ... exiting"
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


################################################################################
##
##  Compute intersection BED file
##
################################################################################

export COMBINED_SAM_FILE_BASE=`echo $TRANSCRIPT_SAM_FILE_BASE | sed -e "s/-transcript/-combined/"`
if [ "$COMBINED_SAM_FILE_BASE" = "$TRANSCRIPT_SAM_FILE_BASE" ]
then
  COMBINED_SAM_FILE_BASE="$TRANSCRIPT_SAM_FILE_BASE-combined"
fi
export COMBINED_SAM_FILE="$COMBINED_SAM_FILE_BASE.$TRANSCRIPT_SAM_FILE_EXT"
export INTERSECTION_BED_FILE="$COMBINED_SAM_FILE_BASE-intersection.bed.gz"

## Call script to compute INTERSECTION_BED_FILE

if [ "$RECOMPUTE" = "TRUE" -o ! -f $BED_DIR/$INTERSECTION_BED_FILE -o ! -f $SAM_DIR/$COMBINED_SAM_FILE ]
then
  ## Options: -r -> recompute, -w -> wait (at the end of the script for completion of all background jobs)
  $UTIL_BIN_DIR/compute-intersection-bed-file.sh -r -w $BED_OPTION_LIST $BED_ARGUMENT_LIST
  if [ $? -ne 0 ]
  then
    echo "ERROR: Problem with $UTIL_BIN_DIR/compute-intersection-bed-file.sh -w $BED_OPTION_LIST $BED_ARGUMENT_LIST ... exiting."
    exit 1
  fi
fi

if [ ! -f $BED_DIR/$INTERSECTION_BED_FILE ]
then
  echo "File $INTERSECTION_BED_FILE not found ... exiting."
  exit 1
fi


################################################################################
##
##  Create a SAM header file
##  We just need the @HD and @PG from the transcript SAM header
##
################################################################################

SAM_HEADER_FILE="$SAM_DIR/header.sam"
if [ ! -f $SAM_HEADER_FILE ]
then
  zcat $TRANSCRIPT_SAM_PATH | head -1000000 | grep "^@" | egrep "@HD|@PG" | sed -e "s/ID:/ID:EQP2.0-/" > $SAM_HEADER_FILE
fi


################################################################################
##
##
##  Compute genomic SAM/BAM file
##
##
################################################################################

JAVA_CLASS_DIR="$JAVA_DIR/classes:$JAVA_DIR"
JAVA="java -oss8M -ss8M -ms6G -mx6G -cp ${JAVA_CLASS_DIR}:${CLASSPATH}"
echo "JAVA=$JAVA"


################################################################################
##
##  Compute simple SAM file
##
################################################################################

GENOME_FILE_OPTION=""

echo "Computing SAM file without correction"
echo "Command pipeline:"  
echo "  zcat $BED_DIR/$INTERSECTION_BED_FILE | \ "
echo "  \$JAVA ComputeGenomeSamFile -c $CHROMOSOME_ID_FILE -s $SAM_DIR/$COMBINED_SAM_FILE -b - -o - $QUIET_OPTION $GENOME_FILE_OPTION | \ "
echo "  cat $SAM_HEADER_FILE $GENOME_SAM_HEADER_FILE - | grep '[A-z0-9]' > $SAM_DIR/$EQP_GENOME_SAM_FILE_BASE-simple.sam"

zcat $BED_DIR/$INTERSECTION_BED_FILE | \
$JAVA ComputeGenomeSamFile -c $CHROMOSOME_ID_FILE -s $SAM_DIR/$COMBINED_SAM_FILE -b - -o - $QUIET_OPTION $GENOME_FILE_OPTION | \
cat $SAM_HEADER_FILE $GENOME_SAM_HEADER_FILE - | grep '[A-z0-9]' > $SAM_DIR/$EQP_GENOME_SAM_FILE_BASE-simple.sam
if [ $? -ne 0 ]
then
  echo "ERROR: Problem with creation of genome SAM file without correction ... exiting."
  OUTPUT_FILE="$SAM_DIR/$EQP_GENOME_SAM_FILE_BASE-simple.sam"
  if [ -f $OUTPUT_FILE ]
  then
    rm $OUTPUT_FILE
  fi
  exit 1
fi
date

## Piped version
if [ "$USE_PIPED_COMPUTATION" = "TRUE" ]
then
  
  echo "Computing final SAM file"
  echo "Command pipeline:"
  echo "  $TOOLS_DIR/samtools view -bS -h $SAM_DIR/$EQP_GENOME_SAM_FILE_BASE-simple.sam | \ "
  echo "  $BEDTOOLS_EXE bamtobed -split -bed12 -i stdin | \ "
  echo "  $BEDTOOLS_EXE getfasta -split -name -fi $GENOME_FILE -bed stdin -fo stdout | \ "
  echo "  $UTIL_LIB_DIR/adjustSamEntries.py -s $SAM_DIR/$EQP_GENOME_SAM_FILE_BASE-simple.sam -f - -o $SAM_DIR/$EQP_GENOME_SAM_FILE_BASE.sam"

  $TOOLS_DIR/samtools view -bS -h $SAM_DIR/$EQP_GENOME_SAM_FILE_BASE-simple.sam | \
  $BEDTOOLS_EXE bamtobed -split -bed12 -i stdin | \
  $BEDTOOLS_EXE getfasta -split -name -fi $GENOME_FILE -bed stdin -fo stdout > $SAM_DIR/$EQP_GENOME_SAM_FILE_BASE-simple.fa
  if [ $? -ne 0 ]
  then
    echo "ERROR: Problem with creation of EQP genome Fasta file ... exiting."
    OUTPUT_FILE="$EQP_GENOME_SAM_DIR/$EQP_GENOME_SAM_FILE_BASE-simple.fa"
    if [ -f $OUTPUT_FILE ]
    then
      rm $OUTPUT_FILE
    fi
    exit 1
  fi
  echo "Fasta file $EQP_GENOME_SAM_FILE_BASE-simple.fa created."
  date

  $UTIL_LIB_DIR/adjustSamEntries.py -s $SAM_DIR/$EQP_GENOME_SAM_FILE_BASE-simple.sam -f $SAM_DIR/$EQP_GENOME_SAM_FILE_BASE-simple.fa \
     -o $EQP_GENOME_SAM_DIR/$EQP_GENOME_SAM_FILE_BASE.sam
  if [ $? -ne 0 ]
  then
    echo "ERROR: Problem with creation of EQP genome SAM file ... exiting."
    OUTPUT_FILE="$EQP_GENOME_SAM_DIR/$EQP_GENOME_SAM_FILE_BASE.sam"
    if [ -f $OUTPUT_FILE ]
    then
      rm $OUTPUT_FILE
    fi
    exit 1
  fi
  date
else

  echo "Converting the genome SAM file to a BAM file"
  date
  $TOOLS_DIR/samtools view -bS -h $SAM_DIR/$EQP_GENOME_SAM_FILE_BASE-simple.sam > $SAM_DIR/$EQP_GENOME_SAM_FILE_BASE-simple.bam
  if [ $? -ne 0 ]
  then
    echo "ERROR: Problem with converstion of $SAM_DIR/$EQP_GENOME_SAM_FILE_BASE.sam to BAM file ... exiting."
    exit 1
  fi
  
  echo "Creating a genomic BED file from the BAM file"
  date
  $BEDTOOLS_EXE bamtobed -split -bed12 -i $SAM_DIR/$EQP_GENOME_SAM_FILE_BASE-simple.bam > $BED_DIR/$EQP_GENOME_SAM_FILE_BASE-simple.bed
  if [ $? -ne 0 ]
  then
    echo "ERROR: Problem with creation of the BED file $SAM_DIR/$EQP_GENOME_SAM_FILE_BASE-simple.bed from the BAM file ... exiting."
    exit 1
  fi

  echo "Creating a genomic Fasta file from the BED file"
  date
  $BEDTOOLS_EXE getfasta -split -name -fi $GENOME_FILE -bed $BED_DIR/$EQP_GENOME_SAM_FILE_BASE-simple.bed -fo $SAM_DIR/$EQP_GENOME_SAM_FILE_BASE-simple.fa
  if [ $? -ne 0 ]
  then
    echo "ERROR: Problem with creation of the Fasta file $SAM_DIR/$EQP_GENOME_SAM_FILE_BASE-simple.fa from the BAM file ... exiting."
    exit 1
  fi
  
  echo "Adjusting alignments of SAM entries"
  date
  $UTIL_LIB_DIR/adjustSamEntries.py -s $SAM_DIR/$EQP_GENOME_SAM_FILE_BASE-simple.sam -f $SAM_DIR/$EQP_GENOME_SAM_FILE_BASE-simple.fa \
     -o $EQP_GENOME_SAM_DIR/$EQP_GENOME_SAM_FILE_BASE.sam
  if [ $? -ne 0 ]
  then
    echo "ERROR: Problem with adjusting SAM entries ... exiting."
    exit 1
  fi
fi


################################################################################
##
##  Clean up
##
################################################################################

if [ "$OUTPUT_PREFIX" != "" ]
then
  if [ -d $OUTPUT_PREFIX ]
  then
    echo "Removing intermediate results in directory $OUTPUT_PREFIX"
    rm -rf $OUTPUT_PREFIX
  fi
else
  echo "Cleaning up"
  rm -f $SAM_DIR/$EQP_GENOME_SAM_FILE_BASE-simple.[bs]am $SAM_DIR/$EQP_GENOME_SAM_FILE_BASE-simple.fa
fi


################################################################################
##
##  Creating BAM file
##
################################################################################

if [ "$CREATE_BAM_FILE" = "TRUE" -o "$EQP_GENOME_SAM_FILE_EXT" = "bam" ]
then

  echo "Converting the adjusted genome SAM file to a BAM file"
  date
  $SAMTOOLS_EXE view -bS -h $EQP_GENOME_SAM_DIR/$EQP_GENOME_SAM_FILE_BASE.sam | \
  $SAMTOOLS_EXE sort -m 5000000000 - $EQP_GENOME_SAM_DIR/$EQP_GENOME_SAM_FILE_BASE
  if [ $? -ne 0 ]
  then
    echo "ERROR: Problem with converstion of $EQP_SAM_FILE_BASE.sam to BAM file ... exiting."
    exit 1
  fi
  
  echo "Indexing the genome BAM file"
  date
  $SAMTOOLS_EXE index $EQP_GENOME_SAM_DIR/$EQP_GENOME_SAM_FILE_BASE.bam
  if [ $? -ne 0 ]
  then
    echo "ERROR: Problem with indexing of $EQP_SAM_FILE_BASE.bam ... exiting."
    exit 1
  fi

  echo "Coordinate-sorted genome BAM file created."
  date
fi


################################################################################
##
##  Compressing SAM file
##
################################################################################

if [ "$EQP_GENOME_SAM_FILE_EXT" = "sam.gz" ]
then
  echo "Compressing $EQP_GENOME_SAM_FILE_BASE.sam"
  date
  gzip -f $EQP_GENOME_SAM_DIR/$EQP_GENOME_SAM_FILE_BASE.sam
  if [ $? -ne 0 ]
  then
    echo "ERROR: Problem with compressing $EQP_GENOME_SAM_FILE_BASE.sam ... exiting."
    exit 1
  fi
  echo "Compression done"
fi

echo "Compute EQP genome alignments successfully completed."
