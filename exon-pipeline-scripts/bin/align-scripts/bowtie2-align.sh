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
#$ -R y
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
      elif [ "$TOOL_NAME" = "bowtie2" ]
      then
        ACTUAL_TOOL_VERSION=`$TOOL_EXE --version | fgrep version | fgrep -v "Compiler" | sed -e "s/.* version \([0-9]*[.][0-9]*[.][0-9]*\).*/\1/"`
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
## Read options
##
################################################################################

PAIRED_END=FALSE
FIRST_FASTQ_FILE=
SECOND_FASTQ_FILE=
ALIGN_GENOME=FALSE
CREATE_BAM_FILE=FALSE
NUM_PROCESSORS=3
MAX_ALIGN_NUM=
GZIP=TRUE
COMMON_PREFIX=
SENSITIVITY_OPTION=""
SCORE_MIN_OPTION="--score-min L,0,-0.1"
SORT="FALSE"
STRAND_SPECIFIC_OPTION=""
MIN_FRAGMENT_LENGTH_OPTION=
while [ "$1" = "-1" -o "$1" = "-2" -o "$1" = "-g" -o "$1" = "-bam" -o "$1" = "-p" -o "$1" = "-m"  -o "$1" = "-nogzip" -o \
        "$1" = "-j" -o "$1" = "-C" -o "$1" = "-clean_up" -o "$1" = "-S" -o "$1" = "-sm" -o "$1" = "-sort" -o "$1" = "-I" -o \
	"$1" = "-norc" -o "$1" = "-nofw" ]
do
  if [ "$1" = "-1" ]
  then
    shift
    FIRST_FASTQ_FILE=$1
    shift
  fi

  if [ "$1" = "-2" ]
  then
    shift
    SECOND_FASTQ_FILE=$1
    shift
    PAIRED_END=TRUE
  fi

  if [ "$1" = "-g" ]
  then
    shift
    ALIGN_GENOME=TRUE
  fi

  if [ "$1" = "-bam" ]
  then
    shift
    CREATE_BAM_FILE=TRUE
  fi

  if [ "$1" = "-nogzip" ]
  then
    shift
    GZIP=FALSE
  fi

  if [ "$1" = "-p" ]
  then
    shift
    NUM_PROCESSORS=$1
    shift
  fi

  if [ "$1" = "-m" ]
  then
    shift
    MAX_ALIGN_NUM=$1
    shift
  fi

  if [ "$1" = "-C" ]
  then
    shift
    COMMON_PREFIX="$1"
    shift
  fi

  if [ "$1" = "-S" ]
  then
    shift
    SENSITIVITY_OPTION="--sensitive"
  fi

  # For backwards compatibility we accept the option -j (join files) and -clean_up
  if [ "$1" = "-j" -o "$1" = "-clean_up" ]
  then
    shift
  fi

  if [ "$1" = "-sm" ]
  then
    shift
    SCORE_MIN_OPTION="--score-min L,0,$1"
    shift
  fi

  if [ "$1" = "-sort" ]
  then
    shift
    SORT="TRUE"
  fi

  if [ "$1" = "-I" ]
  then
    shift
    MIN_FRAGMENT_LENGTH_OPTION="-I $1"
    shift
  fi

  if [ "$1" = "-norc" ]
  then
    shift
    STRAND_SPECIFIC_OPTION="--norc"
  fi

  if [ "$1" = "-nofw" ]
  then
    shift
    STRAND_SPECIFIC_OPTION="--nofw"
  fi

done


################################################################################
##
## Print help
##
################################################################################

if [ "$FIRST_FASTQ_FILE" = "" -o "$1" = "" -o "$2" = "" -o "$3" = "" ]
then
  echo "Usage: $PROG_NAME -1 <fastq file 1> [-2 <fastq file 2>] <options>"
  echo "          <project dir> <sample dir> <bowtie2 index> [<sam file suffix>]"
  echo
  echo "where <options> is"
  echo "  [-g] [-bam] [-p <num processors>] [-m <max num alignments>] [-nogzip]"
  echo
  echo
  echo "-1 STRING: fastq file with the first read of each pair"
  echo "-2 STRING: fastq file with the second read of each pair [optional]"
  echo "project dir: directory where the relevant files are stored"
  echo "sample dir: the directory into which to write the SAM files."
  echo "bowtie2 index: which index to use for bowtie2"
  echo "sam file suffix: suffix appended to the SAM file names (before the .sam"
  echo "  extension) to allow for several SAM files in one directory "
  echo
  echo "Options:"
  echo "-g:    align against genome (uses different bowtie2 options)"
  echo "-bam:  use samtools to merge paired-end and single read alignments"
  echo "-nogzip: do not compress SAM files with gzip"
  echo "-p INT: the number of processors to use for bowtie2 [1]"
  echo "-m INT: max. no of alignments reported by Bowtie2 (only used with option -g"
  echo "        when aligning against the genome) [100]"
  echo "-S STRING: subdir of sam-files dir for storing the alignment files"
  exit
fi


################################################################################
##
## Read arguments
##
################################################################################

PROJECT_DIR=$1
shift

SAMPLE_DIR=$1
shift

BOWTIE2_INDEX=$1
shift

SAM_FILE_SUFFIX=""
if [ "$1" != "" ]
then
  SAM_FILE_SUFFIX="$1"
  shift
fi

if [ "$SAM_FILE_SUFFIX" = "genome" ]
then
  ALIGN_GENOME=TRUE
elif [ "$SAM_FILE_SUFFIX" = "transcript" ]
then
  ALIGN_TRANSCRIPT=TRUE
elif [ "$SAM_FILE_SUFFIX" = "junction" ]
then
  ALIGN_JUNCTION=TRUE
elif [ "$SAM_FILE_SUFFIX" != "" ]
then
  echo "Unknown SAM file suffix $SAM_FILE_SUFFIX ... exiting."
  exit 1
elif [ "$ALIGN_GENOME" = "TRUE" ]
then
  SAM_FILE_SUFFIX="genome"
elif [ "$ALIGN_TRANSCRIPT" = "TRUE" ]
then
  SAM_FILE_SUFFIX="transcript"
elif [ "$ALIGN_JUNCTION" = "TRUE" ]
then
  SAM_FILE_SUFFIX="junction" 
fi

###### Set MAX_ALIGN_NUM ###########################
if [ "$MAX_ALIGN_NUM" = "" ]
then
  if [ "$ALIGN_GENOME" = "TRUE" ]
  then
    MAX_ALIGN_NUM="100"
  elif [ "$ALIGN_TRANSCRIPT" = "TRUE" ]
  then
    MAX_ALIGN_NUM="2000"
  elif [ "$ALIGN_JUNCTION" = "TRUE" ]
  then
    MAX_ALIGN_NUM="400"
  fi
fi

################################################################################
##
## Set and create directories
##
################################################################################

EXON_PIPELINE_DIR=$PROJECT_DIR/exon-pipeline-files
EXON_PIPELINE_SCRIPTS_DIR=$PROJECT_DIR/exon-pipeline-scripts

FASTA_DIR=$EXON_PIPELINE_DIR/fasta-files
BED_DIR=$EXON_PIPELINE_DIR/bed-files

export TOOLS_DIR=$EXON_PIPELINE_SCRIPTS_DIR/tools
checkTool $TOOLS_DIR samtools 0.1.17
SAMTOOLS_EXE=$TOOL_EXE
checkTool $TOOLS_DIR bowtie2 2.0.5
BOWTIE2_EXE=$TOOL_EXE

export PICARD_JAR_DIR=$EXON_PIPELINE_SCRIPTS_DIR/tools/picard

JAVA_CLASS_DIR=$EXON_PIPELINE_SCRIPTS_DIR/java/classes
JAVA="java -oss8M -ss8M -ms16G -mx16G -cp $HOME/java:$HOME/java/classes:${JAVA_CLASS_DIR}:."

if [ ! -d $SAMPLE_DIR ]
then
  mkdir $SAMPLE_DIR
fi

if [ ! -d $SAMPLE_DIR/sam-files ]
then
  mkdir $SAMPLE_DIR/sam-files
fi


################################################################################
##
## Set fastq file names
##
################################################################################

FASTQ_FILE_BASE=`basename $FIRST_FASTQ_FILE | sed -e 's/[.]gz$//' | sed -e 's/[.]fq$//' | sed -e 's/[.]fastq$//' | sed -e 's/[-_]1$//'`
FASTQ_CHUNK=`echo $FASTQ_FILE_BASE | sed -e "s/.*-\(C[0-9][0-9][0-9]*\)/\1/"`

FASTQ_GZIPPED=`echo $FIRST_FASTQ_FILE | grep '.gz$'`
if [[ ! -z "$FASTQ_GZIPPED" ]]
then
  if [ ! -f $FIRST_FASTQ_FILE ]
  then
    FIRST_FASTQ_FILE=`echo $FIRST_FASTQ_FILE | sed -e 's/.gz$//'`    
  fi

  if [ ! -f $SECOND_FASTQ_FILE ]
  then
    SECOND_FASTQ_FILE=`echo $SECOND_FASTQ_FILE | sed -e 's/.gz$//'`
  fi

fi

if [ ! -f "$FIRST_FASTQ_FILE" ]
then
  echo "Cannot find the files $FIRST_FASTQ_FILE ... exiting"
  exit 1
fi

if [ "$PAIRED_END" = "TRUE" ]
then
  if [ ! -f "$SECOND_FASTQ_FILE" ]
  then
    echo "Cannot find the files $SECOND_FASTQ_FILE ... exiting"
    exit 1
  fi
fi


################################################################################
##
## Set bowtie2 options
##
################################################################################

# Align exon reads
# -p: number of threads/processes used
# -reorder: output the reads in the same order as in the input file(s)
# -k: Max. number of reported alignments
# $SENSITIVITY_OPTION: -D 15 -R 2 -L 22 -i S,1,1.15
# $SCORE_MIN_OPTION: --score-min L,0,-0.1
# where
#  -D <int>  Up to <int> consecutive seed extension attempts can "fail" before Bowtie 2 moves on, using the alignments found so far
#  -R <int>  <int> is the maximum number of times Bowtie 2 will "re-seed" reads with repetitive seeds.
#  -L <int>  Sets the length of the seed substrings to align during multiseed alignment. Smaller values make alignment slower but more senstive.
#  -i <func> sets a function governing the interval between seed substrings to use during multiseed alignment. (L = linear, S = square root)
BOWTIE2_OPTIONS="-p $NUM_PROCESSORS --reorder -k $MAX_ALIGN_NUM --phred33 $SENSITIVITY_OPTION $SCORE_MIN_OPTION"

if [ "$ALIGN_GENOME" = "FALSE" ]
then
  # Align exon reads against the transcripts
  # $MIN_FRAGMENT_LENGTH_OPTION: -I <min fragment length> - can be used to force single read alignments for junctions
  # $STRAND_SPECIFIC_OPTION: --norc -> bowtie2 does not align reads to the reverse complement (for single strand specific processing of samples)
  #                          --nofw -> bowtie2 does not align reads to the forward strand (for single strand specific processing of samples)
  #                          Note: this does not make sense for genomic alignments
  BOWTIE2_OPTIONS="$BOWTIE2_OPTIONS $MIN_FRAGMENT_LENGTH_OPTION $STRAND_SPECIFIC_OPTION"
fi

# Only output the header for the first chunk of transcript and junction files. We need it for genome files. This gets in the way
# when we try to sort SAM files with Picard
if [ "$FASTQ_CHUNK" != "$FASTQ_FILE_BASE" -a "$FASTQ_CHUNK" != "C001" -a "$ALIGN_GENOME" = "FALSE" ]
then
  BOWTIE2_OPTIONS="$BOWTIE2_OPTIONS --sam-no-hd"
fi


################################################################################
##
## Call bowtie2
##
################################################################################

SAM_DIR=$SAMPLE_DIR/sam-files
mkdir -p $SAM_DIR

if [ "$PAIRED_END" = "FALSE" ]
then

  SAM_FILE=$SAM_DIR/${FASTQ_FILE_BASE}-bowtie2-${SAM_FILE_SUFFIX}-sr.sam
  echo "Aligning $FIRST_FASTQ_FILE against $BOWTIE2_INDEX giving $SAM_FILE"
  
  CMD="$BOWTIE2_EXE $BOWTIE2_OPTIONS -x $BOWTIE2_INDEX -U $FIRST_FASTQ_FILE"

else

  SAM_FILE=$SAM_DIR/${FASTQ_FILE_BASE}-bowtie2-${SAM_FILE_SUFFIX}-pe.sam
  echo "Aligning $FIRST_FASTQ_FILE and $SECOND_FASTQ_FILE against $BOWTIE2_INDEX giving $SAM_FILE"
  
  # -X: max. insert size
  BOWTIE2_OPTIONS="$BOWTIE2_OPTIONS -X 500"
  CMD="$BOWTIE2_EXE $BOWTIE2_OPTIONS -x $BOWTIE2_INDEX -1 $FIRST_FASTQ_FILE -2 $SECOND_FASTQ_FILE"

fi

BOWTIE_VERSION=`$BOWTIE2_EXE --version | grep "bowtie.*version" | sed -e 's/.* \([0-9a-zA-Z.]*\)$/\1/'`

echo "ALIGNMENT_PROGRAM=bowtie2"
echo "ALIGNMENT_VERSION=$BOWTIE_VERSION"
echo "ALIGNMENT_OPTIONS=$BOWTIE2_OPTIONS"

## Run Bowtie2
date

echo "Bowtie2 call: $CMD | gzip > $SAM_FILE.gz"
$CMD | gzip > $SAM_FILE.gz
if [ $? -ne 0 ]
then
  echo "Bowtie2 alignment failed ... exiting"
  exit 1
fi

# echo "Compressing $SAM_FILE"
# if [[ -z "$COMMON_PREFIX" ]]
# then
#   gzip $SAM_FILE
# else
#   sed -e 's/$COMMON_PREFIX//' $SAM_FILE | gzip > $SAM_FILE.gz
# fi

# if [ $? -ne 0 ]
# then
#   echo "Compression failed ... exiting"
#   exit 1
# fi

echo "Alignment done."
date


################################################################################
##
## Create BAM file
##
################################################################################

if [ "$CREATE_BAM_FILE" = "TRUE" ]
then

  BAM_FILE=`echo $SAM_FILE | sed -e 's/[.]sam$/.bam/'`

  echo "Creating paired-end BAM file $BAM_FILE"
  gunzip $SAM_FILE.gz
  $SAMTOOLS_EXE view -bSh $SAM_FILE -o $BAM_FILE
  gzip $SAM_FILE
  echo "BAM file created."
  date

fi


################################################################################
##
## Sort SAM file by name
##
################################################################################

if [ "$SORT" = "TRUE" ]
then
  SAM_SORT_INPUT_FILE=`echo $SAM_FILE | sed -e 's/[.]sam$/-unsorted.sam/'`
  mv $SAM_FILE.gz $SAM_SORT_INPUT_FILE
  zcat $SAM_SORT_INPUT_FILE | $SAMTOOLS_EXE view -Sb - | $SAMTOOLS_EXE sort -m 3000000000 - $SAM_FILE-sorted
  $SAMTOOLS_EXE view -h $SAM_FILE-sorted.bam | gzip > $SAM_FILE.gz
  rm $SAM_SORT_INPUT_FILE $SAM_FILE-sorted.bam
  date
fi


################################################################################
##
## Decompress SAM file
##
################################################################################

if [ "$GZIP" = "FALSE" ]
then
  echo "Uncompressing $SAM_FILE.gz"
  gunzip -f $SAM_FILE.gz
  date
fi


################################################################################
##
## Compress FASTQ files if necessary
##
################################################################################

if [[ ! -z "$FASTQ_GZIPPED" ]]
then
  FASTQ_GZIPPED=`echo $FIRST_FASTQ_FILE | grep '.gz$'`
  if [[ -z "$FASTQ_GZIPPED" ]]
  then
    echo "Compressing $FIRST_FASTQ_FILE"
    gzip $FIRST_FASTQ_FILE &
  fi

  if [ "$FASTQ_GZIPPED" = "" -a "$PAIRED_END" = "TRUE" ]
  then
    echo "Compressing $SECOND_FASTQ_FILE"
    gzip $SECOND_FASTQ_FILE &
  fi
  wait
fi

echo "Alignment of Fastq files successfully completed"
date

exit 0
