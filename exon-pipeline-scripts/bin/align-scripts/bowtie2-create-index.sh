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

PROG_NAME=`basename $0`
PROG_DIR=`dirname $0`
VERSION=1.0


################################################################################
##
## Read options
##
################################################################################

PRINT_HELP="FALSE"
while [ "$1" = "-h" ]
do
  if [ "$1" = "-h" ]
  then
    shift
    PRINT_HELP="TRUE"
  fi
done


################################################################################
##
## Print help
##
################################################################################

if [ "$3" = "" ]
then
  echo "Usage: $PROG_NAME [-h] <project dir> <Fasta file> <Bowtie2 index>"
  echo
  echo "project dir: directory where the relevant files are stored"
  echo "Fasta file: The Fasta file for which to create the Bowtie2 index."
  echo "Bowtie2 index: The Bowtie2 index"
  echo
  echo "Options:"
  echo "-h: print this message"
  exit
fi


################################################################################
##
## Read arguments
##
################################################################################

PROJECT_DIR=$1
shift

if [ ! -d $PROJECT_DIR ]
then
  echo "Directory $PROJECT_DIR does not exist ... exiting."
  exit 1
fi


BOWTIE2_FASTA_FILE="$1"
shift

if [ ! -f "$BOWTIE2_FASTA_FILE" ]
then
  echo "Fasta file $BOWTIE2_FASTA_FILE for the creation"
  echo "of the Bowtie2 index not found ... exiting"
  exit 1
fi

BOWTIE2_INDEX="$1"
shift

BOWTIE2_INDEX_DIR=`dirname $BOWTIE2_INDEX_DIR`
if [ ! -d "$BOWTIE2_INDEX_DIR" ]
then
  echo "Creating directory $BOWTIE2_INDEX_DIR"
  mkdir -p $BOWTIE2_INDEX_DIR
fi

BOWTIE2_INDEX_BASE=`basename $BOWTIE2_INDEX`
BOWTIE2_INDEX_FILE_NUM=`ls -1 $BOWTIE2_INDEX_DIR | fgrep $BOWTIE2_INDEX_BASE | wc -l`

if [ $BOWTIE2_INDEX_FILE_NUM -ge 6 ]
then
  return 0
fi


################################################################################
##
## Create the Bowtie2 index
##
################################################################################

echo "Bowtie2 index for $BOWTIE2_FASTA_FILE not found ... creating."

BOWTIE2_BUILD_EXE="$PROJECT_DIR/exon-pipeline-scripts/tools/bowtie2-build"

if [ ! -f $BOWTIE2_BUILD_EXE ]
then
  echo "Executable $PROJECT_TOOLS_DIR/bowtie2-build not found."
  echo "Cannot create Bowtie2 index for Fasta file $BOWTIE2_FASTA_FILE."
  echo " ... exiting"
  exit 1
fi

BOWTIE2_INDEX_DATE=`date +%y%m%d-%H%M`

$BOWTIE2_BUILD_EXE $BOWTIE2_FASTA_FILE $BOWTIE2_INDEX
if [ $? -ne 0 ]
then
  echo "Bowtie2 index creation failed for $BOWTIE2_FASTA_FILE ... exiting."
  echo "See file $PROJECT_DIR/log-files/bowtie-index-$BOWTIE2_INDEX_DATE.log"
  exit 1
fi

echo "Bowtie2 index successfully computed."
exit 1

  
