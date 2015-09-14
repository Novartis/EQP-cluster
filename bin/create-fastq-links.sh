#!/bin/sh

PROG_NAME=`basename $0`
PROG_DIR=`dirname $0`
VERSION=1.0

cd $PROG_DIR
PROG_DIR=`pwd`

PROJECT_DIR=`echo $PROG_DIR | sed -e "s;/bin;;"`
BIN_DIR="$PROJECT_DIR/exon-pipeline-scripts/bin"


################################################################################
##
## Print help
##
################################################################################

ANNOTATION_DIR=$PROJECT_DIR/annotation-files
if [ ! -d $ANNOTATION_DIR ]
then
  ANNOTATION_DIR=$PROJECT_DIR
fi
ANNOTATION_FILES=`ls -1 $ANNOTATION_DIR | grep 'Data_Release_info.*txt$' | tr "\n" " "`

FIRST_ARGUMENT="$1"
SINGLE_READ_OPTION=""
if [ "$FIRST_ARGUMENT" = "-sr" ]
then
  SINGLE_READ_OPTION="-S"
  shift
  FIRST_ARGUMENT="$1"
fi

if [ "$ANNOTATION_FILES" = "" -a "$FIRST_ARGUMENT" != "-d" ]
then
  FIRST_ARGUMENT="-h"
fi

if [ "$FIRST_ARGUMENT" = "-h" -o "$FIRST_ARGUMENT" = "" ]
then
  if [ "$ANNOTATION_FILES" != "" ]
  then
    echo "Usage: $PROG_NAME [-sr] [-d] [<annotation files>]"
    echo
    echo "annotation files: the files containing the release information for the project"
    echo
    echo "-d: Use the default annotation file(s):"
    echo $ANNOTATION_FILES | tr " " "\n"
    echo
  else
    echo "Usage: $PROG_NAME [-sr] [-d] <annotation files>"
    echo
    echo "annotation files: the files containing the release information for the project"
    echo
    echo "(There are no default files of the format .*Data_Release_info.* in:"
    echo "$ANNOTATION_DIR)"
    echo
  fi    
  exit
fi


################################################################################
##
## Create samples directory if necessary
##
################################################################################

SAMPLES_DIR=$PROJECT_DIR/samples
if [ ! -d $SAMPLES_DIR ]
then
  mkdir $SAMPLES_DIR
fi


################################################################################
##
## Get the annotation file
##
################################################################################

if [ "$FIRST_ARGUMENT" != "-d" ]
then
  ANNOTATION_FILES=
  while [ "$1" != "" ]
  do
    ANNOTATION_FILES="$ANNOTATION_FILES $1"
    shift
  done
fi

if [ "$ANNOTATION_FILES" = "" ]
then
  echo "No Data_Release_info files found ... exiting."
  exit 1
fi


ANNOTATION_FILES_FULL=""
for ANNOTATION_FILE in $ANNOTATION_FILES
do
  if [ -f $ANNOTATION_DIR/$ANNOTATION_FILE ]
  then
    ANNOTATION_FILES_FULL="$ANNOTATION_FILES_FULL $ANNOTATION_DIR/$ANNOTATION_FILE"
  elif [ -f $ANNOTATION_FILE ]
  then
    ANNOTATION_FILES_FULL="$ANNOTATION_FILES_FULL $ANNOTATION_FILE"
  else
    echo "Neither file $ANNOTATION_FILE nor $ANNOTATION_DIR/$ANNOTATION_FILE not found ... exiting."
    exit 1
  fi
done


################################################################################
##
## Create the file links
##
################################################################################

$BIN_DIR/util-lib/createFileLinks.py $SINGLE_READ_OPTION -a $ANNOTATION_FILES_FULL -s $SAMPLES_DIR
