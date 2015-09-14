#!/bin/ksh

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

PROG_NAME=`basename $0`
PROG_DIR=`dirname $0`
VERSION=2.0


############################################################################
##
##  Functions
##
############################################################################

cleanUp ()
{
  ######################  Remove combined files ##############################
  
  if [ "$COMBINED_FASTQ_FILE_1" != "" -a -f "$COMBINED_FASTQ_FILE_1" -a "$REMOVE" = "TRUE" ]
  then
    echo "Removing files $COMBINED_FASTQ_FILE_1"
    rm $COMBINED_FASTQ_FILE_1
  fi
  
  if [ "$COMBINED_FASTQ_FILE_2" != "" -a -f "$COMBINED_FASTQ_FILE_2" -a "$REMOVE" = "TRUE" ]
  then
    echo "Removing files $COMBINED_FASTQ_FILE_2"
    rm $COMBINED_FASTQ_FILE_2
  fi
  
  ######################  Remove unique files ##############################
  
  if [ "$UNIQUE_FASTQ_FILE_1" != "$FASTQ_FILE_1" -a -f "$UNIQUE_FASTQ_FILE_1" ]
  then
    echo "Remove unique file $UNIQUE_FASTQ_FILE_1"
    rm $UNIQUE_FASTQ_FILE_1
  fi
  
  if [ "$UNIQUE_FASTQ_FILE_2" != "$FASTQ_FILE_2" -a -f "$UNIQUE_FASTQ_FILE_2" ]
  then
    echo "Remove unique file $UNIQUE_FASTQ_FILE_2"
    rm $UNIQUE_FASTQ_FILE_2
  fi
  
}


############################################################################
##
##  Main program
##
############################################################################

############################# Read in options ##################################

FASTQ_FILE_1=
FASTQ_FILE_2=
SPLIT_READ_NUM=25000000
PAIRED_END="FALSE"
REMOVE_DUPLICATES="FALSE"
COMPUTE_DUPLICATES="FALSE"
TALLY_OPTIONS="--with-quality"
while [ "$1" = "-1" -o "$1" = "-2" -o "$1" = "-s" -o "$1" = "-r" -o "$1" = "-remove-duplicates" -o "$1" = "-d" -o "$1" = "-tally" ]
do
  if [ "$1" = "-1" ]
  then
    shift
    FASTQ_FILE_1=$1
    shift
  fi

  if [ "$1" = "-2" ]
  then
    shift
    FASTQ_FILE_2=$1
    shift
    PAIRED_END=TRUE
    TALLY_OPTIONS="$TALLY_OPTIONS --pair-by-offset"
  fi

  if [ "$1" = "-s" ]
  then
    shift
    SPLIT_READ_NUM=$1
    shift
  fi
  
  if [ "$1" = "-r" -o "$1" = "-remove-duplicates" ]
  then
    shift
    REMOVE_DUPLICATES="TRUE"
  fi

  if [ "$1" = "-d" ]
  then
    shift
    COMPUTE_DUPLICATES="TRUE"
  fi

  ## This is actually not needed as the standard record format is anyway
  ## @%I%#%R%n+%#%Q%n
  if [ "$1" = "-tally" ]
  then
    TALLY_OPTIONS="$TALLY_OPTIONS -record-format '@%I%b%X%#%R%n+%#%Q%n'"
  fi

done


################################ Usage #########################################

if [ "$FASTQ_FILE_1" = "" -o "$1" = "" ]
then
  echo "Usage: $PROG_NAME <options> -1 <fastq file 1> [-2 <fastq file 2>]"
  echo "       <project dir> <sample dir> <fastq dir>"
  echo
  echo "where <options> is:"
  echo "  [-s <split read num>] [-r|-remove-duplicates]"
  echo
  echo "-s INT: set the number of reads per file [$SPLIT_READ_NUM]"
  echo "-r|-remove-duplicates: remove duplicate reads with reaper/tally"
  exit
fi


############################ Read in arguments #################################

PROJECT_DIR=$1
shift

SAMPLE_DIR=$1
shift

FASTQ_DIR=$1
shift

if [ -d $FASTQ_DIR ]
then
  rm -r $FASTQ_DIR
fi
echo "Creating directory $FASTQ_DIR"
mkdir $FASTQ_DIR


cd $SAMPLE_DIR

GZIP_SUFFIX=`echo "$FASTQ_FILE_1" | grep '[.]gz[#]\?$'`
if [[ -z "$GZIP_SUFFIX" ]]
then
  CAT="cat"
  GZIP="cat"
else
  CAT="zcat"
  GZIP="gzip"
fi

echo "FASTQ_FILE_1=$FASTQ_FILE_1, GZIP_SUFFIX=$GZIP_SUFFIX, CAT=$CAT, GZIP=$GZIP"
WILDCARD=`echo "$FASTQ_FILE_1" | grep "[%#]"`

SAMPLE_NAME=`basename $SAMPLE_DIR`

COMBINED_FASTQ_FILE_1=
COMBINED_FASTQ_FILE_2=
if [[ ! -z $WILDCARD ]]
then
  COMBINED_FASTQ_FILE_1=$FASTQ_DIR/$SAMPLE_NAME-combined_1.fq.gz
  COMBINED_FASTQ_FILE_2=$FASTQ_DIR/$SAMPLE_NAME-combined_2.fq.gz

  echo "Creating file $COMBINED_FASTQ_FILE_1"
  date
  echo "$CAT `echo $FASTQ_FILE_1 | sed -e 's/%/*/' | tr '#' ' '` | $GZIP > $COMBINED_FASTQ_FILE_1 &"
  ($CAT `echo $FASTQ_FILE_1 | sed -e 's/%/*/' | tr '#' ' '` | $GZIP > $COMBINED_FASTQ_FILE_1; echo "$COMBINED_FASTQ_FILE_1 complete") &
  FASTQ_FILE_1=$COMBINED_FASTQ_FILE_1

  if [ "$PAIRED_END" = "TRUE" ]
  then
    echo "Creating file $COMBINED_FASTQ_FILE_2"
    echo "$CAT `echo $FASTQ_FILE_2 | sed -e 's/%/*/' | tr '#' ' '` | $GZIP > $COMBINED_FASTQ_FILE_2 &"
    ($CAT `echo $FASTQ_FILE_2 | sed -e 's/%/*/' | tr '#' ' '` | $GZIP > $COMBINED_FASTQ_FILE_2; echo "$COMBINED_FASTQ_FILE_2 complete") &
    FASTQ_FILE_2=$COMBINED_FASTQ_FILE_2
  fi
  
  REMOVE="TRUE"
  wait
fi


############################ Remove duplicate #################################

UNIQUE_FASTQ_FILE_1="$FASTQ_FILE_1"
UNIQUE_FASTQ_FILE_2="$FASTQ_FILE_2"
if [ "$REMOVE_DUPLICATES" = "TRUE" -o "$COMPUTE_DUPLICATES" = "TRUE" ]
then
  TOOLS_DIR=$PROJECT_DIR/exon-pipeline-scripts/tools
  
  NUM_LINES_FASTQ_FILE=`$CAT $FASTQ_FILE_1 | wc -l`
  echo "Number of reads before duplicate removal: " $(echo "$NUM_LINES_FASTQ_FILE / 4" | bc)

  echo "Removing duplicate reads"
  UNIQUE_FASTQ_FILE_1=`echo "$FASTQ_FILE_1" | sed -e 's/1[.]f[ast]*q/uniq-1.fastq/'`
  if [ "$PAIRED_END" = "TRUE" ]
  then
    UNIQUE_FASTQ_FILE_2=`echo "$FASTQ_FILE_2" | sed -e 's/2[.]f[ast]*q/uniq-2.fastq/'`
    echo "Tally command: $TOOLS_DIR/reaper/tally -i $FASTQ_FILE_1 -j $FASTQ_FILE_2 -o $UNIQUE_FASTQ_FILE_1 -p $UNIQUE_FASTQ_FILE_2 $TALLY_OPTIONS"
    $TOOLS_DIR/reaper/tally -i $FASTQ_FILE_1 -j $FASTQ_FILE_2 -o $UNIQUE_FASTQ_FILE_1 -p $UNIQUE_FASTQ_FILE_2 $TALLY_OPTIONS
    ls -l $UNIQUE_FASTQ_FILE_1 $UNIQUE_FASTQ_FILE_2
  else
    $TOOLS_DIR/reaper/tally -i $FASTQ_FILE_1 -o $UNIQUE_FASTQ_FILE_1 $TALLY_OPTIONS
  fi
  echo "Tally done"

  echo "Searching for double hyphen"
  echo "$CAT $UNIQUE_FASTQ_FILE_1 | grep '^--$'"
  DOUBLE_HYPHEN=`$CAT $UNIQUE_FASTQ_FILE_1 | grep '^--$'`
  if [ "$DOUBLE_HYPHEN" != "" ]
  then
    echo "Double hyphen found in $UNIQUE_FASTQ_FILE_1 ... exiting."
    exit 1
  fi

  echo "Computing the number of lines"
  echo "$CAT $UNIQUE_FASTQ_FILE_1 | wc -l"
  NUM_LINES_UNIQUE_FASTQ_FILE=`$CAT $UNIQUE_FASTQ_FILE_1 | wc -l`
  echo "Number of reads after duplicate removal: " $(echo "$NUM_LINES_UNIQUE_FASTQ_FILE / 4" | bc)

  if [ "$COMPUTE_DUPLICATES" = "TRUE" ]
  then
    FILTER_FILE="$SAMPLE_DIR/filter-duplicates.tmp"
    echo "Filter file: $FILTER_FILE"
    i=0
    while [ "$i" -lt 10 ]
    do
      if [ "$i" = "0" ]
      then
        echo " 1$i" > $FILTER_FILE
      else
        echo " 1$i" >> $FILTER_FILE
      fi
      
      if [ "$i" != "1" -a "$i" != "0" ]
      then
        echo " $i" >> $FILTER_FILE
      fi
      
      i=`echo $i + 1 | bc`
    done

    DUPLICATES_SAMPLE_DIR=`echo "$SAMPLE_DIR" | sed -e 's/samples-Rem-dup/samples-Dup/'`
    DUPLICATES_FASTQ_FILE_1=`echo "$FASTQ_FILE_1" | sed -e 's/1[.]f[ast]*q/dup-1.fastq/'`

    echo "Creating duplicates file for $UNIQUE_FASTQ_FILE_1"
    #  Remove group separator '--' placed by grep between consecutive hits for -A
    $CAT $UNIQUE_FASTQ_FILE_1 | fgrep -A 3 -f $FILTER_FILE | grep -v '^--$' | $GZIP > $DUPLICATES_SAMPLE_DIR/$DUPLICATES_FASTQ_FILE_1   # $DUPLICATES_ONLY_FASTQ_FILE_1

    DOUBLE_HYPHEN=`$CAT $DUPLICATES_SAMPLE_DIR/$DUPLICATES_FASTQ_FILE_1 | grep '^--$'`
    if [ "$DOUBLE_HYPHEN" != "" ]
    then
      echo "Double hyphen found in $DUPLICATES_SAMPLE_DIR/$DUPLICATES_FASTQ_FILE_1 ... exiting."
      exit 1
    fi
  
    if [ "$PAIRED_END" = "TRUE" ]
    then
      DUPLICATES_FASTQ_FILE_2=`echo "$FASTQ_FILE_2" | sed -e 's/2[.]f[ast]*q/dup-2.fastq/'`
      echo "Creating duplicates file for $UNIQUE_FASTQ_FILE_2"
      #  Remove group separator '--' placed by grep between consecutive hits for -A
      $CAT $UNIQUE_FASTQ_FILE_2 | fgrep -A 3 -f $FILTER_FILE | grep -v '^--$' | $GZIP > $DUPLICATES_SAMPLE_DIR/$DUPLICATES_FASTQ_FILE_2   # $DUPLICATES_ONLY_FASTQ_FILE_2
    fi

    # rm $FILTER_FILE

    if [ "$REMOVE_DUPLICATES" != "TRUE" ]
    then
      cleanUp
      date
      
      echo "Done."
      exit 0
    fi
  fi
fi


########################## Set java directories ################################

PROJECT_JAVA_DIR=$PROJECT_DIR/exon-pipeline-scripts/java

JAVA_CLASS_DIR="$PROJECT_JAVA_DIR/classes:$PROJECT_JAVA_DIR"
JAVA="java -oss8M -ss8M -ms2G -mx2G -cp ${CLASSPATH}:${JAVA_CLASS_DIR}"


#########################  Split fastq files ###################################

date

if [ "$PAIRED_END" = "TRUE" ]
then
  echo "Creating new fastq file(s) for $UNIQUE_FASTQ_FILE_1 and $UNIQUE_FASTQ_FILE_2"
  $JAVA ChangeFastqIdsAndSplit -s $SPLIT_READ_NUM -f $UNIQUE_FASTQ_FILE_1 -F $UNIQUE_FASTQ_FILE_2 -o $FASTQ_DIR/${SAMPLE_NAME}_1.fq.gz \
     -O $FASTQ_DIR/${SAMPLE_NAME}_2.fq.gz
else 
  echo "Creating new fastq file(s) for $UNIQUE_FASTQ_FILE_1"
  $JAVA ChangeFastqIdsAndSplit -s $SPLIT_READ_NUM -S -f $UNIQUE_FASTQ_FILE_1 -o $FASTQ_DIR/${SAMPLE_NAME}_1.fq.gz
fi

date

FASTQ_FILES=`ls $FASTQ_DIR | grep '.*C[0-9][0-9][0-9].*[.]f[ast]*q' | grep -v '[.]gz$'`
if [[ ! -z $FASTQ_FILES ]]
then
  for FASTQ_FILE in $FASTQ_FILES
  do
    gzip $FASTQ_DIR/$FASTQ_FILE &
  done
  wait
fi

######################  Create chunk directories ##############################

FASTQ_FILE_CHUNKS=`ls -1 $FASTQ_DIR | sed -e 's/[.]gz$//' | sed -e 's/[.]fq$//' | sed -e 's/[.]fastq$//' |  sed -e 's/[-_]1$//' | grep 'C[0-9][0-9][0-9]$' | tr "\n" " "`

for FASTQ_FILE_CHUNK in $FASTQ_FILE_CHUNKS
do

  CHUNK=`echo $FASTQ_FILE_CHUNK | sed -e 's/.*\(C[0-9][0-9][0-9]\)$/\1/'`
  CHUNK_DIR=$SAMPLE_DIR/$CHUNK

  if [ ! -d $CHUNK_DIR ]
  then
    echo "Making directory $CHUNK_DIR"
    mkdir $CHUNK_DIR
  fi

done

cleanUp

echo "Done"
date




