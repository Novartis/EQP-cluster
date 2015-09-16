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
    PID_EXIT_STATUS="$?"
  fi
  
  if [ "$LOG_FILE" != "" ]
  then
    if [ -f "$LOG_FILE" ]
    then
      if [ "$PID_EXIT_STATUS" != 0 ]
      then
        cat $LOG_FILE
      fi
      rm  $LOG_FILE
    else
      echo "WARNING: File $LOG_FILE not found."
    fi
  fi

  if [ "$PID_EXIT_STATUS" != "0" ]
  then
    echo "ERROR: Problem with $CMD"
    echo "Exit status: $PID_EXIT_STATUS ... exiting."
    if [ -e $OUTPUT_FILE ]
    then
      if [ -d $OUTPUT_FILE ]
      then
        rm $OUTPUT_FILE/*
      else
        rm $OUTPUT_FILE
      fi
    fi
    EXIT_STATUS=1
  else
    NUM_LINES=0
    if [ -f "$OUTPUT_FILE" ]
    then
      NUM_LINES=`cat $OUTPUT_FILE | head -100 | wc -l`
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

  if [ "$VERSION2_SECOND" -gt "$VERSION1_SECOND" ]
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
  TOOL_NAME="$1"
  TOOL_VERSION="$2"
  
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
    elif [ "$TOOL_NAME" = "bowtie2" -o "$TOOL_NAME" = "bowtie2-build" ]
    then
      ACTUAL_TOOL_VERSION=`$TOOL_EXE --version | fgrep version | fgrep -v "Compiler" | sed -e "s/.* version \([0-9]*[.][0-9]*[.][0-9]*\).*/\1/"`
    else
      echo "Unknown tool $TOOL_NAME"
      return
    fi
    ## echo "ACTUAL_TOOL_VERSION: $ACTUAL_TOOL_VERSION"
    checkVersion $TOOL_NAME $ACTUAL_TOOL_VERSION $TOOL_VERSION
  fi
}


################################################################################
##
##
## Main program
##
##
################################################################################

PROG_DIR=`dirname $0`
PROG_NAME=`basename $0`
VERSION=2.0

cd $PROG_DIR
PROG_DIR=`pwd`

PROJECT_DIR=`echo $PROG_DIR | sed -e 's;/bin[/]*$;;' | sed -e 's;/exon-pipeline-files[/]*$;;'`
# echo "Using project directory: $PROJECT_DIR"

cd $PROJECT_DIR

PROJECT_DATA_DIR="$PROJECT_DIR/exon-pipeline-files"
if [ ! -d $PROJECT_DATA_DIR ]
then
 echo "Directory "
 echo "  $PROJECT_DIR/exon-pipeline-files"
 echo "not found. Please make sure that the download from GitHub is complete and no"
 echo "file or directory is moved ... exiting."
 exit 1
fi

cd exon-pipeline-files
if [ -f exon-pipeline-files.tgz.00 ]
then
  echo "Extracting directory exon-pipeline-files"
  cat exon-pipeline-files.tgz.* | tar xz
  if [ $? -ne 0 ]
  then
    echo "Extraction of directory exon-pipeline-files failed ... exiting."
    exit 1
  fi
  rm exon-pipeline-files.tgz.*
else
  echo "File"
  echo "  exon-pipeline-files.tgz.00"
  echo "found. Please make sure that the download from GitHub is complete and no"
  echo "file or directory is moved ... exiting."
  exit 1
fi

cd ..

FORCE="TRUE"

################################################################################
##
##  Read relevant project information
##
################################################################################

SETUP_FILE=$PROJECT_DATA_DIR/setup.sh
if [ -f $SETUP_FILE ]
then
  source $SETUP_FILE
  ## echo "Copying file $SETUP_FILE for read length $READ_LENGTH to $PROJECT_DIR/bin"
  cp -f $SETUP_FILE $PROJECT_DIR/bin
else
  echo "File $SETUP_FILE not found ... exiting."
  exit 1  
fi


################################################################################
##
##  Checking the versions of samtools and bowtie2
##
################################################################################

checkTool samtools 0.1.17
SAMTOOLS_EXE=$TOOL_EXE
checkTool bowtie2-build 2.0.5
BOWTIE2_BUILD_EXE=$TOOL_EXE


################################################################################
##
##  Check transcript and junction Fasta file
##
################################################################################

FASTA_DIR="$PROJECT_DATA_DIR/fasta-files"
if [ ! -d $FASTA_DIR ]
then
  echo "Subdirectory fasta-files of directory"
  echo "$PROJECT_DATA_DIR not found."
  echo "Please make sure all data are extracted correctly ... exiting"
  exit 1
fi

FILE_BASE=`cat $PROJECT_DIR/bin/setup.sh | fgrep "FILE_BASE" | sed -e 's/export FILE_BASE=//'`
RADIUS=`cat $PROJECT_DIR/bin/setup.sh | fgrep "RADIUS" | sed -e 's/export RADIUS=//'`

JUNCTION_FASTA_FILE=$FASTA_DIR/${FILE_BASE}-${RADIUS}-junction.fasta
if [ ! -f $JUNCTION_FASTA_FILE ]
then
  echo "File ${FILE_BASE}-${RADIUS}-junction.fasta in directory"
  echo "$FASTA_DIR not found."
  echo "Please make sure all data are extracted correctly ... exiting"
  exit 1
fi


################################################################################
##
##  Download genome Fasta file
##
################################################################################

GENOME_DIR="$PROJECT_DATA_DIR/genome-files"
if [ ! -d $GENOME_DIR ]
then
  mkdir $GENOME_DIR
fi
cd $GENOME_DIR

GENOME_FASTA_FILE="$GENOME_DIR/genome.fa"
if [ ! -f $GENOME_FASTA_FILE -a ! -f $GENOME_FASTA_FILE.gz ]
then
  if [ -f $GENOME_DIR/genome.url ]
  then
    GENOME_FASTA_URL=`cat $GENOME_DIR/genome.url`
    
    DOWNLOAD_FILE=`basename $GENOME_FASTA_URL`
    DECOMPRESSED_DOWNLOAD_FILE=`basename $GENOME_FASTA_URL .gz`
    if [ ! -f $DOWNLOAD_FILE -a ! -f $DECOMPRESSED_DOWNLOAD_FILE ]
    then
      echo "Download of genome Fasta file from URL:"
      echo "  $GENOME_FASTA_URL"
      wget $GENOME_FASTA_URL
      if [ $? -ne 0 ]
      then
        echo "Download of genome Fasta file from URL failed. Please download the file"
        echo "manually and move it to the directory:"
        echo "  $GENOME_DIR"
        echo "and try again ... exiting"
        exit 1
      fi
    elif [ -f $DOWNLOAD_FILE ]
    then
      echo "Using file $GENOME_DIR/$DOWNLOAD_FILE as genome Fasta file"
    else
      echo "Using file $GENOME_DIR/$DECOMPRESSED_DOWNLOAD_FILE as genome Fasta file"
    fi
  
    if [ ! -f "$DECOMPRESSED_DOWNLOAD_FILE" -a -f "$DOWNLOAD_FILE" ]
    then
      echo "Decompressing the genome Fasta file $DOWNLOAD_FILE"
      gunzip $DOWNLOAD_FILE
      if [ $? -ne 0 ]
      then
        echo "gunzip of file $DOWNLOAD_FILE failed. Maybe,"
        echo "the download was not complete. Please remove the file, download it"
        echo "manually, and move it to the directory:"
        echo "  $GENOME_DIR"
        echo "and try again ...... exiting"
        exit 1
      fi
    fi
  
    mv $DECOMPRESSED_DOWNLOAD_FILE $GENOME_FASTA_FILE
  else
    echo "Neither file genome.fa nor genome.url found in directory"
    echo "  $GENOME_DIR"
    echo "found. Please make sure that all files are unpacked ... exiting."
    exit 1
  fi
elif [ ! -f $GENOME_FASTA_FILE ]
then
  echo "Decompressing file $GENOME_FASTA_FILE.gz"
  gunzip $GENOME_FASTA_FILE.gz
# else
#   echo "Using $GENOME_FASTA_FILE as genome Fasta file"
fi
  
if [ ! -f $GENOME_FASTA_FILE.fai ]
then
  echo "Creating the Fasta index file for the genome Fasta file"
  $SAMTOOLS_EXE faidx $GENOME_FASTA_FILE
  if [ $? -ne 0 ]
  then
    echo "Could not create Fasta index file for file"
    echo "  $GENOME_FASTA_FILE"
    echo "with samtools ... exiting"
    exit 1
  fi
fi


################################################################################
##
##  Download transcript Fasta file
##
################################################################################

cd $FASTA_DIR

TRANSCRIPT_FASTA_FILE=$FASTA_DIR/${FILE_BASE}.fasta
if [ ! -f $TRANSCRIPT_FASTA_FILE -a ! -f $TRANSCRIPT_FASTA_FILE.gz ]
then
  if [ -f $FASTA_DIR/transcript.url ]
  then
    TRANSCRIPT_FASTA_URL=`cat $FASTA_DIR/${FILE_BASE}.url`
    
    DOWNLOAD_FILE=`basename $TRANSCRIPT_FASTA_URL`
    DECOMPRESSED_DOWNLOAD_FILE=`basename $TRANSCRIPT_FASTA_URL .gz`
    if [ ! -f $DOWNLOAD_FILE -a ! -f $DECOMPRESSED_DOWNLOAD_FILE ]
    then
      echo "Download of transcript Fasta file from URL:"
      echo "  $TRANSCRIPT_FASTA_URL"
      wget $TRANSCRIPT_FASTA_URL
      if [ $? -ne 0 ]
      then
        echo "Download of transcript Fasta file from URL failed. Please download the file"
        echo "manually and move it to the directory:"
        echo "  $FASTA_DIR"
        echo "and try again ... exiting"
        exit 1
      fi
    elif [ -f $DOWNLOAD_FILE ]
    then
      echo "Using file $FASTA_DIR/$DOWNLOAD_FILE as genome Fasta file"
    else
      echo "Using file $FASTA_DIR/$DECOMPRESSED_DOWNLOAD_FILE as genome Fasta file"
    fi
  
    if [ ! -f "$DECOMPRESSED_DOWNLOAD_FILE" -a -f "$DOWNLOAD_FILE" ]
    then
      echo "Decompressing the transcript Fasta file $DOWNLOAD_FILE"
      gunzip $DOWNLOAD_FILE
      if [ $? -ne 0 ]
      then
        echo "gunzip of file $DOWNLOAD_FILE failed. Maybe,"
        echo "the download was not complete. Please remove the file, download and"
        echo "gunzip it manually, and move it to the directory:"
        echo "  $FASTA_DIR"
        echo "and try again ...... exiting"
        exit 1
      fi
    fi
  
    mv $DECOMPRESSED_DOWNLOAD_FILE $TRANSCRIPT_FASTA_FILE
  else
    echo "Neither file ${FILE_BASE}.fasta nor ${FILE_BASE}.url found in directory"
    echo "  $FASTA_DIR"
    echo "found. Please make sure that all files are unpacked ... exiting."
    exit 1
  fi
fi


################################################################################
##
##  Download GTF file
##
################################################################################

GTF_DIR="$PROJECT_DATA_DIR/gtf-files"
if [ ! -d $GTF_DIR ]
then
  mkdir $GTF_DIR
fi
cd $GTF_DIR

GTF_FILE=$GTF_DIR/${FILE_BASE}.gtf
if [ ! -f $GTF_FILE ]
then
  if [ -f $GTF_DIR/gtf.url ]
  then
    GTF_URL=`cat $GTF_DIR/gtf.url`
    
    DOWNLOAD_FILE=`basename $GTF_URL`
    DECOMPRESSED_DOWNLOAD_FILE=`basename $GTF_URL .gz`
    if [ ! -f $DOWNLOAD_FILE -a ! -f $DECOMPRESSED_DOWNLOAD_FILE ]
    then
      echo "Download of GTF file from URL:"
      echo "  $GTF_URL"
      wget $GTF_URL
      if [ $? -ne 0 ]
      then
        echo "Download of GTF file from URL failed. Please download the file"
        echo "manually and move it to the directory:"
        echo "  $GTF_DIR"
        echo "and try again ... exiting"
        exit 1
      fi
    elif [ -f $DOWNLOAD_FILE ]
    then
      echo "Using file $GTF_DIR/$DOWNLOAD_FILE as GTF file"
    else
      echo "Using file $GTF_DIR/$DECOMPRESSED_DOWNLOAD_FILE as GTF file"
    fi
  
    if [ ! -f "$DECOMPRESSED_DOWNLOAD_FILE" -a -f "$DOWNLOAD_FILE" ]
    then
      echo "Decompressing the GTF file $DOWNLOAD_FILE"
      gunzip $DOWNLOAD_FILE
      if [ $? -ne 0 ]
      then
        echo "gunzip of file $DOWNLOAD_FILE failed. Maybe,"
        echo "the download was not complete. Please remove the file, download it"
        echo "manually, and move it to the directory:"
        echo "  $GTF_DIR"
        echo "and try again ...... exiting"
        exit 1
      fi
    fi
  
    mv $DECOMPRESSED_DOWNLOAD_FILE $GTF_FILE
  else
    echo "Neither file ${FILE_BASE}.gtf nor ${FILE_BASE}.url found in directory"
    echo "  $GTF_DIR"
    echo "found. Please make sure that all files are unpacked ... exiting."
    exit 1
  fi
elif [ ! -f $GTF_FILE ]
then
  echo "Decompressing file $GTF_FILE.gz"
  gunzip $GTF_FILE.gz
# else
#   echo "Using $GTF_FILE as genome Fasta file"
fi


################################################################################
##
##  Create Bowtie2 index for transcript, genome, and junction Fasta file
##
################################################################################

##Make Bowtie2 index for transcript Fasta file
TRANSCRIPT_BOWTIE2_INDEX_DIR=$FASTA_DIR/${FILE_BASE}-bowtie2Index
if [ ! -d "$TRANSCRIPT_BOWTIE2_INDEX_DIR" -o "$FORCE" = "TRUE" ]
then
  if [ ! -d "$TRANSCRIPT_BOWTIE2_INDEX_DIR" ]
  then
    mkdir $TRANSCRIPT_BOWTIE2_INDEX_DIR
  fi
  
  echo "Starting to build the Bowtie2 index for the Transcript Fasta file"
  $BOWTIE2_BUILD_EXE $TRANSCRIPT_FASTA_FILE $TRANSCRIPT_BOWTIE2_INDEX_DIR/bowtie2Index > $TRANSCRIPT_BOWTIE2_INDEX_DIR/bowtie2-build.log 2>&1 &
  CREATE_TRANSCRIPT_BOWTIE2_INDEX_PID=$!

fi

##Make Bowtie2 index for junction Fasta file
JUNCTION_BOWTIE2_INDEX_DIR=$FASTA_DIR/${FILE_BASE}-${RADIUS}-junction-bowtie2Index
if [ ! -d "$JUNCTION_BOWTIE2_INDEX_DIR" -o "$FORCE" = "TRUE" ]
then
  if [ ! -d "$JUNCTION_BOWTIE2_INDEX_DIR" ]
  then
    mkdir $JUNCTION_BOWTIE2_INDEX_DIR
  fi
  
  echo "Starting to build the Bowtie2 index for the Junction Fasta file"
  $BOWTIE2_BUILD_EXE $JUNCTION_FASTA_FILE $JUNCTION_BOWTIE2_INDEX_DIR/bowtie2Index > $JUNCTION_BOWTIE2_INDEX_DIR/bowtie2-build.log 2>&1 &
  CREATE_JUNCTION_BOWTIE2_INDEX_PID=$!
fi

##Make Bowtie2 index for genome Fasta file
GENOME_BOWTIE2_INDEX_DIR="$GENOME_DIR/bowtie2Index"
if [ ! -d "$GENOME_BOWTIE2_INDEX_DIR" -o "$FORCE" = "TRUE" ]
then
  if [ ! -d "$GENOME_BOWTIE2_INDEX_DIR" ]
  then
    mkdir $GENOME_BOWTIE2_INDEX_DIR

    cd $GENOME_BOWTIE2_INDEX_DIR
    ln -s $GENOME_FASTA_FILE
  fi
  
  echo "Starting to build the Bowtie2 index for the Genome Fasta file"
  $BOWTIE2_BUILD_EXE $GENOME_FASTA_FILE $GENOME_BOWTIE2_INDEX_DIR/bowtie2Index > $GENOME_BOWTIE2_INDEX_DIR/bowtie2-build.log 2>&1 &
  CREATE_GENOME_BOWTIE2_INDEX_PID=$!
fi


################################################################################
##
##  Wait for processes to finish and check exit status
##
################################################################################

EXIT_STATUS=0
if [ "$CREATE_TRANSCRIPT_BOWTIE2_INDEX_PID" != "" ]
then
  waitPid $CREATE_TRANSCRIPT_BOWTIE2_INDEX_PID $TRANSCRIPT_BOWTIE2_INDEX_DIR "$BOWTIE2_BUILD_EXE $TRANSCRIPT_FASTA_FILE $TRANSCRIPT_BOWTIE2_INDEX_DIR/bowtie2Index" \
    "Waiting for creation of Bowtie2 index for transcript Fasta file to finish" $TRANSCRIPT_BOWTIE2_INDEX_DIR/bowtie2-build.log
fi

if [ "$CREATE_GENOME_BOWTIE2_INDEX_PID" != "" ]
then
  waitPid $CREATE_GENOME_BOWTIE2_INDEX_PID $GENOME_BOWTIE2_INDEX_DIR "$BOWTIE2_BUILD_EXE $GENOME_FASTA_FILE $GENOME_BOWTIE2_INDEX_DIR/bowtie2Index" \
    "Waiting for creation of Bowtie2 index for genome Fasta file to finish" $GENOME_BOWTIE2_INDEX_DIR/bowtie2-build.log
fi

if [ "$CREATE_JUNCTION_BOWTIE2_INDEX_PID" != "" ]
then
  waitPid $CREATE_JUNCTION_BOWTIE2_INDEX_PID $JUNCTION_BOWTIE2_INDEX_DIR "$BOWTIE2_BUILD_EXE $JUNCTION_FASTA_FILE $JUNCTION_BOWTIE2_INDEX_DIR/bowtie2Index" \
    "Waiting for creation of Bowtie2 index for junction Fasta file to finish" $JUNCTION_BOWTIE2_INDEX_DIR/bowtie2-build.log
fi

if [ "$EXIT_STATUS" != "0" ]
then
  echo "Project setup failed."
  exit $EXIT_STATUS
fi

echo "Project setup complete"


