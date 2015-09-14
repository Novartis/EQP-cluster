#!/bin/ksh

################################################################################
##
## Function definitions
##
################################################################################

NUM_PROCESSES=0


################################################################################
##
## submitJob - uses upto seven arguments
##
##  1. mem_token, e.g. 2G
##  2. job name, e.g. SP-SEQC-A-BC1-s_1
##  3. log file
##  4. the command to be executed
##  5. the number of processors to use (if applicable)
##  6. the names of the jobs on whose completion the job to be submitted depends
##     (if applicable)
##  7. a run time limit; default is 4hrs.
## 
################################################################################

submitJob ()
{
  QSUB_MEM_REQUESTED="$1"
  shift
  QSUB_JOB_NAME="$1"
  shift
  QSUB_LOG_FILE="$1"
  shift
  QSUB_COMMAND="$1"
  shift
  
  QSUB_NUM_PROCESSORS=1
  if [ "$1" != "" ]
  then
    QSUB_NUM_PROCESSORS="$1"
    shift
  fi

  QSUB_DEPENDENCY_STRING=""
  if [ "$1" != "" ]
  then
    QSUB_DEPENDENCY_STRING="-hold_jid $1"
    shift
  fi

  QSUB_H_RT="14400"
  if [ "$1" != "" ]
  then
    QSUB_H_RT="$1"
    shift
  fi

  SGE_MEM_FREE_OPTION=
  if [ "$SGE_MEM_FREE" != "" ]
  then
    SGE_MEM_FREE_OPTION="-l $SGE_MEM_FREE=$QSUB_MEM_REQUESTED"
  fi

  H_RT_OPTION=
  if [ "$H_RT" != "" ]
  then
    H_RT_OPTION="-l $H_RT=$QSUB_H_RT"
  fi

  if [ "$USE_QSUB" = "TRUE" ]
  then
    if [ "$QSUB_NUM_PROCESSORS" = "" ]
    then
      qsub -V $SGE_MEM_FREE_OPTION $H_RT_OPTION -N $QSUB_JOB_NAME -o $QSUB_LOG_FILE $QSUB_ADDITIONAL_OPTIONS $QSUB_COMMAND
    else
      qsub -V $SGE_MEM_FREE_OPTION $H_RT_OPTION -N $QSUB_JOB_NAME $QSUB_DEPENDENCY_STRING -o $QSUB_LOG_FILE -pe smp $QSUB_NUM_PROCESSORS \
        $RESERVE_OPTION $QSUB_ADDITIONAL_OPTIONS $QSUB_COMMAND
    fi
  else
    if [ $NUM_PROCESSES_THRESH -lt $QSUB_NUM_PROCESSORS ]
    then
      QSUB_NUM_PROCESSORS=$NUM_PROCESSES_THRESH
    fi
  
    NEW_NUM_PROCESSES=`echo $NUM_PROCESSES + $QSUB_NUM_PROCESSORS | bc`
    if [ $NEW_NUM_PROCESSES -gt $NUM_PROCESSES_THRESH ]
    then
      echo "Waiting for previous processes to finish."
      wait
      NUM_PROCESSES=0
    fi
    echo "Running command:"
    echo $QSUB_COMMAND
    echo "Output in:"
    echo $QSUB_LOG_FILE
    $QSUB_COMMAND > $QSUB_LOG_FILE 2>&1 &
    NUM_PROCESSES=`echo $NUM_PROCESSES + $QSUB_NUM_PROCESSORS | bc`
  fi
  
}


################################################################################
## showStatusAndSetSubmit
################################################################################

showStatusAndSetSubmit ()
{
  RECOMPUTE_OPTION_SAMPLE="$RECOMPUTE_OPTION"
  SUBMIT="TRUE"
  if [ "$USE_QSUB" = "TRUE" ]
  then
    if [ "$SHOW_STATUS" = "TRUE" -o "$RESUBMIT" = "TRUE" ]
    then
      STATUS_LINE=`$UTIL_LIB_DIR/getJobStatus.py -E -s $SAMPLES_DIR -S $SAMPLE -C $CHUNK -A $EQP_ALIGNER -J $JOB_NAME -O $OPERATION`
      # echo "$UTIL_LIB_DIR/getJobStatus.py -E -s $SAMPLES_DIR -S $SAMPLE -C $CHUNK -A $EQP_ALIGNER -J $JOB_NAME -O $OPERATION"
      # exit 0
      STATUS=`echo "$STATUS_LINE" | cut -f 5`
      if [ "$SHOW_STATUS" = "TRUE" ]
      then
        echo "$STATUS_LINE"
      fi
  
      SUBMIT="FALSE"
      if [ "$RESUBMIT" = "TRUE" -a \( "$STATUS" = "failed" -o "$STATUS" = "crashed" -o "$STATUS" = "SGE-error" -o "$STATUS" = "not submitted" \) ]
      then
        SUBMIT="TRUE"
        if [ "$STATUS" = "SGE-error" ]
        then
          qdel $JOB_NAME
        fi
      fi
    fi
    
    JOB_NUMBER=`qstat -j $JOB_NAME 2>&1 | fgrep "job_number" | sed -e "s/job_number: *//"`
    if [ "$JOB_NUMBER" != "" ]
    then
      JOB_STATUS=`qstat | sed -e "s/^ *//" | grep "^$JOB_NUMBER" | sed -e "s/  */ /g" | cut -f 5 -d " "`
      JOB_STATUS_ADDITION=`echo $JOB_STATUS | sed -e "s/[wrht]//g"`
      if [ "$JOB_STATUS" != "" -a "$JOB_STATUS_ADDITION" != "d" -a "$JOB_STATUS_ADDITION" != "E" -a "$SHOW_STATUS" != "TRUE" -a "$FORCE_SUBMIT" != "TRUE" ]
      then
        echo "$OPERATION job $JOB_NAME is already in"
        echo "the job queue with status: $JOB_STATUS ... not submitted"
        SUBMIT="FALSE"
      fi
    fi
  fi
}


################################################################################
## clearLogFiles
################################################################################

clearLogFiles ()
{
  if [ "$KEEP_LOG_FILES" = "TRUE" -a "$OPERATION" = "compute-counts" ]
  then
    LAST_LOG_FILE=`ls -1 $LOG_DIR | grep "^$OPERATION-$EQP_ALIGNER-$CHUNK" | fgrep -v "~" | tail -1`
    echo "Moving $LAST_LOG_FILE to $OPERATION-$LOG_FILE_SUFFIX"
    cat $LOG_DIR/$LAST_LOG_FILE | grep -v "Compute counts sucessfully completed" > $LOG_DIR/$OPERATION-$LOG_FILE_SUFFIX
    rm $LOG_DIR/$LAST_LOG_FILE
  elif [ "$OPERATION" = "merge-genome-alignments" ]
  then
    NUM_LOGS=`ls -1 $LOG_DIR | fgrep $OPERATION-$EQP_ALIGNER | wc -l`
    if [ "$NUM_LOGS" != "0" ]
    then	  
      rm $LOG_DIR/$OPERATION-$EQP_ALIGNER-*.log
    fi
  else
    NUM_LOGS=`ls -1 $LOG_DIR | fgrep $OPERATION-$EQP_ALIGNER-$CHUNK | wc -l`
    if [ "$NUM_LOGS" != "0" ]
    then	  
      rm $LOG_DIR/$OPERATION-$EQP_ALIGNER-$CHUNK-*.log
    fi
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
## Main program
##
################################################################################

PROG_NAME=`basename $0`
PROG_DIR=`dirname $0`
VERSION=2.0

cd $PROG_DIR
PROG_DIR=`pwd`


################################################################################
##
## Read in options
##
################################################################################

SPLIT_READ_NUM_OPTION="-s 25000000"
NUM_PROCESSORS=6
KEEP_LOG_FILES="FALSE"
RESERVE_OPTION=""
SAMPLES_OPT=""
CHUNKS_OPT=""
MEM_REQUESTED="2G"
SHOW_STATUS="FALSE"
RESUBMIT="FALSE"
RECOMPUTE_OPTION="-r"
EQP="FALSE"
UNAMBIGUOUS_OPTION=""
UNWEIGHTED_OPTION=""
OVERLAP_OPTION=""
READ_WEIGHT_OPTION=""
CONTAINMENT_OPTION=""
SAMPLES_SUBDIR="samples"
USE_QSUB="TRUE"
FORCE_SUBMIT="FALSE"
PRINT_HELP="FALSE"
NUM_PROCESSES_THRESH=1
while [ "$1" = "-sp" -o "$1" = "-p" -o "$1" = "-k" -o "$1" = "-R" -o "$1" = "-s" -o "$1" = "-c" -o "$1" = "-A" -o \
        "$1" = "-m" -o "$1" = "-stat" -o "$1" = "-status" -o "$1" = "-resubmit" -o "$1" = "-imd" -o "$1" = "-sr" -o \
	"$1" = "-from-bed" -o "$1" = "-eqp" -o "$1" = "-S" -o "$1" = "--unweighted" -o "$1" = "--unambig" -o \
	"$1" = "-sd" -o "$1" = "-noqsub" -o "$1" = "-O" -o "$1" = "-W" -o "$1" = "-C" -o "$1" = "-M" -o "$1" = "-F" -o \
	"$1" = "-qsub_options" -o "$1" = "-h" -o "$1" = "--help" ]
do
  if [ "$1" = "-sp" ]
  then
    shift
    SPLIT_READ_NUM_OPTION="-s $1"
    shift
    SPLIT="TRUE"
  fi

  if [ "$1" = "-p" ]
  then
    shift
    NUM_PROCESSORS="$1"
    shift
  fi

  if [ "$1" = "-k" ]
  then
    shift
    KEEP_LOG_FILES="TRUE"
  fi

  if [ "$1" = "-R" ]
  then
    shift
    RESERVE_OPTION="-R y"
  fi

  if [ "$1" = "-S" ]
  then
    shift
    STRAND_SPECIFIC="TRUE"
    STRAND_SPECIFIC_DIRECTION="$1"
    shift
  fi
  
  if [ "$1" = "-from-bed" ]
  then
    shift
    RECOMPUTE_OPTION=""
  fi
  
  if [ "$1" = "-s" ]
  then
    shift
    SAMPLES_OPT=`echo $1 | sed -e "s/[:;,|]/ /g"`
    shift
  fi
  
  if [ "$1" = "-c" ]
  then
    shift
    CHUNKS_OPT=`echo $1 | sed -e 's/[,:;|]/\\\\|/g'`
    shift
  fi

  if [ "$1" = "-m" ]
  then
    shift
    MEM_REQUESTED="$1"
    shift
  fi

  if [ "$1" = "-stat" -o "$1" = "-status" ]
  then
    shift
    SHOW_STATUS="TRUE"
  fi

  if [ "$1" = "-resubmit" ]
  then
    shift
    RESUBMIT="TRUE"
  fi

  if [ "$1" = "-sr" ]
  then
    shift
    PAIRED_END="FALSE"
  fi

  if [ "$1" = "-A" ]
  then
    shift
    ALIGNER=` echo "$1" | tr [A-Z] [a-z]`
    shift
  fi
  
  if [ "$1" = "-eqp" ]
  then
    shift
    EQP="TRUE"
  fi
  
  if [ "$1" = "--unweighted" ]
  then
    UNWEIGHTED_OPTION="$1"
    shift
  fi

  if [ "$1" = "--unambig" ]
  then
    UNAMBIGUOUS_OPTION="$1"
    shift
  fi

  if [ "$1" = "-sd" ]
  then
    shift
    SAMPLES_SUBDIR="$1"
    shift
  fi

  if [ "$1" = "-noqsub" ]
  then
    shift
    USE_QSUB="FALSE"
    NUM_PROCESSES_THRESH="$1"
    shift
  fi

  if [ "$1" = "-O" ]
  then
    shift
    OVERLAP_OPTION="-E $1 -J $1"
    shift
  fi

  if [ "$1" = "-W" ]
  then
    shift
    READ_WEIGHT_OPTION="-W $1"
    shift
  fi
  
  if [ "$1" = "-C" ]
  then
    shift
    CONTAINMENT_OPTION="-C"
  fi

  if [ "$1" = "-M" ]
  then
    shift
    INTERSECTION_MODE="$1"
    shift
  fi

  if [ "$1" = "-F" ]
  then
    shift
    FORCE_SUBMIT="TRUE"
  fi

  if [ "$1" = "-qsub_options" ]
  then
    shift
    QSUB_ADDITIONAL_OPTIONS="$1"
    shift
  fi

  if [ "$1" = "-h" -o "$1" = "--help" ]
  then
    shift
    PRINT_HELP="TRUE"
  fi

done


################################################################################
##
## Set SGE_MEM_FREE and H_RT
##
################################################################################

SGE_MEM_FREE="m_mem_free"
H_RT="h_rt"

if [ -f /etc/redhat-release ]
then
  RED_HAT_RELEASE=`cat /etc/redhat-release | sed -e 's/.* release \([0-9]*\)[.][0-9]* .*/\1/'`
  if [ "$RED_HAT_RELEASE" = "5" ]
  then
    SGE_MEM_FREE="mem_token"
    H_RT=""
  fi
fi


################################################################################
##
## Set project dir
##
################################################################################

PROJECT_DIR=`echo $PROG_DIR | sed -e "s;/bin;;"`

## Set ORGANISM, ORGANISM_SHORT, READ_LENGTH, FILE_BASE, GENE_MODEL, and RADIUS etc.
if [ ! -f $PROG_DIR/setup.sh ]
then
  if [ -f $PROJECT_DIR/exon-pipeline-files/setup.sh ]
  then
    cp $PROJECT_DIR/exon-pipeline-files/setup.sh $PROG_DIR
  else
    echo "File setup.sh neither in $PROG_DIR nor in"
    echo "$PROJECT_DIR/exon-pipeline-files found."
    echo "Please make sure that all files are properly unpacked ... exiting"
    exit 1
  fi
fi
source $PROG_DIR/setup.sh
if [ "$OPERATION" = "compute-counts" -o "$OPERATION" = "align-all" -o "$OPERATION" = "align-transcripts" ]
then
  if [ "$STRAND_SPECIFIC" = "TRUE" ]
  then
    echo "Processing data in a strand-specific manner (direction $STRAND_SPECIFIC_DIRECTION). If your data"
    echo "is not strand-specific, then please edit setup.sh to change this setting."
  else
    echo "Processing data in a strand-unspecific manner."
  fi
fi
REFERENCE_SET=`echo $FILE_BASE | sed -e "s/_.*//"`


################################################################################
##
## Print help
##
################################################################################

if [ "$1" = "" -o "$PRINT_HELP" = "TRUE" ]
then
  echo "Usage: $PROG_NAME [<options>] <processing option>"
  echo
  echo "where <options> are given by:"
  echo "[-sp <num reads>] [-s <sample>] [-c <chunk>] [-noqsub <num cores>] [-status]"
  echo "[-resubmit] [-sr] [-S <direction>] [-p <num. cores>] [-m <memory per core>] "
  echo "[-eqp] [-sd <samples sub dir>] [-O <min overlap>] [-W <min read weight>]"
  echo "[-qsub_option <qsub option>]"
  echo
  echo "where <processing option> is (one of):"
  echo "   split|split-and-clean|align-genome|align-transcripts|align-junctions|"
  echo "   align-all|compute-counts|combine-gene-counts|combine-counts|compute-fpkm|"
  echo "   merge-genome-alignments|clean-up"
  echo
  echo "split: just split the Fastq files (remove all chunk Fastq files and directories"
  echo "       see option -sp)"
  echo "align-all: combination of align-genome, align-transcripts, and align-junctions"
  echo
  echo "-sp INT: split fastq into chunks of INT reads and process each chunk"
  echo "           separately (-sp enforces the split even if chunk files already"
  echo "           already exist) - only valid with \"split\" processing option."
  echo "-s STRING: align only the given sample(s) (several sample ids can be"
  echo "    separated by ,;: or |)"
  echo "-c STRING:  align only the given chunks(s) (several chunk ids can be"
  echo "  separated by ,;: or |"
  echo "-noqsub INT: execute the jobs directly and do not submit to SGE via qsub; use"
  echo "   up to INT cores"
  echo "-status: show the status of the jobs for the chosen operation (one of"
  echo "       not submitted, running, waiting, SGE-error, complete, or failed)"
  echo "-resubmit: resubmit the jobs for which the status is not submitted or failed"
  echo "-sr: treat the samples as single read samples"
  echo "-S STRING: treat the read as originating from a strand-specific protocol with"
  echo "   direction STRING: forward (for forward-reverse) or backward (for reverse-forward)"
  echo "-p INT: number of processors to use for alignment [6]"
  echo "-m: memory to use per core (-l $SGE_MEM_FREE=X) [$MEM_REQUESTED]"
  echo "-eqp: use EQP genome alignment file as input or output (for align-genome,"
  echo "      merge-genome-alignments, compute-counts, and combine counts)."
  echo "-sd STRING: Use STRING as the sub directory containing the samples [$SAMPLES_SUBDIR]"
  echo "-O INT: Minimal overlap of a read with an exon and a junction"
  echo "-W FLOAT: Minimal weight of a read; reads with a lower weight are disregarded"
  echo "-qsub_options STRING: STRING contains additional options for qsub"
  echo
  echo "Project dir: $PROJECT_DIR"
  echo "Organism: $ORGANISM, read length: $READ_LENGTH, gene model/org: $FILE_BASE, radius: $RADIUS"
  
  exit
fi

if [ ! -d $PROJECT_DIR/$SAMPLES_SUBDIR ]
then
  echo "Directory $PROJECT_DIR/$SAMPLES_SUBDIR not found ... exiting."
  exit 1
fi


################################################################################
## score min
################################################################################

SCORE_MIN_OPTION="-sm -0.1"

################################################################################
## Aligner type
################################################################################

if [ "$ALIGNER" = "" ]
then
  ALIGNER="bowtie2"
fi

GENOME_ALIGNER="$ALIGNER"

ALIGNER_UNGAPPED_FILTER=`echo $ALIGNER | sed -n '/bowtie2/Ip'`
if [[ ! -z "$ALIGNER_UNGAPPED_FILTER" ]]
then
  ALIGNER_TYPE="ungapped"
  ALIGNMENT_REF="combined"
fi

ALIGNER_GENOMIC_FILTER=`echo $ALIGNER | sed -n '/tophat2\|star\|splicemap\|eqp\|hisat/Ip'`
if [[ ! -z "$ALIGNER_GENOMIC_FILTER" ]]
then
  ALIGNER_TYPE="genomic"
  ALIGNMENT_REF="genome"
fi

if [[ -z "$ALIGNER_GENOMIC_FILTER" && -z "$ALIGNER_UNGAPPED_FILTER" ]]
then
  echo "Unknown aligner: $ALIGNER ... exiting."
  exit 1
fi


################################################################################
##
## Aligner options
##
################################################################################


################################################################################
## strand specific option
################################################################################

ALIGNER_STRAND_SPECIFIC_OPTION=""
COUNTS_STRAND_SPECIFIC_OPTION=""
if [ "$STRAND_SPECIFIC" = "TRUE" ]
then
  if [ "$STRAND_SPECIFIC_DIRECTION" != "forward" -a "$STRAND_SPECIFIC_DIRECTION" != "backward" ]
  then
    echo "Unknown strand specific direction: $STRAND_SPECIFIC_DIRECTION ... exiting"
    exit 1
  fi
  
  if [ "$ALIGNER" = "bowtie2" ]
  then
    if [ "$STRAND_SPECIFIC_DIRECTION" = "forward" ]
    then
      ALIGNER_STRAND_SPECIFIC_OPTION="-norc"
    else
      ALIGNER_STRAND_SPECIFIC_OPTION="-nofw"
    fi
  elif [ "$OPERATION" = "align-transcripts" ]
  then
    echo "Strand-specific alignment option not available for aligner $ALIGNER"
  fi
  
  COUNTS_STRAND_SPECIFIC_OPTION="-s $STRAND_SPECIFIC_DIRECTION"
fi


################################################################################
## paired end option
################################################################################

if [ "$PAIRED_END" != "TRUE" -a "$PAIRED_END" != "FALSE" ]
then
  if [ "$PAIRED_END" != "" ]
  then
    echo "Unknown value for paired-end flag: $PAIRED_END (different from TRUE or FALSE)"
    echo "Assuming paired end reads."
  fi
  PAIRED_END="TRUE"
fi

if [ "$PAIRED_END" = "TRUE" ]
then
  READ_TYPE="pe"
  SINGLE_READ_OPTION=""
else
  READ_TYPE="sr"
  SINGLE_READ_OPTION="-sr"  
fi


################################################################################
## inner mate distance option
################################################################################

INNER_MATE_DISTANCE_OPTION=
if [ "$INNER_MATE_DISTANCE_FILE" != "" -a ! -f "$INNER_MATE_DISTANCE_FILE" ]
then
  INNER_MATE_DISTANCE=$INNER_MATE_DISTANCE_FILE
  INNER_MATE_DISTANCE_OPTION="-r $INNER_MATE_DISTANCE"
  echo "Using uniform inner mate distance $INNER_MATE_DISTANCE for all samples"
fi


################################################################################
##
## Print version
##
################################################################################

if [ "$1" = "-v" -o "$1" = "--version" ]
then
  echo "$PROG_NAME version $VERSION"
  exit 0
fi


################################################################################
##
## Set directories and files
##
################################################################################

PROJECT_BIN_DIR=$PROJECT_DIR/exon-pipeline-scripts/bin
PROCESS_DIR="$PROJECT_BIN_DIR/process-scripts"
UTIL_LIB_DIR="$PROJECT_BIN_DIR/util-lib"
ALIGN_DIR="$PROJECT_BIN_DIR/align-scripts"

PROJECT_TOOLS_DIR=$PROJECT_DIR/exon-pipeline-scripts/tools
checkTool $PROJECT_TOOLS_DIR bedtools 2.24.0
checkTool $PROJECT_TOOLS_DIR samtools 0.1.17
SAMTOOLS_EXE=$TOOL_EXE
checkTool $PROJECT_TOOLS_DIR bowtie2 2.0.5

GENOME_DIR=$PROJECT_DIR/exon-pipeline-files/genome-files
FASTA_DIR=$PROJECT_DIR/exon-pipeline-files/fasta-files
GTF_DIR=$PROJECT_DIR/exon-pipeline-files/gtf-files
MAP_DIR=$PROJECT_DIR/exon-pipeline-files/map-files

GTF_FILE=$GTF_DIR/$FILE_BASE.gtf

export JAVA_DIR=$PROJECT_DIR/exon-pipeline-scripts/java
if [ ! -d $JAVA_DIR ]
then
  echo "Directory $JAVA_DIR does not exist ... exiting"
  exit 1
fi

## Make sure that CLASSPATH is not empty
if [[ -z $CLASSPATH ]]
then
  export CLASSPATH=$HOME
fi

NUM_READS_FILE=$PROJECT_DIR/statistic-files/samples-num-reads.txt

mkdir -p $PROJECT_DIR/statistic-files


################################################################################
##
## Read in operation
##
################################################################################

SPLIT="FALSE"
ALIGN_TRANSCRIPT="FALSE"
ALIGN_JUNCTION="FALSE"
ALIGN_GENOME="FALSE"
COMPUTE_COUNTS="FALSE"

COMBINE_GENE_COUNTS="FALSE"
COMBINE_EXON_COUNTS="FALSE"
COMBINE_JUNCTION_COUNTS="FALSE"

COMPUTE_FPKM="FALSE"
CLEAN_UP="FALSE"

MERGE_GENOME_ALIGNMENTS="FALSE"
ALIGNMENT_PROGRESS="FALSE"

COMPUTE_ALIGNED_READ_NUM="FALSE"
while [ "$1" != "" ]
do

  OPERATION=""
  if [ "$1" = "split" ]
  then
    OPERATION="$1"
    SPLIT="TRUE"
  fi

  if [ "$1" = "split-and-clean" ]
  then
    OPERATION="split"
    SPLIT="TRUE"
    CLEAN_UP="TRUE"
  fi

  if [ "$1" = "align-genome" ]
  then
    OPERATION="$1"
    ALIGN_GENOME="TRUE"
  fi

  if [ "$1" = "align-transcripts" ]
  then
    OPERATION="$1"
    ALIGN_TRANSCRIPT="TRUE"
  fi

  if [ "$1" = "align-junctions" ]
  then
    OPERATION="$1"
    ALIGN_JUNCTION="TRUE"
  fi

  if [ "$1" = "align-all" ]
  then
    OPERATION="$1"
    ALIGN_TRANSCRIPT="TRUE"
    ALIGN_GENOME="TRUE"
    ALIGN_JUNCTION="TRUE"
  fi

  if [ "$1" = "compute-counts" ]
  then
    OPERATION="$1"
    COMPUTE_COUNTS="TRUE"
  fi
    
  if [ "$1" = "combine-counts" ]
  then
    OPERATION="$1"
    COMBINE_GENE_COUNTS="TRUE"
    COMBINE_EXON_COUNTS="TRUE"
    COMBINE_JUNCTION_COUNTS="TRUE"
    COUNT_TYPE="cnt"
  fi
  
  if [ "$1" = "combine-gene-counts" ]
  then
    OPERATION="$1"
    COMBINE_GENE_COUNTS="TRUE"
    COUNT_TYPE="cnt"
  fi

  if [ "$1" = "combine-exon-counts" ]
  then
    OPERATION="$1"
    COMBINE_EXON_COUNTS="TRUE"
    COUNT_TYPE="cnt"
  fi

  if [ "$1" = "combine-gene-exon-counts" ]
  then
    OPERATION="$1"
    COMBINE_GENE_COUNTS="TRUE"
    COMBINE_EXON_COUNTS="TRUE"
    COUNT_TYPE="cnt"
  fi
 
  if [ "$1" = "combine-junction-counts" ]
  then
    OPERATION="$1"
    COMBINE_JUNCTION_COUNTS="TRUE"
    COUNT_TYPE="cnt"
  fi
  
  if [ "$1" = "compute-fpkm" ]
  then
    OPERATION="$1"
    COMPUTE_FPKM="TRUE"
  fi

  if [ "$1" = "clean-up" ]
  then
    OPERATION="$1"
    CLEAN_UP="TRUE"
  fi
    
  if [ "$1" = "merge-genome-alignments" ]
  then
    OPERATION="$1"
    MERGE_GENOME_ALIGNMENTS="TRUE"
  fi
  
  if [ "$1" = "compute-aligned-read-num" ]
  then
    OPERATION="$1"
    COMPUTE_ALIGNED_READ_NUM="TRUE"
  fi

  if [ "$1" = "uniquely-aligned-reads" ]
  then
    OPERATION="$1"
    UNIQUELY_ALIGNED_READS="TRUE"
  fi
  
  MINUS_PREFIX=`echo $1 | sed -e "s/^-//"`
  if [ "$MINUS_PREFIX" != "$1" ]
  then
    echo "Unknown option: $1 ... exiting"
    exit 1
  elif [ "$OPERATION" = "" ]
  then
    echo "Unknown operation: $1 ... exiting"
    exit 1
  fi
  
  shift
  
done


################################################################################
##
## Check that the operation and the quantifier are correct
##
################################################################################

if [ "$OPERATION" = "" ]
then
 echo "No operation given ... exiting."
 exit 1
fi

if [ "$OPERATION" = "combine-counts" -o "$OPERATION" = "combine-gene-counts"  -o "$OPERATION" = "combine-exon-counts" -o "$OPERATION" = "combine-junction-counts" ]
then
 if [ "$SHOW_STATUS" = "TRUE" ]
 then
   echo "The option -status cannot be combined with operation $OPERATION"
   exit 1
 fi
fi



################################################################################
##
## Set aligner variables
##
################################################################################

ALIGN_SCRIPT=$ALIGN_DIR/$ALIGNER-align.sh
if [ "$ALIGNER" = "bowtie2" ]
then
  ALIGN_SCRIPT="$ALIGN_SCRIPT $SCORE_MIN_OPTION"
fi

GENOME_ALIGN_SCRIPT=$ALIGN_DIR/$GENOME_ALIGNER-align.sh

if [ "$ALIGNER" != "$GENOME_ALIGNER" ]
then
  EQP_ALIGNER="$ALIGNER-$GENOME_ALIGNER"
  # In compute-counts.sh the name of the combined SAM file and the counts file is based
  # on the transcript alignment file which has only the aligner in its name.
  if [ "$COMBINE_GENE_COUNTS" = "TRUE" -o "$COMBINE_EXON_COUNTS" = "TRUE" -o "$COMBINE_JUNCTION_COUNTS" = "TRUE" ]
  then
    EQP_ALIGNER="$ALIGNER"
  fi
  if [ "$OPERATION" = "align-genome" -o "$OPERATION" = "align-all" ]
  then
    echo "It is not possible to use combined aligner $EQP_ALIGNER for genome"
    echo "alignments ... exiting."
    exit 1
  fi
else
  EQP_ALIGNER="$ALIGNER"
fi

if [ "$EQP" = "TRUE" ]
then
  if [ "$ALIGNER_TYPE" = "genomic" ]
  then
    echo "Option -eqp ignored for aligner $ALIGNER"
    EQP = "FALSE"
  else
    EQP_ALIGNER="$EQP_ALIGNER-eqp"
    ALIGNER_TYPE="genomic"
    ALIGNMENT_REF="genome"
  fi
fi


################################################################################
##
## Prepare genome alignment
##
################################################################################

if [ "$GENOME_ALIGNER" = "bowtie2" -o "$GENOME_ALIGNER" = "tophat2" ]
then
  GENOME_INDEX_DIR=$GENOME_DIR/bowtie2Index
  GENOME_INDEX=$GENOME_INDEX_DIR/bowtie2Index
elif [ "$GENOME_ALIGNER" = "hisat" ]
then
  GENOME_INDEX_DIR=$GENOME_DIR/hisatIndex
  GENOME_INDEX=$GENOME_INDEX_DIR/hisatIndex
else
  GENOME_INDEX_DIR=$GENOME_DIR/${GENOME_ALIGNER}Index
  GENOME_INDEX=$GENOME_DIR/${GENOME_ALIGNER}Index
fi

if [ "$ALIGN_GENOME" = "TRUE" -a ! -d "$GENOME_INDEX_DIR" -a ! -h "$GENOME_INDEX_DIR" ]
then
  echo "Directory $GENOME_INDEX_DIR not found"
  echo "Please create a ${GENOME_ALIGNER} index for the genome Fasta file in the above directory."
  echo "... exiting"
  exit 1
fi


################################################################################
##
## Prepare transcript and junction alignment
##
################################################################################

if [ "$ALIGN_TRANSCRIPT" = "TRUE" -o "$ALIGN_JUNCTION" = "TRUE" ]
then
  if [ "$ALIGNER" = "bowtie2" ]
  then
    TRANSCRIPT_INDEX_DIR=$FASTA_DIR/${FILE_BASE}-bowtie2Index
    JUNCTION_INDEX_DIR=$FASTA_DIR/${FILE_BASE}-${RADIUS}-junction-bowtie2Index
  
    TRANSCRIPT_INDEX=$TRANSCRIPT_INDEX_DIR/bowtie2Index
    JUNCTION_INDEX=$JUNCTION_INDEX_DIR/bowtie2Index
  else
    TRANSCRIPT_INDEX_DIR=$FASTA_DIR/${FILE_BASE}-${ALIGNER}Index
    JUNCTION_INDEX_DIR=$FASTA_DIR/${FILE_BASE}-${RADIUS}-junction-${ALIGNER}Index
    
    TRANSCRIPT_INDEX=$TRANSCRIPT_INDEX_DIR
    JUNCTION_INDEX=$JUNCTION_INDEX_DIR
  fi
  
  if [ ! -d $TRANSCRIPT_INDEX_DIR -a ! -h $TRANSCRIPT_INDEX_DIR ]
  then
    echo "Directory $TRANSCRIPT_INDEX_DIR not found"
    echo "Please create a ${ALIGNER} index for the transcript Fasta file in the above directory."
    echo "... exiting"
    exit 1
  fi
  
  if [ ! -d $JUNCTION_INDEX_DIR -a ! -h $JUNCTION_INDEX_DIR ]
  then
    echo "Directory $JUNCTION_INDEX_DIR not found"
    echo "Please create a ${ALIGNER} index for the junction Fasta file in the above directory."
    echo "... exiting"
    exit 1
  fi
fi

################################################################################
##
## Prepare file variables for combine counts
##
################################################################################

SAMPLE_GENE_COUNT_FILES=""
SAMPLE_EXON_COUNT_FILES=""
SAMPLE_JUNCTION_COUNT_FILES=""
SAMPLE_JUNCTION_SAM_COUNT_FILES=""

SAMPLE_BED_JUNCTION_COUNT_FILES=""
SAMPLE_SAM_JUNCTION_COUNT_FILES=""
SAMPLE_EQP_JUNCTION_COUNT_FILES=""


################################################################################
##
## Loop over samples
##
################################################################################

SAMPLES_DIR=$PROJECT_DIR/$SAMPLES_SUBDIR
if [[ -z "$SAMPLES_OPT" ]]
then
  SAMPLES=`ls -1d $SAMPLES_DIR/*/ | sed -e "s;$SAMPLES_DIR/;;" | sed -e 's;/$;;' | fgrep -v bin | fgrep -v files | tr '\n' ' '`
else
  SAMPLES=$SAMPLES_OPT
fi

if [ "$SHOW_STATUS" = "FALSE" ]
then
  echo "Samples: $SAMPLES"
fi

FIRST_SAMPLE=
FIRST_CHUNK=
I=0
for SAMPLE in $SAMPLES
do

  if [ "$FIRST_SAMPLE" = "" ]
  then
    FIRST_SAMPLE=$SAMPLE
  fi

  if [ "$SHOW_STATUS" = "FALSE" ]
  then
    echo "Sample $SAMPLE"
  else
    if [ "$FIRST_SAMPLE" = "$SAMPLE" ]
    then
      if [ "$SPLIT" == "FALSE" ]
      then
        echo "Sample/Chunk/Operation/Job name/Status" | tr "/" "\t"
      else
        echo "Sample/Operation/Job name/Status" | tr "/" "\t"
      fi
    fi
  fi
  
  SAMPLE_DIR=$SAMPLES_DIR/$SAMPLE
  if [ ! -d $SAMPLE_DIR ]
  then
     echo "Directory $SAMPLE_DIR not found ... exiting"
     exit 1
  fi

  FASTQ_DIR=$SAMPLE_DIR/fastq-files
  if [ ! -d $FASTQ_DIR ]
  then
    echo "Making directory $FASTQ_DIR"
    mkdir $FASTQ_DIR
  fi

  LOG_DIR=$SAMPLE_DIR/log-files
  if [ ! -d $LOG_DIR ]
  then
    echo "Creating directory $LOG_DIR"
    mkdir $LOG_DIR
  fi

  SAMPLE_RANDOM_FILE="$SAMPLE_DIR/.random-number.txt"
  if [ -f "$SAMPLE_RANDOM_FILE" ]
  then
    SAMPLE_RANDOM_NUMBER=`cat $SAMPLE_RANDOM_FILE | head -1`
  else
    SAMPLE_RANDOM_NUMBER=`awk -v min=0 -v max=10000 'BEGIN{srand(); print int(min+rand()*(max-min+1))}'`
    echo $SAMPLE_RANDOM_NUMBER > $SAMPLE_RANDOM_FILE
  fi

  JOB_NAME=SP-$SAMPLE-$SAMPLE_RANDOM_NUMBER
  SPLIT_IN_JOB_QUEUE=`qstat -j $JOB_NAME 2>&1 | tr '\n' ' ' | fgrep -v "Following jobs do not exist" | fgrep -v "not found"`
  if [ "$SPLIT_IN_JOB_QUEUE" != "" -a "$FORCE_SUBMIT" != "TRUE" ]
  then
    echo "Split job $JOB_NAME is still in job queue ... skipping sample $SAMPLE"
    continue
  fi

  CHUNKS=`ls -1 $SAMPLE_DIR | grep '^C[0-9][0-9][0-9]$' | tr '\n' ' '`
  if [[ -z "$CHUNKS" && "$SPLIT" = "FALSE" ]]
  then
    echo "Please split up the Fastq file(s) of sample $SAMPLE into chunks by calling:"
    echo
    echo "    process-samples.sh -s $SAMPLE split"
    echo
    echo "otherwise the sample $SAMPLE cannot be processed further ... skipping sample $SAMPLE"
    continue
  fi

  if [[ -z "$CHUNKS" || "$SPLIT" = "TRUE" ]]
  then
    if [ "$PAIRED_END" = "TRUE" ]
    then
      FASTQ_FILES_1=`ls -1 $SAMPLE_DIR | grep "${SAMPLE}[-_]1.f[ast]*q[.gz]*"`
      FASTQ_FILES_2=`ls -1 $SAMPLE_DIR | grep "${SAMPLE}[-_]2.f[ast]*q[.gz]*"`
    
      if [ "$FASTQ_FILES_1" = "" ]
      then
    
        FASTQ_FILES_1=`ls -1 $SAMPLE_DIR | grep "${SAMPLE}_[1-9]_1.f[ast]*q[.gz]*" | tr '\n' ' ' | sed -e 's/ /#/g'`
        FASTQ_FILES_2=`ls -1 $SAMPLE_DIR | grep "${SAMPLE}_[1-9]_2.f[ast]*q[.gz]*" | tr '\n' ' ' | sed -e 's/ /#/g'`    
    
        if [[ -z $FASTQ_FILES_1 ]]
        then
          echo "Paired-end mode: No Fastq files for first read (_1) in $SAMPLE_DIR found ... exiting."
          exit 1
        fi
      fi

      if [ "$FASTQ_FILES_2" = "" ]
      then
        if [[ -z $FASTQ_FILES_2 ]]
        then
          echo "Paired-end mode: No Fastq files for second read (_2) in $SAMPLE_DIR found ... exiting."
          exit 1
        fi
      fi
    
      FASTQ_FILES_1_NUM=`echo $FASTQ_FILES_1 | tr '#' '\n' | wc -l`
      FASTQ_FILES_2_NUM=`echo $FASTQ_FILES_2 | tr '#' '\n' | wc -l`
    
      if [ "$FASTQ_FILES_1_NUM" != "$FASTQ_FILES_2_NUM" ]
      then
        echo "Differing number of Fastq for read 1 vs read 2 found: $FASTQ_FILES_1_NUM vs $FASTQ_FILES_2_NUM ... exiting"
        exit 1
      fi
      FASTQ_FILES_2_OPTION="-2 $FASTQ_FILES_2"
    else
      FASTQ_FILES_1=`ls -1 $SAMPLE_DIR | grep "${SAMPLE}.*[.]f[ast]*q[.gz]*" | fgrep -v "combined" | tr '\n' '#'`
      if [[ -z $FASTQ_FILES_1 ]]
      then
        echo "Single read mode: No Fastq files in $SAMPLE_DIR found ... exiting."
        exit 1
      fi
    
      FASTQ_FILES_1_NUM=`echo $FASTQ_FILES_1 | tr '#' '\n' | wc -l`
      FASTQ_FILES_2_OPTION=""
    fi
    
    #############################################################################
    ##
    ## Submit split job and create Fastq chunks
    ##
    #############################################################################

    # Note that even without compression and splitting of the files we want to create new fastq files
    # to obtain new identifiers

    DATE=`date +%y%m%d-%H%M`

    OPERATION="split"
    JOB_NAME=SP-$SAMPLE
    CHUNK=none
    showStatusAndSetSubmit
    
    if [ "$SUBMIT" = "TRUE" ]
    then
      if [ "$CLEAN_UP" = "TRUE" ]
      then
        echo "Removing Fastq files and chunk directories for sample $SAMPLE"
        rm -f $FASTQ_DIR/*-C[0-9][0-9][0-9]*.fq.gz
        rm -rf $SAMPLE_DIR/C[0-9][0-9][0-9]
      fi

      NUM_LOGS=`ls -1 $LOG_DIR | fgrep split- | wc -l`
      if [ "$NUM_LOGS" != "0" ]
      then	  
        rm $LOG_DIR/split-*.log
      fi
      
      JOB_NAME=SP-$SAMPLE
      CMD="$PROCESS_DIR/split-fastq-file.sh $SPLIT_READ_NUM_OPTION -1 $FASTQ_FILES_1 $FASTQ_FILES_2_OPTION $PROJECT_DIR $SAMPLE_DIR $FASTQ_DIR"
      
      submitJob 3G $JOB_NAME $LOG_DIR/split-$DATE.log "$CMD"
    fi
    continue
  fi

  if [ "$COMBINE_GENE_COUNTS" = "TRUE" -o "$COMBINE_EXON_COUNTS" = "TRUE" -o "$COMBINE_JUNCTION_COUNTS" = "TRUE" ]
  then
    if [ ! -d $SAMPLE_DIR/count-files ]
    then
      echo "Making directory $SAMPLE_DIR/count-files"
      mkdir $SAMPLE_DIR/count-files
    fi
    
    SAMPLE_GENE_COUNT_FILE=$SAMPLE_DIR/count-files/$SAMPLE-$EQP_ALIGNER-$READ_TYPE-gene.$COUNT_TYPE
    SAMPLE_EXON_COUNT_FILE=$SAMPLE_DIR/count-files/$SAMPLE-$EQP_ALIGNER-$READ_TYPE-exon.$COUNT_TYPE
    SAMPLE_JUNCTION_COUNT_FILE=$SAMPLE_DIR/count-files/$SAMPLE-$EQP_ALIGNER-$READ_TYPE-junction.$COUNT_TYPE
    SAMPLE_JUNCTION_SAM_COUNT_FILE=$SAMPLE_DIR/count-files/$SAMPLE-$EQP_ALIGNER-$READ_TYPE-junction-sam.$COUNT_TYPE
    
    CHUNK_GENE_COUNT_FILES=""
    CHUNK_EXON_COUNT_FILES=""
    CHUNK_JUNCTION_COUNT_FILES=""
    CHUNK_JUNCTION_SAM_COUNT_FILES=""
  fi

  if [ "$MERGE_GENOME_ALIGNMENTS" = "TRUE" -a "$ALIGNER_TYPE" = "genomic" ]
  then
    GENOME_SAM_CHUNK_FILES=""
  fi
  

  ################################################################################
  ##
  ## Loop over chunks
  ##
  ################################################################################

  if [ "$CHUNKS_OPT" != "" ]
  then
    CHUNKS=`echo $CHUNKS | tr ' ' '\n' | grep $CHUNKS_OPT | tr '\n' ' '`
  fi

  for CHUNK in $CHUNKS
  do

    #echo "Chunk $CHUNK"
    FASTQ_FILE_CHUNK_1=`ls -1 $FASTQ_DIR | fgrep "$CHUNK" | grep -v '2.fq.gz$'`
    FASTQ_FILE_CHUNK_2=`echo $FASTQ_FILE_CHUNK_1 | sed -e 's/\(-C[0-9][0-9][0-9]\)_1[.]fq/\1_2.fq/'` 

    DATE=`date +%y%m%d-%H%M`
    CHUNK_DIR=$SAMPLE_DIR/$CHUNK
    if [ ! -d "$CHUNK_DIR" ]
    then
      echo "Directory $CHUNK_DIR not found ... exiting."
      exit 1
    fi

    if [ "$FIRST_CHUNK" = "" ]
    then
      FIRST_CHUNK=$CHUNK
    fi

    CHUNK_RANDOM_FILE="$CHUNK_DIR/.random-number.txt"
    if [ -f "$CHUNK_RANDOM_FILE" ]
    then
      CHUNK_RANDOM_NUMBER=`cat $CHUNK_RANDOM_FILE | head -1`
    else
      CHUNK_RANDOM_NUMBER=`awk -v min=0 -v max=10000 'BEGIN{srand(); print int(min+rand()*(max-min+1))}'`
      echo $CHUNK_RANDOM_NUMBER > $CHUNK_RANDOM_FILE
    fi

    JOB_NAME_SUFFIX=$SAMPLE-$CHUNK-$ALIGNER-$REFERENCE_SET-$CHUNK_RANDOM_NUMBER
    FILE_PREFIX=${SAMPLE}-${CHUNK}-$ALIGNER
    GENOME_FILE_PREFIX=${SAMPLE}-${CHUNK}-$GENOME_ALIGNER
    LOG_FILE_SUFFIX=${ALIGNER}-$CHUNK-$DATE.log

    EQP_JOB_NAME_SUFFIX=$SAMPLE-$CHUNK-$EQP_ALIGNER-$REFERENCE_SET-$CHUNK_RANDOM_NUMBER
    EQP_FILE_PREFIX=${SAMPLE}-${CHUNK}-$EQP_ALIGNER
    EQP_LOG_FILE_SUFFIX=${EQP_ALIGNER}-$CHUNK-$DATE.log
	
    FAST_ALIGN_OPTION=""

    FASTQ_FILE_ALIGN_1=$FASTQ_DIR/$FASTQ_FILE_CHUNK_1
    FASTQ_FILE_ALIGN_2_OPTION=""
    if [ "$PAIRED_END" = "TRUE" ]
    then
      FASTQ_FILE_ALIGN_2=$FASTQ_DIR/$FASTQ_FILE_CHUNK_2
      FASTQ_FILE_ALIGN_2_OPTION="-2 $FASTQ_FILE_ALIGN_2"
    fi

    ################################################################################
    ## Transcript alignments
    ################################################################################
  
    JOB_NAME=TA-$JOB_NAME_SUFFIX
    if [ "$ALIGN_TRANSCRIPT" = "TRUE" ]
    then
      OPERATION="align-transcripts"
      if  [ "$ALIGNER_TYPE" = "ungapped" ]
      then
	showStatusAndSetSubmit
	
	if [ "$SUBMIT" = "TRUE" ]
	then
	  clearLogFiles

          CMD="$ALIGN_SCRIPT $ALIGNER_STRAND_SPECIFIC_OPTION -1 $FASTQ_FILE_ALIGN_1 $FASTQ_FILE_ALIGN_2_OPTION \
                 -p $NUM_PROCESSORS $PROJECT_DIR $CHUNK_DIR $TRANSCRIPT_INDEX transcript"

	  submitJob $MEM_REQUESTED $JOB_NAME $LOG_DIR/$OPERATION-$LOG_FILE_SUFFIX "$CMD" $NUM_PROCESSORS            
	fi
      else
	echo "Operation $OPERATION is not available for the alignment program $ALIGNER"
      fi
    fi
    DEPENDENCY_STRING="$JOB_NAME"
    

    ################################################################################
    ## Junction alignments: Note that since junctions are contatenations of
    ## genome plus strand sequences we cannot use a strand-specific option
    ## for alignment
    ################################################################################
  
    JOB_NAME=JA-$JOB_NAME_SUFFIX
    if [ "$ALIGN_JUNCTION" = "TRUE" ]
    then
      OPERATION="align-junctions"
      if  [ "$ALIGNER_TYPE" = "ungapped" ]
      then
	showStatusAndSetSubmit
	
	if [ "$SUBMIT" = "TRUE" ]
	then
	  clearLogFiles

	  CMD="$ALIGN_SCRIPT -1 $FASTQ_FILE_ALIGN_1 $FASTQ_FILE_ALIGN_2_OPTION -p $NUM_PROCESSORS $PROJECT_DIR $CHUNK_DIR $JUNCTION_INDEX junction"

          submitJob $MEM_REQUESTED $JOB_NAME $LOG_DIR/$OPERATION-$LOG_FILE_SUFFIX "$CMD" $NUM_PROCESSORS "dummy" "12:00:00"
        fi
      else
	echo "Operation $OPERATION is not available for the alignment program $ALIGNER"
      fi
    fi
    DEPENDENCY_STRING="$DEPENDENCY_STRING,$JOB_NAME"
    

    ################################################################################
    ## Genome alignments
    ################################################################################
  
    JOB_NAME=GA-$JOB_NAME_SUFFIX
    if [ "$ALIGN_GENOME" = "TRUE" -a "$EQP" = "FALSE" -a "$ALIGNER" != "tophat2" -a "$ALIGNER" != "hisat" ]
    then
      OPERATION="align-genome"
      showStatusAndSetSubmit
	
      if [ "$SUBMIT" = "TRUE" ]
      then
	clearLogFiles

	if [ "$ALIGNER" = "star" ]
	then
	  MEM_REQUESTED=5G
	fi

	CMD="$GENOME_ALIGN_SCRIPT -1 $FASTQ_FILE_ALIGN_1 $FASTQ_FILE_ALIGN_2_OPTION -p $NUM_PROCESSORS $PROJECT_DIR $CHUNK_DIR $GENOME_INDEX genome"
	
	submitJob $MEM_REQUESTED $JOB_NAME $LOG_DIR/$OPERATION-$LOG_FILE_SUFFIX "$CMD" $NUM_PROCESSORS
      fi
    fi
    DEPENDENCY_STRING="$DEPENDENCY_STRING,$JOB_NAME"
       

    ################################################################################
    ## Compute counts
    ## Note that since compute-counts.sh and compute-genomic-alignments.sh both
    ## call compute-intersection-bed-file.sh we have to make sure that they
    ## are not both run at the same time (and since $EQP_JOB_NAME_SUFFIX is not
    ## guaranteed to contain the <aligner>-eqp suffix, we need to enforce that we use
    ## the correct name by using $JOB_NAME_SUFFIX-eqp instead).
    ##
    ## Options:
    ##  -R FLOAT: Minimal weight of a read; reads with a lower weight are
    ##     disregarded [0.01]
    ##  -u: do not use read weights and compute unweighted read counts only
    ##      (the read weight file is not computed and the minimal weight of a read
    ##      is set to -1)
    ################################################################################
  
    JOB_NAME=CC-$EQP_JOB_NAME_SUFFIX
    if [ "$COMPUTE_COUNTS" = "TRUE" -a "$ALIGNER_TYPE" = "ungapped" -a "$EQP" = "FALSE" ]
    then
      OPERATION="compute-counts"	
      showStatusAndSetSubmit
	
      if [ "$SUBMIT" = "TRUE" ]
      then
	clearLogFiles

	MEMORY_UNIT_IN_GB=6G
	NUM_PROCESSORS_CC=4
	
	MAX_MEMORY_IN_GB=`echo $MEMORY_UNIT_IN_GB | sed -e 's/G$//' | sed -e "s/^/$NUM_PROCESSORS_CC * /" | bc | sed -e 's/$/G/'`

	CMD="$PROCESS_DIR/compute-counts.sh -a -nj $COUNTS_STRAND_SPECIFIC_OPTION $UNAMBIGUOUS_OPTION $RECOMPUTE_OPTION $CONTAINMENT_OPTION \
	       $SINGLE_READ_OPTION $UNWEIGHTED_OPTION $READ_WEIGHT_OPTION $OVERLAP_OPTION -M $MAX_MEMORY_IN_GB $PROJECT_DIR \
               $CHUNK_DIR/sam-files/$FILE_PREFIX-transcript-$READ_TYPE.sam.gz \
               $CHUNK_DIR/sam-files/$FILE_PREFIX-junction-$READ_TYPE.sam.gz \
               $CHUNK_DIR/sam-files/$GENOME_FILE_PREFIX-genome-$READ_TYPE.sam.gz"

	submitJob $MEMORY_UNIT_IN_GB $JOB_NAME $LOG_DIR/$OPERATION-$EQP_LOG_FILE_SUFFIX "$CMD" $NUM_PROCESSORS_CC \
	   "$DEPENDENCY_STRING,GA-$JOB_NAME_SUFFIX-eqp" "12:00:00"
      fi
    fi


    ################################################################################
    ## Align to the genome with EQP
    ################################################################################
  
    JOB_NAME=GA-$EQP_JOB_NAME_SUFFIX
    if [ "$ALIGN_GENOME" = "TRUE" -a "$EQP" = "TRUE" ]
    then
      OPERATION="align-genome"
      showStatusAndSetSubmit
	
      if [ "$SUBMIT" = "TRUE" ]
      then
	clearLogFiles

	CMD="$PROCESS_DIR/compute-genomic-alignments.sh $PROJECT_DIR \
               $CHUNK_DIR/sam-files/$FILE_PREFIX-transcript-$READ_TYPE.sam.gz \
               $CHUNK_DIR/sam-files/$FILE_PREFIX-junction-$READ_TYPE.sam.gz \
               $CHUNK_DIR/sam-files/$GENOME_FILE_PREFIX-genome-$READ_TYPE.sam.gz \
	       $CHUNK_DIR/sam-files/$EQP_FILE_PREFIX-genome-$READ_TYPE.sam.gz"
	
	submitJob $MEM_REQUESTED $JOB_NAME $LOG_DIR/$OPERATION-$EQP_LOG_FILE_SUFFIX "$CMD" 3 "$DEPENDENCY_STRING,CC-$EQP_JOB_NAME_SUFFIX" "12:00:00"
      fi
    fi
    DEPENDENCY_STRING="$DEPENDENCY_STRING,$JOB_NAME"


    ################################################################################
    ## Collect chunk count files
    ################################################################################
  
    if [ "$COMBINE_GENE_COUNTS" = "TRUE" ]
    then
      if [ -f $CHUNK_DIR/count-files/$EQP_FILE_PREFIX-$ALIGNMENT_REF-$READ_TYPE-gene.$COUNT_TYPE ]
      then
        CHUNK_GENE_COUNT_FILES="$CHUNK_GENE_COUNT_FILES $CHUNK_DIR/count-files/$EQP_FILE_PREFIX-$ALIGNMENT_REF-$READ_TYPE-gene.$COUNT_TYPE"
      else
        echo "Gene count file $CHUNK_DIR/count-files/$EQP_FILE_PREFIX-$ALIGNMENT_REF-$READ_TYPE-gene.$COUNT_TYPE not found ... skipping."
      fi
    fi
      
    if [ "$COMBINE_EXON_COUNTS" = "TRUE" ]
    then
      if [ -f $CHUNK_DIR/count-files/$EQP_FILE_PREFIX-$ALIGNMENT_REF-$READ_TYPE-exon.$COUNT_TYPE ]
      then
        CHUNK_EXON_COUNT_FILES="$CHUNK_EXON_COUNT_FILES $CHUNK_DIR/count-files/$EQP_FILE_PREFIX-$ALIGNMENT_REF-$READ_TYPE-exon.$COUNT_TYPE"
      else
        echo "Exon count file $CHUNK_DIR/count-files/$EQP_FILE_PREFIX-$ALIGNMENT_REF-$READ_TYPE-exon.$COUNT_TYPE not found ... skipping."
      fi
    fi
      
    if [ "$COMBINE_JUNCTION_COUNTS" = "TRUE" ]
    then
      if [ -f $CHUNK_DIR/count-files/$EQP_FILE_PREFIX-$ALIGNMENT_REF-$READ_TYPE-junction.$COUNT_TYPE ]
      then
        CHUNK_JUNCTION_COUNT_FILES="$CHUNK_JUNCTION_COUNT_FILES $CHUNK_DIR/count-files/$EQP_FILE_PREFIX-$ALIGNMENT_REF-$READ_TYPE-junction.$COUNT_TYPE"
      else
        echo "Junction count file $CHUNK_DIR/count-files/$EQP_FILE_PREFIX-$ALIGNMENT_REF-$READ_TYPE-junction.$COUNT_TYPE not found ... skipping."
      fi
      
      ## CHUNK_JUNCTION_SAM_COUNT_FILES="$CHUNK_JUNCTION_COUNT_FILES $CHUNK_DIR/count-files/$EQP_FILE_PREFIX-$ALIGNMENT_REF-$READ_TYPE-junction-sam.$COUNT_TYPE"
    fi


    ################################################################################
    ## Merge genome alignments
    ################################################################################
    
    if [ "$MERGE_GENOME_ALIGNMENTS" = "TRUE" -a "$ALIGNER_TYPE" = "genomic" ]
    then
      GENOME_SAM_CHUNK_FILES="$GENOME_SAM_CHUNK_FILES $CHUNK_DIR/sam-files/$EQP_FILE_PREFIX-genome-$READ_TYPE.sam.gz"
    fi
    

    ################################################################################
    ## Compute uniquely aligned reads
    ################################################################################
    
    if [ "$UNIQUELY_ALIGNED_READS" = "TRUE" ]
    then
      if [ "$GENOME_ALIGNER" != "bowtie2" ]
      then
        CHUNK_WEIGHT_FILE="$CHUNK_DIR/weight-files/$EQP_FILE_PREFIX-genome-$READ_TYPE.wgt"
      else
        CHUNK_WEIGHT_FILE="$CHUNK_DIR/weight-files/$EQP_FILE_PREFIX-final.wgt"
      fi
      
      NUM_UNIQUELY_ALIGNED_READS=`cut -f 2 $CHUNK_WEIGHT_FILE | sed -e "s/\(.*\)/#\1#/" | fgrep "#1#" | wc -l`
      NUM_ALIGNED_READS=`cut -f 2 $CHUNK_WEIGHT_FILE | sed -e "s/\(.*\)/#\1#/" | fgrep -v "#0#" | wc -l`
      UNIQUELY_ALIGNED_FRACTION=`echo "scale=2; $NUM_UNIQUELY_ALIGNED_READS * 100 / $NUM_ALIGNED_READS" | bc`
      echo "$SAMPLE $CHUNK $ALIGNER $UNIQUELY_ALIGNED_FRACTION $NUM_UNIQUELY_ALIGNED_READS $NUM_ALIGNED_READS" | tr ' ' '\t'
    fi
              
    ################################################################################
    ## Clean up
    ################################################################################

    if [ "$CLEAN_UP" = "TRUE" ]
    then
      rm $CHUNK_DIR/sam-files/$EQP_FILE_PREFIX-trans-junction-$READ_TYPE.sam.gz
      rm $CHUNK_DIR/sam-files/$EQP_FILE_PREFIX-genome-$READ_TYPE.sam.gz
      rm $CHUNK_DIR/bed-files/*
      rm $CHUNK_DIR/weight-files/*
      if [ -f $CHUNK_DIR/QC/Picard ]
      then
	rm $CHUNK_DIR/QC/Picard/*.[bs]am
      fi
    fi
    
    ################################################################################
    ## Compute aligned read num
    ################################################################################

    if [ "$COMPUTE_ALIGNED_READ_NUM" = "TRUE" ]
    then
      if [ ! -f $CHUNK_DIR/sam-files/$EQP_FILE_PREFIX-genome-$READ_TYPE.sam.gz ]
      then
        echo "File $CHUNK_DIR/sam-files/$EQP_FILE_PREFIX-genome-$READ_TYPE.sam.gz not found ... exiting."
	exit 1
      fi
      
      if [ ! -f $PROJECT_DIR/statistic-files/$EQP_ALIGNER-aligned-read-num.txt ]
      then
        echo -n > $PROJECT_DIR/statistic-files/$EQP_ALIGNER-aligned-read-num.txt
      elif [ "$FORCE_SUBMIT" = "TRUE" -a "$SAMPLE" = "$FIRST_SAMPLE" -a "$CHUNK" = "$FIRST_CHUNK" ]
      then
        echo -n > $PROJECT_DIR/statistic-files/$EQP_ALIGNER-aligned-read-num.txt
      fi

      CHUNK_ALIGNED_READ_NUM_LINE=`grep "^$SAMPLE" $PROJECT_DIR/statistic-files/$EQP_ALIGNER-aligned-read-num.txt | awk "/$SAMPLE\t/" | awk "/\t$CHUNK\t/"`

      if [ "$CHUNK_ALIGNED_READ_NUM_LINE" = "" ]
      then
        echo "$SAMPLE - $CHUNK"
        date
        CHUNK_ALIGNED_READ_NUM=`zcat $CHUNK_DIR/sam-files/$EQP_FILE_PREFIX-genome-$READ_TYPE.sam.gz | $SAMTOOLS_EXE view -S -F 4 - | \
          cut -f 1 | sort -u -S 6G | wc -l`
        if [ $? -ne 0 ]
        then
          echo "Computation of aligned read num failed ... exiting"
          exit 1
        fi
        echo "$SAMPLE#$CHUNK#$CHUNK_ALIGNED_READ_NUM" | tr '#' '\t' >> $PROJECT_DIR/statistic-files/$EQP_ALIGNER-aligned-read-num.txt
      fi
    fi

    
    ################################################################################
    ## Print last lines
    ################################################################################

    DATE_FILTER=`date +"%b %d %H"`
    if [ "$ALIGNMENT_PROGRESS" = "TRUE" ]
    then
      echo "$SAMPLE - $CHUNK"
      if [ -f $CHUNK_DIR/sam-files/${SAMPLE}-${CHUNK}-transcript-$READ_TYPE.sam ]
	then
	  echo "Transcript"
	  tail -1 $CHUNK_DIR/sam-files/${SAMPLE}-${CHUNK}-transcript-$READ_TYPE.sam
	  echo
	elif [ -f $CHUNK_DIR/sam-files/${SAMPLE}-${CHUNK}-transcript-$READ_TYPE.sam.gz ]
	then
	  echo "Transcript"
	  ls -l $CHUNK_DIR/sam-files/${SAMPLE}-${CHUNK}-transcript-$READ_TYPE.sam.gz | fgrep "$DATE_FILTER"
	fi

      if [ -f $CHUNK_DIR/sam-files/${SAMPLE}-${CHUNK}-genome-$READ_TYPE.sam ]
	then
	  echo "Genome"
	  tail -1 $CHUNK_DIR/sam-files/${SAMPLE}-${CHUNK}-genome-$READ_TYPE.sam
	  echo
	elif [ -f $CHUNK_DIR/sam-files/${SAMPLE}-${CHUNK}-genome-$READ_TYPE.sam.gz ]
	then
	  echo "Genome"
	  ls -l $CHUNK_DIR/sam-files/${SAMPLE}-${CHUNK}-genome-$READ_TYPE.sam.gz | fgrep "$DATE_FILTER"
	fi

	if [ "TRUE" = "FALSE" ]
	then
	if [ -f $CHUNK_DIR/sam-files/${SAMPLE}-${CHUNK}-junction-$READ_TYPE.sam ]
	then
	  tail -1 $CHUNK_DIR/sam-files/${SAMPLE}-${CHUNK}-junction-$READ_TYPE.sam
	  echo
	elif [ -f $CHUNK_DIR/sam-files/${SAMPLE}-${CHUNK}-junction-$READ_TYPE.sam.gz ]
	then
	  ls -l $CHUNK_DIR/sam-files/${SAMPLE}-${CHUNK}-junction-$READ_TYPE.sam.gz  | fgrep "$DATE_FILTER"
	fi
	fi

	echo

    fi
  done

  ## Change suffixes to sample level
  JOB_NAME_SUFFIX=$SAMPLE-$ALIGNER
  LOG_FILE_SUFFIX=$ALIGNER-$DATE.log

  ################################################################################
  ## Merge genomic alignments
  ################################################################################
  
  JOB_NAME=MG-$EQP_JOB_NAME_SUFFIX
  if [ "$MERGE_GENOME_ALIGNMENTS" = "TRUE" ]
  then
    if [ "$ALIGNER_TYPE" = "genomic" ]
    then
      OPERATION="merge-genome-alignments"
      CHUNK="none"
      showStatusAndSetSubmit
  	
      if [ "$SUBMIT" = "TRUE" ]
      then
  	clearLogFiles
  
  	if [ ! -d $SAMPLE_DIR/sam-files ]
  	then
  	  mkdir $SAMPLE_DIR/sam-files
  	fi

	CMD="$PROCESS_DIR/merge-bam.sh $PROJECT_DIR $SAMPLE_DIR/sam-files/$SAMPLE-$EQP_ALIGNER-genome.bam $GENOME_SAM_CHUNK_FILES"
  	
        submitJob $MEM_REQUESTED $JOB_NAME $LOG_DIR/$OPERATION-$LOG_FILE_SUFFIX "$CMD" 3 "$DEPENDENCY_STRING"
      fi
    else
      echo "merge-genome-alignments can only be called by splice-aware aligners and aligner $ALIGNER is not recognized as splice-aware."
      exit 1
    fi
  fi
  
  
  ################################################################################
  ## Combine chunk count files and collect sample count files
  ################################################################################

  if [ "$COMBINE_GENE_COUNTS" = "TRUE" ]
  then
    if [ "$CHUNK_GENE_COUNT_FILES" != "" ]
    then
      $UTIL_LIB_DIR/combineCounts.py -A -i $CHUNK_GENE_COUNT_FILES -o $SAMPLE_GENE_COUNT_FILE
      SAMPLE_GENE_COUNT_FILES="$SAMPLE_GENE_COUNT_FILES $SAMPLE_GENE_COUNT_FILE"
    fi
  fi
  
  if [ "$COMBINE_EXON_COUNTS" = "TRUE" ]
  then
    if [ "$CHUNK_EXON_COUNT_FILES" != "" ]
    then
      $UTIL_LIB_DIR/combineCounts.py -A -i $CHUNK_EXON_COUNT_FILES -o $SAMPLE_EXON_COUNT_FILE
      SAMPLE_EXON_COUNT_FILES="$SAMPLE_EXON_COUNT_FILES $SAMPLE_EXON_COUNT_FILE"
    fi
  fi
    
  if [ "$COMBINE_JUNCTION_COUNTS" = "TRUE" ]
  then
    if [ "$CHUNK_JUNCTION_COUNT_FILES" != "" ]
    then
      $UTIL_LIB_DIR/combineCounts.py -A -i $CHUNK_JUNCTION_COUNT_FILES -o $SAMPLE_JUNCTION_COUNT_FILE
      SAMPLE_JUNCTION_COUNT_FILES="$SAMPLE_JUNCTION_COUNT_FILES $SAMPLE_JUNCTION_COUNT_FILE"
      NUM_COMBINED_LINES=`cat $SAMPLE_JUNCTION_COUNT_FILE | wc -l`
      if [ "$NUM_COMBINED_LINES" = "0" ]
      then
        echo "No output generated for sample $SAMPLE"
	echo "Call: $UTIL_LIB_DIR/combineCounts.py -A -i $CHUNK_JUNCTION_COUNT_FILES -o $SAMPLE_JUNCTION_COUNT_FILE"
	exit 1
      fi
    fi
    ## $UTIL_LIB_DIR/combineCounts.py -A -i $CHUNK_JUNCTION_SAM_COUNT_FILES -o $SAMPLE_JUNCTION_SAM_COUNT_FILE
    ## SAMPLE_JUNCTION_SAM_COUNT_FILES="$SAMPLE_JUNCTION_COUNT_FILES $SAMPLE_JUNCTION_SAM_COUNT_FILE"
  fi

  if [ "$CLEAN_UP" = "TRUE" ]
  then
    $UTIL_LIB_DIR/empty-files.sh $SAMPLE_DIR/fastq-files/*
  fi
  
done

if [ "$USE_QSUB" = "FALSE" ]
then
  echo "Waiting for remaining processes to finish."
  wait
fi

################################################################################
## Combine sample count files
################################################################################

if [ "$COMBINE_GENE_COUNTS" = "TRUE" -o "$COMBINE_EXON_COUNTS" = "TRUE" -o "$COMBINE_JUNCTION_COUNTS" = "TRUE" ]
then
  PROJECT_COUNT_DIR=$SAMPLES_DIR/count-files
  if [ ! -d $PROJECT_COUNT_DIR ]
  then
    echo "Creating directory $PROJECT_COUNT_DIR"
    mkdir $PROJECT_COUNT_DIR
  fi
  PROJECT=`basename $PROJECT_DIR`
fi

if [ "$COMBINE_GENE_COUNTS" = "TRUE" ]
then
  if [ "$SAMPLE_GENE_COUNT_FILES" != "" ]
  then
    ## Combine gene count files
    PROJECT_GENE_COUNT_FILE=$PROJECT_COUNT_DIR/$PROJECT-$EQP_ALIGNER-$READ_TYPE-gene.$COUNT_TYPE
    $UTIL_LIB_DIR/combineCounts.py -i $SAMPLE_GENE_COUNT_FILES -o $PROJECT_GENE_COUNT_FILE -s $EQP_ALIGNER-$READ_TYPE-gene.$COUNT_TYPE
  fi
fi

if [ "$COMBINE_EXON_COUNTS" = "TRUE" ]
then
  if [ "$SAMPLE_EXON_COUNT_FILES" != "" ]
  then
    ## Combine exon count files
    PROJECT_EXON_COUNT_FILE=$PROJECT_COUNT_DIR/$PROJECT-$EQP_ALIGNER-$READ_TYPE-exon.$COUNT_TYPE
    $UTIL_LIB_DIR/combineCounts.py -i $SAMPLE_EXON_COUNT_FILES -o $PROJECT_EXON_COUNT_FILE -s $EQP_ALIGNER-$READ_TYPE-exon.$COUNT_TYPE
  fi
fi

if [ "$COMBINE_JUNCTION_COUNTS" = "TRUE" ]
then
  if [ "$SAMPLE_JUNCTION_COUNT_FILES" != "" ]
  then
    if [ ! -f $MAP_DIR/${FILE_BASE}_junction_exon.map.gz ]
    then
      echo "File $MAP_DIR/${FILE_BASE}_junction_exon.map.gz not found ... exiting."
      exit 1
    fi
    ## Combine junction count files
    PROJECT_JUNCTION_COUNT_FILE=$PROJECT_COUNT_DIR/$PROJECT-$EQP_ALIGNER-$READ_TYPE-junction.$COUNT_TYPE
    $UTIL_LIB_DIR/combineCounts.py -i $SAMPLE_JUNCTION_COUNT_FILES -o $PROJECT_JUNCTION_COUNT_FILE \
      -s $EQP_ALIGNER-$READ_TYPE-junction.$COUNT_TYPE # -I $MAP_DIR/${FILE_BASE}_junction_exon.map.gz 
  fi
fi

## Compute FPKM
if [ "$COMPUTE_FPKM" = "TRUE" ]
then
  PROJECT_GENE_FPKM_FILE=`echo $PROJECT_GENE_COUNT_FILE | sed -e 's/[.]cnt$/.fpkm/'`
  $UTIL_LIB_DIR/computeFpkm.py -c $PROJECT_GENE_COUNT_FILE -g $GTF_FILE -n $ALIGNED_READS_NUM_FILE -o $PROJECT_GENE_FPKM_FILE \
       -l $PROJECT_COUNT_DIR/EntrezGene-gene-lengths.txt
fi



