#!/usr/bin/env python

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
##  Script to extract the status of a job from the cluster queue and the
##  log files
##
################################################################################

import sys
import argparse
import re
import os.path
import subprocess

## execfile(os.path.join(os.environ["HOME"], "ngs", "pipelines", "exon-pipeline", "bin", "util-lib", "getJobStatus.py"))

################################################################################
##
## Define command line options
##
################################################################################

parser = argparse.ArgumentParser(description='Extract pipeline meta informations')
parser.add_argument('-s', dest="samplesDir", metavar="<Samples sub directory>", default="samples", help='Full path to sub directory of <project dir> containing the samples.')
parser.add_argument('-S', dest="sample", metavar="<Sample>", help='Name of the sample which is considered.')
parser.add_argument('-C', dest="chunk", default = "none", metavar="<Chunk>", help='Name of the chunk which is considered.')
parser.add_argument('-A', dest="aligner", default = "bowtie2", metavar="<Aligner>", help='Name of the aligner used [bowtie2].')
parser.add_argument('-J', dest="jobName", metavar="<Job name>", help='Name of the job.')
parser.add_argument('-O', dest="operation", metavar="<Operation>", help='Operation to be considered: one of align-genome, align-transcripts, align-junctions, compute-counts')
parser.add_argument('-E', dest="noError", action='store_true', help='Do not allow for any error messages in the log-files.')
parser.add_argument('-w', dest="warningsOn", action='store_true', help='Show warning messages.')


if len(sys.argv) > 1:
  args = parser.parse_args(sys.argv[1:])
else:
  inputArgs = []
  inputArgs.append("-s")
  inputArgs.append(os.path.join(os.environ["HOME"], "ngs", "RNA-seq", "2014", "WNT_YAPC_NG81", "samples"))
  inputArgs.append("-S")
  inputArgs.append("YAPC_DMSO_24h_15")
  inputArgs.append("-C")
  inputArgs.append("C004")
  inputArgs.append("-A")
  inputArgs.append("bowtie2")
  inputArgs.append("-O")
  inputArgs.append("compute-counts")
  inputArgs.append("-J")
  inputArgs.append("CC-YAPC_DMSO_24h_15-C004-bowtie2")
  inputArgs.append("-E")
  
  args = parser.parse_args(inputArgs)

samplesDir = args.samplesDir
sample     = args.sample
chunk      = args.chunk
aligner    = args.aligner
operation  = args.operation
jobName    = args.jobName
noError    = args.noError
warningsOn = args.warningsOn

if len(sys.argv) <= 1:
  print "Using default parameter:"
  print "samplesDir: " + samplesDir
  print "sample    : " + sample    
  print "chunk     : " + chunk     
  print "operation : " + operation 
  print "jobName   : " + jobName
  print "noError   : " + str(noError)


################################################################################
##
## Get the status for split jobs
##
################################################################################

def getSplitStatus (logFileList, fullSampleDir, status):
  
  logDate = {}
  for logFile in logFileList:
    if logFile != "":
      date  = re.sub(".*-([0-9]*-[0-9]*)[.]log[~]*$", "\\1", logFile)
      if date != logFile:
        if date > logDate:
          logDate = date
      else:
        raise Exception ("Log file " + logFile + " does not have right format to extract the date: .*-<yymmdd>-<hhmm>.log")

  try:
    logFilename = os.path.join(fullSampleDir, "log-files", logFile)
    # print "Reading file: " + logFilename
    logFile = open(logFilename)
    if status == "running":
      splitComplete = True
    else:
      splitComplete = False

    fastqFileRead = False
    for line in logFile:
      if "Error" in line or "ERROR" in line or "exception" in line.lower() or "exiting" in line.lower() or line.strip() == "null":
        splitComplete = False
        break
      if line.strip () == "Fastq input complete.":
        fastqFileRead = True
      if fastqFileRead and line.strip() == "Done":
        splitComplete = True
  
    logFile.close()
  
  except IOError, e:
    print logFile + " cannot be opened ... skipping"
    print "Error: " + str(e[0]) + " " + str(e[1])
  
  return  splitComplete


################################################################################
##
## Get the percent aligned reads and the alignment status
##
################################################################################

def getBowtie2AlignmentStatus (logFileList, fullSampleDir, chunks, status):
  
  chunkLogFile = {}
  chunkLogDate = {}
  for logFile in logFileList:
    if logFile != "":
      chunk = re.sub(".*-(C[0-9]*)-.*", "\\1", logFile)
      date  = re.sub(".*-([0-9]*-[0-9]*)[.]log[~]*$", "\\1", logFile)
      if chunk in chunkLogDate:
        if date > chunkLogDate[chunk]:
          chunkLogDate[chunk] = date
          chunkLogFile[chunk] = logFile
      else:
        chunkLogDate[chunk] = date
        chunkLogFile[chunk] = logFile

  numReads          = {}
  percentAligned    = {}
  numAlignedReads   = {}
  alignmentComplete = {}
  for chunk in chunks:
    if chunk in chunkLogFile:
      try:
        logFilename = os.path.join(fullSampleDir, "log-files", chunkLogFile[chunk])
        # print "Reading file: " + logFilename
        logFile = open(logFilename)
        if status[chunk] == "running":
          alignmentComplete[chunk] = True
        else:
          alignmentComplete[chunk] = False
        numLines = 0
        for line in logFile:
          # print "Line: " + line
          if "error" in line.lower() or "exception" in line.lower() or "exiting" in line.lower() or line.strip() == "null":
            alignmentComplete[chunk] = False
            break
          if "reads; of these:" in line:
            numReads[chunk] = int(re.sub("reads; of these:", "", line).strip())
          if "% overall alignment rate" in line:
            percentAligned[chunk] = float (re.sub("% overall alignment rate", "", line).strip())
          if "NUMBER_MAPPED_READS=" in line:
            numAlignedReads[chunk] = int(re.sub(".*NUMBER_MAPPED_READS=", "", line).strip())
          if line.strip()[:len("Alignment done")] == "Alignment done" or "successfully completed" in line:
            alignmentComplete[chunk] = True
          # print "alignmentComplete["+chunk+"]: " + str(alignmentComplete[chunk])
          numLines = +1
          
        if numLines == 0:
          print >> sys.stderr, "Warning: File " + chunkLogFile[chunk] + " is empty."
        logFile.close()
      
      except IOError, e:
        raise Exception (chunkLogFile[chunk] + " cannot be opened ... skipping\n" + \
                         "Unix error number: " + str(e[0]) + " and message: " + str(e[1]))
    
    if chunk in status and status[chunk] != "running" and not chunk in alignmentComplete:
      alignmentComplete[chunk] = False
    if chunk in numReads and chunk in percentAligned and not chunk in numAlignedReads:
      numAlignedReads[chunk] = int(numReads[chunk] * percentAligned[chunk])
    if not chunk in numReads:
      numReads[chunk] = ""
    if not chunk in percentAligned:
      percentAligned[chunk] = ""
    if not chunk in numAlignedReads:
      numAlignedReads[chunk] = ""
  
  return  numReads, percentAligned, alignmentComplete


################################################################################
##
## Get the Tophat2 alignment exit status counts
##
################################################################################

def getTophat2AlignmentStatus (logFileList, fullSampleDir, chunks, status):
  
  chunkLogFile = {}
  chunkLogDate = {}
  for logFile in logFileList:
    if logFile != "":
      chunk = re.sub(".*-(C[0-9]*)-.*", "\\1", logFile)
      date  = re.sub(".*-([0-9]*-[0-9]*)[.]log[~]*$", "\\1", logFile)
      if chunk in chunkLogDate:
        if date > chunkLogDate[chunk]:
          chunkLogDate[chunk] = date
          chunkLogFile[chunk] = logFile
      else:
        chunkLogDate[chunk] = date
        chunkLogFile[chunk] = logFile
  
  alignmentComplete  = {}
  for chunk in chunks:
    if chunk in chunkLogFile:
      try:
        logFilename = os.path.join(fullSampleDir, "log-files", chunkLogFile[chunk])
        # print "Reading file: " + logFilename
        logFile = open(logFilename)
        if status[chunk] == "running":
          alignmentComplete[chunk] = True
        else:
          alignmentComplete[chunk] = False
        for line in logFile:
          if "error" in line.lower() or "exception" in line.lower() or "exiting" in line.lower() or line.strip() == "null":
            alignmentComplete[chunk] = False
            break
          if line.strip() == "Tophat2 alignments successfully completed." or line.strip() == "Tophat2 alignment sucessfully completed.":
            alignmentComplete[chunk] = True
      
        logFile.close()
      
      except IOError, e:
        raise Exception (chunkLogFile[chunk] + " cannot be opened ... skipping\n" + \
                         "Unix error number: " + str(e[0]) + " and message: " + str(e[1]))
    
    if not chunk in alignmentComplete:
      alignmentComplete[chunk] = False
  
  return alignmentComplete


################################################################################
##
## Get the number of aligned and expressed reads as well as the compute counts
## exit status
##
################################################################################

def getComputeCountsStatus (logFileList, fullSampleDir, chunks, prefixList, noError, status):
  
  chunkLogFile = {}
  chunkLogDate = {}
  for logFile in logFileList:
    # print logFile
    if logFile != "":
      chunkFound = False
      for prefix in prefixList:
        chunk = re.sub("^" + prefix + ".*-(C[0-9]*)-.*", "\\1", logFile)
        date  = re.sub("^" + prefix + ".*-([0-9]*-[0-9]*)[.]log[~]*$", "\\1", logFile)
        if chunk != logFile and date != logFile:
          chunkFound = True
          if chunk in chunkLogDate:
            if date > chunkLogDate[chunk]:
              chunkLogDate[chunk] = date
              chunkLogFile[chunk] = logFile
          else:
            chunkLogDate[chunk] = date
            chunkLogFile[chunk] = logFile
          
      if not chunkFound:
        raise Exception ("The filename " + logFile + " not of the form " + "^" + "|".join(prefixList) + "-(C[0-9]*).*-([0-9]*-[0-9]*)[.]log[~]*")

  geneCountsStarted       = {}
  geneCountsCompleted     = {}
  junctionCountsStarted   = {}
  junctionCountsCompleted = {}
  numAlignedReads         = {}
  numExpressedReads       = {}
  operationComplete       = {}
  errorOccured            = {}
  for chunk in chunks:
    if chunk in chunkLogFile:
      try:
        logFilename = os.path.join(fullSampleDir, "log-files", chunkLogFile[chunk])
        ## print >> sys.stderr, "Reading file: " + logFilename
        logFile = open(logFilename)

        if status[chunk] == "running":
          operationComplete[chunk] = True
        else:
          operationComplete[chunk] = False
        errorOccured[chunk] = False
        numLines = 0
        for line in logFile:
          if "NUMBER_MAPPED_READS=" in line:
            numAlignedReads[chunk] = int(re.sub(".*NUMBER_MAPPED_READS=", "", line).strip())
          if "NUMBER_EXPRESSED_READS=" in line:
            numAlignedReads[chunk] = int(re.sub(".*NUMBER_EXPRESSED_READS=", "", line).strip())
          if "ERROR" in line.upper() or "Exception" in line or "exiting" in line or line.strip() == "null":
            errorOccured[chunk] = True
          if "ComputeGeneCountsSam -R" in line:
            junctionCountsStarted[chunk] = True
          if chunk in junctionCountsStarted and "Done" in line:
            junctionCountsCompleted[chunk] = True
          if "ComputeExonCounts -g -w" in line or "ComputeCounts -g -w" in line:
            geneCountsStarted[chunk] = True
          if chunk in geneCountsStarted and "Done" in line:
            geneCountsCompleted[chunk] = True
          if "successfully completed" in line:
            operationComplete[chunk] = True
          numLines = +1
          
        if numLines == 0:
          print >> sys.stderr, "Warning: File " + chunkLogFile[chunk] + " is empty."
        logFile.close()
      
      except IOError, e:
        raise Exception (chunkLogFile[chunk] + " cannot be opened ... skipping\n" + \
                         "Unix error number: " + str(e[0]) + " and message: " + str(e[1]))
    
    if not chunk in operationComplete:
      operationComplete[chunk] = False
    if not chunk in numAlignedReads:
      numAlignedReads[chunk] = ""
    if not chunk in numExpressedReads:
      numExpressedReads[chunk] = ""
        
    if noError and chunk in errorOccured and errorOccured[chunk]:
      operationComplete[chunk] = False
  
  return  numAlignedReads, numExpressedReads, operationComplete


################################################################################
##
## Get the exit status
##
################################################################################

def getJobExitStatus (logFileList, fullSampleDir, chunks, prefix, noError, status):

  chunkLogFile = {}
  chunkLogDate = {}
  for logFile in logFileList:
    ## print logFile
    if logFile != "":
      chunk = re.sub("^" + prefix + ".*-(C[0-9]*)-.*", "\\1", logFile)
      date  = re.sub("^" + prefix + ".*-([0-9]*-[0-9]*)[.]log[~]*$", "\\1", logFile)
      if chunk in chunkLogDate:
        if date > chunkLogDate[chunk]:
          chunkLogDate[chunk] = date
          chunkLogFile[chunk] = logFile
      else:
        chunkLogDate[chunk] = date
        chunkLogFile[chunk] = logFile

  operationComplete  = {}
  errorOccured       = {}
  for chunk in chunks:
    if chunk in chunkLogFile:
      try:
        logFilename = os.path.join(fullSampleDir, "log-files", chunkLogFile[chunk])
        ## print "Reading file: " + logFilename
        logFile = open(logFilename)
      except IOError, e:
        raise Exception(chunkLogFile[chunk] + " cannot be opened ... skipping \n" + \
                        "Unix error code and message: " + str(e[0]) + " " + str(e[1]))
        
      if status[chunk] == "running":
        operationComplete[chunk] = True
      else:
        operationComplete[chunk] = False
      errorOccured[chunk] = False
      numLines = 0
      for line in logFile:
        # print >> sys.stderr, line
        if "ERROR" in line.upper () or "Exception" in line or "exiting" in line or line.strip() == "null":
          errorOccured[chunk] = True
        if "successfully completed" in line.strip():
          operationComplete[chunk] = True
        # print >> sys.stderr, "operationComplete[" + chunk + "]: "+ str(operationComplete[chunk]) + ", errorOccured[" + chunk + "]: "+ str(errorOccured[chunk])
        numLines = +1
          
      if numLines == 0:
        print >> sys.stderr, "Warning: File " + chunkLogFile[chunk] + " is empty."
      logFile.close()
    
    if not chunk in operationComplete:
      operationComplete[chunk] = False


    # print >> sys.stderr, "errorOccured[" + chunk + "]: "+ str(errorOccured[chunk])
    if noError and chunk in errorOccured and errorOccured[chunk]:
      operationComplete[chunk] = False
  
  return operationComplete


################################################################################
##
## Read log file and extract fields
##
################################################################################

sampleDir = os.path.join(samplesDir, sample)
if not os.path.isdir(sampleDir):
  raise Exception ("ERROR: directory " + sampleDir + " not found.")
chunks = [chunk]

try: 
  cmd = " ".join(["qstat -r"])
  qstatInf = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True).stdout.read().strip().split("\n")
except IOError, e:
  raise Exception ("Could not execute qstat -r\n" +  "Unix error code and message: " + str(e[0]) + " " + str(e[1]))

jobStatus = {}
for chunk in chunks:
  jobStatus[chunk] = "not submitted"
  for line in qstatInf:
    fields = line.split()
    if not "job-ID" in line and len(fields) > 5:
     if jobName.startswith(fields[2]):
        ## print ", ".join([fields[2], fields[4]])
        if fields[4] == 'r':
          currentJobStatus = "running"
        elif "E" in fields[4]:
          currentJobStatus = "SGE-error"
        elif "d" in fields[4].lower():
          currentJobStatus = "in deletion"
        elif "h" in fields[4].lower():
          currentJobStatus = "on hold"
        elif "w" in fields[4].lower() or "t" in fields[4].lower():
          currentJobStatus = "waiting"
        else:
          currentJobStatus = "not submitted"
    elif "Full jobname:" in line:
      ## print line + ", " + str(jobName == fields[-1])
      if fields[-1] == jobName:
        jobStatus[chunk] = currentJobStatus
    ## print jobStatus[chunk]
  
  if jobStatus[chunk] == "not submitted" or jobStatus[chunk] == "in deletion" or jobStatus[chunk] == "running":
    
    if operation in ["split"]:
      logFilePrefix = operation
    else:
      logFilePrefix = operation + "-" + aligner

    chunkRegExp="-" + chunk
    if operation in ["split", "merge-genome-alignments"]:
      chunkRegExp=""

    cmd = " ".join(["ls -1", os.path.join(sampleDir, "log-files"), "| grep '^" + logFilePrefix + chunkRegExp + "-'"])
    logFileList = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True).stdout.read().strip().split("\n")
    
    if len(logFileList) > 1 or logFileList[0] != "":
      if operation == "split":
        #Process split log files
        operationComplete = {}
        operationComplete[chunk] = getSplitStatus (logFileList, sampleDir, jobStatus[chunk])
          
      elif operation in ["align-genome", "align-transcripts", "align-junctions"]:
        #Process genome log files    
        numReads, percentAligned, operationComplete = {}, {}, {}
        if aligner == "tophat2" and operation == "align-genome":
          operationComplete = getTophat2AlignmentStatus (logFileList, sampleDir, chunks, jobStatus)
        else:
          numReads, percentAligned, operationComplete = getBowtie2AlignmentStatus (logFileList, sampleDir, chunks, jobStatus)
        
      elif operation == "compute-counts" or operation == "compute-junction-counts":
        #Process compute counts log files
        countNumReads, countExpressedReads, operationComplete = {}, {}, {}
        countNumReads, countExpressedReads, operationComplete = getComputeCountsStatus (logFileList, sampleDir, chunks, \
                                                                                        ["counts", operation + "-" + aligner], noError, jobStatus)    
      else:
        if warningsOn and not operation in ["align-genome-eqp", "merge-genome-alignments"]:
          print >> sys.stderr, "Unknown operation: " + operation 
        operationComplete = {}
        operationComplete = getJobExitStatus (logFileList, sampleDir, chunks, operation, noError, jobStatus)
    
      if chunk in operationComplete:
        if operationComplete[chunk]:
          if jobStatus[chunk] != "running":
            jobStatus[chunk] = "complete"
        else:
          if jobStatus[chunk] != "running":
            jobStatus[chunk] = "failed"
          else:
            jobStatus[chunk] = "crashed"

  if operation in ["split", "merge-genome-alignments"]:
    print "\t".join([sample, operation, jobName, jobStatus[chunk]])
  else:
    print "\t".join([sample, chunk, operation, jobName, jobStatus[chunk]])
