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
##  Script to extract the meta information for all the samples in a project directory 
##
################################################################################

import sys
import argparse
import re
import os.path
import subprocess

## execfile(os.path.join(os.environ["HOME"], "ngs", "pipelines", "exon-pipeline", "bin", "util-lib", "extractPipelineMetaData.py"))

################################################################################
##
## Define command line options
##
################################################################################

parser = argparse.ArgumentParser(description='Extract pipeline meta informations')
parser.add_argument('-p', dest="projectDir", metavar="<Project directory>",  help='Directory containing all files of one project')
parser.add_argument('-d', dest="dataReleaseInfoDir", default="", metavar="<Data release info directory>", help='Directory containing the data release info files.')
parser.add_argument('-f', dest="dataReleaseInfoFile", default="", metavar="<Data release info file>", help='Data release info file.')
parser.add_argument('-s', dest="samplesDir", metavar="<Samples sub directory>", default="samples", help='Sub directory of <project dir> containing the samples [samples]')
parser.add_argument('-A', dest="aligner", metavar="STRING", default="bowtie2", help='The aligner used for alignment [bowtie2]')
parser.add_argument('-Q', dest="quantifier", metavar="STRING", default="eqp", help='The quantification program used [eqp]')
parser.add_argument('-I', dest="intersectionMode", metavar="STRING", default="union", help='The intersection mode used for htseq [union]')

if len(sys.argv) > 1:
  args = parser.parse_args(sys.argv[1:])
else:
  defaultProjectDir = os.path.join(os.environ["HOME"], "ngs", "RNA-seq-test", "SEQC-NG00008.1-EQP2.0-Ensembl")
  inputArgs = []
  inputArgs.append("-p")
  inputArgs.append(defaultProjectDir)
  inputArgs.append("-s")
  inputArgs.append("samples-unique")
  inputArgs.append("-f")
  inputArgs.append(os.path.join(defaultProjectDir, "annotation-files", "IL_NGS_11_067_120104_SN889_0095_AD0902ACXX_Data_Release_info_GA258.finished.txt"))
  inputArgs.append("-A")
  inputArgs.append("star")
  #inputArgs.append("-Q")
  #inputArgs.append("htseq")
  #inputArgs.append("-I")
  #inputArgs.append("intersection-nonempty")

  args = parser.parse_args(inputArgs)

projectDir          = args.projectDir
dataReleaseInfoDir  = args.dataReleaseInfoDir
dataReleaseInfoFile = args.dataReleaseInfoFile
samplesDir          = args.samplesDir
aligner             = args.aligner
quantifier          = args.quantifier
intersectionMode    = args.intersectionMode

print "Project dir: " + projectDir
print "Data release info file: " + dataReleaseInfoFile
print "Aligner: " + aligner

if dataReleaseInfoDir == "" and dataReleaseInfoFile == "":
  dataReleaseInfoDir = os.path.join (projectDir, "annotation-files")


################################################################################
##
## Compute the number of reads
##
################################################################################

def computeNumReads (sampleName, sampleRootDir):
  sampleDir = os.path.join (sampleRootDir, sampleName, "fastq-files")
  cmd = " ".join(["ls -1", sampleDir])
  sampleFileList = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True).stdout.read().strip().split("\n")
  
  numReads = 0
  for sampleFile in sampleFileList:
    if "_1.fq.gz" in sampleFile or "_1.fastq.gz" in sampleFile:
      print >> sys.stderr, "Counting the number of lines of " + sampleFile
      cmd = " ".join(["zcat", os.path.join(sampleDir, sampleFile), "| wc -l"])
      numReads += int(subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True).stdout.read().strip().split("\n")[0]) / 4
        
  return numReads
     

################################################################################
##
## Get the alignment parameter (program, version, and options)
##
################################################################################

def getAlignmentParameter (logFileList, fullLogDir):
  
  chunkLogFile = {}
  chunkLogDate = {}
  for logFile in logFileList:
    if logFile != "":
      chunk = re.sub(".*-(C[0-9]*)-.*", "\\1", logFile)
      date  = re.sub(".*-([0-9]*-[0-9]*)[.]log$", "\\1", logFile)
      if chunk in chunkLogDate:
        if date > chunkLogDate[chunk]:
          chunkLogDate[chunk] = date
          chunkLogFile[chunk] = logFile
      else:
        chunkLogDate[chunk] = date
        chunkLogFile[chunk] = logFile

  alignmentPrograms = []
  alignmentVersions = []
  alignmentOptions  = []
  for chunk in chunkLogFile:
    try:
      logFilename = os.path.join(fullLogDir, chunkLogFile[chunk])
      # print "Reading file: " + logFilename
      logFile = open(logFilename)
      for line in logFile:
        if line.startswith("ALIGNMENT_PROGRAM"):
          alignmentPrograms.append(re.sub("ALIGNMENT_PROGRAM" + "=", "", line.strip()))
        if line.startswith("ALIGNMENT_VERSION"):
          alignmentVersions.append(re.sub("ALIGNMENT_VERSION" + "=", "", line.strip()))
        if line.startswith("ALIGNMENT_OPTIONS"):
          alignmentOptions.append(re.sub("ALIGNMENT_OPTIONS" + "=", "", line.strip()).replace(" --sam-no-hd", ""))
      
      logFile.close()
      
    except IOError, e:
      raise Exception ("File " + chunkLogFile[chunk] + " not found ... skipping\n" + 
                       "Unix error: " + str(e[0]) + " " + str(e[1]))
  
  return alignmentPrograms, alignmentVersions, alignmentOptions


################################################################################
##
## Get the htseq parameter (program, version, and options)
##
################################################################################

def getHtSeqParameter (logFileList, fullSampleDir, sample, htSeqKeyWords, dexSeqKeyWords):

  chunkLogFile = {}
  chunkLogDate = {}
  for logFilename in logFileList:
    if logFilename != "":
      chunk = re.sub(".*-(C[0-9]*)-.*", "\\1", logFilename)
      date  = re.sub(".*-([0-9]*-[0-9]*)[.]log$", "\\1", logFilename)
      if chunk in chunkLogDate:
        if date > chunkLogDate[chunk]:
          chunkLogDate[chunk] = date
          chunkLogFile[chunk] = logFilename
      else:
        chunkLogDate[chunk] = date
        chunkLogFile[chunk] = logFilename
      
  htseqSummaryStarted  = False
  dexseqSummaryStarted = False
  dexseqStarted        = False
  modeFound            = False
  
  countMode   = []
  htSeqCounts = {}
  for keyWord in htSeqKeyWords:
    htSeqCounts[keyWord] = 0

  dexSeqCounts = {}
  for keyWord in dexSeqKeyWords:
    dexSeqCounts[keyWord] = 0
  
  chunkIndex = 0
  for chunk in chunkLogFile:
    try:
      logFilename = os.path.join(fullSampleDir, "log-files", chunkLogFile[chunk])
      # print "Reading file: " + logFilename
      logFile = open (logFilename)
    except IOError, e:
      raise Exception("File " + logFilename + " not found ... exiting.\n" + 
                      "Unix error message: " + str(e[1]) + " and code: " + str(e[0]))

    transcriptJunctionWeightsStarted = False
    for line in logFile:
      if not modeFound and "htseq-count -m" in line:
        modeFound = True
        index = line.index("htseq-count -m")
        line = line[(index + len("htseq-count -m") + 1):]
        # print "Count mode: " + line.split()[0]
        countMode.append(line.split()[0])
      elif line.startswith("HTSeq summary:"):
        if not dexseqStarted:
          htseqSummaryStarted = True
        else:
          dexseqSummaryStarted = True
      elif htseqSummaryStarted and line.startswith("__"):
        keyWord, value = line.split("\t")
        if keyWord in htSeqKeyWords:
          htSeqCounts[keyWord] += int(value)
        else:
          raise Exception ("Unknown HtSeq keyword " + keyWord + " in file " + logFilename)
      elif modeFound and htseqSummaryStarted and "dexseq_count.py" in line:
        htseqSummaryStarted = False
        dexseqStarted       = True
      elif dexseqSummaryStarted and line.startswith("__"):
        keyWord, value = line.split("\t")
        if keyWord in dexSeqKeyWords:
          dexSeqCounts[keyWord] += int(value)
        else:
          raise Exception ("Unknown Dexseq keyword " + keyWord + " in file " + logFilename)
    
    logFile.close()
    
  return htSeqCounts, dexSeqCounts




################################################################################
##
## Get the quantification parameter (program, version, and options)
##
################################################################################

def getQuantificationParameter (logFileList, fullSampleDir, sample):
  
  chunkLogFile = {}
  chunkLogDate = {}
  for logFilename in logFileList:
    if logFilename != "":
      chunk = re.sub(".*-(C[0-9]*)-.*", "\\1", logFilename)
      date  = re.sub(".*-([0-9]*-[0-9]*)[.]log$", "\\1", logFilename)
      if chunk in chunkLogDate:
        if date > chunkLogDate[chunk]:
          chunkLogDate[chunk] = date
          chunkLogFile[chunk] = logFilename
      else:
        chunkLogDate[chunk] = date
        chunkLogFile[chunk] = logFilename
  
  readWeightThreshold  = []
  exonOverlap          = []
  junctionOverlap      = []
  numberMappedReads    = []
  numberExpressedReads = []
  
  numberExpressedReadsTransJunc = []
  numberExpressedReadsCount     = []

  chunkIndex = 0
  for chunk in chunkLogFile:
    numExpressed = 0
    try:
      logFilename = os.path.join(fullSampleDir, "log-files", chunkLogFile[chunk])
      # print "Reading file: " + logFilename
      logFile = open (logFilename)
    except IOError, e:
      raise Exception(chunkLogFile[chunk] + " not found ... exiting.\n" + 
                      "Unix error message: " + str(e[1]) + " and code: " + str(e[0]))
    
    transcriptJunctionWeightsStarted = False
    for line in logFile:
      if line.startswith("READ_WEIGHT_THRESHOLD"):
        readWeightThreshold.append(re.sub("READ_WEIGHT_THRESHOLD" + "=", "", line.strip()))
      if line.startswith("EXON_OVERLAP"):
        exonOverlap.append(re.sub("EXON_OVERLAP" + "=", "", line.strip()))
      if line.startswith("JUNCTION_OVERLAP"):
        junctionOverlap.append(re.sub("JUNCTION_OVERLAP" + "=", "", line.strip()))
      if line.startswith("NUMBER_MAPPED_READS"):
        numberMappedReads.append(int(re.sub(".*NUMBER_MAPPED_READS" + "=", "", line.strip())))
      if line.startswith("NUMBER_EXPRESSED_READS"):
        numberExpressedReads.append(int(re.sub(".*NUMBER_EXPRESSED_READS" + "=", "", line.strip())))
      if line.startswith("Creating transcript junction SAM weights"):
        transcriptJunctionWeightsStarted = True
      if transcriptJunctionWeightsStarted and "fragments are mapped" in line:
        transcriptJunctionWeightsStarted = False
        numberExpressedReadsTransJunc.append(int(line.split(" ")[0]))
    
    logFile.close()

    if len(numberExpressedReads) < chunkIndex + 1 and len(numberExpressedReadsTransJunc) == chunkIndex + 1:
      print >> sys.stderr, "NUM_EXPRESSED_READS field not found in " + chunkLogFile[chunk]
      ## numberExpressedReads.append(numberExpressedReadsTransJunc[chunkIndex])
    ## print "Number of expressed reads for " + logFilename + ": " + str(numberExpressedReads)
    if len(numberExpressedReads) < chunkIndex + 1:
      chunkGeneCountFilename = os.path.join(fullSampleDir, chunk, "count-files", "-".join([sample, chunk, "combined-mixed-gene.cnt"]))
      print >> sys.stderr, "Reading " + chunkGeneCountFilename
      try:
        chunkGeneCountFile = open (chunkGeneCountFilename)
      except IOError, e:
        raise Exception(chunkGeneCountFilename + " not found ... exiting.\n" + 
                        "Unix error message: " + str(e[1]) + " and code: " + str(e[0]))

      sumCounts = 0
      for line in chunkGeneCountFile:
        sumCounts += float(line.strip().split("\t")[1])
      numberExpressedReads.append(int(sumCounts))

    chunkIndex += 1

  return readWeightThreshold, exonOverlap, junctionOverlap, numberMappedReads, numberExpressedReads


################################################################################
##
## Read SMF-Id sample name mapping file
##
################################################################################

def readDataReleaseInfoDir (dataReleaseInfoDir, dataReleaseInfoFileOption):
  
  sampleIdMap   = {}
  numReadsTotal = {}
  
  try:
    if dataReleaseInfoDir != "":
      print "Using files in directory: " + dataReleaseInfoDir
      dataReleaseInfoFilenames = [filename for filename in os.listdir (dataReleaseInfoDir) if "Data_Release_info" in filename and filename.endswith("txt")]
    else:
      dataReleaseInfoFilenames = [dataReleaseInfoFileOption]
    
    for dataReleaseInfoFilename in dataReleaseInfoFilenames:
      if "Data_Release_info" in dataReleaseInfoFilename or dataReleaseInfoFilename == dataReleaseInfoFileOption:
        print "Reading Data Release info file: " + dataReleaseInfoFilename
        dataReleaseInfoFile = open(os.path.join(dataReleaseInfoDir, dataReleaseInfoFilename))
        
        i = 0
        for line in dataReleaseInfoFile:
          if i == 0:
            header = line.lower().strip().split("\t")
            sampleNameIndex = -1
            if "sample name" in header:
              sampleNameIndex = header.index("sample name")
            smfIdIndex = -1
            if "novartis tracking id" in header:
              smfIdIndex = header.index ("novartis tracking id")
            numReadsTotalIndex = -1
            if "number of reads" in header:
              numReadsTotalIndex = header.index ("number of reads")
            if sampleNameIndex == -1 and smfIdIndex == -1:
              raise Exception ("ERROR: Neither Sample name nor Novartis Tracking Id columns not found.")
            if sampleNameIndex == -1:
              sampleNameIndex == smfIdIndex
            if numReadsTotalIndex == -1:
              print >> sys.stderr, "Number of reads column not found ... calculating the number of reads."
          else:
            fields = line.split("\t")
            sampleName = fields[sampleNameIndex]
            if smfIdIndex != -1:
              sampleIdMap[sampleName] = fields[smfIdIndex]
            else:
              sampleIdMap[sampleName] = "unknown"
            if not sampleName in numReadsTotal:
              numReadsTotal[sampleName] = 0
            if numReadsTotalIndex != -1:
              numReads = int (fields[numReadsTotalIndex])
              numReadsTotal[sampleName] = numReadsTotal[sampleName] + numReads
              ## print >> sys.stderr, "Adding " + str(numReads) + " to " + sampleName + " yielding a total of: " + str(numReadsTotal[sampleName])
            else:
              numReadsTotal[sampleName] = -1
          i = i + 1
        
        dataReleaseInfoFile.close()
        
  except IOError, e:
    raise Exception ("File " + dataReleaseInfoFilename + " not found ... skipping\n" + 
                     "Unix error: " + str(e[0]) + " " + str(e[1]))
  
  return sampleIdMap, numReadsTotal


################################################################################
##
## Read pipeline meta data
##
################################################################################

try:
  pipelineMetaDataFilename = os.path.join(projectDir, "exon-pipeline-files", "metadata-files", "pipeline-metadata.txt")
  print "Reading file: " + pipelineMetaDataFilename
  pipelineMetaDataFile = open(pipelineMetaDataFilename)

  keys = []
  values = []
  for line in pipelineMetaDataFile:
    key, value = line.strip().split("=")
    if value == "EQP":
      value = "EQP1.0"
    keys.append(key)
    values.append(value)
    
except IOError, e:
  raise Exception ("File " + pipelineMetaDataFilename + " not found ... skipping\n" + 
                   "Unix error: " + str(e[0]) + " " + str(e[1]))



################################################################################
##
## Read sample directories
##
################################################################################
  
sampleIdMap, numReadsTotal = readDataReleaseInfoDir (dataReleaseInfoDir, dataReleaseInfoFile)

sampleRootDir = os.path.join(projectDir, samplesDir)
cmd = " ".join(["ls -1", sampleRootDir])
sampleDirList = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True).stdout.read()


################################################################################
##
## Read log files and extract fields
##
################################################################################

try:
  sampleMetaDataFilename = os.path.join(sampleRootDir, "sample-metadata.txt")
  print "Writing meta data to file: " + sampleMetaDataFilename
  sampleMetaDataFile = open(sampleMetaDataFilename, 'w')
  
  headings = ["Sample name", "SMF-ID", "GENOME_ALIGNMENT_TOOL", "GENOME_ALIGNMENT_VERSION", "GENOME_ALIGNMENT_PARAMS", "TRANSCRIPT_ALIGNMENT_TOOL", "TRANSCRIPT_ALIGNMENT_VERSION", \
              "TRANSCRIPT_ALIGNMENT_PARAMS", "JUNCTION_ALIGNMENT_TOOL", "JUNCTION_ALIGNMENT_VERSION", "JUNCTION_ALIGNMENT_PARAMS", "READ_WEIGHT_THRESHOLD", \
              "EXON_OVERLAP", "JUNCTION_OVERLAP", "NUM_MAPPED_READS"] + keys
  heading = "\t".join(headings)
  print >> sampleMetaDataFile, heading

  htSeqKeyWords = ["__no_feature", "__ambiguous", "__too_low_aQual", "__not_aligned", "__alignment_not_unique", "__alignment_not_unique_expr", "__counted"]
  dexSeqKeyWords = ['__empty', '__ambiguous', '__lowaqual', '__notaligned', '__ambiguous_readpair_position', "__alignment_not_unique", "__counted"]

  numberMappedReads    = {}
  numberExpressedReads = {}
  htSeqCounts = {}
  dexSeqCounts = {}
  for sampleName in sampleDirList.split("\n"):
    sampleDir = os.path.join(sampleRootDir, sampleName)
    if os.path.isdir(sampleDir) and not "files" in sampleName and not sampleName == "" and os.path.isdir(os.path.join(sampleDir, "log-files")):
      print "Processing sample " + sampleName

      if quantifier == "eqp":
      
        #Process genome log files
        cmd = " ".join(["ls -1", os.path.join(sampleDir, "log-files"), "| grep 'align-genome-" + aligner + "-C'"])
        genomeLogFileList = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True).stdout.read().strip().split("\n")

        genomeAlignmentPrograms = []
        if len(genomeLogFileList) > 1 or genomeLogFileList[0] != "":
          genomeAlignmentPrograms, genomeAlignmentVersions, genomeAlignmentOptions = getAlignmentParameter (genomeLogFileList, os.path.join(sampleDir, "log-files"))
          
          if len(genomeAlignmentPrograms) > 1:
            genomeAlignmentPrograms, genomeAlignmentVersions, genomeAlignmentOptions = \
                zip(*[list (y) for y in [x for x in set(zip(genomeAlignmentPrograms, genomeAlignmentVersions, genomeAlignmentOptions))]])
          if len(genomeAlignmentPrograms) > 1:
            print >> sys.stderr, "Different genome alignment options used for different chunks of the same sample ... exiting"
            print >> sys.stderr, "Genome alignment program: " + "; ".join(genomeAlignmentPrograms)
            print >> sys.stderr, "Genome alignment version: " + "; ".join(genomeAlignmentVersions)
            print >> sys.stderr, "Genome alignment options: " + "; ".join(genomeAlignmentOptions)
        if len(genomeAlignmentPrograms) == 0:
          genomeAlignmentPrograms, genomeAlignmentVersions, genomeAlignmentOptions = [["undefined"], ["undefined"], ["undefined"]]
    
        #Process transcript log files
        cmd = " ".join(["ls -1", os.path.join(sampleDir, "log-files"), "| grep 'align-transcripts-" + aligner + "-C'"])
        transcriptLogFileList = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True).stdout.read().strip().split("\n")
        
        if len(transcriptLogFileList) > 1 or transcriptLogFileList[0] != "":
          transcriptAlignmentPrograms, transcriptAlignmentVersions, transcriptAlignmentOptions = \
                                       getAlignmentParameter (transcriptLogFileList, os.path.join(sampleDir, "log-files"))
          
          if len(transcriptAlignmentPrograms) > 1:
            transcriptAlignmentPrograms, transcriptAlignmentVersions, transcriptAlignmentOptions = \
                zip(*[list (y) for y in [x for x in set(zip(transcriptAlignmentPrograms, transcriptAlignmentVersions, transcriptAlignmentOptions))]])
          if len(transcriptAlignmentPrograms) > 1:
            print >> sys.stderr, "Different transcript alignment options used for different chunks of the same sample ... exiting"
            print >> sys.stderr, "Transcript alignment program: " + "; ".join(transcriptAlignmentPrograms)
            print >> sys.stderr, "Transcript alignment version: " + "; ".join(transcriptAlignmentVersions)
            print >> sys.stderr, "Transcript alignment options: " + "; ".join(transcriptAlignmentOptions)
        else:
          if aligner == "bowtie2":
            print >> sys.stderr, "No transcript alignment log file found for sample " + sampleName
          transcriptAlignmentPrograms, transcriptAlignmentVersions, transcriptAlignmentOptions = [["undefined"], ["undefined"], ["undefined"]]
    
        #Process junction log files
        cmd = " ".join(["ls -1", os.path.join(sampleDir, "log-files"), "| grep 'align-junctions-" + aligner + "-C'"])
        junctionLogFileList = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True).stdout.read().strip().split("\n")
        
        if len(junctionLogFileList) > 1 or junctionLogFileList[0] != "":
          junctionAlignmentPrograms, junctionAlignmentVersions, junctionAlignmentOptions = \
                                     getAlignmentParameter (junctionLogFileList, os.path.join(sampleDir, "log-files"))
          
          if len(junctionAlignmentPrograms) > 1:
            junctionAlignmentPrograms, junctionAlignmentVersions, junctionAlignmentOptions = \
              zip(*[list (y) for y in [x for x in set(zip(junctionAlignmentPrograms, junctionAlignmentVersions, junctionAlignmentOptions))]])
          if len(junctionAlignmentPrograms) > 1:
            print >> sys.stderr, "Different junction alignment options used for different chunks of the same sample ... exiting"
            print >> sys.stderr, "Junction alignment program: " + "; ".join(junctionAlignmentPrograms)
            print >> sys.stderr, "Junction alignment version: " + "; ".join(junctionAlignmentVersions)
            print >> sys.stderr, "Junction alignment options: " + "; ".join(junctionAlignmentOptions)
        else:
          if aligner == "bowtie2":
            print >> sys.stderr, "No junction alignment log file found for sample " + sampleName
          junctionAlignmentPrograms, junctionAlignmentVersions, junctionAlignmentOptions = [["undefined"], ["undefined"], ["undefined"]]
    
        #Process count log files
        countLogFilePrefix = "compute-counts"
          
        cmd = " ".join(["ls -1", os.path.join(sampleDir, "log-files"), "| grep " + "-".join(["'^" + countLogFilePrefix, aligner, "C'"])])
        countLogFileList = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True).stdout.read().strip().split("\n")
        
        readWeightThreshold = []
        exonOverlap = []
        junctionOverlap = []
        if len(countLogFileList) > 1 or countLogFileList[0] != "":
          readWeightThreshold, exonOverlap, junctionOverlap, numberMappedReadsChunks, numberExpressedReadsChunks = \
             getQuantificationParameter (countLogFileList, sampleDir, sampleName)
          
          if len(readWeightThreshold) > 1:
            readWeightThreshold, exonOverlap, junctionOverlap = \
              zip(*[list (y) for y in [x for x in set(zip(readWeightThreshold, exonOverlap, junctionOverlap))]])
          if len(readWeightThreshold) > 1:
            print >> sys.stderr, "Different count options used for different chunks of the same sample ... exiting"
            print >> sys.stderr, "Read weight threshold: " + "; ".join(readWeightThreshold)
            print >> sys.stderr, "Exon count:            " + "; ".join(exonOverlap)
            print >> sys.stderr, "Junction count:        " + "; ".join(junctionOverlap)
            
          numberMappedReads[sampleName]    = sum(map(int, numberMappedReadsChunks))
          numberExpressedReads[sampleName] = sum(map(int, numberExpressedReadsChunks))
          
        numSamples = len(genomeAlignmentPrograms)
        if len(readWeightThreshold) == 0:
          readWeightThreshold = [-1] * numSamples
        if len(exonOverlap) == 0:
          exonOverlap = [-1] * numSamples
        if len(junctionOverlap) == 0:
          junctionOverlap = [-1] * numSamples
          
        parameters = list(zip(*map(list, [genomeAlignmentPrograms, genomeAlignmentVersions, genomeAlignmentOptions, transcriptAlignmentPrograms, \
                                          transcriptAlignmentVersions, transcriptAlignmentOptions, junctionAlignmentPrograms, junctionAlignmentVersions, \
                                          junctionAlignmentOptions, readWeightThreshold, exonOverlap, junctionOverlap]))[0])
  
        sampleId = sampleName
        if sampleName in sampleIdMap:
          sampleId = sampleIdMap[sampleName]
        
        print >> sampleMetaDataFile, "\t".join([sampleName, sampleId] + parameters + [str(numberMappedReads[sampleName])] + values)
        
      elif quantifier == "htseq":          
        countLogFilePrefix = quantifier + "-" + intersectionMode
        cmd = " ".join(["ls -1", os.path.join(sampleDir, "log-files"), "| grep " + "-".join(["'^" + countLogFilePrefix, aligner, "C'"])])
        countLogFileList = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True).stdout.read().strip().split("\n")

        htSeqCountsSample, dexSeqCountsSample = getHtSeqParameter (countLogFileList, sampleDir, sampleName, htSeqKeyWords, dexSeqKeyWords)

        htSeqCounts[sampleName]= htSeqCountsSample
        if sum([dexSeqCountsSample[a] for a in dexSeqCountsSample]) > 0:
          dexSeqCounts[sampleName] = dexSeqCountsSample
        numberExpressedReads[sampleName] = sum([htSeqCountsSample[keyword] for keyword in ["__counted", "__ambiguous", "__alignment_not_unique_expr"]])
        numberMappedReads[sampleName]    = sum([htSeqCountsSample[keyword] for keyword in set(htSeqKeyWords) - set(["__not_aligned"])])
        readWeightThreshold, exonOverlap, junctionOverlap = [["-1"], ["-1"], ["-1"]]
        
      else:
        readWeightThreshold, exonOverlap, junctionOverlap = [["-1"], ["-1"], ["-1"]]
        numberMappedReads[sampleName]    = 0
        numberExpressedReads[sampleName] = 0
  
  sampleMetaDataFile.close()
      
except IOError, e:
  raise Exception ("File " + sampleMetaDataFilename + " not found ... skipping\n" + 
                   "Unix error: " + str(e[0]) + " " + str(e[1]))



################################################################################
##
## Print sample alignment statistics
##
################################################################################

try:
  sampleAlignmentStatisticsFilename = os.path.join(projectDir, "statistic-files", "-".join([samplesDir, aligner, quantifier, "alignment-statistics.txt"]))
  print "Writing alignment statistics to file: " + sampleAlignmentStatisticsFilename
  sampleAlignmentStatisticsFile = open(sampleAlignmentStatisticsFilename, 'w')
  
  headings = ["Sample name", "SMF-ID", "Number reads", "Number mapped reads", "% mapped reads", "Number expressed reads", "% expressed reads" ]
  print >> sampleAlignmentStatisticsFile, "\t".join(headings)
  
  for sampleName in sorted(numberMappedReads.keys()):
    smfId = "unknown"
    if sampleName in sampleIdMap:
      smfId = sampleIdMap[sampleName]
    
    if sampleName in numReadsTotal:
      numReads = numReadsTotal[sampleName]
    else:
      print >> sys.stderr, "Computing the number of reads for sample " + sampleName
      numReads = computeNumReads (sampleName, sampleRootDir)

    if numReads == 0:
      raise Exception ("No number of reads found for sample " + sampleName)
    
    print >> sampleAlignmentStatisticsFile, "\t".join (map (str, [sampleName, smfId, numReads, numberMappedReads[sampleName], \
                                                            int(100.0 * numberMappedReads[sampleName] / numReads), numberExpressedReads[sampleName], \
                                                            int(100.0 * numberExpressedReads[sampleName] / numReads)]))
    
  
  sampleAlignmentStatisticsFile.close()
      
except IOError, e:
  raise Exception ("File " + sampleAlignmentStatisticsFilename + " not found ... skipping\n" + 
                   "Unix error: " + str(e[0]) + " " + str(e[1]))


################################################################################
##
## Print htseq statistics
##
################################################################################

if len(htSeqCounts) > 0:
  
  try:
    sampleHtseqGenestatisticsFilename = os.path.join(projectDir, "statistic-files", "-".join([aligner, quantifier, intersectionMode, "gene-statistics.txt"]))
    sampleHtseqGenestatisticsFile = open(sampleHtseqGenestatisticsFilename, 'w')
  except IOError, e:
    raise Exception (sampleHtseqGenestatisticsFilename + " not found ... skipping\n" + \
                     "Unix error: " + str(e[0]) + " " + str(e[1]))
  
  print "Writing htseq gene statistics to file: " + sampleHtseqGenestatisticsFilename
  
  headings = ["Sample name"] + htSeqKeyWords
  print >> sampleHtseqGenestatisticsFile, "\t".join(headings)
  
  for sampleName in sorted(htSeqCounts.keys()):
    print >> sampleHtseqGenestatisticsFile, "\t".join (map (str, [sampleName] + [htSeqCounts[sampleName][keyword] for keyword in htSeqKeyWords]))
    
  sampleHtseqGenestatisticsFile.close()


################################################################################
##
## Print dexseq statistics
##
################################################################################

if len(dexSeqCounts) > 0:
  
  try:
    sampleDexseqExonstatisticsFilename = os.path.join(projectDir, "statistic-files", "-".join([aligner, quantifier, intersectionMode, "exon-statistics.txt"]))
    sampleDexseqExonstatisticsFile = open(sampleDexseqExonstatisticsFilename, 'w')
  except IOError, e:
    raise Exception (sampleDexseqExonstatisticsFilename + " not found ... skipping\n" + \
                     "Unix error: " + str(e[0]) + " " + str(e[1]))
  
  print "Writing dexseq exon statistics to file: " + sampleDexseqExonstatisticsFilename
  
  headings = ["Sample name"] + dexSeqKeyWords
  print >> sampleDexseqExonstatisticsFile, "\t".join(headings)
  
  for sampleName in sorted(dexSeqCounts.keys()):
    print >> sampleDexseqExonstatisticsFile, "\t".join (map (str, [sampleName] + [dexSeqCounts[sampleName][keyword] for keyword in dexSeqKeyWords]))
    
  sampleDexseqExonstatisticsFile.close()
      
