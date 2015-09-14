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
##  Script to convert gene counts to FPKM
##
################################################################################

import sys
import argparse
import re
import os.path
import numpy

## execfile(os.path.join(os.environ["HOME"], "ngs", "pipelines", "exon-pipeline", "bin", "util-lib", "computeFpkm.py"))


################################################################################
##
## Define command line options
##
################################################################################

parser = argparse.ArgumentParser(description='Combine the count files from different')
parser.add_argument('-d', type=int, default=0, dest="debugLevel", metavar="INT", help='debug level [0, no debug output]')
parser.add_argument('-c', dest="countFile", metavar="STRING", help='File with the counts in a column based format')
parser.add_argument('-n', dest="readNumFile", default="", metavar="STRING",
                    help='File the number of aligned reads for each sample (if third column exists, it is aggregated by first column)')
parser.add_argument('-g', dest="gtfFile", default="", metavar="STRING",  help='GTF file with the exon coordinates')
parser.add_argument('-o', dest="fpkmFile", metavar="STRING", help='File with the FPKM values (output)')
parser.add_argument('-l', dest="geneLengthFile", metavar="STRING", help='File with the lengths of genes')

parser.add_argument('-m', dest="computeMean", action="store_true", help='Write the mean of the control samples (see option -C)')
parser.add_argument('-C', dest="controlSamples", default="", metavar="STRING", help='Control samples (":" separated); if omitted all samples are considered to be control samples')


if len(sys.argv) > 1:
  args = parser.parse_args(sys.argv[1:])
else:
  project = "IPC-298"
  projectDir = os.path.join(os.environ["HOME"], "ngs", "RNA-seq", "2014", "Cell_lines", "human", project)
  inputArgs = []
  inputArgs.append("-c")
  inputArgs.append(os.path.join(projectDir, "samples", "count-files", project + "-mixed-gene.cnt"))
  ## inputArgs.append("-n")
  ## inputArgs.append(os.path.join(projectDir, "statistic-files", "bowtie2-aligned-read-num.txt"))
  inputArgs.append("-g")
  #inputArgs.append(os.path.join(projectDir, "exon-pipeline-files", "gtf-files", "ensembl_rna_hs.gtf"))
  inputArgs.append("-o")
  #inputArgs.append(os.path.join(projectDir, "samples", "count-files", project + "-mixed-gene.fpkm"))
  inputArgs.append(os.path.join(os.environ["HOME"], project + "-mixed-gene.fpkm"))
  inputArgs.append("-l")
  inputArgs.append(os.path.join(projectDir, "exon-pipeline-files", "gtf-files", "ensembl_rna_hs-gene-lengths.txt"))

  
  args = parser.parse_args(inputArgs)
  print "Using default arguments"

debugLevel         = args.debugLevel
countFilename      = args.countFile
readNumFilename    = args.readNumFile
gtfFilename        = args.gtfFile
fpkmFilename       = args.fpkmFile
geneLengthFilename = args.geneLengthFile

computeMean    = args.computeMean

controlSamples = args.controlSamples


################################################################################
##
## readCountFile
##
################################################################################

def readCountFile (countFilename):
  counts = {}
  countValues = ["dummy"]
  colHeaders = [""]
  
  try:
    countFile = open(countFilename)
    print >> sys.stderr, "Opening file: " + countFilename
    for line in countFile:
      line = line.rstrip ()
      countValues = line.split("\t")
      countObjectId = countValues[0]
      if countObjectId != "Id":
        counts[countObjectId] = [float(v) for v in countValues[1:]]
      else:
        colHeaders = countValues[1:]
  
    countFile.close()
  
  except IOError, e:
    print countFilename + " cannot be opened ... skipping"
    print "Error: " + str(e[0]) + " " + str(e[1])
  
  return [counts, colHeaders]


################################################################################
##
## readNumReadsFile
##
################################################################################

def readNumReadsFile (numReadsFilename):
  numReads = {}
  
  try:
    numReadsFile = open(numReadsFilename)
    print >> sys.stderr, "Opening file: " + numReadsFilename
    for line in numReadsFile:
      lineFields = line.rstrip ().split("\t")
      if len(lineFields) >= 3:
        sample, chunk, alignedReads = lineFields
      else:
        sample, alignedReads = lineFields
      if not "Sample" in sample:
        if sample in numReads:
          numReads[sample] = numReads[sample] + int(alignedReads)
        else:
          numReads[sample] = int(alignedReads)
  
    numReadsFile.close()
  
  except IOError, e:
    print numReadsFilename + " cannot be opened ... skipping"
    print "Error: " + str(e[0]) + " " + str(e[1])
  
  return numReads


################################################################################
##
## readCountFile
##
################################################################################

def readGeneLengthFile (geneLengthFilename):
  
  geneLengths = {}
  
  try:
    geneLengthFile = open(geneLengthFilename)
    print >> sys.stderr, "Opening file: " + geneLengthFilename
  except IOError, e:
    raise Exception(geneLengthFilename + " cannot be opened ... skipping\n" + \
                    "Unix error code: " + str(e[0]) + ", message: " + str(e[1]))
                    
  for line in geneLengthFile:
    line = line.rstrip ()
    geneLengthValues = line.split("\t")
    geneLengthObjectId = geneLengthValues[0]
    if geneLengthObjectId != "Gene Id":
      geneLengths[geneLengthObjectId] = int(geneLengthValues[1])

  geneLengthFile.close()
  
  return geneLengths


################################################################################
##
## readGtfFile
##
################################################################################

def readGtfFile (gtfFilename):
  countObjects = {}
 
  try:
    gtfFile = open(gtfFilename)
    print >> sys.stderr, "Opening file: " + gtfFilename
  except IOError, e:
    raise Exception(gtfFilename + " cannot be opened ... skipping\n" + \
                    "Unix error code: " + str(e[0]) + " and message " + str(e[1]))

  lineNum = 0
  for line in gtfFile:
    line = line.rstrip ()
    gtfEntries = line.split("\t")
    
    chromosome = gtfEntries[0]
    gtfType    = gtfEntries[2]
    interval   = [int(gtfEntries[3]), int(gtfEntries[4])]
    
    annotationEntries = gtfEntries[8].split(";")
    for annotationEntry in annotationEntries:
      if annotationEntry != "":
        annotationType, annotationValue = annotationEntry.strip().split(" ")
        if annotationType == "gene_id" or annotationType == "transcript_id":
          countObjectId = annotationValue.strip("\"")

          if countObjectId in countObjects:
           if chromosome in countObjects[countObjectId]:
             countObjects[countObjectId][chromosome].append(interval)
           else:
             countObjects[countObjectId][chromosome] = [interval]
          else:
            countObjects[countObjectId] = {}
            countObjects[countObjectId][chromosome] = [interval]
            
    lineNum += 1
    if lineNum % 100000 == 0:
      sys.stderr.write(".")
      sys.stderr.flush()
  
  gtfFile.close()

  if lineNum > 100000:
    sys.stderr.write("\n")
    sys.stderr.flush()

  return countObjects


################################################################################
##
## computeLength
##
################################################################################

def computeLength (intervalList):
  
  intervalList.sort(key=lambda tup: tup[0])
  oldInterval = [-1,-1]
  combinedList = []
  for interval in intervalList:
    if not oldInterval[1] < interval[0]:
      oldInterval = [oldInterval[0], max(oldInterval[1], interval[1])]
    else:
      combinedList.append(oldInterval)
      oldInterval = interval
  
  combinedList.append(oldInterval)
  
  length = 0
  for interval in combinedList[1:]:
    length = length + interval[1] - interval[0] + 1
  
  return length


################################################################################
##
## computeCountObjectLengths
##
################################################################################

def computeCountObjectLengths (countObjects):
  countObjectLengths = {}
  
  for countObjectId in countObjects:
    countObject = countObjects[countObjectId]
    
    maxCountObjectLength = 0
    for chromosome in countObject:
      countObjectLength = computeLength(countObject[chromosome])
      if countObjectLength > maxCountObjectLength:
        maxCountObjectLength = countObjectLength
    
    countObjectLengths[countObjectId] = maxCountObjectLength
  
  return countObjectLengths

################################################################################
##
## Write FPKM file
##
################################################################################

def outputFpkmFile (fpkmFilename, computeMean, controlSamples, samples, fpkm):
  try:
    fpkmFile = open(fpkmFilename, 'w')
  except IOError, e:
    raise Exception ("File " + fpkmFilename + " cannot be opened ... skipping\n" + \
                     "Unix error code: " + str(e[0]) + " and message: " + str(e[1]))
    
  if computeMean:
    if controlSamples != "":
      controlSampleInd = [samples.index(controlSample) for controlSample in controlSamples.split(":")]
    else:
      controlSampleInd = range(len(samples))
    
    print "Writing mean of control sample FPKM values to " + fpkmFilename
    print >> fpkmFile, "\t".join(["Gene Id"] + ["Mean control FPKM"])
    for countObjectId in sorted(fpkm.iterkeys()):
      print >> fpkmFile, countObjectId + "\t" + str(numpy.mean([fpkm[countObjectId][i] for i in controlSampleInd]))
      
  else:
    print "Writing FPKM values to " + fpkmFilename
    print >> fpkmFile, "\t".join(["Gene Id"] + samples)
    for countObjectId in sorted(fpkm.iterkeys()):
      print >> fpkmFile, countObjectId + "\t" + "\t".join([str(v) for v in fpkm[countObjectId]])
  
  fpkmFile.close()


################################################################################
##
## Main program
##
################################################################################

## Read counts
counts, samples = readCountFile(countFilename)

## Read gene lengths
countObjectLengths = readGeneLengthFile (geneLengthFilename)

if readNumFilename != "":
  numReads = readNumReadsFile(readNumFilename)
else:
  print >> sys.stderr, "Using the sum of expression values per sample for normalization"
  numReads = {}
  for i in range(len(samples)):
    numReads[samples[i]] = sum([c[i] for c in counts.values()])


## Read GTF file if provided
if gtfFilename != "":
  gtfCountObjects = readGtfFile(gtfFilename)
  gtfCountObjectLengths = computeCountObjectLengths (gtfCountObjects)
  
  ## Add missing gene lengths
  lengthsAdded = 0
  for countObjectId in counts:
    if not countObjectId in countObjectLengths:
      if countObjectId in gtfCountObjectLengths:
        lengthsAdded += 1
        countObjectLengths[countObjectId] = gtfCountObjectLengths[countObjectId]
        
  if lengthsAdded > 0:
    print >> sys.stderr, "Lengths for " + str(lengthsAdded) + " genes/transcripts added."

## Compute FPKM
fpkm = {}
genesWithoutLength = set([])
transcriptsWithoutLength = set([])
for countObjectId in counts:
  if not countObjectId in countObjectLengths:
    if countObjectId.startswith("NM_") or countObjectId.startswith("XM_") or countObjectId.startswith("NR_") or countObjectId.startswith("XR_") or \
       countObjectId.startswith("ENST") or countObjectId.startswith("ENSMUST") or countObjectId.startswith("ENSRNOT"):
      transcriptsWithoutLength.add(countObjectId)
    else:
      genesWithoutLength.add(countObjectId)
  else:
    countObjectLength = countObjectLengths[countObjectId]
    
    fpkm[countObjectId] = [c / (countObjectLength / 1000.0) for c in counts[countObjectId]]
    
    for i in range(len(samples)):
      numReadsSample = numReads[samples[i]]
      fpkm[countObjectId][i] = fpkm[countObjectId][i] / (numReadsSample / (1000.0 * 1000.0))

if len(genesWithoutLength) > 0:
  print >> sys.stderr, "There are " + str(len(genesWithoutLength)) + " genes without length."

if len(transcriptsWithoutLength) > 0:
  print >> sys.stderr, "There are " + str(len(transcriptsWithoutLength)) + " transcripts without length."


outputFpkmFile (fpkmFilename, computeMean, controlSamples, samples, fpkm)

