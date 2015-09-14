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
##  Combines all overlapping exons of a gene into one BED entry
##
################################################################################

import sys
import argparse
import re
import os

## execfile(os.path.join(os.environ["HOME"], "ngs", "pipelines", "exon-pipeline", "bin", "process-scripts", "util-scripts", "combineGenomeExonBedEntries.py"))

################################################################################
##
## Define command line options
##
################################################################################

parser = argparse.ArgumentParser(description='Combines all overlapping exons of a gene into one BED entry.')
parser.add_argument('-d', default=0, dest="debugLevel", metavar="INT", help='debug level [0, no debug output] ')
parser.add_argument('-g', dest='gtfFile', metavar="STRING", help='GTF input file')
parser.add_argument('-b', dest='bedFile', metavar="STRING", help='BED file with combined exon genome BED entries')
parser.add_argument('-m', dest='mapFile', metavar="STRING", help='Combined exon to gene id map file')

if len(sys.argv) > 1:
  args = parser.parse_args(sys.argv[1:])
else:
  projectDir = os.path.join(os.environ["HOME"], "ngs", "RNA-seq-test", "SEQC-NG00008.1-EQP2.0-Ensembl", "exon-pipeline-files")
  inputArgs = []
  inputArgs.append('-g')
  inputArgs.append(os.path.join(projectDir, "gtf-files", "ensembl_rna_hs.gtf"))
  inputArgs.append('-b')
  inputArgs.append(os.path.join(projectDir, "bed-files", "ensembl_rna_hs_genome_exons_combined.bed"))
  inputArgs.append('-m')
  inputArgs.append(os.path.join(projectDir, "map-files", "ensembl_rna_hs_exon_gene_combined.map"))
  args = parser.parse_args(inputArgs)

debugLevel  = args.debugLevel
gtfFilename = args.gtfFile
bedFilename = args.bedFile
mapFilename = args.mapFile


################################################################################
##
## Check for overlap of an interval with a list of intervals
##
################################################################################
      
def overlaps (interval1, interval2):
  
  if (interval1[0] <= interval2[0] and interval2[0] <= interval1[1]) or (interval2[0] <= interval1[0] and interval1[0] <= interval2[1]):
    return True
  
  return False


################################################################################
##
## union
##
################################################################################

def union (interval1, interval2):

  return [min (interval1[0], interval2[0]), max (interval1[1], interval2[1])]


################################################################################
##
## readGtfFile
##
## Use the exonNumber fields of each transcript to infer how many alignments a
## transcript has. This is captured by the variable alignmentNumber.
##
## Each of the three transcript output variables has an entry for each each transcript
## and each alignmentNumber of each transcript.
##
################################################################################

def readGtfFile (gtfFilename):
  
  transcriptExonList         = {}
  transcriptChromosomeStrand = {}

  geneTranscriptMap     = {}
  
  transcriptAlignmentExonNumbers = {}

  lineNumChunkSize = 50000
 
  try:
    gtfFile = open(gtfFilename)
  except IOError, e:
    raise Exception(gtfFilename + " cannot be opened ... skipping\n" + 
                    "Unix error code and message: " + str(e[0]) + " " + str(e[1]))

  print >> sys.stderr, "Opening file: " + gtfFilename
  lineNum = 0
  oldTranscriptId = ""
  alignmentNumber = 1

  for line in gtfFile:
    line = line.rstrip ()
    gtfEntries = line.split("\t")
    
    chromosome  = gtfEntries[0]
    source      = gtfEntries[1]
    featureType = gtfEntries[2]
    interval    = [int(gtfEntries[3]), int(gtfEntries[4])]
    strand      = gtfEntries[6]
    frame       = gtfEntries[7]
    
    score = 0
    if gtfEntries[5] != ".":
      score = float(gtfEntries[5])
    
    exonId      = "/".join(map(str, [chromosome] + interval + [strand]))
    
    if featureType == "exon":
      annotationEntries = gtfEntries[8].split(";")
      geneId = ""
      transcriptId = ""
      exonNumber = -1
      for annotationEntry in annotationEntries:
        if annotationEntry != "":
          annotationType, annotationValue, annotationEmpty = annotationEntry.strip().split('"')
          if annotationEmpty != "":
            raise Exception ("Unknown format of annotation entry: " + annotationEntry)

          annotationType = annotationType.strip()
          if annotationType == "gene_id":
            geneId = annotationValue
          elif annotationType == "transcript_id":
            transcriptId = annotationValue
          elif annotationType == "exon_number":
            exonNumber = int(annotationValue)
            
      if oldTranscriptId != transcriptId:
        alignmentNumber = 1
        oldTranscriptId = transcriptId
        
      if exonNumber == -1:
        raise Exception ("WARNING: no exon number for exon " + exonId + " of transcript " + transcriptId)
               
      if not transcriptId in transcriptExonList:
        transcriptExonList            [transcriptId] = {}
        transcriptChromosomeStrand    [transcriptId] = {}
        transcriptAlignmentExonNumbers[transcriptId] = {}

      ## Increment the alignment number if the current exon number already exists in the current alignment
      ## Note that GTF entries of other transcripts may occur between two alignments of the same transcript
      while alignmentNumber in transcriptAlignmentExonNumbers[transcriptId] and \
             exonNumber in transcriptAlignmentExonNumbers[transcriptId][alignmentNumber]:
        alignmentNumber += 1

      if not alignmentNumber in transcriptAlignmentExonNumbers[transcriptId]:
        transcriptAlignmentExonNumbers[transcriptId][alignmentNumber] = [exonNumber]
      else:
        transcriptAlignmentExonNumbers[transcriptId][alignmentNumber].append(exonNumber)

      if not alignmentNumber in transcriptChromosomeStrand[transcriptId]:
        #print "transcriptChromosomeStrand of " + transcriptId + "." + str(alignmentNumber) + " set to " + chromosome + strand
        transcriptChromosomeStrand[transcriptId][alignmentNumber] = chromosome + "/" + strand
      elif transcriptChromosomeStrand[transcriptId][alignmentNumber] != chromosome + "/" + strand:
        print >> sys.stderr, "WARNING: Exon number " + str(exonNumber) + " of transcript " + transcriptId + " on chromosome/strand " + chromosome + strand + \
              " is assigned to alignment " + str(alignmentNumber) + " on chromosome/strand " + transcriptChromosomeStrand[transcriptId][alignmentNumber]
      
      if alignmentNumber in transcriptExonList[transcriptId]:
        transcriptExonList[transcriptId][alignmentNumber].append(interval)
      else:
        transcriptExonList[transcriptId][alignmentNumber] = [interval]
        
      if exonNumber in transcriptExonList[transcriptId][alignmentNumber]:
        print >> sys.stderr ("Exon number: " + str(exonNumber) + " already stored for alignment " + str(alignmentNumber) + " of " + transcriptId)
        sys.exit(1)
      
      if geneId != "" and transcriptId != "":
        if geneId in geneTranscriptMap:
          if not transcriptId in geneTranscriptMap[geneId]:
            geneTranscriptMap[geneId].append(transcriptId)
        else:
          geneTranscriptMap[geneId] = [transcriptId]
    
    lineNum = lineNum + 1
    if lineNum % lineNumChunkSize == 0:
      sys.stdout.write(".")
      sys.stdout.flush()
  
  if lineNum > lineNumChunkSize:
    sys.stdout.write("\n")
  
  gtfFile.close()
  return transcriptChromosomeStrand, transcriptExonList, geneTranscriptMap


################################################################################
##
## computeIntervalListLength
##
################################################################################

def combineIntervals (intervalList):
  
  intervalList.sort(key=lambda tup: tup[0])
  oldInterval = [-1,-1]
  combinedList = []
  for interval in sorted(intervalList):
    if oldInterval[1] + 1 < interval[0]:
      combinedList.append(oldInterval)
      oldInterval = interval
    else:
      oldInterval = [oldInterval[0], max(oldInterval[1], interval[1])]
  
  combinedList.append(oldInterval)
    
  return combinedList[1:]


################################################################################
##
## Main program
##
## Read GTF file and compute transcript length lists
##
################################################################################

#raise Exception ("stop")

transcriptChromosomeStrand, transcriptExonList, geneTranscriptMap = readGtfFile(gtfFilename)


################################################################################
##
## For each gene and each chromosome/strand combination that the gene maps to
## compute the combined exon list and output it to the BED file as well as exonId
## to gene mappings
##
################################################################################

try:
  bedFile = open(bedFilename, 'w')
except IOError, e:
  raise Exception(bedFilename + " cannot be opened ... skipping\n" + 
                  "Unix error code and message: " + str(e[0]) + " " + str(e[1]))
print >> sys.stderr, "Writing combined BED entries to file " + bedFilename

try:
  mapFile = open(mapFilename, 'w')
except IOError, e:
  raise Exception(mapFilename + " cannot be opened ... skipping\n" + 
                  "Unix error code and message: " + str(e[0]) + " " + str(e[1]))
print >> sys.stderr, "Writing combined exon to gene mappings to file " + mapFilename


geneLengths = {}
for geneId in geneTranscriptMap:
  geneTranscriptList = geneTranscriptMap[geneId]
  
  ## Create exon lists for each chromosome/strand that a transcript of the gene aligns to.
  geneExonList = {}
  for transcriptId in geneTranscriptMap[geneId]:
    for alignmentNumber in transcriptChromosomeStrand[transcriptId]:
      if not transcriptChromosomeStrand[transcriptId][alignmentNumber] in geneExonList:
        geneExonList[transcriptChromosomeStrand[transcriptId][alignmentNumber]] = []
      geneExonList[transcriptChromosomeStrand[transcriptId][alignmentNumber]].extend(transcriptExonList[transcriptId][alignmentNumber])
  
  for chromosomeStrand in geneExonList:
    geneExonList[chromosomeStrand] = combineIntervals(geneExonList[chromosomeStrand])
    chromosome, strand = chromosomeStrand.split("/")
    for exonInterval in geneExonList[chromosomeStrand]:
      exonStart, exonEnd = exonInterval
      exonId = "/".join(map(str, [chromosome, exonStart, exonEnd, strand]))
      print >> bedFile, "\t".join(map(str, [chromosome, exonStart - 1, exonEnd, exonId, 0, strand]))
      print >> mapFile, "\t".join(map(str, [exonId, geneId]))

bedFile.close()
mapFile.close()
      
