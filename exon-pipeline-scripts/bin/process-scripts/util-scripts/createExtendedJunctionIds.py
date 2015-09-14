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
##  Script to create main junction ids based on the junction extensions.
##  Needed to identify possible additional main junction ids.
##
################################################################################

import sys
import argparse
import os
import re
import os.path

## execfile(os.path.join(os.environ["HOME"], "ngs", "pipelines", "exon-pipeline", "bin", "util-lib", "createExtendedJunctionIds.py"))

################################################################################
##
## Define command line options
##
################################################################################

parser = argparse.ArgumentParser(description='Create extended junction ids')
parser.add_argument('-d', type=int, default=0, dest="debugLevel", metavar="INT", help='debug level [0, no debug output]')
parser.add_argument('-j', dest="junctionIdFile", metavar="<Junction Id file>", help='File with junction ids')
parser.add_argument('-o', dest="outputFile", metavar="<Output file>", help='File with mapping of internal exon ids to external junction ids')
parser.add_argument('-w', dest="printWarnings", action='store_true', help='Flag to enable the printing of warnings')


if len(sys.argv) > 1:
  args = parser.parse_args(sys.argv[1:])
else:
  projectDir = os.path.join(os.environ["HOME"], "ngs", "RNA-seq-test", "EQP2.0-speed-test", "samples", "SEQC-test", "C001", "sam-files")
  
  inputArgs = []
  inputArgs.append("-j")
  inputArgs.append(os.path.join(projectDir, "SEQC-test-C001-bowtie2-junction-pe-reference.ids"))
  inputArgs.append("-o")
  inputArgs.append(os.path.join(projectDir, "SEQC-test-C001-bowtie2-junction-pe-reference.ids.new"))
  inputArgs.append("-w")
  
  print "No arguments given - using files in " + projectDir
  args = parser.parse_args(inputArgs)

debugLevel         = args.debugLevel
junctionIdFilename = args.junctionIdFile
outputFilename     = args.outputFile
printWarnings      = args.printWarnings


################################################################################
##
##  Main program
##
################################################################################
  
try:
  if junctionIdFilename == "-":
    junctionIdFile = sys.stdin
    print >> sys.stderr, "Reading junction ids from stdin"
  else:
    junctionIdFile = open(junctionIdFilename)
    print >> sys.stderr, "Reading junction ids from file: " + junctionIdFilename
except IOError:
  raise Exception("Could not open file " + junctionIdFilename)

try:
  if outputFilename == "-":
    outputFile = sys.stdout
    print >> sys.stderr, "Writing junction ids to stdout"
  else:
    outputFile = open(outputFilename, 'w')
    print >> sys.stderr, "Writing junction ids to file: " + outputFilename
except IOError:
  raise Exception("Could not open file " + outputFilename)
  
lineNumber = 0
junctionIdSet = set([])
for line in junctionIdFile:
  if "-" in line:
    junctionId, leftExtension, rightExtension = line.rstrip ().split (":")
    if leftExtension != "":
      leftExtensionIds = leftExtension.split("-")
    else:
      leftExtensionIds = []
      
    if rightExtension != "":
      rightExtensionIds = rightExtension.split("-")
    else:
      rightExtensionIds = []
      
    geneId, fillWord, junctionId1, junctionId2 = junctionId.split("-")
  
    junctionIds = leftExtensionIds + [junctionId1, junctionId2] + rightExtensionIds
  
    for i in range(len(junctionIds) - 1):
      curJunctionId = "-".join([geneId, fillWord, junctionIds[i], junctionIds[i+1]])
      if not curJunctionId in junctionIdSet:
        print >> outputFile, curJunctionId
        junctionIdSet.add(curJunctionId)
       
    lineNumber = lineNumber + 1
    if lineNumber % 50000 == 0:
      sys.stderr.write(".")
      sys.stderr.flush()

if junctionIdFilename != "-":
  junctionIdFile.close()
  
if outputFilename != "-":
  outputFile.close()

if lineNumber > 50000:
  sys.stderr.write("\n")
  

