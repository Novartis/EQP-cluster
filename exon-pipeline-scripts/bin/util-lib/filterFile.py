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
##  Script to filter a file consisting of fields (default: tab-separated) for a
##  the line where at least one of the fields is contained in set of keys created
##  from a filter file
##
################################################################################

import sys
import argparse
import os
import gzip

## execfile(os.path.join(os.environ["HOME"], "ngs", "pipelines", "exon-pipeline", "bin", "util-lib", "filterFile.py"))


################################################################################
##
## Define command line options
##
################################################################################

parser = argparse.ArgumentParser(description='Extract pipeline meta informations')
parser.add_argument('-d', type=int, default=0, dest="debugLevel", metavar="INT", help='debug level [0, no debug output]')
parser.add_argument('-i', dest="inputFile", default="-", metavar="STRING",  help='File to be filtered')
parser.add_argument('-I', dest="index", default="", metavar="INT",  help='List of indices of the fields of the input file to be used for filtering (colon separated)')
parser.add_argument('-S', dest="splitString", default="\t", metavar="STRING",  help='String to be used to split the input into fields')
parser.add_argument('-s', dest="separatorString", default="", metavar="STRING",  help='String to be used to split the filter entries and the key field entries.')
parser.add_argument('-f', dest="filterFile", metavar="STRING",  help='File to be used for filtering')
parser.add_argument('-F', dest="filterIndex", default="", metavar="INT",  help='List of indices of the field of the filter file to be used for filtering (colon separated)')
parser.add_argument('-v', dest="invertFilter", action="store_true", default=False, \
                    help='Flag to invert the filter: output only lines that do not match any of the entries in the filter file')
parser.add_argument('-o', dest="outputFile", default="-", metavar="STRING",  help='Filtered file')

if len(sys.argv) > 1:
  args = parser.parse_args(sys.argv[1:])
else:
  inputArgs = []
  projectDir = os.path.join(os.environ["HOME"], "ngs", "RNA-seq-test", "Sim-NG00008.1-EQP2.0-Ensembl")
  inputArgs.append("-f")
  inputArgs.append(os.path.join(projectDir, "samples-010x", "S1.1", "C001", "sam-files", "S1.1-C001-bowtie2-junction-pe-external-junction.ids"))
  inputArgs.append("-i")
  inputArgs.append(os.path.join(projectDir, "exon-pipeline-files", "map-files", "ensembl_rna_hs_equivalent_junctions.map"))
  inputArgs.append("-o")
  inputArgs.append(os.path.join(projectDir, "samples-010x", "S1.1", "C001", "sam-files", "S1.1-C001-bowtie2-junction-pe-equivalent-junction.ids.new"))
  args = parser.parse_args(inputArgs)

debugLevel        = args.debugLevel
inputFilename     = args.inputFile
outputFilename    = args.outputFile
filterFilename    = args.filterFile
invertFilter      = args.invertFilter
indexString       = args.index
splitString       = args.splitString
separatorString   = args.separatorString
filterIndexString = args.filterIndex


################################################################################
##
## Main program
##
################################################################################

index = [0]
if indexString != "":
  index = [int(i) for i in args.index.split(":")]

filterIndex = [0]
if filterIndexString != "":
  filterIndex = [int(i) for i in args.filterIndex.split(":")]

if len(index) != len(filterIndex):
  raise Exception ("Filter indices do not have the same length: " + indexString + " - " + str(len(index)) + " vs " + filterIndexString + " - " + str(len(filterIndex)))


################################################################################
##
## open a file
##
################################################################################

def openFile (filename, mode="r"):
  try:
    if filename.endswith(".gz"):
      fileLink = gzip.open(filename, mode + "b")
    elif filename != "-":
      fileLink = open(filename, mode)
    elif "r" in mode:
      fileLink = sys.stdin
    elif "w" in mode:
      fileLink =sys.stdout
    else:
      raise Exception("Unknown mode for using stdin or stdout: " + mode)

    return fileLink
      
  except IOError, e:
    raise Exception (filename + " cannot be opened ... exiting\n" +
                     "Unix error code and message: " + str(e[0]) + " " + str(e[1]))
  

################################################################################
##
## Reading filter file
##
################################################################################

try:
  filterFile = openFile(filterFilename)
except IOError, e:
  raise Exception (filterFilename + " cannot be opened ... exiting\n" + 
                   "Unix error code and message on opening: " + str(e[0]) + " " + str(e[1]))
print >> sys.stderr, "Reading file " + filterFilename


filterEntries = set([])
lineNum = 0
for line in filterFile:
  lineFields = line.strip().split(splitString)
  if separatorString != "":
    lineFields = [lineField.split(separatorString)[0] for lineField in lineFields]
  filterEntries.add(splitString.join([lineFields[i] for i in filterIndex]))
  
  lineNum += 1
  if lineNum % 500000 == 0:
    sys.stderr.write(".")
    sys.stderr.flush()
    
if lineNum > 500000:
  sys.stderr.write("\n")
  sys.stderr.flush ()

print >> sys.stderr, str(lineNum) + " filter entries read."
filterFile.close()


################################################################################
##
## Reading and filtering the input file
##
################################################################################

try:
  if inputFilename == "" or inputFilename == "-":
    inputFile = sys.stdin
    print >> sys.stderr, "Reading from stdin"
  else:
    inputFile = openFile(inputFilename)
    print >> sys.stderr, "Reading file " + inputFilename
except IOError, e:
  raise Exception (inputFilename + " cannot be opened ... exiting\n" + 
                   "Unix error code and message on opening: " + str(e[0]) + " " + str(e[1]))

try:
  if outputFilename == "" or outputFilename == "-":
    outputFile = sys.stdout
    print >> sys.stderr, "Writing to stdout"
  else:
    outputFile = openFile(outputFilename, 'w')
    print >> sys.stderr,  "Writing to file " + outputFilename
except IOError, e:
  raise Exception (outputFilename + " cannot be opened ... exiting\n" + 
                   "Unix error code and message on opening: " + str(e[0]) + " " + str(e[1]))

lineNum = 0
numSelectedInputEntries = 0
for line in inputFile:
  # print >> sys.stderr, line + " " + str(idIndex)
  line = line.strip()
  lineFields = line.strip().split(splitString)
  if separatorString != "":
    lineFields = [lineField.split(separatorString)[0] for lineField in lineFields]

  lineField  = splitString.join([lineFields[i] for i in index])
  fieldFound = lineField in filterEntries

  if (fieldFound and not invertFilter) or (not fieldFound and invertFilter):
    print >> outputFile, line
    numSelectedInputEntries += 1
    
  lineNum += 1
  if lineNum % 500000 == 0:
    sys.stderr.write(".")
    sys.stderr.flush()

if lineNum > 500000:
  sys.stderr.write("\n")
  sys.stderr.flush ()
      
print >> sys.stderr, str(numSelectedInputEntries) + " Input entries of " + str(lineNum) + " selected."
inputFile.close()
outputFile.close ()

