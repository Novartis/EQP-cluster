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

## execfile(os.path.join(os.environ["HOME"], "ngs", "pipelines", "exon-pipeline", "bin", "util-lib", "combineCountsTest.py"))

################################################################################
##
##  Script to convert a gene exon file to a gene rank map file
##
################################################################################

import sys
import argparse
import re
import os.path
import subprocess
import array
import numpy as np
import time
import gzip

################################################################################
##
## Define command line options
##
################################################################################

parser = argparse.ArgumentParser(description='Combine the count files from different')
parser.add_argument('-d', type=int, default=0, dest="debugLevel", metavar="INT", help='debug level [0, no debug output]')
parser.add_argument('-n', dest="outputZeroes", action="store_false", default=True, help='output only rows with at least one non-zero count')
parser.add_argument('-A', dest="addUp", action="store_true",  default=False, help='Add counts instead of storing them as separate columns')
parser.add_argument('-M', dest="multiColumn", action="store_true",  default=False, help='Flag to indicate whether the files to be combined contain multiple columns')
parser.add_argument('-s', dest="fileSuffix", default="", metavar="<Input file suffix>", help='Suffix to remove from file names to obtain column headings (see option -c)')
parser.add_argument('-i', dest="inputFiles", metavar="<Input count file>", nargs="+", help='One or more files with counts')
parser.add_argument('-f', dest="filterFile", default="", metavar="<Filer file>", help='File containing the set of strings needed to filter out non-zero rows')
parser.add_argument('-F', dest="useFilter", default=False, action="store_true", help='Flag to indicate whether to filter input files for non-zero rows.')
parser.add_argument('-I', dest="countObjectIdFile", default="", metavar="<Id file>", help='Optional file containing the complete set of ids.')
parser.add_argument('-o', dest="combinedCountFile", metavar="<Combined count file>", help='File with combined counts')


if len(sys.argv) > 1:
  args = parser.parse_args(sys.argv[1:])
else:
  projectDir = os.path.join (os.environ["HOME"], "ngs", "RNA-seq-test", "SEQC-NG00007.0-EQP2.0-Ensembl")
  fileSuffix = "combined-mixed-gene.cnt"
  fileSuffix = "mixed-junction-sam.cnt"
  
  inputArgs = []
  inputArgs.append("-A")
  inputArgs.append("-s")
  inputArgs.append(fileSuffix)
  inputArgs.append("-i")
  inputArgs.append(os.path.join (projectDir, "samples", "SEQC-A-BC01-s_1", "count-files", "SEQC-A-BC01-s_1-EQP-genome-junction.cnt"))
  inputArgs.append(os.path.join (projectDir, "samples", "SEQC-A-BC01-s_3", "count-files", "SEQC-A-BC01-s_3-EQP-genome-junction.cnt"))
  inputArgs.append("-o")
  inputArgs.append(os.path.join (projectDir, "samples", "count-files", "SEQC-NG00007.0-EQP2.0-Ensembl-EQP-genome-junction-test.cnt"))

  args = parser.parse_args(inputArgs)

debugLevel            = args.debugLevel
outputZeroes          = args.outputZeroes
addUp                 = args.addUp
multiColumn           = args.multiColumn
fileSuffix            = args.fileSuffix
inputFiles            = args.inputFiles
filterFilename        = args.filterFile
useFilter             = args.useFilter
countObjectIdFilename = args.countObjectIdFile
combinedCountFile     = args.combinedCountFile

if len(sys.argv) <= 1:
  print >> sys.stderr, "Using default arguments"

lineNumChunkSize = 500 * 1000


################################################################################
##
## createFilterFile: Create a file to filter for non-zero counts
##
################################################################################

def createFilterFile (filterFilename):
  
  print >> sys.stderr, "Creating filter file" + filterFilename
  try:
    filterFile = open (filterFilename, 'w')
  except IOError, e:
    raise Exception(filterFilename + " cannot be opened for writing... exiting" + \
                    "Error: " + str(e[0]) + " " + str(e[1]))
  
  print >> filterFile, "Id\t" 
  for i in range (10):
    print >> filterFile, "\t0.0" + str(i)
  
  for i in range (9):
    print >> filterFile, "\t0." + str(i + 1)
  
  for i in range (9):
    print >> filterFile, "\t" + str(i + 1)
  
  filterFile.close ()


################################################################################
##
## getCounts
##
################################################################################

def getCounts (inputFilename, filterFilename, useFilter):
  
  colHeaders = [""]

  if useFilter and filterFilename != "":
    if not os.path.isfile (filterFilename):
      createFilterFile (filterFilename)

    print >> sys.stderr, "Reading filtered file: " + inputFilename
    cmd = ["fgrep", "-f", filterFilename, inputFilename]
    try:
      inputCounts = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).communicate()[0].split("\n")
    except OSError, e:
      print >> sys.stderr, 'Execution of command: "' + " ".join(cmd) + '" failed.'
      raise Exception ("ERROR: " + str(e))
  else:
    print >> sys.stderr, "Reading complete file: " + inputFilename
    try:
      inputCounts = open (inputFilename)
    except IOError, e:
      raise Exception(inputFilename + " cannot be opened for reading... exiting" + \
                      "Error: " + str(e[0]) + " " + str(e[1]))

  lineNum           = 0
  maxLenCountValues = 0
  count = {}
  for line in inputCounts:
    line = line.rstrip ()
    countValues = line.split("\t")
    maxLenCountValues = max([maxLenCountValues, len(countValues) - 1])
    countObjectId = countValues[0]
    if countObjectId != "":
      if countObjectId != "Id":
        count[countObjectId] = array.array ('d', [float(v) for v in countValues[1:]])
      else:
        colHeaders = countValues[1:]
      
    lineNum = lineNum + 1
    if lineNum % lineNumChunkSize == 0:
      sys.stderr.write(".")
      sys.stderr.flush()
    
  if lineNum > lineNumChunkSize:
    sys.stderr.write("\n")

  if not useFilter:
    inputCounts.close()

  return [maxLenCountValues, colHeaders, count]


################################################################################
##
## Main program: Reading the input files
##
################################################################################

combinedCounts = {}
colHeaders = []

if useFilter:
  if filterFilename == "":
    filterFileDir  = os.path.dirname(inputFiles[0])
    filterFilename = os.path.join(filterFileDir, "." + inputFiles[0].split(".")[0] + "-filter.tmp")
  
inputFilenum = 0
for inputFile in inputFiles:
  # start = time.clock()    
  length, colHeadersFile, count = getCounts (inputFile, filterFilename, useFilter)
  # print >> sys.stderr, "Elapsed time: " + str(time.clock() - start) + "s"

  if length > 0:
    if len(colHeadersFile) == 1 and colHeadersFile[0] == "":
      inputFileBase = re.sub(fileSuffix, '', os.path.basename(inputFile)).rstrip(" -_")
      colHeaders.append(inputFileBase)
    else:
      colHeaders = colHeaders + colHeadersFile
  for countObjectId in count:
    if addUp:
      if countObjectId in combinedCounts:
         if len(combinedCounts[countObjectId]) == len(count[countObjectId]):
           countsArray = [combinedCounts[countObjectId].tolist (), count[countObjectId].tolist ()]
           combinedCounts[countObjectId] = array.array("d", [sum(c) for c in zip(*countsArray)])
         else:
           raise Exception ("Length of combined counts " + str(len(combinedCounts[countObjectId])) + " does not match the length of the added counts " + \
             str(len(count[countObjectId])) + " for " + countObjectId)
      elif outputZeroes or count[countObjectId] != 0:
         combinedCounts[countObjectId] = count[countObjectId]
    else:
      if countObjectId in combinedCounts:
         combinedCounts[countObjectId] = combinedCounts[countObjectId] + count[countObjectId]
      elif not multiColumn:
        if outputZeroes or count[countObjectId] != 0:
           combinedCounts[countObjectId] = array.array("d", [0.0] * inputFilenum) + count[countObjectId]
      else:
        raise Exception ("Count object " + countObjectId + " occurs in file " + inputFile + " but not in the other files and option -M is set.")

  if len(colHeaders) != len(combinedCounts[countObjectId]):
    print >> sys.stderr, "Differing number of columns for header " + str(len(colHeaders)) + " and counts " + str(len(combinedCounts[countObjectId]))

  filteredCountObjects = [c for c in combinedCounts if not c in count]
  if len(filteredCountObjects) > 0 and multiColumn:
    raise Exception ("There are " + str(len(filteredCountObjects)) + " count objects missing in file " + inputFile + " and option -M is set.")

  if not addUp:
    for countObjectId in [c for c in combinedCounts if not c in count]:
      combinedCounts[countObjectId].append(0)
      
  inputFilenum += 1

  ## If we do not have multicolumn files, we can use the filter for the files after the
  ## first file.
  if not multiColumn:
    useFilter = True


if os.path.isfile (filterFilename):
  print >> sys.stderr, "Removing file " + filterFilename
  cmd = "rm " + filterFilename
  try:
    retcode = subprocess.call(cmd, shell=True)
  except OSError, e:
    print >>sys.stderr, 'Execution of command: "' + cmd + '" failed:', e
  if retcode < 0:
    print >>sys.stderr, 'Command: "' + cmd + '" was terminated by signal', -retcode
  elif retcode > 0:
    print >>sys.stderr, 'Command: "' + cmd + '" returned ', retcode


################################################################################
##
## Writing the output
##
################################################################################

try:
  outputFile = open(combinedCountFile, 'w')
except IOError, e:
  raise Exception(combinedCountFile + " cannot be opened for writing... exiting\n" + "ERROR: " + str(e[0]) + " " + str(e[1]))

print >> sys.stderr, "Writing combined counts to " + combinedCountFile
if not addUp:
  print >> outputFile, "Id" + "\t" + "\t".join(colHeaders)
lineNum = 0

if countObjectIdFilename != "":
  try:
    print >> sys.stderr, "Reading count objects ids from file " + countObjectIdFilename
    if countObjectIdFilename.endswith(".gz"):
      countObjectIdFile = gzip.open(countObjectIdFilename)
    else:
      countObjectIdFile = open(countObjectIdFilename)
  except IOError, e:
    raise Exception(countObjectIdFilename + " cannot be opened for reading... exiting\n" + "ERROR: " + str(e[0]) + " " + str(e[1]))

  columnNum = len(combinedCounts[combinedCounts.keys()[0]])
  for line in countObjectIdFile:
    countObjectId = line.split("\t")[0]
    if countObjectId in combinedCounts:
      print >> outputFile, countObjectId + "\t" + "\t".join([str(v) for v in combinedCounts[countObjectId]])
    else:
      print >> outputFile, countObjectId + "\t" + "\t".join(["0.0"] * columnNum)

else:
  for countObjectId in sorted(combinedCounts.iterkeys()):
    print >> outputFile, countObjectId + "\t" + "\t".join([str(v) for v in combinedCounts[countObjectId]])
    
    lineNum = lineNum + 1
    if lineNum % lineNumChunkSize == 0:
      sys.stderr.write(".")
      sys.stderr.flush()
    
if lineNum > lineNumChunkSize:
  sys.stderr.write("\n")
  
outputFile.close()
print >> sys.stderr, "combineCounts done."
