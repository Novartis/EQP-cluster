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
##  Script to create symbolic links from a data release info sheet.
##
################################################################################

import sys
import argparse
import re
import os.path
import subprocess

## execfile(os.path.join(os.environ["HOME"], "ngs", "pipelines", "exon-pipeline", "bin", "util-lib", "createFileLinks.py"))

################################################################################
##
## Define command line options
##
################################################################################

parser = argparse.ArgumentParser(description='Combine the count files from different')
parser.add_argument('-d', type=int, default=0, dest="debugLevel", metavar="INT", help='debug level [0, no debug output]')
parser.add_argument('-a', dest="annotationFiles", metavar="<Annotation files>", nargs="+", help='Files with the sample to fastq file mappings')
parser.add_argument('-s', dest="samplesDir", metavar="<Samples directory>", help='Directory containing all the samples')
parser.add_argument('-i', dest="insertSizeFile", default="", metavar="<Insert size file>", help='Directory containing all the samples')
parser.add_argument('-S', dest="singleReadData", action="store_true", default=False, help='Flag to indicate if the data are single read data.')

if len(sys.argv) >= 2:
  args = parser.parse_args(sys.argv[1:])
else:
  projectDir = os.path.join(os.environ["HOME"], "ngs", "RNA-seq", "2014", "SEQC-strand-specific")
  projectDir = os.path.join(os.environ["HOME"], "ngs", "RNA-seq-test", "eyedevtc")
  inputArgs = []
  inputArgs.append("-a")
  inputArgs.append(os.path.join(projectDir, "annotation-files", "IL_NGS_11_067_120104_SN889_0095_AD0902ACXX_Data_Release_info_GA258.finished.txt"))
  inputArgs.append(os.path.join(projectDir, "annotation-files", "sample_files.txt"))
  inputArgs.append("-s")
  inputArgs.append(os.path.join(projectDir, "samples-test"))
  inputArgs.append("-S")
  args = parser.parse_args(inputArgs)

debugLevel         = args.debugLevel
annotationFiles    = args.annotationFiles
samplesDir         = args.samplesDir
insertSizeFilename = args.insertSizeFile
singleReadData     = args.singleReadData


################################################################################
##
## Extract information for samples from annotation files
##
################################################################################

sampleNames = []
sampleFields = {}

for annotationFilename in annotationFiles:
  try:
    annotationFile = open(annotationFilename)
    if insertSizeFilename != "":
      insertSizeFile = open(insertSizeFilename, 'w')
    
    print >> sys.stderr, "Opening file: " + annotationFilename
    for line in annotationFile:
      line = line.rstrip ()
      if line != "":
        annotationFields = line.split("\t")
        header = [f.lower() for f in annotationFields]
        if "sample name" in header:
           sampleNameIndex  = header.index("sample name")
           fastqPathIndex   = header.index("fastq path")
           fastqNameIndex   = header.index("fastq name")
           if "library insert size (bp)" in header:
             insertSizeIndex = header.index("library insert size (bp)")
           else:
             print >> sys.stderr, "WARNING: field \"library insert size (bp)\" not found in file " + annotationFilename
             insertSizeIndex = -1
        else:
          sampleName = annotationFields[sampleNameIndex].replace(" ", "-").replace("--", "-")
          if "&" in sampleName or "|" in sampleName or "*" in sampleName or "(" in sampleName or ")" in sampleName or \
             "$" in sampleName or "!" in sampleName or "@" in sampleName or "#" in sampleName or "%" in sampleName or \
             '"' in sampleName or "'" in sampleName or ";" in sampleName or "/" in sampleName or ">" in sampleName or \
             "<" in sampleName or "?" in sampleName:
            raise Exception ("Sample name contains forbidden characters (&|*()$!@#%\"';/><?): " + sampleName)
          fastqDir   = annotationFields[fastqPathIndex].rstrip("/")
          if not singleReadData:
            fastqFile1 = annotationFields[fastqNameIndex].replace("2.fastq", "1.fastq")
            fastqFile2 = annotationFields[fastqNameIndex].replace("1.fastq", "2.fastq")
          else:
            fastqFile1 = annotationFields[fastqNameIndex]
            fastqFile2 = None
          
          if insertSizeIndex != -1 and insertSizeFilename != "":
            insertSize = annotationFields[insertSizeIndex]
            print >> insertSizeFile, "\t".join([sampleName, insertSize])
          else:
            insertSize = -1
          
          if not sampleName in sampleFields:
            sampleFields[sampleName] = [[insertSize, fastqDir, fastqFile1, fastqFile2]]
          else:
            sampleFields[sampleName].append([insertSize, fastqDir, fastqFile1, fastqFile2])
    
    annotationFile.close()
  
  except IOError, e:
    print annotationFilename + " cannot be opened ... skipping"
    print "Error: " + str(e[0]) + " " + str(e[1])
    
  
################################################################################
##
## Create file links for samples
##
################################################################################

gzipSuffix=".gz"
for sampleName in sampleFields:
  numFastqFiles = len(sampleFields[sampleName])
  print >> sys.stderr, "Creating " + str(numFastqFiles) + " links for sample " + sampleName
  sampleIndex = 1
  for sampleField in sampleFields[sampleName]:
    if numFastqFiles > 1:
      sampleFileIndex = "_" + str(sampleIndex)
    else:
      sampleFileIndex = ""
    
    insertSize, fastqDir, fastqFile1, fastqFile2 = sampleField

    if not gzipSuffix in fastqFile1:
      print >> sys.stderr, "File " + fastqFile1 + " is not gzipped ... skipping"
    else:
      if not sampleName in sampleNames:
        print >> sys.stderr, "mkdir -p " + samplesDir + "/" + sampleName
        subprocess.call("mkdir -p " + samplesDir + "/" + sampleName, shell=True)
        sampleNames.append(sampleName)

      if os.path.exists(fastqDir + "/" + fastqFile1):
        fastqLink1 = samplesDir + "/" + sampleName + "/" + sampleName + sampleFileIndex + "_1.fastq" + gzipSuffix
        if not os.path.exists(fastqLink1):
          cmd1 = "ln -s " + fastqDir + "/" + fastqFile1 + " " + fastqLink1
          print >> sys.stderr, cmd1
          subprocess.call(cmd1, shell=True)
        else:
          print >> sys.stderr, fastqLink1 + " exists ... skipping"
      else:
        raise Exception ("File " + fastqDir + "/" + fastqFile1 + " not found ... exiting")

      if not singleReadData:
        if os.path.exists(fastqDir + "/" + fastqFile2):
          fastqLink2 = samplesDir + "/" + sampleName + "/" + sampleName + sampleFileIndex + "_2.fastq" + gzipSuffix
          if not os.path.exists(fastqLink2):
            cmd2 = "ln -s " + fastqDir + "/" + fastqFile2 + " " + fastqLink2
            print >> sys.stderr, cmd2
            subprocess.call(cmd2, shell=True)
          else:
            print >> sys.stderr, fastqLink2 + " exists ... skipping"
        else:
          raise Exception ("File " + fastqDir + "/" + fastqFile2 + " not found ... exiting")
        
    sampleIndex = sampleIndex + 1
 
