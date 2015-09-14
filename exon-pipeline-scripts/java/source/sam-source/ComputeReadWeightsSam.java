/**File: ComputeReadWeightsSam.java 

Original Author: Sven Schuierer
Date: 19/12/2011

Copyright 2015 Novartis Institutes for BioMedical Research
Inc.Licensed under the Apache License, Version 2.0 (the "License"); you
may not use this file except in compliance with the License. You may
obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing,
software distributed under the License is distributed on an "AS IS"
BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
implied. See the License for the specific language governing
permissions and limitations under the License.
*/


import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;


/***********************************************************************************
 *
 * 
 *                           Class ComputeReadWeightsSam
 *
 *  
 ***********************************************************************************/


public class ComputeReadWeightsSam {


  private static int debugLevel = 0;


  /***********************************************************************************/

  private static HashSet<String> readChromosomeNames (BufferedReader chromosomeReader) throws IOException {

    HashSet<String> chromosomeNames = new HashSet<String> (200);
   
    String line = chromosomeReader.readLine();
    while (line != null) {
      StringTokenizer st = new StringTokenizer (line, "\t");
      if (st.hasMoreTokens ()) {
	chromosomeNames.add(st.nextToken());
      }
      line = chromosomeReader.readLine();
    }

    return chromosomeNames;
  }
  

  /***********************************************************************************/

  private static void printHelp () {
    System.out.println("ComputeReadWeightsSam.java\n" +                                              
      "   -- Script to compute the number of occurences and the edit distance\n" +
      "      sum for a SAM file\n" +
      "\n" +
      "USAGE: ComputeReadWeightsSam [-c <chromosome file name>] [-s <sam file>]\n" +
      "          [-o <output file>]\n" +
      "\n" +
      "Reads a SAM file from <sam file> and outputs a the number of\n" +
      "alignments and the edit distance for each read\n" +
      "<sam file>: the input SAM file (- for STDIN) [default: -].\n" +
      "<output file>: the output SAM file without header (- for STDOUT) [default: -].\n" +
      "<chromosome file name>: tab-separated file with the chromosome ids in the\n" +
      "   first column\n");
  }

  
  /***********************************************************************************/

  public static void main (String [] args) {

    int     debugLevel = 0;

    boolean pairedEndOnly = false;
    boolean printLines = false;
    String  outputFilename = "-";
    String  samFilename = "-";
    String  chromosomeFilename = "";

    String  idCutOffString = "";
    String  commonPrefix = "";

    boolean warningsOn = false;
   
    Getopt g = new Getopt("ComputeReadWeightsSam.java", args, "c:d:Mo:pr:s:SRwh");
    
    int c;
    String arg = "";
    
    c = g.getopt();
    
    while (c  != -1) {
      switch(c) {
      case 'c':
	chromosomeFilename = g.getOptarg();
	break;
      case 'd':
	UtilLib.setDebugLevel(Integer.parseInt(g.getOptarg()));
	break;
      case 'M':
	FragmentEntry.setOutputNumSoftMaskedBases ();
	break;
      case 'o':
	outputFilename = g.getOptarg();
	break;
      case 's':
	samFilename = g.getOptarg();
	break;
      case 'R':
	SamProcessorFragmentEntry.setReadsMode ();
	break;
      case 'S':
	FragmentEntry.setOutputNumSplicedFragments ();
	break;
      case 'p':
	pairedEndOnly = true;
	break;
      case 'w':
	warningsOn = true;
	break;
      case 'h':
	printHelp();
	System.exit(0);
	break;
      default:
	System.err.print("Error: getopt() returned " + c + "\n");
      }
      c = g.getopt();
    }

    if (debugLevel >= 2) {
      System.out.println("Options read.");
    }

    UtilLib.setIdCutOffString (idCutOffString);
    UtilLib.setCommonPrefix   (commonPrefix);

    String line = "";
    try {

      BufferedReader reader      = UtilLib.getBufferedReader (samFilename);
      PrintWriter    printWriter = UtilLib.getPrintWriter    (outputFilename);

      /*   A SamReader object reads a SAM file line by line. It constructs a SamRecord object for each
	   line that does not start with a "@". SamRecords with the same fragment name are collected.
	   Once the fragment name changes, the collection of SamRecords is sent to/processed by the
	   processSamRecords method of the samProcessor object. */
      SamReader samReader = new SamReader (reader, samFilename.equals("-")?"std in":samFilename);

      /* The reads are processed by fragmentName and only if the fragmentName changes, then a fragmentEntry is
	 output. For each paired-end read a mate is identified. Either the SamRecord of read and mate (for paired-end
	 alignments) or of the read alone are then added to a FragmentEntry. */
      SamProcessorFragmentEntry samProcessor = new SamProcessorFragmentEntry (printWriter);
      samProcessor.setWarningsOn(warningsOn);

      if (! chromosomeFilename.equals("")) {
	BufferedReader chromosomeReader = UtilLib.getBufferedReader (chromosomeFilename);
	HashSet<String> chromosomeIds = readChromosomeNames (chromosomeReader);
	samProcessor.setReferenceSequenceIdSet (chromosomeIds);
      }

      samReader.readSamFile (samProcessor);

      /* Process last fragment entry */
      samProcessor.outputFragmentEntry ();

      reader.close ();
      printWriter.close ();

    }
    catch (Exception e) {
      
      System.err.println ("Problem in line: " + line + ": " + e==null?"No error message":e.getMessage());
      System.exit (1);
      
    }
  }
}

