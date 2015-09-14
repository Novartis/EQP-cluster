/**File: ComputeGenomicBedEntries.java 

Original Author: Sven Schuierer
Date: 19/07/2012

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
 *                           Class ComputeGenomicBedEntries
 *
 *  
 ***********************************************************************************/


public class ComputeGenomicBedEntries {

  private static int debugLevel = 0;


  /***********************************************************************************/

  private static void printHelp () {
    System.out.println("ComputeGenomicBedEntries\n" +                                              
    "USAGE: ComputeGenomicBedEntries -b <intersect. bed file> -o <outputFile>\n" +
    "\n" +
    " -b STRING: intersect. bed file - the bed file containing the intersection of\n" +
    "     the exons transcript intervals and the reads mapped to the transcripts.\n" +
    "     (- for STDIN) [-]\n" +
    " -o STRING: output file - the file to which the output is written (- for STDOUT)\n" + 
    "    [-]\n" +
    "\n" +
    "Example: ComputeGenomicBedEntries -b ~/projects/exons/projects/SEQC/samples/SEQC-A-BC01/s_1/bed-files/SEQC-A-BC01-s_1-combined-intersection.bed.gz\n" +
    "            -o ~/projects/exons/projects/SEQC/samples/SEQC-A-BC01/s_1/bed-files/SEQC-A-BC01-s_1-combined-intersection-genomic.bed.gz");

  }
                                     

  /***********************************************************************************/

  public static void main (String [] args) {

    String intersectionFilename = "-";
    String outputFilename = "-";

    String oldFragmentId = "";

    WeightObject weightObject = null;
    final int countUnit = 2 * 1250 * 1000;
    
    Getopt g = new Getopt("ComputeGenomicBedEntries.java", args, "b:d:gm:o:O:pr:R:w:Wh");
    
    int c;
    String arg = "";
    
    c = g.getopt();
    
    while (c  != -1) {
      switch(c) {
      case 'b':
	intersectionFilename = g.getOptarg();
	break;
      case 'd':
	UtilLib.setDebugLevel (Integer.parseInt(g.getOptarg()));
	break;
      case 'o':
	outputFilename = g.getOptarg();
	break;
      case 'W':
	UtilLib.setWarningsOn (true);
	break;	
      case 'h':
	printHelp();
	System.exit(0);
	break;
      default:
	System.out.print("Error: getopt() returned " + c + "\n");
      }
      c = g.getopt();
    }

    debugLevel = UtilLib.getDebugLevel ();
    String line = "";
    try {
      
      BufferedReader reader       = UtilLib.getBufferedReader (intersectionFilename);
      PrintWriter    outputWriter = UtilLib.getPrintWriter    (outputFilename);

      line = reader.readLine();
      int lineNumber = 1;
      System.err.println("Reading bed file: " + intersectionFilename + " (. = " + countUnit + " lines)");

      while (line != null) {
	
	BedRecord bedRecord = new BedRecord (line);

	if (bedRecord.getOverlap () > 0 && ! bedRecord.isTranscriptExonAlignment()) {
	  
	  /* We consider fragments as weight objects. If there are two single-read alignments of one fragment
	     against a count object, then we count this as one. */
	  String  fragmentId = bedRecord.getFragmentId ();
	  	  
	  if (weightObject == null) {
	    
	    weightObject = new WeightObject (bedRecord);
	    
	  } else if (! fragmentId.equals(oldFragmentId)) {
	    
	    /* Finish old weight object */
	    weightObject.printGenomicBedEntries(outputWriter);	  
	    
	    /* Create new weight object */
	    weightObject = new WeightObject (bedRecord);
	    
	  } else {

	    weightObject.addBedRecord (bedRecord);
	    
	  }
	  
	  oldFragmentId = fragmentId;
	}
	  
	if (lineNumber % countUnit == 0) {
	  System.err.print(".");
	}
	lineNumber++;
	  
	
	line = reader.readLine();

	if (debugLevel >= 1) {
	  System.out.println ("Line: " + line);
	}
	
      }

      if (debugLevel >= 1) {
	System.out.println ("Processing last fragment");
      }

      /* process last fragment */
      if (weightObject != null) {

	weightObject.printGenomicBedEntries (outputWriter);	  

      }

      if (lineNumber >= countUnit) {
	System.err.println();
      }
      
      reader.close();
      outputWriter.close ();

    }
    catch (Exception e) {
      System.out.println ("Problem in line: " + line + ": " + e==null?"No error message":e.getMessage());
    }
  }
}

