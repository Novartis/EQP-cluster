/**File: ComputeReadWeights.java 

Original Author: Sven Schuierer
Date: 09/12/2011

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
 *                           Class ComputeReadWeights
 *
 *  
 ***********************************************************************************/


public class ComputeReadWeights {

  private static int debugLevel = 0;

  /***********************************************************************************/

   private static void printHelp () {
    System.out.println("ComputeReadWeights.java\n" +                                              
    "USAGE: ComputeReadWeights [-w] -b <intersection bed file> -o <outputFile>\n" +
    "\n" +
    " intersection bed file: the bed file containing the intersection of the\n" +
    "                        exons transcript intervals and the reads mapped\n" +
    "                        to the transcripts (- for STDIN) [-].\n" +
    " output file: the file to which the output is written (- for STDOUT) [-].\n" +
    " -w: display warning messages.\n" +
    "\n" +
    "Example: ComputeReadWeights sma-sample1-refseq-mm9.bed 51 sma-sample1-refseq-mm9.wgt\n");


  }
                                     

  /***********************************************************************************/

  public static void main (String [] args) {

    int numDeletions  = 0;
    int numInsertions = 0;

    String separator = "#";
    
    String intersectionBedFilename = "-";
    String outputFilename = "-";

    boolean printLines = false;
    
    Getopt g = new Getopt("ComputeReadWeights.java", args, "b:d:Go:t:Wh");
    
    int c;
    String arg = "";
    
    c = g.getopt();
    
    while (c  != -1) {
      switch(c) {
      case 'b':
	intersectionBedFilename = g.getOptarg();
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

      BufferedReader reader       = UtilLib.getBufferedReader (intersectionBedFilename);
      PrintWriter    weightWriter = UtilLib.getPrintWriter (outputFilename);

      if (debugLevel >= 2) {
	System.out.println("Options read.");
      }
      
      line = reader.readLine();

      int lineNumber = 1;
      
      String oldFragmentId = "";

      WeightObjectAlignment weightObjectAlignment = null;
      WeightObject          weightObject = null;
      
      int countUnit = 5 * 1000 * 1000;
      while (line != null) {
	
	BedRecord bedRecord  = new BedRecord (line);

	if (bedRecord.getOverlap () > 0) {
	
	  /* We consider fragments as weight objects. If there are two single-read alignments of one fragment,
	     then we add the edit distance of both and count this as two alignments of the weight object (not
	     as the alignment of a different weight object or even the same aligment). */
	  String  fragmentId = bedRecord.getFragmentId ();

	  if (printLines) {
	    System.out.println (bedRecord.toString());
	  }

	  if (weightObject == null) {
	    
	    weightObject = new WeightObject (bedRecord);
	    
	  } else if (! fragmentId.equals(oldFragmentId)) {

	    /* Finish the old weight object */
	    weightWriter.println(weightObject.toWeightString ());

	    if (debugLevel >= 2) {
	      System.err.println("BED entries:");
	      System.err.println(weightObject.toBedString ());
	    }

	    /* Create a new weight object */
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
	
      }

      /* process last weightObject */
      if (weightObject != null) {
	if (UtilLib.warningsOn ()) {
	  System.err.println ("\nPrinting last weight object: " + weightObject);
	}
	weightWriter.println(weightObject.toWeightString ());
	if (debugLevel >= 2) {
	  System.err.println("BED entries:");
	  System.err.println(weightObject.toBedString ());
	}
      } 

      if (lineNumber >= countUnit) {
	System.err.println();
      }

      weightWriter.close ();
      reader.close();

    }
    catch (Exception e) {
      System.out.println ("Problem in line: " + line + ": " + e==null?"No error message":e.getMessage());
      System.exit (1);
    }
  }
}

