/**File: ConvertSamBed.java 

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
 *                           Class ConvertSamBed
 *
 *  
 ***********************************************************************************/


public class ConvertSamBed {


  private static int countUnit = 5 * 1000 * 1000;

   
  /***********************************************************************************
   * 
   *                           loadBedFile
   *
   ***********************************************************************************/

  public static HashSetTable<String, String> loadBedFile (String transcriptBedFilename) throws IOException {

    if (transcriptBedFilename == "") {
      return (null);
    }
    
    BufferedReader reader = UtilLib.getBufferedReader(transcriptBedFilename);

    String line = reader.readLine();
    HashSetTable<String, String> transcriptVersionIds = new  HashSetTable<String, String> (500000);
    String transcriptId = "";
    int lineNumber = 1;
    while (line != null) {
      StringTokenizer st = new StringTokenizer (line, "\t");
      if (st.hasMoreTokens()) {
	String transcriptVersionId = st.nextToken ();
	transcriptId = transcriptVersionId;

	/* Remove invert strand suffix from transcript id */
	String invertStrandSuffix = UtilLib.getInvertStrandSuffix ();
	if (transcriptId.endsWith (invertStrandSuffix)) {
	  transcriptId = transcriptId.substring(0, transcriptId.length() - invertStrandSuffix.length());
	}
	
	/* Remove version from transcript id */
	int lastDotIndex = transcriptId.lastIndexOf('.');
	if (lastDotIndex != -1 && transcriptId.length() - lastDotIndex <= 4) {
	  transcriptId = transcriptId.substring(0, lastDotIndex);
	}

	transcriptVersionIds.putValue (transcriptId, transcriptVersionId);
	
      }


      if (lineNumber % countUnit == 0 && false) {
	System.err.print(".");
      }

      line = reader.readLine();
      lineNumber++;
      
    }

    if (lineNumber > countUnit) {
      System.err.println();
    }

    return (transcriptVersionIds);
    
  }
   

  /***********************************************************************************/

   private static void printHelp () {
    System.err.println("ConvertSamBed.java\n" +                                              
      "   -- Script to convert SAM files to BED files\n" +
      "\n" +
    "USAGE: ConvertSamBed [-n] [-N] [-M <mapping file>] [-T <transcript bed file>]\n" +
    "       [-c <read id cut-off string>] [-C <common prefix>] [-s <sam file>]\n" +
    "       [-S] [-o <output bed file>] [-e <edit distance threshold>]\n" +
    "\n" +
    "-n: create new identifiers of the form F0000001\n" +
    "-N: do not splice reads + every alignment is one contiguous interval\n" +
    "-M STRING: mapping file: a file with mappings sam fragment ids to bed fragment\n" +
    "    ids\n" +
    "-T STRING: transcript bed file - the bed file for the exon to transcript \n" +
    "  mappings (used to identify transcript ids with version numbers (.1, .2, etc)\n" +
    "  for which the SAM alignments are duplicated)\n" +
    "-c STRING: read id cut-off string - The part of a read id before the first\n" +
    "   occurence of this string is cut off for the bed file read ids. [<empty>]\n" +
    "-C STRING: common prefix - Common prefix of all read ids which is removed.\n" +
    "   [<empty>]\n" +
    "-s STRING: sam file - file with the sam alignments (- for stdin) [-]\n" +
    "-S: convert strand-specific, i.e. invert the strand of the second read of a\n" +
    "    pair.\n" +
    "-o STRING: output bed file - file with the bed records (- for stdout) [-]\n" + 
    "-e INT: edit distance threshold - convert only SAM records with an edit distance\n" + 
    "     of at most INT\n" + 
    "-p: convert only primary alignments\n" +    "\n" +
    "Reads a SAM file from STDIN or from <sam file> and output a BED file with\n" +
    "converted fragment ids (of the form F0000001) to STDOUT\n");
  }
                                     

  /***********************************************************************************/

  public static void main (String [] args) {

    int     debugLevel = 0;
    
    String  mappingFilename = "";
    String  transcriptBedFilename = "";
    String  samFilename = "-";
    String  bedFilename = "-";

    String  idCutOffString = "";
    String  commonPrefix = "";
       
    boolean newIdentifiers = false;
    boolean allowIncompletePairs = false;
    boolean noSplice = false;
    boolean warningsOn = false;

    int    alignmentScoreThreshold = Integer.MIN_VALUE;
    int    editDistanceThreshold   = Integer.MAX_VALUE;

    String strandSpecificDirection = "none";

    boolean primaryAlignmentsOnly = false;
      
    boolean printLines = false;
    
    Getopt g = new Getopt("ConvertSamBed.java", args, "a:c:C:d:e:nNo:pPs:S:T:wh");
    
    int c;
    String arg = "";
    
    c = g.getopt();
    
    while (c  != -1) {
      switch(c) {
      case 'c':
	idCutOffString = g.getOptarg();
	newIdentifiers = true;
	break;
      case 'a':
	alignmentScoreThreshold = (new Integer (g.getOptarg())).intValue ();
	break;	
      case 'C':
	commonPrefix = g.getOptarg();
	newIdentifiers = true;
	break;
      case 'd':
	arg = g.getOptarg();
	debugLevel = (new Integer(arg)).intValue();
	UtilLib.setDebugLevel (debugLevel);
	break;
      case 'e':
	arg = g.getOptarg();
	editDistanceThreshold = (new Integer(arg)).intValue();
	break;
      case 'M':
	mappingFilename = g.getOptarg();
	break;
      case 'n':
	newIdentifiers = true;
	break;
      case 'N':
	noSplice = true;
	break;	
      case 'o':
	bedFilename = g.getOptarg();
	break;
      case 'p':
	primaryAlignmentsOnly = true;
	break;
      case 'P':
	printLines = true;
	break;
      case 's':
	samFilename = g.getOptarg();
	break;
      case 'S':
	strandSpecificDirection = g.getOptarg();
	break;
      case 'T':
	transcriptBedFilename = g.getOptarg();
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
    
    SamRecord.init ();
    SamRecord.setNoSplice (noSplice);

    UtilLib.setIdCutOffString (idCutOffString);
    UtilLib.setCommonPrefix   (commonPrefix);

    if (debugLevel >= 2) {
      System.err.println("Options read.");
    }

    String line = "";
    try {
      
      BufferedReader reader        = UtilLib.getBufferedReader (samFilename);
      PrintWriter    mappingWriter = UtilLib.getPrintWriter    (mappingFilename);
      PrintWriter    bedWriter     = UtilLib.getPrintWriter    (bedFilename);

      HashSetTable<String, String> transcriptVersionIds = null;
      if (! "".equals(transcriptBedFilename)) {
	System.err.println ("Loading file " + transcriptBedFilename);
	transcriptVersionIds = loadBedFile (transcriptBedFilename);
      }

      SamProcessorBed samProcessorBed =
	new SamProcessorBed (bedWriter, newIdentifiers, mappingWriter, transcriptVersionIds, alignmentScoreThreshold,
			     editDistanceThreshold, strandSpecificDirection, primaryAlignmentsOnly);
      
      SamReader samReader = new SamReader (reader, samFilename.equals("-")?"std in":samFilename);
      samReader.readSamFile (samProcessorBed);

      bedWriter.close ();

      if (mappingWriter != null) {
	mappingWriter.close ();
      }

      System.err.println ("NUMBER_MAPPED_READS=" + samReader.getNumMappedReads());
      System.err.flush();


    }
    catch (Exception e) {
      System.err.println ("Problem in line: " + line + ": " + e==null?"No error message":e.getMessage());
      System.exit (1);
    }
  }
}

