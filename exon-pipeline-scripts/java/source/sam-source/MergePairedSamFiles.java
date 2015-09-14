/**File: MergePairedSamFiles.java 

Original Author: Sven Schuierer
Date: 18/01/2012

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
 *                           Class MergePairedSamFiles
 *
 *  
 ***********************************************************************************/


public class MergePairedSamFiles {

  private static int debugLevel = 0;

  private static String  samFilename1 = "-";
  private static String  samFilename2 = "-";
  
  private static int lineNumber1 = 0;
  private static int lineNumber2 = 0;
  private static int lineNumber  = 0;
  private static final int countUnit = 5 * 1000 * 1000;

  private static PrintWriter outputWriter = null;

  
  /***********************************************************************************/

  public static SamRecord getNextFragment (BufferedReader reader, String filename, String fragmentName,
					   SamRecord samRecord, boolean printSamRecord, String suffix) throws Exception {

    if (samRecord == null) {
      System.err.println ("WARNING: Call to get next fragment with null samRecord");
      return null;
    }

    if (debugLevel >= 1) {
      System.out.println ("fragmentName: " + fragmentName + ", samRecord fragment name: " + samRecord.getFragmentName() + ", suffix: " + suffix);
    }

    String line = "";
    while (line != null && samRecord.getFragmentName().equals(fragmentName)) {

      if (printSamRecord) {
	if (samRecord.isMapped () && samRecord.hasMate ()) {
	  outputWriter.println (samRecord.toString ());
	} else {
	  if (debugLevel >= 2) {
	    System.out.println ("Writing: " + samRecord.toString (suffix));
	  }
	  outputWriter.println (samRecord.toString (suffix));
	}
      }
      outputWriter.flush();
      
      line = reader.readLine();

      if (debugLevel >= 2) {
	System.out.println ("line: " + line);
      }
      
      if (line == null || line.trim() == "" || line.trim().length() < 50) {
	return null;
      }
      
      if (line.startsWith("@")) {
	throw new Exception ("Encountered line: \n" + line + "\n" + "in record region of " + reader.toString ());
      }
      
      if (debugLevel >= 2) {
	System.out.println ("samRecord fragment name: " + samRecord.getFragmentName());
      }

      if (filename.equals (samFilename1)) {
	lineNumber1++;
	samRecord = new SamRecord (line, lineNumber1, filename);	      
      } else {
	lineNumber2++;
	samRecord = new SamRecord (line, lineNumber2, filename);	      
      }

      lineNumber++;
      if (lineNumber % countUnit == 0) {
	System.err.print(".");
      }

      if (debugLevel >= 2 && lineNumber >= 100 * 1000) {
	System.exit(0);
      }
      
    }

    return (samRecord);

  }




  /***********************************************************************************/

  private static void printHelp () {
    System.out.println("MergePairedSamFiles.java\n" +                                              
     "   -- Script to select the best alignment of two SAM files for the same set of reads.\n" +
     "\n" +
     "USAGE: java MergePairedSamFiles [-b] [-c <read id cut-off string>]\n" +
    "          [-C <common prefix>] [-1 <sam file 1>] [-2 <sam file 2>]\n" +
     "	       [-o <output file>]\n" + 
     "\n" +
     "Read two SAM files and merge them. Preferentially, alignments of the first\n" +
     "file are kept. The reads of the SAM records in both files must have the same\n" +
     "order (even though some records may be missing in the second file.)\n" +
     "\n" +
     "-b: Keep aligned reads from both files. The sequence of reads must be exactly the\n" +
     "  same in both files if specified (though the number of alignments can differ). A\n" +
     "  suffix of \"/1\" and \"/2\" is attached to the identifiers of sam file 1 and\n" +
     "  sam file 2, resp., if not present.\n" +
     "<sam file 1>: the first SAM file (- for STDIN). SAM records are kept of aligned\n" +
     "   reads are kept. [default: -]\n" +
     "<sam file 2>: the second SAM file (- for STDIN). SAM records are kept of aligned\n" +
     "   if they are not aligned in the first file or -b is specified. [default: -].\n" +
     "<output file>: the output SAM file without header (- for STDOUT) [default: -].\n" +
     "-c STRING: read id cut-off string - The part of a read id before the first\n" +
     "   occurence of this string is cut off for the bed file read ids. Ignored\n" +
     "   if -k is specified. [<empty>]\n" +
     "-C STRING: common prefix - Common prefix of all read ids which is removed.\n" +
     "   Ignored if -k is specified. [<empty>]\n");
  }
                                    

  /***********************************************************************************/

  public static void main (String [] args) {

    boolean printLines = false;
    
    String  outputFilename = "-";
    boolean outputBoth = false;

    String  idCutOffString = "";
    String  commonPrefix = "";
    
    Getopt g = new Getopt("FilterSamFile.java", args, "1:2:c:C:d:bo:h");
    
    int c;
    String arg = "";
    
    c = g.getopt();
    
    while (c  != -1) {
      switch(c) {
      case '1':
	samFilename1 = g.getOptarg();
	break;
      case '2':
	samFilename2 = g.getOptarg();
	break;
      case 'b':
	outputBoth = true;
	break;
      case 'c':
	idCutOffString = g.getOptarg();
	break;
      case 'C':
	commonPrefix = g.getOptarg();
	break;
      case 'd':
	debugLevel = Integer.parseInt(g.getOptarg());
	break;
      case 'o':
	outputFilename = g.getOptarg();
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

    if (debugLevel >= 2) {
      System.out.println("Options read.");
    }

    UtilLib.setIdCutOffString (idCutOffString);
    UtilLib.setCommonPrefix   (commonPrefix);

    String line1 = "";
    String line2 = "";

    String fragmentName1 = null;	
    String fragmentName2 = null;

    SamRecord samRecord1 = null;
    SamRecord samRecord2 = null;

    String suffix1 = "";
    String suffix2 = "";

    if (outputBoth) {
      suffix1 = "/1";
      suffix2 = "/2";
    }

    try {

      if (samFilename1.equals ("-") && samFilename2.equals ("-")) {
	throw new IOException ("One of the input SAM must be different from '-' (option -1 and -2).");
      }

      if (samFilename1.equals ("-")) {
	System.err.println ("Reading SAM file 1 from std in.");
      }

      if (samFilename2.equals ("-")) {
	System.err.println ("Reading SAM file 2 from std in.");
      }

      if (debugLevel >= 2) {
	System.out.println("File name 1: " + samFilename1 + ", file name 2: " + samFilename2);
      }

      BufferedReader reader1   = UtilLib.getBufferedReader (samFilename1);
      BufferedReader reader2   = UtilLib.getBufferedReader (samFilename2);

      outputWriter = UtilLib.getPrintWriter (outputFilename);

      line1 = reader1.readLine();
      while (line1 != null && line1.startsWith("@")) {
	outputWriter.println (line1);
	line1 = reader1.readLine();
	lineNumber1++;
      }

      line2 = reader2.readLine();
      while (line2 != null && line2.startsWith("@")) {	
	line2 = reader2.readLine();
	lineNumber2++;
      }

      if (line1 != null && line2 != null) {
	
	samRecord1 = new SamRecord (line1);
	samRecord2 = new SamRecord (line2);

	int lastLineOutputFile2 = 0;
	while (samRecord1 != null && samRecord2 != null) {
	
	  fragmentName1 = samRecord1.getFragmentName();	
	  fragmentName2 = samRecord2.getFragmentName();

	  if (fragmentName1.indexOf(":1:120:18760:21434#0") != -1 && false) {
	    debugLevel = 2;
	  }

	  if ((debugLevel >= 2 || lineNumber - lastLineOutputFile2 > 200000) && false) {
	    System.out.println("fragmentName1: " + fragmentName1 + ", fragmentName2: " + fragmentName2);
	  }

	  if (fragmentName1 == null || fragmentName2 == null) {
	      throw new Exception ("ERROR: One fragment name is null: " + fragmentName1 + " vs " + fragmentName2);
	    }

	  if (outputBoth) {
	  
	    if (! fragmentName1.equals(fragmentName2)) {
	      throw new Exception ("Fragment name difference: " + fragmentName1 + " vs " + fragmentName2);
	    }
	  
	    samRecord1 = getNextFragment (reader1, samFilename1, fragmentName1, samRecord1, true, suffix1);
	    samRecord2 = getNextFragment (reader2, samFilename2, fragmentName2, samRecord2, true, suffix2);
	  
	  } else {

	    if (samRecord1.isMapped ()) {

	      samRecord1 = getNextFragment (reader1, samFilename1, fragmentName1, samRecord1, true, suffix1);
	      if (fragmentName1.equals(fragmentName2)) {
		samRecord2 = getNextFragment (reader2, samFilename2, fragmentName2, samRecord2, false, suffix2);
	      }
	      
	    } else {
	    
	      if (fragmentName1.equals(fragmentName2)) {
		samRecord1 = getNextFragment (reader1, samFilename1, fragmentName1, samRecord1, false, suffix1);
		if (debugLevel >= 2) {
		  System.out.println("Outputting records of " + samFilename2 + " for " + fragmentName2);
		}
		samRecord2 = getNextFragment (reader2, samFilename2, fragmentName2, samRecord2, true, suffix2);
		lastLineOutputFile2 = lineNumber;
	      } else {
		samRecord1 = getNextFragment (reader1, samFilename1, fragmentName1, samRecord1, true, suffix1);
	      }

	    }
	    
	  }
	}
      }

      System.err.println();
      System.err.println("Done with main loop.");
      System.err.flush();

      /* Process potentially remaining entries */
      if (samRecord1 != null) {
	System.err.println();
	System.err.println("Processing remaining entries of " + samFilename1 + " starting with: " + samRecord1.getQueryName());
	System.err.flush();

	while (samRecord1 != null) {
	  fragmentName1 = samRecord1.getFragmentName();	
	  samRecord1 = getNextFragment (reader1, samFilename1, fragmentName1, samRecord1, true, suffix1);
	}
      }
      
      if (samRecord2 != null) {
	System.err.println();
	System.err.println("Processing remaining entries of " + samFilename2 + " starting with: " + samRecord2.getQueryName());
	System.err.flush();

	while (samRecord2 != null) {
	  fragmentName2 = samRecord2.getFragmentName();	
	  samRecord2 = getNextFragment (reader2, samFilename2, fragmentName2, samRecord2, true, suffix2);
	}
      }
      

      System.err.println();

      reader1.close();
      reader2.close();
      
      outputWriter.close ();

    }
    catch (Exception e) {
      System.err.println ("Problem in line: " + line1 + ": " + e==null?"No error message":e.getMessage());
    }
  }
}

