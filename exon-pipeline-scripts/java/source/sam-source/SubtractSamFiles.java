/**File: SubtractSamFiles.java 

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
 *                           Class SubtractSamFiles
 *
 *  
 ***********************************************************************************/


public class SubtractSamFiles {

  private static int debugLevel = 0;
  private static boolean warningsOn = false;


  /***********************************************************************************/

  private static void printHelp () {
    System.out.println("SubtractSamFiles.java\n" +                                              
     "USAGE: java SubtractSamFiles [-s <slack constant>] [-t <distance threshold>]\n" +
     "       -1 <sam file 1> -2 <sam file 2> -w <read weight file1>\n" +
     "       -W <read weight file2> -o <output file> [-O <output weight file>]\n" +
     "\n" +
     "Reads two SAM files and selects the better alignment for each read based on the \n" +
     "edit distance (see <sam file 1> and <sam file 2> for more details). The reads\n" +
     "should be in (roughly) the same order in both files. Small differences due to\n" +
     "multi-threaded execution of the alignment are compensated for by storing the\n" +
     "SAM records of the minor file (the one whose index is different from\n" +
     "<main file index>) in a hash table.\n" +
     "\n" +
     "<sam file 1>: the first SAM file (- for STDIN). A SAM record is kept if the\n" +
     "   av. edit distance of the alignment is at most the av. edit distance of the\n" +
     "   alignment for the read in the second file minus a slack constant (0.25, see\n" +
     "   option -s). Paired-end alignments are kept if there are only single read\n" +
     "   alignments in the second file (independent of the edit distance). [deflt: -]\n" +
     "<sam file 2>: the secound SAM file (- for STDIN).\n" +
     "<output file>: the output SAM file without header (- for STDOUT) [default: -].\n" +
     "<output weight file>: file with the weights for the SAM entries that are\n" +
     "              reported.\n" +
     "<read weight file1>: the read weight file for sam file 1.\n" +
     "<read weight file2>: the read weight file for sam file 2.\n" +
     "-s DOUBLE: slack constant for the comparison of the av. edit distance of the SAM\n" +
     "           records; twice for paired-end reads [0.25].\n" + 
     "-t INT: edit distance threshold for aligned reads. Alignments with an edit\n" +
     "        distance more than the threshold are not reported (-1 for no\n" +
     "        threshold [default: -1]\n");
  }
                                    

  /***********************************************************************************/

  public static void main (String [] args) {

    boolean printLines = false;

    String samFilename1 = "-";
    String samFilename2 = "-";
    
    BufferedReader samReader1 = null;
    BufferedReader samReader2 = null;
    
    String  readWeightFilename1 = "";
    String  readWeightFilename2 = "";
    String  outputFilename = "-";
    String  outputWeightFilename = "";
    boolean pairedEndOnly = false;
    double  distanceSlack = 0.25;
    String  readIdCutOffString = ":";
    int     distanceThreshold = -1;

    final int countUnit = 5 * 1000 * 1000;

    Getopt g = new Getopt("SubtractSamFiles", args, "1:2:d:o:O:s:t:w:W:h");
    
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
      case 'd':
	debugLevel = Integer.parseInt (g.getOptarg());
	UtilLib.setDebugLevel (debugLevel);
	CombineSamFiles.setDebugLevel (debugLevel);
	break;
      case 'o':
	outputFilename = g.getOptarg();
	break;
      case 'O':
	outputWeightFilename = g.getOptarg();
	break;
      case 's':
	distanceSlack = new Double(g.getOptarg()).doubleValue ();
	break;
      case 't':
	distanceThreshold = Integer.parseInt (g.getOptarg());
	break;
      case 'w':
	readWeightFilename1 = g.getOptarg();
	break;
      case 'W':
	readWeightFilename2 = g.getOptarg();
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


    String line = "";
    try {

      if (samFilename1.equals ("-") && samFilename2.equals ("-")) {
	throw new IOException ("ERROR: One of the input SAM must be different from '-' (option -1 and -2).");
      }

      if (readWeightFilename1.equals ("")) {
	throw new IOException ("ERROR: Read weight file 1 must be specified (option -w).");
      }

      if (readWeightFilename2.equals ("")) {
	throw new IOException ("ERROR: Read weight file 2 must be specified (option -W).");
      }


      if (debugLevel >= 2) {
	System.out.println("File name 1: " + samFilename1 + ", file name 2: " + samFilename2);
      }

      PrintWriter outputWriter       = UtilLib.getPrintWriter (outputFilename);
      PrintWriter outputWeightWriter = null;
      if (outputWeightFilename != "") {
	System.err.println ("Writing combined read weights to " + outputWeightFilename);
	outputWeightWriter = UtilLib.getPrintWriter (outputWeightFilename);
      }


      samReader1 = UtilLib.getBufferedReader (samFilename1);
      CombineSamFiles.setSamReader1 (samReader1);
      samReader2 = UtilLib.getBufferedReader (samFilename2);

      BufferedReader weightReader1 = UtilLib.getBufferedReader (readWeightFilename1);
      BufferedReader weightReader2 = UtilLib.getBufferedReader (readWeightFilename2);

      FragmentEntry fragmentEntry1 = null;
      FragmentEntry fragmentEntry2 = null;

      Vector<SamRecord> samRecords1 = null;
      Vector<SamRecord> samRecords2 = null;
      
      samRecords1 = CombineSamFiles.getSamRecords (samReader1);
      samRecords2 = CombineSamFiles.getSamRecords (samReader2);
      
      String fragmentName    = null;
      String oldFragmentName = "";
	
      int numSelected1 = 0;
      int numTotal1    = 0;
      int numLines     = 0;
      while (samRecords1 != null && samRecords2 != null) {

	if (debugLevel >= 2) {
	  System.out.println ("SamRecords1: " + samRecords1);
	  System.out.println ("SamRecords2: " + samRecords2);
	}

	if (samRecords1 != null) {
	  fragmentName = samRecords1.firstElement().getFragmentName();
	} else if (samRecords2 != null) {
	  fragmentName = samRecords2.firstElement().getFragmentName();
	} else {
	  throw new IOException ("ERROR: samRecords1 and samRecords2 are null.");
	}

	if (oldFragmentName.equals(fragmentName)) {
	  throw new IOException ("ERROR: Duplicate fragment name: old " + oldFragmentName + " vs new " + fragmentName);
	}
	
	if (debugLevel >= 2) {
	  System.out.println("Fragment name: " + fragmentName);
	}
	
	boolean isMapped1 = CombineSamFiles.isMapped (samRecords1);
	boolean isMapped2 = CombineSamFiles.isMapped (samRecords2);

	fragmentEntry1 = null;
	fragmentEntry2 = null;

	if (isMapped1) {
	  if (! isMapped2) {
	    numTotal1++;
	    numSelected1++;
	    // System.out.println("2. fragmentName: " + fragmentName);
	    fragmentEntry1 = CombineSamFiles.getFragmentEntry (weightReader1, fragmentName);
	    CombineSamFiles.printSamRecords (outputWriter, fragmentName, samRecords1, fragmentEntry1, distanceThreshold, outputWeightWriter);
	  } else {
	    
	    // System.out.println("4. fragmentName: " + fragmentName);
	    fragmentEntry1 = CombineSamFiles.getFragmentEntry (weightReader1, fragmentName);
	    numTotal1++;
	    fragmentEntry2 = CombineSamFiles.getFragmentEntry (weightReader2, fragmentName);

	    if (fragmentEntry1 != null) {
	      if (fragmentEntry2 == null) {		
		CombineSamFiles.printSamRecords (outputWriter, fragmentName, samRecords1, fragmentEntry1, distanceThreshold, outputWeightWriter);
		numSelected1++;
		
	      } else {

		int compMultiplicity = fragmentEntry1.compareMultiplicity (fragmentEntry2);

		/* Note that compMultiplicity is non-zero if and only if the multiplicity of exactly one of the two fragmentEntries
		   is two. */
		if (compMultiplicity == 1) {
		  numSelected1++;
		  CombineSamFiles.printSamRecords (outputWriter, fragmentName, samRecords1, fragmentEntry1, distanceThreshold, outputWeightWriter);
		} else if (compMultiplicity == 0) {

		  if (fragmentEntry1.getEditDistance () < 0) {
		    throw new IOException ("ERROR: No edit distance for " + fragmentName + " in " + readWeightFilename1 + " given.");
		  }

		  if (fragmentEntry2.getEditDistance () < 0) {
		    throw new IOException ("ERROR: No edit distance for " + fragmentName + " in " + readWeightFilename2 + " given.");
		  }

		  /* Note that the first samRecord in printSamRecords is the preferred one whereas we need the second to be
		     the preferred on in SubtractSamFiles */
		  int compDist = fragmentEntry2.compareEditDistance (fragmentEntry1, distanceSlack);

		  /* If compDist == 1, then fragmentEntry2 wins, if compDist == -1, then fragmentEntry1 wins */
		  if (compDist == -1) {
		    CombineSamFiles.printSamRecords (outputWriter, fragmentName, samRecords1, fragmentEntry1, distanceThreshold, outputWeightWriter);
		    numSelected1++;      
		  }
		}
	      }
	    } else {
	      throw new IOException ("No weight information found for fragment" + fragmentName);
	    }
	  }
	}
		      
	oldFragmentName = fragmentName;

	// System.out.println("Getting samRecords1 after " + fragmentName + " (w/o comparison)");
	samRecords1 = CombineSamFiles.getSamRecords (samReader1, samRecords1);
	// System.out.println("Getting samRecords2 after " + fragmentName + " (w/o comparison)");
	samRecords2 = CombineSamFiles.getSamRecords (samReader2, samRecords2);

	numLines++;
	if (numLines % countUnit == 0) {
	  System.err.print(".");
	  System.err.flush();
	}

	int samRecordCompValue = CombineSamFiles.getCompValue(samRecords1, samRecords2);
	while (samRecordCompValue != 0) {
	  
	  oldFragmentName = fragmentName;
	  
	  if (samRecordCompValue > 0) {
	    fragmentName = samRecords2.firstElement().getFragmentName();
	    samRecords2 = CombineSamFiles.getSamRecords (samReader2, samRecords2);
	    
	  } else {
	    
	    numTotal1++;
	    numSelected1++;

	    fragmentName = samRecords1.firstElement().getFragmentName();
	    FragmentEntry fragmentEntry = null;
	    if (CombineSamFiles.isMapped (samRecords1)) {
	      fragmentEntry = CombineSamFiles.getFragmentEntry (weightReader1, fragmentName);
	    }
	    CombineSamFiles.printSamRecords (outputWriter, fragmentName, samRecords1, fragmentEntry, distanceThreshold, outputWeightWriter);

	    // System.out.println("Getting samRecords1 after " + fragmentName);
	    samRecords1 = CombineSamFiles.getSamRecords (samReader1, samRecords1);
	  }

	  // System.out.println("Comparing samRecords1 and samRecords2.");
	  samRecordCompValue = CombineSamFiles.getCompValue(samRecords1, samRecords2);

	}

	numLines++;
	if (numLines % countUnit == 0) {
	  System.err.print(".");
	  System.err.flush();
	}
	
      }

      if (samRecords1 != null) {
	throw new IOException ("ERROR: First file: " + samFilename1 + " not completely exhausted: " + samRecords1.size() + " elements still remaining.\n " +
			       "Starting with: " + samRecords1);
      }

      if (samRecords2 != null) {
	throw new IOException ("ERROR: Second file: " + samFilename2 + " not completely exhausted: " + samRecords2.size() + " elements still remaining.\n " +
			       "Starting with: " + samRecords2.firstElement ());
      }


      if (numLines >= countUnit) {
	System.err.println();
      }

      System.err.println(numSelected1 + " records of " + numTotal1 + " aligned records selected from file " + samFilename1);

      samReader1.close();
      samReader2.close();

      weightReader1.close();
      weightReader2.close();
      
      outputWriter.close ();

      if (outputWeightWriter != null) {
	outputWeightWriter.close();
      }

    }
    catch (IOException e) {
      System.err.println ("Problem in line: " + line + ": " + e==null?"No error message":e.getMessage());
    }
  }
}

