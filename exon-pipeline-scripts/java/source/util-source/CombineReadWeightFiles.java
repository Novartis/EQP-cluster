/**File: CombineReadWeightFiles.java 

Original Author: Sven Schuierer
Date: 20/01/2012

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
 *                           Class FragmentEntryFactory
 *
 *  
 ***********************************************************************************/


class FragmentEntryFactory {

  private int     debugLevel = 0;
  private String  readWeightFilename = null;
  private String  line               = "";
  private String  fragmentEntryId    = "";
  private boolean fragmentEntrySet   = false;
  
  private FragmentEntry  fragmentEntry = null;
  private BufferedReader weightReader  = null;

  protected static String specialFragmentName = "";

  FragmentEntryFactory (String readWeightFilename) throws IOException {
    this.readWeightFilename = readWeightFilename;
    if (readWeightFilename != null && readWeightFilename != "") {
      weightReader = UtilLib.getBufferedReader (readWeightFilename);
      line         = weightReader.readLine ();
      fragmentEntry = new FragmentEntry ();
    } else {
      line = null;
    }
  }

  public FragmentEntry getNextFragmentEntry () throws IOException {

    fragmentEntry = null;
    
    if (line != null) {
      fragmentEntry = new FragmentEntry (line);
      line = weightReader.readLine ();
    }

    fragmentEntrySet = true;

    return fragmentEntry;
    
  }

  public FragmentEntry getNextFragmentEntry (FragmentEntry queryFragmentEntry) throws IOException {

    if (queryFragmentEntry == null) {
      return null;
    }

    if (weightReader == null) {
      return null;
    }

    String queryFragmentEntryId = queryFragmentEntry.getFragmentName ();
    while ((! fragmentEntrySet || queryFragmentEntry != null) && fragmentEntry != null && fragmentEntryId.compareTo(queryFragmentEntryId) < 0) {
      fragmentEntry = getNextFragmentEntry ();
      
      if (debugLevel >= 1 && queryFragmentEntryId.equals(specialFragmentName)) {
	System.err.println ("fragmentEntry: " + fragmentEntry);
      }

      if (fragmentEntry != null) {
	fragmentEntryId = fragmentEntry.getFragmentName ();
      }
    }

    if (debugLevel >= 2) {
      if (fragmentEntry == null || ! fragmentEntryId.equals(queryFragmentEntryId)) {
	System.err.println ("No fragment " + queryFragmentEntryId + " found in file " + readWeightFilename);
      }
    }
    
    return fragmentEntry;
    
  }

  public void close () throws IOException {
    if (weightReader != null) {
      weightReader.close ();
    }
  }
}


/***********************************************************************************
 *
 * 
 *                           Class CombineReadWeightFiles
 *
 *  
 ***********************************************************************************/


public class CombineReadWeightFiles {

  private static int debugLevel = 0; 
  private static int countUnit = 5 * 1000 * 1000;
  private static double distanceSlack = 0.25;

  private static int fragmentEntryNum1 = 0;
  private static int fragmentEntryNum2 = 0;
  
  private static String readWeightFilename1 = "-";
  private static String readWeightFilename2 = "";

  private static String specialFragmentName = "";

  /***********************************************************************************/
  
  private static int getCompValue (FragmentEntry f1, FragmentEntry f2) throws IOException {

    String s1 = null;
    String s2 = null;

    if (f1 != null) {
      s1 = f1.getFragmentName ();
    }

    if (f2 != null) {
      s2 = f2.getFragmentName ();
    }

    if ((s1 == null || s1.equals("")) && (s2 == null || s2.equals(""))) {
      return 0;
    }

    if (s1 == null) {
      return 1;
    }

    if (s2 == null) {
      return -1;
    }

    return s1.compareTo(s2);
    
  }


  /***********************************************************************************/
  
  private static FragmentEntry getRelevantFragmentEntry (FragmentEntry fragmentEntry1, FragmentEntry fragmentEntry2, FragmentEntry auxFragmentEntry,
							 int compValue, double distanceSlack, boolean useNumAlignments, boolean useEditDistance) throws IOException {

    FragmentEntry returnFragmentEntry = null;

    if (fragmentEntry1 == null && fragmentEntry2 == null) {
      return returnFragmentEntry;
    }

    if (compValue < 0) {
      returnFragmentEntry = fragmentEntry1;
    } else if (compValue > 0) {
      returnFragmentEntry = fragmentEntry2;
    } else /* fragmentEntry1.getFragmentName() == fragmentEntry2.getFragmentName () */ {

      if (fragmentEntry1 == null || fragmentEntry2 == null) {
	throw new IOException ("Comp value: " + compValue + " and fragmentEntry1: " + fragmentEntry1 + " and fragmentEntry2: " + fragmentEntry2);
      }

      int compMultiplicity = 0;
      if (auxFragmentEntry != null) {
	compMultiplicity = auxFragmentEntry.compareMultiplicity (fragmentEntry2);
      }

      if (debugLevel >= 1 && (fragmentEntry1.getFragmentName().equals(specialFragmentName) || fragmentEntry2.getFragmentName().equals(specialFragmentName))) {
	System.err.println ("Aux fragment entry: " + auxFragmentEntry);
	System.err.println ("compMultiplicity: " + compMultiplicity);
      }
      
      if (compMultiplicity == 1) {
	returnFragmentEntry = fragmentEntry1;
      } else if (compMultiplicity == -1) {
	returnFragmentEntry = fragmentEntry2;
      } else {
	
	if (useNumAlignments) {
	  if (fragmentEntry1.getNumAlignments () >= fragmentEntry2.getNumAlignments ()) {	  
	    returnFragmentEntry = fragmentEntry1;
	  } else {
	    returnFragmentEntry = fragmentEntry2;
	  }
	}
	
	if (useEditDistance) {
	  if (fragmentEntry1.getEditDistance () < 0) {
	    throw new IOException ("ERROR: No edit distance for " + fragmentEntry1 + " in " + readWeightFilename1 + " given.");
	  }

	  if (fragmentEntry2.getEditDistance () < 0) {
	    throw new IOException ("ERROR: No edit distance for " + fragmentEntry2 + " in " + readWeightFilename2 + " given.");
	  }

	  if (fragmentEntry1.compareEditDistance (fragmentEntry2, distanceSlack) >= 0) {
	    returnFragmentEntry = fragmentEntry1;
	  } else {
	    returnFragmentEntry = fragmentEntry2;
	  }
	}
      }
    }

    if (returnFragmentEntry == fragmentEntry1) {
      if (debugLevel >= 1 && (fragmentEntry1.getFragmentName().equals(specialFragmentName) || fragmentEntry2.getFragmentName().equals(specialFragmentName))) {
	System.err.println ("Outputting fragment 1");
      }
      fragmentEntryNum1++;
    } else {
      if (debugLevel >= 1 && (fragmentEntry1.getFragmentName().equals(specialFragmentName) || fragmentEntry2.getFragmentName().equals(specialFragmentName))) {
	System.err.println ("Outputting fragment 2");
      }
      fragmentEntryNum2++;
    }

    return returnFragmentEntry;
    
  }


  /***********************************************************************************/

  private static void printHelp () {
    System.out.println("CombineReadWeightFiles.java\n" +                                              
     "   -- Program to combine two weight files. The entries in the new weight files\n" +
     "      will be selected either by the number of alignments (default, the more,\n" +
     "      the better) or by the edit distance (option -e). The priority of chosing\n" +
     "      weight file entries is according to the following categories:\n" +
     "      o paired-end\n" +
     "      o single read both reads mapped\n" +
     "      o single read\n" +
     "\n" +
     "USAGE: java CombineReadWeightFiles [-C <common prefix>] -1 <read weight file1>\n" +
     "        -2 <read weight file2> -o <output file>\n" +
     "\n" +
     "Reads two weight files and combines them into one. The fastq file is used to\n" +
     "obtain the sequence of the reads; the sequence of read ids in both weight files\n" +
     "needs to be a subsequence of the reads in the fastq file.\n" +
     "\n" +
     "-C STRING: common prefix - Common prefix of all read ids which is removed.\n" +
     "   [<empty>]\n" +
     "-e: use the edit distance instead of the number of alignments for chosing the\n" +
     "    fragmentEntry of the combined file. There is a slight bias towards entries\n" +
     "    from <read weight file1>. If two entries fall into the same category (see\n" +
     "     above), then the second entry is chosen either if\n" +
     "    o the edit distance is no larger than for the first entry and the number of\n" +
     "      alignments is larger than for the first entry or (if this does not hold)\n" +
     "    o the edit distance is smaller than the edit distance of the first entry\n" +
      "     minus " + distanceSlack + " per read (This slack value favors entries of the first file.)\n" +
     "-1 STRING: read weight file1 - the first read weight file (- for STDIN) [-].\n" +
     "-2 STRING: read weight file2 - the second read weight file (- for STDIN) [-].\n" +
     "-m INT: weights in file1 are capped at this value; no capping if negative [-1].\n" +
     "-M INT: weights in file2 are capped at this value; no capping if negative [-1].\n" +
     "-o STRING: output file - the read weight file with the maximum of the alignment\n" +
     "   numbers (- for STDOUT) [-].\n");
  }
                                    

  /***********************************************************************************/

  public static void main (String [] args) {

    boolean printLines = false;

    String  auxReadWeightFilename = "";
    String  outputFilename        = "-";
    boolean compareNumAlignmentsOnly = false;
    String  commonPrefix = "";

    boolean useNumAlignments = true;
    boolean useEditDistance  = false;

    int max1 = -1;
    int max2 = -1;
    
    final int countUnit = 5 * 1000 * 1000;

    Getopt g = new Getopt("CombineReadWeightFiles.java", args, "1:2:a:d:C:eF:m:M:o:h");
    
    int c;
    String arg = "";
    
    c = g.getopt();
    
    while (c  != -1) {
      switch(c) {
      case '1':
	readWeightFilename1 = g.getOptarg();
	break;
      case '2':
	readWeightFilename2 = g.getOptarg();
	break;
      case 'a':
	auxReadWeightFilename = g.getOptarg();
	break;
      case 'd':
	debugLevel = Integer.parseInt(g.getOptarg());
	break;
      case 'C':
	commonPrefix = g.getOptarg();
	break;
      case 'e':
	System.err.println("Using edit distance to combine read weight files.");
	useEditDistance = true;
	useNumAlignments = false;
	break;
      case 'F':
	specialFragmentName = g.getOptarg();
	FragmentEntryFactory.specialFragmentName = specialFragmentName;
	break;
      case 'o':
	outputFilename = g.getOptarg();
	break;
      case 'm':
	max1 = Integer.parseInt(g.getOptarg());
	break;
      case 'M':
	max2 = Integer.parseInt(g.getOptarg());
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

    if (commonPrefix != "") {
      UtilLib.setCommonPrefix (commonPrefix);
    }

    
    try {

      if ((readWeightFilename1.equals("-") && readWeightFilename2.equals("-"))) {
	throw new IOException ("Only one of <read weight filename 1> and <read weight filename 2> can be set to std in.");
      }

      FragmentEntryFactory fragmentEntryFactory1   = new FragmentEntryFactory (readWeightFilename1);
      FragmentEntryFactory fragmentEntryFactory2   = new FragmentEntryFactory (readWeightFilename2);
      FragmentEntryFactory auxFragmentEntryFactory = new FragmentEntryFactory (auxReadWeightFilename);
      
      PrintWriter writer = UtilLib.getPrintWriter (outputFilename);

      
      /***********************************************************************************
       *
       * Reading the entries of both weight files
       *
       ***********************************************************************************/
 
      FragmentEntry fragmentEntry1   = fragmentEntryFactory1.getNextFragmentEntry ();
      FragmentEntry fragmentEntry2   = fragmentEntryFactory2.getNextFragmentEntry ();
      FragmentEntry auxFragmentEntry = auxFragmentEntryFactory.getNextFragmentEntry (fragmentEntry1);

      int lineNumber = 1;
      FragmentEntry oldFragmentEntry1 = null;
      FragmentEntry oldFragmentEntry2 = null;

      System.err.println("Combining read weight files:");
      System.err.println(readWeightFilename1 + " and");
      System.err.println(readWeightFilename2);
      
      try {
	boolean specialFragmentNameFound = false;
	while (fragmentEntry1 != null || fragmentEntry2 != null) {

	  if (fragmentEntry1 != null) {
	    fragmentEntry1.capNumAlignments (max1);
	  }
	  
	  if (fragmentEntry2 != null) {
	    fragmentEntry2.capNumAlignments (max2);
	  }

	  int compValue = getCompValue (fragmentEntry1, fragmentEntry2);
	  	  
	  if (debugLevel >= 1 && (fragmentEntry1.getFragmentName().equals(specialFragmentName) || fragmentEntry2.getFragmentName().equals(specialFragmentName))) {
	    System.err.println ("Fragment entry 1: " + fragmentEntry1);
	    System.err.println ("Fragment entry 2: " + fragmentEntry2);
	    System.err.println ("Comp value: " + compValue);
	    specialFragmentNameFound = true;
	  }

	  FragmentEntry outputEntry =
	    getRelevantFragmentEntry (fragmentEntry1, fragmentEntry2, auxFragmentEntry, compValue, distanceSlack, useNumAlignments, useEditDistance);
	  writer.println (outputEntry.toPrintString ());
	  writer.flush();

	  if (compValue <= 0) {
	    oldFragmentEntry1 = fragmentEntry1;
	    fragmentEntry1    = fragmentEntryFactory1.getNextFragmentEntry ();
	    auxFragmentEntry  = auxFragmentEntryFactory.getNextFragmentEntry (fragmentEntry1);
	    if (debugLevel >= 1 && fragmentEntry1 != null && fragmentEntry1.getFragmentName().equals(specialFragmentName)) {
	      System.err.println ("Aux fragment entry obtained: " + auxFragmentEntry);
	    }

	    if (fragmentEntry1 != null && getCompValue(oldFragmentEntry1, fragmentEntry1) >= 0) {
	      throw new IOException ("oldFramentId1: " + oldFragmentEntry1.getFragmentName () + " is not smaller than the new fragmentId1: " +
				     fragmentEntry1.getFragmentName ());
	    }
	  }

	  if (compValue >= 0) {
	    oldFragmentEntry2 = fragmentEntry2;
	    fragmentEntry2    = fragmentEntryFactory2.getNextFragmentEntry ();
	    if (debugLevel >= 1 && fragmentEntry2 != null && fragmentEntry2.getFragmentName().equals(specialFragmentName)) {
	      System.err.println (specialFragmentName + " found for fragmentEntry2");
	    }

	    
	    if (fragmentEntry2 != null && getCompValue(oldFragmentEntry2, fragmentEntry2) >= 0) {
	      throw new IOException ("oldFramentId2: " + oldFragmentEntry2.getFragmentName () + " is not smaller than the new fragmentId2: " +
				     fragmentEntry2.getFragmentName ());
	    }
	  }

	  if (lineNumber % countUnit == 0) {
	    System.err.print (".");
	    System.err.flush();
	  }
	  lineNumber++;
	  
	}
	
	fragmentEntryFactory1.close ();
	fragmentEntryFactory2.close ();
	auxFragmentEntryFactory.close ();
	writer.close ();
	
      }
      catch (IOException e) {
	System.err.println (e==null?"Null message":e.getMessage());
      }
      
      if (lineNumber >= countUnit) {
	System.err.println ();
      }
	

      System.err.println (fragmentEntryNum1 + " read weights are taken from " + readWeightFilename1 + " and\n" +
			  fragmentEntryNum2 + " read weights are taken from " + readWeightFilename2 + "\n" +
			  "Total: " + (fragmentEntryNum1 + fragmentEntryNum2)); 

    }
    catch (IOException e) {
      System.err.println (e==null?"No error message":"At end: " + e.getMessage());
    }
  }
}

