/**File: CombineSamFiles.java 

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
 *                           Class CombineSamFiles
 *
 *  
 ***********************************************************************************/


public class CombineSamFiles {

  private static int debugLevel = 0;
  private static boolean warningsOn = false;

  private static String samFilename1 = "-";
  private static String samFilename2 = "";

  private static int samCounter1 = 0;
  private static int samCounter2 = 0;

  private static String samLine1 = null;
  private static String samLine2 = null;
  
  private static SamRecord samRecord1 = null;
  private static SamRecord samRecord2 = null;

  private static BufferedReader samReader1 = null;
  private static BufferedReader samReader2 = null;

  private static PrintWriter outputWeightWriter = null;

  private static int numOutput = 0;
  private static int numFragmentEntries = 0;

  public static void setDebugLevel (int value) {
    debugLevel = value;
  }


  public static void setSamReader1 (BufferedReader samReader) {
    samReader1 = samReader;
  }

  public static String selectSamLine (BufferedReader samReader) {
    return (samReader.equals(samReader1)?samLine1:samLine2);
  }

  public static SamRecord selectSamRecord (BufferedReader samReader) {
    return (samReader.equals(samReader1)?samRecord1:samRecord2);    
  }

  public static String selectSamFilename (BufferedReader samReader) {
    return (samReader.equals(samReader1)?samFilename1:samFilename2);    
  }

  public static int selectSamCounter (BufferedReader samReader) {
    return (samReader.equals(samReader1)?samCounter1:samCounter2);    
  }

  public static void setSamLineAndRecord (BufferedReader samReader, String s, SamRecord r, int c) {
    if (samReader.equals(samReader1)) {
      samLine1    = s;
      samRecord1  = r;
      samCounter1 = c;
    } else {
      samLine2    = s;
      samRecord2  = r;
      samCounter2 = c;
    }
  }
  
  
  /***********************************************************************************
   * 
   *                           getSamRecords
   *
   ***********************************************************************************/
  
  public static Vector<SamRecord> getSamRecords (BufferedReader samReader) throws IOException {

    Vector<SamRecord> samRecords = null;

    /* We need to take the last line read into account since
       we read one line beyond the lines that we needed to read previously */
    String    samLine     = selectSamLine     (samReader);
    SamRecord samRecord   = selectSamRecord   (samReader);
    String    samFilename = selectSamFilename (samReader);
    int       counter     = selectSamCounter  (samReader);

    if (samLine == null) {
      
      samLine = samReader.readLine ();
      if (samLine == null) {
	return null;
      }
      
      while (samLine.startsWith ("@")) {
	samLine = samReader.readLine ();
	counter++;
      }
      samRecord = new SamRecord (samLine, counter, samFilename);
      // System.err.println("SAM record: " + samRecord.getFragmentName() + " read.");
    }

    if (debugLevel >= 2) {
      System.err.println("Adding: " + samRecord);
    }
    
    if (samLine != null) {

      samRecords = new Vector<SamRecord> ();
      
      String fragmentName = samRecord.getFragmentName ();	
      while (samLine != null && samRecord.getFragmentName ().equals(fragmentName)) {

	if (debugLevel >= 2) {
	  System.err.println("Adding: " + samRecord);
	}
	
	samRecords.add(samRecord);
	samLine = samReader.readLine ();
	counter++;

	if (debugLevel >= 2) {
	  System.err.println("Reading: " + samLine);
	}
	if (samLine != null) {
	  samRecord = new SamRecord (samLine, counter, samFilename);
	  // System.err.println("SAM record: " + samRecord.getFragmentName() + " read.");
	} else {
	  samRecord = null;
	}
	
      }
    }

    /* save the last line */
    setSamLineAndRecord (samReader, samLine, samRecord, counter);

    return samRecords;
    
  }


  /***********************************************************************************
   * 
   *                           getSamRecords
   *
   ***********************************************************************************/
  
  public static Vector<SamRecord> getSamRecords (BufferedReader samReader, Vector<SamRecord> oldSamRecords) throws IOException {

    // System.err.println("First element of oldSamRecords: " + oldSamRecords.firstElement().getFragmentName());

    Vector<SamRecord> samRecords = getSamRecords (samReader);

    if (samRecords != null) {
      if (oldSamRecords == null) {
	throw new IOException ("ERROR: Old sam records null for " + samRecords.firstElement());
      }

      if (samRecords.firstElement() == null || oldSamRecords.firstElement() == null) {
	throw new IOException ("ERROR: Null elements in SAM record collections: " + samRecords + " and " + oldSamRecords);
      }

      if (getCompValue(samRecords, oldSamRecords) < 0) {
	throw new IOException ("ERROR: SAM records out of order " + samRecords.firstElement().getFragmentName() + " and " + oldSamRecords.firstElement().getFragmentName());	
      }
      
    }

    return samRecords;
    
  }


  /***********************************************************************************
   * 
   *                           printSamRecordsSingleRead
   *
   ***********************************************************************************/

  public static void printSamRecordsSingleRead (PrintWriter outputWriter, String fragmentName, Vector<SamRecord> samRecords,
						FragmentEntry fragmentEntry, int distanceThreshold, int readIndex,
						PrintWriter outputWeightWriter) throws IOException {

    if (samRecords == null || outputWriter == null) {
      return;
    }

    boolean samRecordPrinted = false;
    for (int i = 0; i < samRecords.size(); i++) {
      SamRecord samRecord = samRecords.get(i);
      if (samRecord.getReadIndex () == readIndex &&
	  (distanceThreshold == -1 || samRecord.getEditDistance () <= distanceThreshold)) {
	samRecordPrinted = true;
	outputWriter.println (samRecord.toString ());
      }
    }

    if (outputWeightWriter != null && fragmentEntry != null) {
      outputWeightWriter.println(fragmentEntry.toPrintString ());
      numOutput++;
    } 
    
  }


  /***********************************************************************************
   * 
   *                           printSamRecords
   *
   ***********************************************************************************/

  public static void printSamRecords (PrintWriter outputWriter, String fragmentName, Vector<SamRecord> samRecords,
				      FragmentEntry fragmentEntry, int distanceThreshold, PrintWriter outputWeightWriter) throws IOException {

    if (samRecords == null || outputWriter == null) {
      return;
    }

    boolean samRecordPrinted1 = false;
    boolean samRecordPrinted2 = false;

    for (int i = 0; i < samRecords.size(); i++) {
      
      SamRecord samRecord = samRecords.get(i);
      if (samRecord.hasMate () && i < samRecords.size() - 1) {
	i++;
	SamRecord mateSamRecord = samRecords.get(i);
	if (distanceThreshold == -1 || (samRecord.getEditDistance () + mateSamRecord.getEditDistance ()) / 2 <= distanceThreshold) {
	  samRecordPrinted1 = true;
	  samRecordPrinted2 = true;
	  outputWriter.println (samRecord.toString ());
	  outputWriter.println (mateSamRecord.toString ());
	}
      } else {
	if (distanceThreshold == -1 || samRecord.getEditDistance () <= distanceThreshold) {
	  if (samRecord.getReadIndex () == 1) {
	    samRecordPrinted1 = true;
	  } else if (samRecord.getReadIndex () == 2) {
	    samRecordPrinted2 = true;
	  } else {
	    throw new IOException ("ERROR: Sam record not paired end and readIndex != 1 and 2: " + samRecord.toString ());
	  }
	  outputWriter.println (samRecord.toString ());
	}
      }
    }

    if (outputWeightWriter != null && fragmentEntry != null) {
      outputWeightWriter.println(fragmentEntry.toPrintString ());
      numOutput++;
    }

  }


  /***********************************************************************************
   * 
   *                           printSamRecords
   *
   ***********************************************************************************/

  public static void printSamRecords (PrintWriter outputWriter, String fragmentName, Vector<SamRecord> samRecords1,
				       Vector<SamRecord> samRecords2, FragmentEntry fragmentEntry1, FragmentEntry fragmentEntry2,
				       int distanceThreshold, String samFilename1, String samFilename2, PrintWriter outputWeightWriter) throws IOException {

    if (samRecords1 == null || outputWriter == null) {
      return;
    }

    if (samFilename1.equals(samFilename2)) {
      throw new IOException (samFilename1 + " equals samFilename1 and samFilename2");
    }
    
    if (fragmentEntry1.hasMate ()) {
      printSamRecords (outputWriter, fragmentName, samRecords1, fragmentEntry1, distanceThreshold, outputWeightWriter);
      return;
    }

    if (fragmentEntry2.hasMate ()) {
      printSamRecords (outputWriter, fragmentName, samRecords2, fragmentEntry2, distanceThreshold, outputWeightWriter);
      return;
    }

    String readIndices1 = getReadIndexSet (samRecords1);

    if (readIndices1.indexOf ("1") >= 0 && readIndices1.indexOf ("2") >= 0) {
      printSamRecords (outputWriter, fragmentName, samRecords1, fragmentEntry1, distanceThreshold, outputWeightWriter);
      return;
    }

    String readIndices2 = getReadIndexSet (samRecords2);
    
    if (readIndices1.indexOf ("1") >= 0) {
      if (fragmentEntry2 != null && fragmentEntry1.getNumAlignments() < fragmentEntry2.getNumAlignments() && readIndices2.indexOf ("2") >= 0 && warningsOn) {
	// System.err.println (fragmentName + ": additional (and more than for read 1 from " + samFilename1 + ") SAM records for read 2 from " + samFilename2);
	// System.err.println (fragmentEntry1.toString () + " vs " + fragmentEntry2.toString());
      }

      if (fragmentEntry2 != null && readIndices2.indexOf ("2") >= 0) {
	fragmentEntry1.add(fragmentEntry2);
      }

      printSamRecordsSingleRead (outputWriter, fragmentName, samRecords1, fragmentEntry1, distanceThreshold, 1, outputWeightWriter);
      printSamRecordsSingleRead (outputWriter, fragmentName, samRecords2, null, distanceThreshold, 2, outputWeightWriter);
      return;
    }

    if (readIndices1.indexOf ("2") >= 0) {
      if (fragmentEntry2 != null && fragmentEntry1.getNumAlignments() < fragmentEntry2.getNumAlignments() && readIndices2.indexOf ("1") >= 0 && warningsOn) {
	// System.err.println (fragmentName + ": additional (and more than for read 2 for " + samFilename1 + ") SAM records for read 1 from " + samFilename2);
	// System.err.println (fragmentEntry1.toString () + " vs " + fragmentEntry2.toString());
      }

      if (fragmentEntry2 != null && readIndices2.indexOf ("1") >= 0) {
	fragmentEntry1.add(fragmentEntry2);
      }

      printSamRecordsSingleRead (outputWriter, fragmentName, samRecords2, null, distanceThreshold, 1, outputWeightWriter);
      printSamRecordsSingleRead (outputWriter, fragmentName, samRecords1, fragmentEntry1, distanceThreshold, 2, outputWeightWriter);
      return;
    }

    throw new IOException ("ERROR: unmapped SAM records for " + fragmentName + "\n" +
			   fragmentEntry1.toString () + " vs " + fragmentEntry2.toString());
    
  }
  

  /***********************************************************************************
   * 
   *                           getFragmentEntry
   *
   * FragmentEntries are read from weightReader until a matching one is found. More than
   * one entry per fragmentName is allowed, the entries are added then.
   *
   ***********************************************************************************/

  public static FragmentEntry getFragmentEntry (BufferedReader weightReader, String fragmentName) throws IOException {

    if (fragmentName == null) {
      return null;
    }
      
    String line = weightReader.readLine ();
    if (line == null) {
      return null;
    }
    
    FragmentEntry fragmentEntry   = new FragmentEntry (line);
    String fragmentNameFromReader = fragmentEntry.getFragmentName ();
    while (line != null && fragmentName.compareTo(fragmentNameFromReader) > 0) {

      line = weightReader.readLine ();
      if (line != null) {
	fragmentEntry = new FragmentEntry (line);
	fragmentNameFromReader = fragmentEntry.getFragmentName ();
      }
	  
    }

    if (! fragmentName.equals(fragmentNameFromReader)) {
      throw new IOException ("ERROR: No fragment entry corresponding to fragment " + fragmentName + " found - fragment name found: " + fragmentNameFromReader);
    }

    return fragmentEntry;
    
  }

    
  /***********************************************************************************
   * 
   *                           getReadIndexSet
   *
   ***********************************************************************************/

  public static String getReadIndexSet  (Vector<SamRecord> samRecords) {

    if (samRecords == null) {
      return "";
    }

    String readIndexSet = "";

    for (SamRecord samRecord: samRecords) {

      if (samRecord.isMapped ()) {
	if (samRecord.getReadIndex () == 0) {
	  return "12";
	} else if (samRecord.getReadIndex () == 1 && readIndexSet.indexOf("1") == -1) {
	  readIndexSet = readIndexSet + "1";
	} else if (samRecord.getReadIndex () == 2 && readIndexSet.indexOf("2") == -1) {
	  readIndexSet = readIndexSet + "2";
	}

	if (readIndexSet.length () == 2) {
	  return readIndexSet;
	}
      }
      
    }

    return readIndexSet;
    
  }

  
  /***********************************************************************************
   * 
   *                           isMapped
   *
   ***********************************************************************************/

  public static boolean isMapped (Vector<SamRecord> samRecords) {

    if (samRecords == null) {
      return false;
    }

    for (SamRecord samRecord: samRecords) {

      if (debugLevel >= 2) {
	System.err.println ("Checking: " + samRecord);
      }

      if (samRecord.isMapped ()) {
	if (debugLevel >= 2) {
	  System.err.println ("Returning true");
	}
	return true;	
      }
    }

    if (debugLevel >= 2) {
      System.err.println ("Returning false");
    }
    return false;
    
  }


    
  /***********************************************************************************
   * 
   *                           getCompValue
   *
   ***********************************************************************************/

  public static int getCompValue (Vector<SamRecord> samRecords1, Vector<SamRecord> samRecords2) throws IOException {

    if (samRecords1 == null && samRecords2 == null) {
      return 0;
    }
    
    if (samRecords1 != null && samRecords2 == null) {
      return 1;
    }

    if (samRecords1 == null && samRecords2 != null) {
      return -1;
    }


    if (samRecords1.firstElement() == null) {
      throw new IOException ("Sam records set 1: first element is null.");
    }

    if (samRecords2.firstElement() == null) {
      throw new IOException ("Sam records set 2: first element is null.");
    }

    // System.err.println("Comparing " + samRecords1.firstElement().getFragmentName() + " and " + samRecords2.firstElement().getFragmentName() + " and returning: " +
    //		          samRecords1.firstElement().getFragmentName().compareTo(samRecords2.firstElement().getFragmentName()));

    if (samRecords1.firstElement().getFragmentName() == null) {
      throw new IOException ("Problem with reading first SAM file.");
    }
    
    if (samRecords2.firstElement().getFragmentName() == null) {
      throw new IOException ("Problem with reading second SAM file.");
    }
    
    return samRecords1.firstElement().getFragmentName().compareTo(samRecords2.firstElement().getFragmentName());
    
  }



  /***********************************************************************************/

  private static void printHelp () {
    System.err.println("CombineSamFiles.java\n" +                                              
     "   -- Script to select the best alignment of two SAM files for the same set of reads.\n" +
     "\n" +
     "USAGE: java CombineSamFiles [-s <slack constant>] [-t <distance threshold>]\n" +
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
     "   alignment for the read in the second file plus a slack constant (0.25, see\n" +
     "   option -s). Paired-end alignments are kept if there are only single read\n" +
     "   alignments in the second file (independent of the edit distance). [deflt: -]\n" +
     "<sam file 2>: the second SAM file (- for STDIN). A SAM record is kept if the\n" +
     "   av. edit distance of the alignment is at most the av. edit distance of the\n" +
     "   alignment for the read in the second file minus a slack constant (0.25, see\n" +
     "   option -s). Paired-end alignments are kept if there are only single read\n" +
     "   alignments in the second file (independent of the edit distance). [deflt: -]\n" +
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

    String  readWeightFilename1 = "";
    String  readWeightFilename2 = "";
    String  outputFilename = "-";
    String  outputWeightFilename = "";
    boolean pairedEndOnly = false;
    double  distanceSlack = 0.25;
    String  readIdCutOffString = ":";
    int     distanceThreshold = -1;

    boolean outputNumExpressedReads = false;

    final int countUnit = 5 * 1000 * 1000;

    Getopt g = new Getopt("CombineSamFiles", args, "1:2:ed:o:O:s:t:w:W:h");
    
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
	break;
      case 'e':
	outputNumExpressedReads = true;
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
	System.err.print("Error: getopt() returned " + c + "\n");
      }
      c = g.getopt();
    }

    if (debugLevel >= 2) {
      System.err.println("Options read.");
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
	System.err.println("File name 1: " + samFilename1 + ", file name 2: " + samFilename2);
      }

      PrintWriter outputWriter       = UtilLib.getPrintWriter (outputFilename);
      PrintWriter outputWeightWriter = null;
      if (outputWeightFilename != "") {
	System.err.println ("Writing combined read weights to " + outputWeightFilename);
	outputWeightWriter = UtilLib.getPrintWriter (outputWeightFilename);
      }

      samReader1 = UtilLib.getBufferedReader (samFilename1);
      samReader2 = UtilLib.getBufferedReader (samFilename2);

      BufferedReader weightReader1 = UtilLib.getBufferedReader (readWeightFilename1);
      BufferedReader weightReader2 = UtilLib.getBufferedReader (readWeightFilename2);

      FragmentEntry fragmentEntry1 = null;
      FragmentEntry fragmentEntry2 = null;

      Vector<SamRecord> samRecords1 = null;
      Vector<SamRecord> samRecords2 = null;
      
      samRecords1 = getSamRecords (samReader1);
      samRecords2 = getSamRecords (samReader2);
      
      String fragmentName    = null;
      String oldFragmentName = "";
	
      int numSelectedMappedFrag1 = 0;
      int numSelected1           = 0;
      int numTotal1              = 0;
      int numSelectedMappedFrag2 = 0;
      int numSelected2           = 0;
      int numTotal2              = 0;
      int numLines               = 0;
      
      while (samRecords1 != null && samRecords2 != null) {

	if (debugLevel >= 2) {
	  System.err.println ("SamRecords1: " + samRecords1);
	  System.err.println ("SamRecords2: " + samRecords2);
	}

       String fragmentName1 = null;
       String fragmentName2 = null;
	if (samRecords1 != null && samRecords2 != null) {
	  fragmentName1 = samRecords1.firstElement().getFragmentName();
	  fragmentName2 = samRecords2.firstElement().getFragmentName();
	  if (fragmentName1.compareTo(fragmentName2) < 0) {
	    fragmentName = fragmentName1;
	  } else {
	    fragmentName = fragmentName2;
	  }
	} else if (samRecords1 != null) {
	  fragmentName1 = samRecords1.firstElement().getFragmentName();
	  fragmentName  = fragmentName1;
	} else if (samRecords2 != null) {
	  fragmentName2 = samRecords2.firstElement().getFragmentName();
	  fragmentName  = fragmentName2;
	} else {
	  throw new IOException ("ERROR: samRecords1 and samRecords2 are null.");
	}

	if (oldFragmentName.equals(fragmentName)) {
	  throw new IOException ("ERROR: Duplicate fragment name: old " + oldFragmentName + " vs new " + fragmentName);
	}
	
	if (debugLevel >= 2) {
	  System.err.println("Fragment name: " + fragmentName);
	}
	
	boolean isMapped1 = isMapped (samRecords1) && fragmentName.equals(fragmentName1);
	boolean isMapped2 = isMapped (samRecords2) && fragmentName.equals(fragmentName2);

	fragmentEntry1 = null;
	fragmentEntry2 = null;

	if (! isMapped1 && ! isMapped2) {
	  if (samRecords1 != null && fragmentName.equals(fragmentName1)) {
	    numTotal1++;
	    numSelected1++;
	    // System.err.println("1. fragmentName: " + fragmentName);
	    printSamRecords (outputWriter, fragmentName, samRecords1, fragmentEntry1, distanceThreshold, outputWeightWriter);
	  } else if (samRecords2 != null && fragmentName.equals(fragmentName2)) {
	    numTotal2++;
	    numSelected2++;
	  // System.err.println("1. fragmentName: " + fragmentName);
	    printSamRecords (outputWriter, fragmentName, samRecords2, fragmentEntry2, distanceThreshold, outputWeightWriter);
	  }
	  
	} else if (! isMapped2) {
	  
	  if (! fragmentName.equals (oldFragmentName)) {
	    numSelectedMappedFrag1++;
	  }
	  
	  numTotal1++;
	  numSelected1++;
	  // System.err.println("2. fragmentName: " + fragmentName);
	  fragmentEntry1 = getFragmentEntry (weightReader1, fragmentName);
	  printSamRecords (outputWriter, fragmentName, samRecords1, fragmentEntry1, distanceThreshold, outputWeightWriter);
	 
	  
	} else if (! isMapped1) {
	  
	  if (! fragmentName.equals (oldFragmentName)) {
	    numSelectedMappedFrag2++;
	  }
	  
	  numTotal2++;
	  numSelected2++;
	  // System.err.println("3. fragmentName: " + fragmentName);
	  fragmentEntry2 = getFragmentEntry (weightReader2, fragmentName);
	  printSamRecords (outputWriter, fragmentName, samRecords2, fragmentEntry2, distanceThreshold, outputWeightWriter);
	  
	} else {
	    
	  // System.err.println("4. fragmentName: " + fragmentName);
	  fragmentEntry1 = getFragmentEntry (weightReader1, fragmentName);
	  numTotal1++;
	  fragmentEntry2 = getFragmentEntry (weightReader2, fragmentName);
	  numTotal2++;

	  if (debugLevel >= 1) {
	    System.err.println ("fragment entry 1: " + fragmentEntry1);
	    System.err.println ("fragment entry 2: " + fragmentEntry2);
	  }

	  if (fragmentEntry2 == null) {
	    
	    printSamRecords (outputWriter, fragmentName, samRecords1, fragmentEntry1, distanceThreshold, outputWeightWriter);
	    numSelected1++;	      
	    
	  } else if (fragmentEntry1 == null) {
	    
	    printSamRecords (outputWriter, fragmentName, samRecords2, fragmentEntry2, distanceThreshold, outputWeightWriter);
	    numSelected2++;	      
	    
	  } else {

	    int compMultiplicity = fragmentEntry1.compareMultiplicity (fragmentEntry2);

	    /* Note that compMultiplicity is non-zero if and only if the multiplicity of exactly one of the two fragmentEntries
	       is two. */
	    if (compMultiplicity == 1) {
	      if (! fragmentName.equals (oldFragmentName)) {
		numSelectedMappedFrag1++;
	      }
	      numSelected1++;
	      printSamRecords (outputWriter, fragmentName, samRecords1, fragmentEntry1, distanceThreshold, outputWeightWriter);
	    } else if (compMultiplicity == -1) {
	      if (! fragmentName.equals (oldFragmentName)) {
		numSelectedMappedFrag2++;
	      }
	      numSelected2++;
	      printSamRecords (outputWriter, fragmentName, samRecords2, fragmentEntry2, distanceThreshold, outputWeightWriter);
	    } else {	    

	      if (fragmentEntry1.getEditDistance () < 0) {
		throw new IOException ("ERROR: No edit distance for " + fragmentName + " in " + readWeightFilename1 + " given.");
	      }

	      if (fragmentEntry2.getEditDistance () < 0) {
		throw new IOException ("ERROR: No edit distance for " + fragmentName + " in " + readWeightFilename2 + " given.");
	      }


	      int compDist = fragmentEntry1.compareEditDistance (fragmentEntry2, distanceSlack);

	      if (compDist == 1) {	      
		printSamRecords (outputWriter, fragmentName, samRecords1, samRecords2, fragmentEntry1, fragmentEntry2,
				 distanceThreshold, samFilename1, samFilename2, outputWeightWriter);
		if (! fragmentName.equals (oldFragmentName)) {
		  numSelectedMappedFrag1++;
		}
		numSelected1++;	      
	      } else {
		printSamRecords (outputWriter, fragmentName, samRecords2, samRecords1, fragmentEntry2, fragmentEntry1,
				 distanceThreshold, samFilename2, samFilename1, outputWeightWriter);
		if (! fragmentName.equals (oldFragmentName)) {
		  numSelectedMappedFrag2++;
		}
		numSelected2++;	      
	      }
	    }	    
	  }
	}
		      
	oldFragmentName = fragmentName;

	if (debugLevel >= 2) {
	  System.err.println("Getting samRecords1 after " + fragmentName + " (w/o comparison)");
	}
	samRecords1 = getSamRecords (samReader1, samRecords1);
	if (debugLevel >= 2) {
	  System.err.println("Getting samRecords2 after " + fragmentName + " (w/o comparison)");
	}
	samRecords2 = getSamRecords (samReader2, samRecords2);
	if (debugLevel >= 2) {
	  System.err.println("samRecords2 read: " + samRecords2);
	}

	numLines++;
	if (numLines % countUnit == 0) {
	  System.err.print(".");
	  System.err.flush();
	}

	if (debugLevel >= 2) {
	  System.err.println("Computing comparison value");
	}
	int samRecordCompValue = getCompValue(samRecords1, samRecords2);
	if (debugLevel >= 2) {
	  System.err.println("Done: " + samRecordCompValue);
	}
	while (samRecordCompValue != 0 && samRecords1 != null && samRecords2 != null) {
	  
	  oldFragmentName = fragmentName;
	  
	  if (samRecordCompValue > 0) {
	    
	    numTotal2++;
	    numSelected2++;

	    fragmentName     = samRecords2.firstElement().getFragmentName();
	    boolean isMapped = isMapped (samRecords2);
	    FragmentEntry fragmentEntry = null;
	    if (isMapped (samRecords2)) {
	      if (! fragmentName.equals (oldFragmentName)) {
		numSelectedMappedFrag2++;
	      }
	      fragmentEntry = getFragmentEntry (weightReader2, fragmentName);
	    }	    
	    printSamRecords (outputWriter, fragmentName, samRecords2, fragmentEntry, distanceThreshold, outputWeightWriter);

	    if (debugLevel >= 2) {
	      System.err.println("Getting samRecords2 after " + fragmentName);
	    }
	    samRecords2 = getSamRecords (samReader2, samRecords2);
	    
	  } else {
	    
	    numTotal1++;
	    numSelected1++;

	    fragmentName = samRecords1.firstElement().getFragmentName();
	    FragmentEntry fragmentEntry = null;
	    if (isMapped (samRecords1)) {
	      if (! fragmentName.equals (oldFragmentName)) {
		numSelectedMappedFrag1++;
	      }
	      fragmentEntry = getFragmentEntry (weightReader1, fragmentName);
	    }
	    printSamRecords (outputWriter, fragmentName, samRecords1, fragmentEntry, distanceThreshold, outputWeightWriter);

	    if (debugLevel >= 2) {
	      System.err.println("Getting samRecords1 after " + fragmentName);
	    }
	    samRecords1 = getSamRecords (samReader1, samRecords1);
	  }

	  if (debugLevel >= 2) {
	    System.err.println("Comparing samRecords1 and samRecords2.");
	  }
	  samRecordCompValue = getCompValue(samRecords1, samRecords2);

	}

	numLines++;
	if (numLines % countUnit == 0) {
	  System.err.print(".");
	  System.err.flush();
	}
	
      }

      if (samRecords1 != null || samRecords2 != null) {
	Vector<SamRecord> samRecords   = samRecords1;
	BufferedReader    samReader    = samReader1;
	BufferedReader    weightReader = weightReader1;
	if (samRecords2 != null) {
	  samRecords   = samRecords1;
	  samReader    = samReader2;
	  weightReader = weightReader2;
	}
	
	while (samRecords != null) {
	  oldFragmentName = fragmentName;
	  	      
	  fragmentName = samRecords.firstElement().getFragmentName();
	  FragmentEntry fragmentEntry = null;
	  if (isMapped (samRecords)) {
	    if (! fragmentName.equals (oldFragmentName) && samRecords1 != null) {
	      numSelectedMappedFrag1++;
	    }
	    fragmentEntry = getFragmentEntry (weightReader, fragmentName);
	  }
	  printSamRecords (outputWriter, fragmentName, samRecords, fragmentEntry, distanceThreshold, outputWeightWriter);
	  samRecords = getSamRecords (samReader, samRecords);

	  if (samRecords1 != null) {
	    numTotal1++;
	    numSelected1++;
	  } else {
	    numTotal2++;
	    numSelected2++;
	  }
	}
      }


      if (numLines >= countUnit) {
	System.err.println();
      }

      System.err.println(numSelected1 + " records of " + numTotal1 + " records selected from file " + samFilename1);
      System.err.println(numSelected2 + " records of " + numTotal2 + " records selected from file " + samFilename2);
      
      if (outputWeightWriter != null) {
	System.err.println("Total number of fragments output: " + numOutput);
      }

      if (outputNumExpressedReads) {
	System.err.println("NUMBER_EXPRESSED_READS=" + numSelectedMappedFrag1);
      }


      samReader1.close();
      samReader2.close();

      weightReader1.close();
      weightReader2.close();
      
      outputWriter.close ();

      if (outputWeightWriter != null) {
	outputWeightWriter.close();
      }

    }
    catch (Exception e) {
      System.err.println ("Problem in line: " + line + ": " + e==null?"No error message":e.getMessage());
      System.exit (1);
    }
  }
}

