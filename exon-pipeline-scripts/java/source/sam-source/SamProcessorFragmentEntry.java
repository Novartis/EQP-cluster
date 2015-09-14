/**File: SamProcessorEntryTable.java 

Original Author: Sven Schuierer
Date: 03/01/2012

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


/***********************************************************************************
 *
 *                              Class SamProcessorFragmentEntry
 *
 ***********************************************************************************/

class SamProcessorFragmentEntry implements SamProcessor {

  private static int     debugLevel = UtilLib.getDebugLevel (); 
  private static boolean warningsOn = false;

  public static void setDebugLevel (int value) {
    debugLevel = value;
  }

  public static void setWarningsOn (boolean value) {
    warningsOn = value;
  }

  private static HashSet<String> referenceSequenceIdSet = null;
  public static void setReferenceSequenceIdSet (HashSet<String> s) {
    referenceSequenceIdSet = s;
  }

  private static boolean readsMode = false;
  public  static void setReadsMode () {
    readsMode = true;
  }
  
  /***********************************************************************************
   *
   *                      Object variables and methods   
   *
   ***********************************************************************************/

  private FragmentEntry fragmentEntry = null;
  private String        fragmentName  = "";
  private PrintWriter   outputWriter  = null;

  SamProcessorFragmentEntry (PrintWriter outputWriter) {
    debugLevel = UtilLib.getDebugLevel ();
    this.outputWriter = outputWriter;
  }

  public void init (SamRecord samRecord) {
  }

  /***********************************************************************************
   *
   *                        outputFragmentEntry
   *
   ***********************************************************************************/

  public void outputFragmentEntry () throws IOException {

    if (fragmentEntry != null) {
      outputWriter.println(fragmentEntry.toPrintString ());
    }

  }

  
  /***********************************************************************************
   *
   *                           processSamRecords
   *      The reads are processed by fragmentName and only if the fragmentName
   *      changes a fragmentEntry is output. For each paired-end read a mate
   *      is identified. Either the SamRecords of read and mate (for paired-end
   *      alignments) or of the read alone are then added to the fragmentEntry.
   *
   *      Note that fragmentName is an object variable and it is in principle possible
   *      to process several vectors which contain the alignments of one fragment
   *      if they are processed consecutively.
   *
   ***********************************************************************************/

  public void processSamRecords (Vector<SamRecord> samRecords) throws Exception {

    if (debugLevel >= 2) {
      System.err.println ("Number of SAM records: " + samRecords.size ());
      System.err.println("Processing sam records: " + samRecords);
    }

    if (samRecords.size() == 0) {
      return;
    }

    if (debugLevel >= 1) {
      warningsOn = true;
    }

    if (readsMode) {
      /* Note that samRecords.compareTo(SamRecord s) is based on the queryName */
      Collections.sort(samRecords);
    }


    int sumEditDistance = 0;
    int numFragments = 0;

    HashSet<Integer>  processedIndices   = new HashSet<Integer>  (2 * samRecords.size());
    Vector<SamRecord> unprocessedRecords = new Vector<SamRecord> (samRecords.size());

    SamRecord samRecord = null;
    String curFragmentName = "";
    for (int i = 0; i < samRecords.size(); i++) {

      Integer iInt = new Integer(i);
      if (! processedIndices.contains(iInt)) {

	samRecord = samRecords.get(i);	
	if (! readsMode) {
	  curFragmentName = samRecord.getFragmentName ();
	} else {
	  curFragmentName = samRecord.getReadId ();
	}

	if (debugLevel >= 2) {
	  System.err.println("debugLevel: " + debugLevel + ", curFragmentName: " + curFragmentName + ", readsMode: " + readsMode);
	}

	if (! fragmentName.equals(curFragmentName)) {
	  /* We only allow the old fragmentName to differ from the current fragment name if the old fragment name is from the previous
	     call to processSamRecords or if we are in readsMode. */
	  if (i == 0 || readsMode) {
	    if (fragmentName != "") {
	      /* Output current fragmentEntry: Note that we need to be able to output the last fragment from ComputeReadWeigths
		 and, hence, we use a parameter less object method to print the current fragmentEntry */
	      outputFragmentEntry ();
	    }

	    if (i > 0 && readsMode && warningsOn) {
	      System.err.println("WARNING: Two different fragments (" + fragmentName + " and " + curFragmentName +
				 " in SAM record collection:\n" + samRecords.toString());
	    }

	    fragmentName  = curFragmentName;
	    fragmentEntry = new FragmentEntry ();
	    
	  } else {
	    throw new Exception ("ERROR: Two different fragments (" + fragmentName + " and " + curFragmentName +
				 " in SAM record collection:\n" + samRecords + "\n");
	  }
	}

	if (fragmentEntry == null) {
	  throw new Exception ("ERROR: fragmentEntry is null.");
	}

	boolean   mateFound = false;
	SamRecord mateSamRecord = null;
	if (samRecord.hasMate() && ! readsMode) {

	  mateSamRecord = samRecord.findMate (samRecords, i, processedIndices);

	  if (debugLevel >= 2) {
	    System.err.println ("SamRecord: " + samRecord);
	    System.err.println ("Mate SamRecord: " + mateSamRecord);
	    System.err.println ();
	  }
 
	  if (mateSamRecord != null) {

	    processedIndices.add(iInt);
	    if (referenceSequenceIdSet == null || referenceSequenceIdSet.contains(samRecord.getReferenceName ())) {
	      fragmentEntry.add (samRecord, mateSamRecord);
	    }

	  } else {
	    unprocessedRecords.add(samRecord);
	    if (warningsOn) {
	      throw new Exception ("WARNING: " + samRecord + " is paired end and no mate is found.");
	    }
	  }
	} else {

	  if (debugLevel >= 2) {
	    System.err.println ("Single end SamRecord: " + samRecord);
	    System.err.println ();
	  }
	  
	  processedIndices.add(iInt);
	  if (referenceSequenceIdSet == null || referenceSequenceIdSet.contains(samRecord.getReferenceName ())) {
	    fragmentEntry.add (samRecord, readsMode);
	  } 
	}
      }
    }

    if (debugLevel >= 2) {
      System.err.println ("Number of paired end alignments for fragment " + fragmentEntry.getFragmentName () + ": " + fragmentEntry.getNumAlignmentsPe ());
      System.err.println ("Number of single read alignments for fragment " + fragmentEntry.getFragmentName () + ": " + fragmentEntry.getNumAlignmentsSr ());
    }

  }

}
