/**File: SamProcessorSpliced.java 

Original Author: Sven Schuierer
Date: 25/02/2015

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
 *                              Class SamProcessorSpliced
 *
 ***********************************************************************************/

class SamProcessorSpliced implements SamProcessor {

  private static int     debugLevel = 0; 
  private static boolean warningsOn = false;

  public static void setDebugLevel (int value) {
    debugLevel = value;
  }

  public static void setWarningsOn (boolean value) {
    warningsOn = value;
  }
  
  /***********************************************************************************
   *
   *                      Object variables and methods   
   *
   ***********************************************************************************/

  private PrintWriter   outputWriter  = null;

  SamProcessorSpliced (PrintWriter outputWriter) {
    debugLevel = UtilLib.getDebugLevel ();
    this.outputWriter = outputWriter;
  }

  public void init (SamRecord samRecord) {
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


    int sumEditDistance = 0;
    int numFragments = 0;

    HashSet<Integer>  processedIndices   = new HashSet<Integer>  (2 * samRecords.size());
    Vector<SamRecord> unprocessedRecords = new Vector<SamRecord> (samRecords.size());

    SamRecord samRecord    = null;
    for (int i = 0; i < samRecords.size(); i++) {

      Integer iInt = new Integer(i);
      if (! processedIndices.contains(iInt)) {

	samRecord = samRecords.get(i);
	boolean   fragmentIsSpliced = samRecord.isSpliced ();
	
	boolean   mateFound = false;
	SamRecord mateSamRecord = null;
	if (samRecord.hasMate()) {

	  mateSamRecord = samRecord.findMate (samRecords, i, processedIndices);
	  if (mateSamRecord != null) {
	    fragmentIsSpliced = fragmentIsSpliced | mateSamRecord.isSpliced();
	  }
	}

	if (fragmentIsSpliced) {
	  outputWriter.println(samRecord.toString ());
	  if (mateSamRecord != null) {
	    outputWriter.println(mateSamRecord.toString ());
	  }
	}
	  
	processedIndices.add(iInt);
      }
    }
  }
}
