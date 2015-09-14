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
 *                              Class SamProcessorBed
 *
 ***********************************************************************************/

public class SamProcessorBed implements SamProcessor {

  private static final int debugLevel = 0;
  

  /***********************************************************************************
   *
   *                     Object variables and methods
   *
   ***********************************************************************************/

  private HashSetTable<String, String> transcriptVersionIds = null;

  private String  readId = "";
  private String  fragmentId = "";
  private String  oldFragmentId = "";
  private Counter alignmentCounter;
  private Counter fragmentCounter;
  
  private boolean     warningsOn     = false;
  private boolean     newIdentifiers = false;
  private PrintWriter mappingWriter  = null;
  private PrintWriter bedWriter      = null;

  private int alignmentScoreThreshold = Integer.MIN_VALUE;
  private int editDistanceThreshold   = Integer.MAX_VALUE;

  private String strandSpecificDirection = "none";

  private boolean primaryAlignmentsOnly;

  public SamProcessorBed () {
    fragmentCounter  = new Counter(9);
    alignmentCounter = new Counter (5);
  }


  public SamProcessorBed (PrintWriter bedWriter, HashSetTable<String, String> transcriptVersionIds) {

    this ();
    this.newIdentifiers = false;

  }

  public SamProcessorBed (PrintWriter bedWriter, PrintWriter mappingWriter, HashSetTable<String, String> transcriptVersionIds) {

    this ();
    this.mappingWriter        = mappingWriter;
    this.transcriptVersionIds = transcriptVersionIds;
    this.newIdentifiers       = false;

  }


  public SamProcessorBed (PrintWriter bedWriter, boolean newIdentifiers, HashSetTable<String, String> transcriptVersionIds) {

    this ();
    this.newIdentifiers       = newIdentifiers;
    this.transcriptVersionIds = transcriptVersionIds;

  }


  public SamProcessorBed (PrintWriter bedWriter, boolean newIdentifiers, PrintWriter mappingWriter,
			  HashSetTable<String, String> transcriptVersionIds, int alignmentScoreThreshold,
			  int editDistanceThreshold) {

    this ();
    
    this.bedWriter            = bedWriter;
    this.newIdentifiers       = newIdentifiers;
    this.mappingWriter        = mappingWriter;
    this.transcriptVersionIds = transcriptVersionIds;

    this.alignmentScoreThreshold = alignmentScoreThreshold;
    this.editDistanceThreshold   = editDistanceThreshold;
    
  }

  public SamProcessorBed (PrintWriter bedWriter, boolean newIdentifiers, PrintWriter mappingWriter,
			  HashSetTable<String, String> transcriptVersionIds, int alignmentScoreThreshold,
			  int editDistanceThreshold, String strandSpecificDirection, boolean primaryAlignmentsOnly) {

    this (bedWriter, newIdentifiers, mappingWriter, transcriptVersionIds, alignmentScoreThreshold, editDistanceThreshold);
    
    this.strandSpecificDirection = strandSpecificDirection;
    this.primaryAlignmentsOnly   = primaryAlignmentsOnly;
    
  }


  /***********************************************************************************
   *
   *                     init
   *
   * Note that the fragment id defined here is used in processSamRecords to invoke
   * samRecord.printBedRecords where it is the basis of the identifiers of the bed
   * entries for the sam record.
   *
   ***********************************************************************************/
  
  public void init (SamRecord samRecord) throws IOException {

    readId     = samRecord.getQueryName();
    fragmentId = samRecord.getFragmentName();

    if (! oldFragmentId.equals(fragmentId)) {
      alignmentCounter = new Counter (5);
    }

    oldFragmentId = fragmentId;
    
  }


  /***********************************************************************************
   *
   *                           processSamRecords   
   *
   ***********************************************************************************/

  public void processSamRecords (Vector<SamRecord> samRecords) throws IOException {

    if (transcriptVersionIds == null && warningsOn) {
      System.err.println ("WARNING: transcriptVersionIds not set.");
    }

    HashSet<Integer> processedIndices = new HashSet<Integer> (2 * samRecords.size());

    if (debugLevel >= 2) {
      System.out.println("Processing " + samRecords.size() + " sam records.");
    }

    int readIndex = 0;
    for (int i = 0; i < samRecords.size(); i++) {

      Integer iInt = Integer.valueOf (i);
      if (! processedIndices.contains(iInt)) {

	processedIndices.add(iInt);

	SamRecord samRecord = samRecords.get(i);
	if (! primaryAlignmentsOnly || samRecord.isPrimary ()) {
	  
	  SamRecord mateSamRecord = null;
	  if (samRecord.hasMate ()) {	    
	    mateSamRecord = samRecord.findMate(samRecords, i, processedIndices);
	    if (primaryAlignmentsOnly && ! mateSamRecord.isPrimary ()) {
	      mateSamRecord = null;
	    }
	    
	  }

	  if (debugLevel >= 2) {
	    System.out.println("Printing BED entry for " + samRecord + " and " + mateSamRecord);
	  }
	  
	  samRecord.printBedRecords (mateSamRecord, fragmentId, alignmentCounter, transcriptVersionIds, alignmentScoreThreshold,
				     editDistanceThreshold, strandSpecificDirection, bedWriter);
	  alignmentCounter.inc();
	}

      }

    }

    if (debugLevel >= 2) {
      System.out.println("Sam records processed.");
    }

    alignmentCounter.dec();
    samRecords.clear();

  }


}
