/**File WeightObject.java 

Original Author: Sven Schuierer
Date: 20/12/2011

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
 *                              Class WeightObject
 *
 *  The weightObject class represents the aligned fragments or reads.
 *
 *  One weight object has one or more weight object alignments associated to it. The
 *  alignments are given by the bed records for the read or fragment. A new weight object
 *  alignment is created for each new alignment id (-ANNNN) of the bed records. It
 *  is assumed that the bed records for one alignment ids are consecutive.
 *
 *  Once all the bed records are collected the genomic alignment strings can be computed.
 *  The weight object alignments with the same genomic alignment string are unified.
 *  In order to do this, first ambiguous genomic mappings (which may result from length
 *  differences of genomic and transcript exon coordinates) are adjusted to match
 *  other genomic mappings if possible. Then the exons associated to weight object 
 *  alignments with the same genomic alignment strings are combined as it may happen
 *  that a read or fragment maps to the same genomic location via a different exon
 *  from a different transcript or even gene and the association to this exon or gene
 *  is lost if the weight object alignment is discarded.
 *
 ***********************************************************************************/

class WeightObject {

  private static int debugLevel = UtilLib.getDebugLevel ();

  public static void setDebugLevel (int value) {
    debugLevel = value;
  }

  private String readAlignmentId = "";
  private String readId = "";
  private String fragmentId = "";
  private String curAlignmentId = "";
  private String curAlignmentBaseId = "";
  private double weight = 1;

  private Vector<WeightObjectAlignment>  weightObjectAlignments   = new Vector<WeightObjectAlignment>  (20);
  private HashSet<WeightObjectAlignment> weightObjectAlignmentSet = new HashSet<WeightObjectAlignment> (20);
  private WeightObjectAlignment lastWeightObjectAlignment = null;

  
  /***********************************************************************************
   *
   *                         Constructors
   *
   ***********************************************************************************/


  public WeightObject () {

    debugLevel = UtilLib.getDebugLevel ();
    
  }

  public WeightObject (BedRecord bedRecord) throws Exception {

    this ();
          
    readAlignmentId = bedRecord.getReadAlignmentId ();
    fragmentId      = bedRecord.getFragmentId ();
    
    int    readIndex = bedRecord.getReadIndex ();
    String readAlignmentIdEnd = readAlignmentId.substring(readAlignmentId.length()-2, readAlignmentId.length());
    
    if (readAlignmentIdEnd.equals("S1") || readAlignmentIdEnd.equals("S2")) {
      readId = fragmentId + "/S" + readIndex;
    } else {
      readId = fragmentId + "/P" + readIndex;
    }

    addWeightObjectAlignment (bedRecord);
    
  }

  
  /***********************************************************************************
   *
   *                         addWeightObjectAlignment
   *
   ***********************************************************************************/


  public void addWeightObjectAlignment (BedRecord bedRecord) throws Exception {

    if (debugLevel >= 2) {
      System.out.println("Creating new weightObjectAlignment for" + bedRecord.getReadAlignmentId () +
			 ", number of weight object alignments: " + weightObjectAlignments.size());
    }
 
    /* alignmentBaseId is the read id (<fragment id>/S[12]) for single-read alignments and the fragment id
       for paired-end alignments. It is part of the genomic alignment as we want to distinguish between the
       genomic alignments of the first and second read of a fragment even if they have the same genomic
       cooridinate. */

    curAlignmentId     = bedRecord.getAlignmentId ();
    curAlignmentBaseId = bedRecord.getAlignmentBaseId ();
    
    lastWeightObjectAlignment = new WeightObjectAlignment (curAlignmentId, curAlignmentBaseId, fragmentId);

    weightObjectAlignments.add(lastWeightObjectAlignment);
    lastWeightObjectAlignment.addBedRecord (bedRecord);
    if (debugLevel >= 1) {
      System.out.println("Adding bed record " + bedRecord + " to " + lastWeightObjectAlignment.getAlignmentId ());
      lastWeightObjectAlignment.printExonSet();
    }
  }


  /***********************************************************************************
   *
   *                         addBedRecord
   *
   * Checks whether the bedRecord belongs to the same or a new alignment and creates
   * a new alignment if necessary. Sometimes it may happen that one weight object (i.e.
   * fragment) has single-read alignments for the first and second read. In this case
   * we need to create a new alignment even if the alignment id remains the same
   * (e.g. A0001/S1 vs A0001/S2).
   *
   ***********************************************************************************/
  
  public void addBedRecord (BedRecord bedRecord) throws Exception {

    if (debugLevel >= 2) {
      System.out.println("weightObject curAlignmentId: " + curAlignmentId +              " and curAlignmentBaseId: " + curAlignmentBaseId);
      System.out.println("bedRecord    alignmentId:    " + bedRecord.getAlignmentId () + " and alignmentBaseId:    " + bedRecord.getAlignmentBaseId ());      
    }
    
    if (bedRecord.getAlignmentId ().equals(curAlignmentId) && bedRecord.getAlignmentBaseId ().equals(curAlignmentBaseId)) {
      if (debugLevel >= 1) {
	System.out.println("Adding bed record " + bedRecord + " to " + lastWeightObjectAlignment.getAlignmentId ());
	lastWeightObjectAlignment.printExonSet();
      }      
      lastWeightObjectAlignment.addBedRecord (bedRecord);
    } else {
      if (debugLevel >= 1) {
	System.out.println("bedRecords for " + lastWeightObjectAlignment.getAlignmentId ());
	System.out.println("bedRecords1: " + lastWeightObjectAlignment.getBedRecords1 ());
	System.out.println("bedRecords2: " + lastWeightObjectAlignment.getBedRecords2 ());
      }
      addWeightObjectAlignment (bedRecord);
    }
  }


  /***********************************************************************************
   *
   *                         adjustGenomicAlignments
   *
   * Adjust the genomic coordinates of ambiguous genomic alignments
   *
   ***********************************************************************************/
  
  public void adjustGenomicAlignments () throws IOException {

    if (debugLevel >= 2) {
      System.out.println("Adjust genomic aligments");
    }

    for (WeightObjectAlignment weightObjectAlignment: weightObjectAlignments) {
      weightObjectAlignment.initializeGenomicAlignmentString ();
      if (debugLevel >= 2) {
	System.out.println("Genomic alignment string before: " + weightObjectAlignment.getAlignmentId () + ": " + weightObjectAlignment.getGenomicAlignmentString ());
      }
    }

    for (WeightObjectAlignment weightObjectAlignment: weightObjectAlignments) {
      /* We need to recompute the genomic alignment strings since although we try to reflect
	 the changes of adjustGenomicAlignmentString in the genomicAlignmentString we sometimes
	 need to change single genomic intervals of weightAlignments whose genomicAlignmentStrings
	 were already adjusted. */
      weightObjectAlignment.adjustGenomicAlignmentString (weightObjectAlignments);
      if (debugLevel >= 2) {
	System.out.println("Genomic alignment string after " + weightObjectAlignment.getAlignmentId () + ": " + weightObjectAlignment.getGenomicAlignmentString ());
      }

    }
  }

  

  /***********************************************************************************
   *
   *                         mergeWeightAlignments
   *
   * Combine the exons of weight object alignments with the same genomic alignment string
   *
   ***********************************************************************************/
  
  public void mergeWeightAlignments () throws IOException {

    Hashtable<String, WeightObjectAlignment>  weightObjectAlignmentTable = new Hashtable<String, WeightObjectAlignment> ();

    /* Note that getGenomicAlignmentString is called as via equals when we add weightObjectAlignment to
       the hash set weigthObjectAlignmentSet ! */
    for (WeightObjectAlignment weightObjectAlignment: weightObjectAlignments) {

      String alignmentString = weightObjectAlignment.getGenomicAlignmentString ();

      if (debugLevel >= 2) {
	System.out.println("Processing weightObjectAlignment: " + weightObjectAlignment);
      }

      if (! alignmentString.equals ("")) {
		
	WeightObjectAlignment weightObjectAlignmentFromTable = weightObjectAlignmentTable.get(alignmentString);
	if (weightObjectAlignmentFromTable == null) {
	  weightObjectAlignmentTable.put (alignmentString, weightObjectAlignment);	  
	  weightObjectAlignmentSet.add (weightObjectAlignment);
	} else {
	  if (debugLevel >= 2) {
	    System.out.println("Combining exons of " + weightObjectAlignment + " with " + weightObjectAlignmentFromTable);
	  }
	  
	  weightObjectAlignmentFromTable.combineExons(weightObjectAlignment);
	}
      }
    }

    weightObjectAlignments.clear ();

    if (debugLevel >= 1) {
      System.out.println("weightObjectAlignmentSet: " + weightObjectAlignmentSet);
    }
    
  }


  /***********************************************************************************
   * 
   *                          printGenomicBedEntries
   *
   ***********************************************************************************/

  public void printGenomicBedEntries (PrintWriter outputWriter) throws IOException {

    /* Compute genomic alignment strings and adjust them if necessary */
    adjustGenomicAlignments ();
    
    /* Combine the exons of weight object alignments with the same genomic alignment string */
    mergeWeightAlignments ();

    Counter counter = new Counter (5);
    for (WeightObjectAlignment weightObjectAlignment: getWeightObjectAlignmentSet()) {

      if (! weightObjectAlignment.getGenomicAlignmentString ().equals(weightObjectAlignment.getAlignmentBaseId() + ":")) {

	GenomicAlignment genomicAlignment1 = weightObjectAlignment.getGenomicAlignment1();
	GenomicAlignment genomicAlignment2 = weightObjectAlignment.getGenomicAlignment2();

	boolean pairedEnd = genomicAlignment1 != null && genomicAlignment2 != null;

	String bedEntryId = fragmentId + "-A" + counter.toString() + "-C/" + (pairedEnd?"P":"S");

	if (genomicAlignment1 != null) {
	  Vector<String> genomicBedEntries1 = genomicAlignment1.getGenomicBedEntryString (bedEntryId + "1");
	  for (String genomicBedEntryString1: genomicBedEntries1) {
	    outputWriter.println(genomicBedEntryString1);
	  }
	}
      
	if (genomicAlignment2 != null) {
	  Vector<String> genomicBedEntries2 = genomicAlignment2.getGenomicBedEntryString (bedEntryId + "2");
	  for (String genomicBedEntryString2: genomicBedEntries2) {
	    outputWriter.println(genomicBedEntryString2);
	  }
	}

	counter.inc();

      }
    }
  }

  /***********************************************************************************
   * 
   * getExonOverlapOnTranscript
   *
   ***********************************************************************************/

  public int getExonOverlapOnTranscript (GenomicInterval genomicInterval, GenomicInterval oldGenomicInterval) throws IOException {

    if (oldGenomicInterval == null || genomicInterval == null) {
      return 0;
    }

    if (genomicInterval.getChrEnd () - genomicInterval.getChrStart () + 1 != genomicInterval.getBedRecord().getOverlap () &&
	! genomicInterval.getBedRecord().getReferenceSequenceId ().equals(genomicInterval.getBedRecord().getExon().getChromosome())) {
      throw new IOException ("Genomic interval where BED overlap " + genomicInterval.getBedRecord().getOverlap () + " is different from length: " +
			     genomicInterval.getChrEnd () + " - " + genomicInterval.getChrStart () + " + 1");
    }

    int exonTranscriptOverlap = 0;
    if (oldGenomicInterval != null) {
      if (genomicInterval.getBedRecord().getExonReferenceStart () >= oldGenomicInterval.getBedRecord().getExonReferenceStart ()) {
	exonTranscriptOverlap = oldGenomicInterval.getBedRecord().getExonReferenceEnd () - genomicInterval.getBedRecord().getExonReferenceStart () + 1;
      } else {
	exonTranscriptOverlap = genomicInterval.getBedRecord().getExonReferenceEnd () - oldGenomicInterval.getBedRecord().getExonReferenceStart () + 1;
      }
    }

    int minBedOverlap = Math.min (oldGenomicInterval.getBedRecord(). getOverlap (), genomicInterval.getBedRecord().getOverlap ());

    if (debugLevel >= 2) {
      System.out.println("genomicIntervalLength = (" + genomicInterval.getBedRecord().getOverlap () + " - Math.min(" +
			 Math.max(0, exonTranscriptOverlap) + ", " + minBedOverlap + ")");
      System.out.flush();
    }
    
    return Math.min (minBedOverlap, Math.max(0, exonTranscriptOverlap));
    
  }
  
  /***********************************************************************************
   * 
   * getGenomicIntervalLength
   *
   ***********************************************************************************/

  public int getGenomicIntervalLength (GenomicInterval genomicInterval, GenomicInterval oldGenomicInterval) throws IOException {

    int genomicIntervalLength = 0;
    if (oldGenomicInterval == null) {
      genomicIntervalLength = genomicInterval.getChrEnd () - genomicInterval.getChrStart () + 1;
      if (debugLevel >= 2) {
	System.out.println("Computing genomicIntervalLength from genomic interval: " + genomicIntervalLength);
	System.out.flush();
      }
    } else {
      /* The length of the genomic interval which we want to subtract from the CIAGR operations is the length of the
	 BED overlap minus the length of the overlap of the consecutive exons on the transcript (here marked with "="):

	                    
	    oldGenomicInterval	 genomicInterval
	         |------------|----------------------|-----...  read
		          |====----------------------|  exon 2
	  |-------------------|                         exon 1

       Note that the read may end before the end of exon 1 (or start after the beginning of exon 2) in which case the
       additional exon overlap is limited by the minimum of the BED overlap of the read with exon 1 and exon 2. */

      genomicIntervalLength =
	genomicInterval.getBedRecord().getOverlap () - getExonOverlapOnTranscript (genomicInterval, oldGenomicInterval);


    }

    return genomicIntervalLength;

  }
  
  
  /***********************************************************************************
   * 
   * cleanGenomicIntervals: Remove overlaps that are covered by overlaps with
   * other exons and transcriptExons
   *
   ***********************************************************************************/

  public Vector<GenomicInterval> cleanGenomicIntervals (Vector<GenomicInterval> genomicIntervals) throws IOException {

    if (debugLevel >= 2) {
      System.out.println ("Cleaning interval list: " + genomicIntervals);
    }
    

    Vector<GenomicInterval> cleansedGenomicIntervals = new Vector<GenomicInterval> (genomicIntervals.size ());

    GenomicInterval oldGenomicInterval = null;
    int oldGenomicIntervalLength = -1;
    for (GenomicInterval genomicInterval : genomicIntervals) {
      int genomicIntervalLength   = genomicInterval.getChrEnd () - genomicInterval.getChrStart () + 1;
      int exonOverlapOnTranscript = getExonOverlapOnTranscript (genomicInterval, oldGenomicInterval);
      
      if (debugLevel >= 2) {
	System.out.println ("Genomic interval: " + genomicInterval + ", length: " + genomicIntervalLength + ", old length: " + oldGenomicIntervalLength +
			    ", overlap: " + exonOverlapOnTranscript);
      }
       
      if (exonOverlapOnTranscript >= oldGenomicIntervalLength) {
	/* Do not add old genomic interval */
	oldGenomicInterval       = genomicInterval;
	oldGenomicIntervalLength = genomicIntervalLength;
      } else if (exonOverlapOnTranscript >= genomicIntervalLength) {
	/* Keep oldGenomicInterval and skip genomicInterval */
      } else {
	if (debugLevel >= 2) {
	  System.out.println ("Adding genomic interval: " + oldGenomicInterval);
	}
	if (! UtilLib.isTranscriptExon(oldGenomicInterval.getChromosome ()) || true) {
	  cleansedGenomicIntervals.add (oldGenomicInterval);
	}
	oldGenomicInterval       = genomicInterval;
	oldGenomicIntervalLength = genomicIntervalLength;
      }
    }

    if (! UtilLib.isTranscriptExon(oldGenomicInterval.getChromosome ()) || true) {
      cleansedGenomicIntervals.add (oldGenomicInterval);
    }

    return cleansedGenomicIntervals;
  }


  /***********************************************************************************
   * 
   * getGenomicSamRecord
   *
   ***********************************************************************************/

  public SamRecord getGenomicSamRecord (SamRecord samRecord, Vector<GenomicInterval> originalGenomicIntervals, Vector<GenomicInterval> mateGenomicIntervals,
					String referenceOrientation, HashSet<String> chromosomeIds, boolean isPrimaryAlignment) throws IOException {

    Vector<GenomicInterval> genomicIntervals = cleanGenomicIntervals (originalGenomicIntervals);

    if (debugLevel >= 2 ) {
      System.out.println("\n\ngetGenomicSamRcord");
      System.out.println("samRecord: " + samRecord);
      System.out.println("originalGenomicIntervals: " + originalGenomicIntervals);
      System.out.println("genomicIntervals: " + genomicIntervals);
      System.out.println("mateGenomicIntervals: " + mateGenomicIntervals);
      System.out.flush();
    }

    int i = 0;
    GenomicInterval genomicInterval = genomicIntervals.get(i++);
    String samReference = genomicInterval.getChromosome ();
    int    samPosition  = genomicInterval.getChrStart ();
    String samStrand    = genomicInterval.getStrand ();
    while ((genomicInterval == null || UtilLib.isTranscriptExon(samReference)) && i < genomicIntervals.size ()) {
      genomicInterval = genomicIntervals.get(i++);
      samPosition  = genomicInterval.getChrStart ();
      samReference = genomicInterval.getChromosome ();
      samStrand    = genomicInterval.getStrand ();
    } 
    
    if (debugLevel >= 2 || false) {
      System.out.println("samRecord: " + samRecord);
      System.out.flush();
    }


    String mateSamReference = "*";
    int    mateSamPosition  = 0;
    String mateSamStrand    = "+";
    String insertSize       = "0";
    if (mateGenomicIntervals != null) {
      int j = 0;
      GenomicInterval mateGenomicInterval = null;
      while ((mateGenomicInterval == null || UtilLib.isTranscriptExon(mateSamReference)) && j < mateGenomicIntervals.size ()) {
	mateGenomicInterval = mateGenomicIntervals.get(j++);
	mateSamPosition  = mateGenomicInterval.getChrStart ();
	mateSamReference = mateGenomicInterval.getChromosome ();
	mateSamStrand    = mateGenomicInterval.getStrand ();
	insertSize       = samRecord.getInsertSize ();
      }
    }

    CiagrString ciagrString = new CiagrString (samRecord.getCiagrString (), samRecord.getSequence(), samRecord.getQualityString());
    if (referenceOrientation.equals("-")) {
      if (debugLevel >= 2 || false) {
	System.out.println("Inverting ciagrString: " + ciagrString);
      }
      ciagrString.invert ();
    }
    
    CiagrString newCiagrString = new CiagrString ();

    int genomicIntervalIndex = 0;
    int leftReadPosition     = -1;
    int rightReadPosition    = 0;
    int curExonEnd           = -1;

    GenomicInterval oldGenomicInterval  = null;
    boolean inverted = false;
    while (genomicIntervalIndex < genomicIntervals.size ()) {

      genomicInterval = genomicIntervals.get (genomicIntervalIndex++);

      if (referenceOrientation.equals("*")) {
	String bedReferenceOrientation = genomicInterval.getBedRecord().getReferenceGenomicAlignmentOrientation ();
	if ((bedReferenceOrientation.equals("-") && ! inverted) || (bedReferenceOrientation.equals("+") && inverted)) {
	  System.out.println("Inverting ciagrString: " + ciagrString);
	  ciagrString.invert ();
	  inverted = ! inverted;
	}
      }

      int genomicIntervalLength = genomicInterval.getChrEnd () - genomicInterval.getChrStart () + 1;
      if (debugLevel >= 2 || false) {
	System.out.println ("Genomic interval: " + genomicInterval);
	System.out.println ("leftReadPosition: " + leftReadPosition + ", rightReadPosition: " +  rightReadPosition + ", genomic interval length: " +
			    genomicIntervalLength  + ", exon transcript overlap: " + getExonOverlapOnTranscript (genomicInterval, oldGenomicInterval));
	System.out.println ("Genomic interval start: " + genomicInterval.getChrStart () + ", current exon end: " + curExonEnd);
	System.out.flush();
      }
      
      int exonOverlapOnTranscript = getExonOverlapOnTranscript (genomicInterval, oldGenomicInterval);
      leftReadPosition  = rightReadPosition + 1 - exonOverlapOnTranscript;
      rightReadPosition = leftReadPosition + genomicIntervalLength - 1;

      if (debugLevel >= 2 || false) {
	System.out.println ("New leftReadPosition: " + leftReadPosition + ", new rightReadPosition: " +  rightReadPosition);
	System.out.flush();
      }

      if (genomicInterval.getChromosome ().equals(samReference) && curExonEnd >= 0) {
	/* Note that the intron length is just the length between the end of the last exon and the beginning
	   of the current exon - independent of the exon overlap on the transcript. */
	int intronLength = genomicInterval.getChrStart () - curExonEnd - 1;
	newCiagrString.addField (intronLength, "N");
      }
      boolean isAtAlignmentMargin = genomicIntervalIndex == 1 || genomicIntervalIndex == genomicIntervals.size ();
      ciagrString.transferInterval (leftReadPosition, rightReadPosition, newCiagrString, UtilLib.isTranscriptExon(genomicInterval.getChromosome()), isAtAlignmentMargin);

      oldGenomicInterval = genomicInterval;
      if (oldGenomicInterval.getChromosome ().equals(samReference)) {
	curExonEnd = oldGenomicInterval.getChrEnd ();
      }

      if (referenceOrientation.equals("*") && genomicInterval.getBedRecord().getReferenceGenomicAlignmentOrientation ().equals("-")) {
	System.out.println("Inverting ciagrString: " + ciagrString);
	ciagrString.invert ();
      }
    }

    if (debugLevel >= 2) {
      System.err.println("samRecord.isPairedInSequencing (): " + samRecord.isPairedInSequencing ());
    }
    
    int samRecordFlag = samRecord.getFlag ();
    samRecordFlag = UtilLib.setBit (samRecordFlag, 1, samRecord.hasMate () && mateGenomicIntervals != null);
    samRecordFlag = UtilLib.setBit (samRecordFlag, 3, samRecord.isPairedInSequencing () && mateGenomicIntervals == null);
    samRecordFlag = UtilLib.setBit (samRecordFlag, 4, samStrand.equals("-"));
    samRecordFlag = UtilLib.setBit (samRecordFlag, 5, samRecord.hasMate () && mateSamStrand.equals("-"));
    samRecordFlag = UtilLib.setBit (samRecordFlag, 8, ! isPrimaryAlignment);

    int mapQuality = samRecord.getMapQuality ();
    if (samRecord.isMapped()) {
      mapQuality = Math.max(mapQuality, 10);
    } else {
      mapQuality = 0;
    }
    
    SamRecord genomicSamRecord = new SamRecord (samRecord.getQueryName (), samRecord.getFragmentName (), samRecord.getFlag (), samReference, samPosition,
						mapQuality, newCiagrString.toString(), mateSamReference, mateSamPosition, insertSize,
						newCiagrString.getSequence (), newCiagrString.getQualityString(), samRecord.getOptionalParameters());

    /* Set the bits of the SAM flag. Take some of the SAM flag from the original flag: 0 is for paired-end sequencing, 2 whether the SAM record
       is mapped, 6 and 7 for first and second read of the pair - all of these do not change */
    genomicSamRecord.setHasMate (samRecord.hasMate () && mateGenomicIntervals != null);
    genomicSamRecord.setMateIsUnmapped (genomicSamRecord.isPairedInSequencing () && mateGenomicIntervals == null);
    genomicSamRecord.setIsOnReverseStrand (samStrand.equals("-"));
    genomicSamRecord.setMateIsOnReverseStrand (genomicSamRecord.hasMate() && mateSamStrand.equals("-"));
    genomicSamRecord.setIsSecondaryAlignment (! isPrimaryAlignment);

    if (genomicSamRecord.getFlag () != samRecordFlag) {
      throw new IOException ("Different flag values: " + genomicSamRecord.getFlag () + " vs " + samRecordFlag);
    }

    if (debugLevel >= 2) {
      System.err.println("samRecord: " + samRecord);
      System.err.println("genomicSamRecord: " + genomicSamRecord);
    }

    
    return genomicSamRecord;

  }

  
  /***********************************************************************************
   * 
   *                          computeGenomeSamPairs
   *
   ***********************************************************************************/

  public Vector <SamRecordPair> computeGenomeSamPairs (Vector<SamRecord> samRecords, HashSet<String> chromosomeIds) throws IOException {

    if (debugLevel >= 2 || false) {
      System.out.println("Print genomic entries.");
    }

    /* Compute genomic alignment strings and adjust them if necessary */
    adjustGenomicAlignments ();
    
    /* Combine the exons of weight object alignments with the same genomic alignment string */
    mergeWeightAlignments ();

    Hashtable<String, SamRecord> samRecordIndex = new Hashtable<String, SamRecord> (2 * samRecords.size());
    boolean pairedEnd = false;
    if (debugLevel >= 2 || false) {
      System.out.println("SAM records: " + samRecords);
    }

    for (SamRecord samRecord: samRecords) {
      String indexString = samRecord.getReferenceName () + "/" + samRecord.getPosition () + "/" + (samRecord.isSecondRead()?"2":"1");
      if (debugLevel >= 2) {
	System.err.println (indexString + ": Putting SAM record: " + samRecord);
      }
      samRecordIndex.put(indexString, samRecord);
      if (! pairedEnd && samRecord.isPairedInSequencing()) {
	pairedEnd = true;
      }
    }

    Counter counter = new Counter (5);
    SamRecord samRecord1 = null;
    SamRecord samRecord2 = null;
    HashSet<SamRecordPair> genomicSamPairSet    = new HashSet<SamRecordPair> (10);
    Vector <SamRecordPair> genomicSamPairVector = new Vector <SamRecordPair> (10);

    boolean isPrimaryAlignment = true;
    WeightObjectAlignment oldWeightObjectAlignment = null;
    for (WeightObjectAlignment weightObjectAlignment: getWeightObjectAlignmentSet()) {

      SamRecordPair samRecordPair = new SamRecordPair ();

      if (debugLevel >= 2  || false) {
	System.out.println("weightObjectAlignment: " + weightObjectAlignment + ", oldWeightObjectAlignment: " + oldWeightObjectAlignment);
	System.out.println("weightObjectAlignment - genomic alignment string " + weightObjectAlignment.getGenomicAlignmentString ()  +
			   ", alignment base id: " + weightObjectAlignment.getAlignmentBaseId() + ":");
	System.out.flush();
      }

      if (! weightObjectAlignment.getGenomicAlignmentString ().equals(weightObjectAlignment.getAlignmentBaseId() + ":") &&
	  (oldWeightObjectAlignment == null || ! weightObjectAlignment.getGenomicAlignmentString ().equals(oldWeightObjectAlignment.getGenomicAlignmentString ()))) {

	if (debugLevel >= 2) {
	  System.out.println("Computing genomic SAM records for: " + weightObjectAlignment);
	  System.out.println ("Genomic alignment string for " + weightObjectAlignment + ": " + weightObjectAlignment.getGenomicAlignmentString ());
	  System.out.println ("Old genomic alignment string for " + oldWeightObjectAlignment + ": " +
			      ((oldWeightObjectAlignment == null)?"":oldWeightObjectAlignment.getGenomicAlignmentString ()));
	  if (oldWeightObjectAlignment != null) {
	    System.out.println (weightObjectAlignment.getGenomicAlignmentString () + ".equals(" + oldWeightObjectAlignment.getGenomicAlignmentString () + "): " +
				weightObjectAlignment.getGenomicAlignmentString ().equals(oldWeightObjectAlignment.getGenomicAlignmentString ()));
	  }
	}

	GenomicAlignment genomicAlignment1 = weightObjectAlignment.getGenomicAlignment1();
	GenomicAlignment genomicAlignment2 = weightObjectAlignment.getGenomicAlignment2();
	String referenceName = weightObjectAlignment.getReferenceId ();

	String bedEntryId = fragmentId + "-A" + counter.toString() + "-C/" + (pairedEnd?"P":"S");

	Vector<String> genomicBedEntries1 = null;
	Vector<GenomicInterval> genomicIntervals1 = null;
	if (genomicAlignment1 != null && ! genomicAlignment1.isTranscriptExonAlignment ()) {
	  int readAlignStart1 = genomicAlignment1.getReadAlignStart ();
	  samRecord1          = samRecordIndex.get(referenceName + "/" + readAlignStart1 + "/1");

	  if (samRecord1 == null) {
	    throw new IOException ("No SAM record for " + bedEntryId + "1 found.");
	  }
	  genomicIntervals1  = new Vector<GenomicInterval> (genomicAlignment1.getGenomicIntervals());
	}

	Vector<String> genomicBedEntries2 = null;
	Vector<GenomicInterval> genomicIntervals2 = null;
	if (genomicAlignment2 != null  && ! genomicAlignment2.isTranscriptExonAlignment ()) {
	  int readAlignStart2 = genomicAlignment2.getReadAlignStart ();
	  samRecord2          = samRecordIndex.get(referenceName + "/" + readAlignStart2 + "/2");

	  if (samRecord2 == null) {
	    throw new IOException ("No SAM record for " + bedEntryId + "2 found.");
	  }
	  genomicIntervals2  = new Vector<GenomicInterval> (genomicAlignment2.getGenomicIntervals());
	}

	if (debugLevel >= 2 || false) {
	  System.out.println("Sam record 1: " + samRecord1);
	  System.out.println("Sam record 2: " + samRecord2);
	  System.out.println("genomicBedEntries1 for " + bedEntryId + ": " + genomicBedEntries1 + "(" + genomicAlignment1 +")");
	  System.out.println("genomicBedEntries2: " + genomicBedEntries2 + "(" + genomicAlignment2 +")");
	  if (genomicAlignment1 != null) {
	    System.out.println("genomicBedRecords1: " + genomicAlignment1.getBedRecords());
	  }	  
	  if (genomicAlignment2 != null) {
	    System.out.println("genomicBedRecords2: " + genomicAlignment2.getBedRecords());
	  }
	  System.out.flush ();
	}
	
	if (genomicAlignment1 != null) {
	  if (genomicIntervals1 != null) {
	    String referenceOrientation1 = genomicAlignment1.getReferenceGenomicAlignmentOrientation ();
	    SamRecord genomicSamRecord1 =
	      getGenomicSamRecord (samRecord1, genomicIntervals1, genomicIntervals2, referenceOrientation1, chromosomeIds, isPrimaryAlignment);

	    if (! UtilLib.isTranscriptExon(genomicSamRecord1.getReferenceName ())) {
	      samRecordPair.addSamRecord (genomicSamRecord1);
	    } else {
	      System.out.println ("Reference name for genomicAlignment1: " + genomicSamRecord1.getReferenceName ());
	    }
	  }
	}
      
	if (debugLevel >= 2) {
	  System.out.println ("genomicAlignment2: " + genomicAlignment2);
	}
	    
	if (genomicAlignment2 != null) {
	  if (debugLevel >= 2) {
	    System.out.println ("genomicIntervals2: " + genomicIntervals2);
	  }
	  if (genomicIntervals2 != null) {
	    String referenceOrientation2 = genomicAlignment2.getReferenceGenomicAlignmentOrientation ();
	    SamRecord genomicSamRecord2 =
	      getGenomicSamRecord (samRecord2, genomicIntervals2, genomicIntervals1, referenceOrientation2, chromosomeIds, isPrimaryAlignment);
	    	    
	    if (! UtilLib.isTranscriptExon(genomicSamRecord2.getReferenceName ())) {
	      samRecordPair.addSamRecord (genomicSamRecord2);
	    } else {
	      System.out.println ("Reference name for genomicAlignment2: " + genomicSamRecord2.getReferenceName ());
	    }
	  }
	}

	/* Check for duplicate genomic SAM records: due to differing tr-exons reads may align to only one genomic location yet may
	   have a number of different genomicAlignmentStrings; this is only a problem if the count mode is set to gene as only in this
	   case tr-exons contribute to the genomicAlignmentString and the count mode is set to gene for the creation of the genomic
	   alignments but not for the computation of the read weights, for instance. */
	if (! genomicSamPairSet.contains (samRecordPair)) {
	  genomicSamPairSet.add(samRecordPair);
	  genomicSamPairVector.add(samRecordPair);
	}

	counter.inc();

      }

      isPrimaryAlignment = false;

      oldWeightObjectAlignment = weightObjectAlignment;
    }

    return genomicSamPairVector;
    
  }

  
  /***********************************************************************************
   *
   *                         getCountObjectIds
   *
   ***********************************************************************************/


  public HashSet<String> getCountObjectIds (int overlapThreshold, Hashtable<String, Integer> countObjectLengthTable) throws IOException {

    HashSet<String> countObjectIds = new HashSet<String> ();
    
    for (WeightObjectAlignment weightObjectAlignment: weightObjectAlignments) {
      countObjectIds.addAll (weightObjectAlignment.getCountObjectIds (overlapThreshold, countObjectLengthTable));
    }

    return countObjectIds;
    
  }


  /***********************************************************************************
   *
   *                         printExonSet
   *
   ***********************************************************************************/


  public void printExonSet () {

    System.out.println("Exon set of " + this + " (number of weightObjectAlignments: " + weightObjectAlignments.size() + ")");

    /* Excatly one of the two sets weightObjectAlignments and weightObjectAlignmentSet is non-empty */
    for (WeightObjectAlignment weightObjectAlignment: weightObjectAlignments) {
      weightObjectAlignment.printExonSet ();
    }

    for (WeightObjectAlignment weightObjectAlignment: weightObjectAlignmentSet) {
      weightObjectAlignment.printExonSet ();
    }
    
  }


    
  /***********************************************************************************
   *
   *                         Basic methods
   *
   ***********************************************************************************/

  public int getNumGenomicAlignments () throws IOException {

    if (debugLevel >= 2) {
      System.out.println("weight object alignments:");
      for (WeightObjectAlignment weightObjectAlignment: weightObjectAlignmentSet) {
	System.out.println(weightObjectAlignment);	
      }
    }
        
    adjustGenomicAlignments ();

    if (debugLevel >= 2) {
      System.err.println ("Printing exon set after adjusting genomic alignments.");
      printExonSet ();
    }

    mergeWeightAlignments ();

    if (debugLevel >= 2) {
      System.err.println ("Printing exon set after merging genomic alignments.");
      printExonSet ();
    }
    
    int leftSize  = 0;
    int rightSize = 0;
    
    boolean isPairedEnd                = false;
    boolean hasTranscriptExonAlignment = false;
    for (WeightObjectAlignment weightObjectAlignment: weightObjectAlignmentSet) {

      if (! weightObjectAlignment.isTranscriptExonAlignment ()) {
	if (leftSize + rightSize > 0 && weightObjectAlignment.isPairedEnd () != isPairedEnd && ! weightObjectAlignment.fromGenomicBedRecords ()) {
	  System.err.println ("WARNING: Change of alignment mode to " + (weightObjectAlignment.isPairedEnd ()?"paired-end":"single-read") +
			      " for fragment " + fragmentId + " in aligment " + weightObjectAlignment);
	}
	
	isPairedEnd = weightObjectAlignment.isPairedEnd ();
      
	String alignmentString = weightObjectAlignment.getGenomicAlignmentString ();
	
	if (alignmentString != "") {
	  if (alignmentString.indexOf(":r-") > 0) {
	    rightSize++;
	  } else {
	    leftSize++;
	  }
	}

	if (debugLevel >= 2) {
	  System.err.println("Alignment string: " + alignmentString + ", leftSize: " + leftSize + ", rightSize: " + rightSize);
	}

      } else {
	hasTranscriptExonAlignment = true;
      }
    } 

    /* Consistent with FragmentEntry we define the weight of a read a the maximum of the number of
       alignments of the first and second read */
    return Math.max(leftSize, rightSize) + (hasTranscriptExonAlignment?1:0);
    
  }

  public String getReadAlignmentId () {
    return readAlignmentId;
  }


  public String getFragmentId () {
    return fragmentId;
  }

  public HashSet<WeightObjectAlignment> getWeightObjectAlignmentSet () {
    return weightObjectAlignmentSet;
  }

    
  public boolean equals (Object o) {
    
    WeightObject w = (WeightObject) o;
    return fragmentId.equals(w.getFragmentId());
    
  }

  public String toString () {
    
    try {
      if (false) {
	throw new IOException ("");
      }
      
      return fragmentId;
    }
    catch (IOException e) {
      System.out.println ("ERROR: " + (e==null?"Null message":e.getMessage()));
      System.exit (1);
    }

    return "";
  }


  public String toWeightString () throws IOException {
    return fragmentId + "\t" + getNumGenomicAlignments ();
  }
  

  public String toBedString () throws IOException {

    adjustGenomicAlignments ();
    mergeWeightAlignments ();

    String bedString = "";
    for (WeightObjectAlignment weightObjectAlignment: weightObjectAlignmentSet) {
      Vector<BedEntry> bedEntries = weightObjectAlignment.getBedEntries ();
      for (BedEntry bedEntry: bedEntries) {
	bedEntry.setName (fragmentId);
	bedString = bedString + bedEntry.toString () + "\n";
      }
    }

    return bedString;
    
  }
  
}
