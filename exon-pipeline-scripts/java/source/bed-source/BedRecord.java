/**File: BedRecord.java 

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
 *                              Class BedRecord
 *
 * The class BedRecord represents one line of the BED intersection file created by
 * intersect bed. It can be sorted according to referenceSequenceId, read index, read start
 * and end, and exon start and end.
 *
 * A bed record can be converted to a genomic interval.
 *
 ***********************************************************************************/

class BedRecord implements Comparable<BedRecord> {
  
  private static int debugLevel = UtilLib.getDebugLevel ();

  private static Pattern versionPattern = Pattern.compile("[.][0-9]+$");
  private static HashSet<String> incorrectStartCoordinates = new HashSet<String> (10000);

  public static void setDebugLevel (int value) {
    debugLevel = value;
  }
  
  /***********************************************************************************
   *
   *                       Object variables
   *
   ***********************************************************************************/

  /* From the input bed record */
  private String referenceSequenceId = "";
  private int    readAlignStart  = -1;
  private int    readAlignEnd    = -1;
  private String bedRecordReadId = "";    //read id and alignment id as in the bed file (e.g. <fragment id>-A0001-L<len>-I<#ins>-D<#del>-[SEBC]/P1)
  private String readAlignmentId = "";    //read id and alignment id as in the bed file (e.g. <fragment id>-A0001/P1)
  private String readId          = "";    //<fragment id>/[SP][12]  
  private String alignmentBaseId = "";    //read id for single read alignments and fragment id for paired-end
  private double alignScore      = 0;
  private String alignStrand     = "";

  private int fieldIndex = 0;
  private int maxFieldIndex = 0;

  private int exonReferenceStart  = -1;
  private int exonReferenceEnd    = -1;
  private String exonId              = "";
  private Exon   exon                = null;
  private double exonReferenceScore  = -1;
  private String exonReferenceStrand = "";
  private int    overlap             = 0;

  private String genomicIntervalStrand = "";

  /* Computed */
  private int adaptedExonReferenceStart = -1;
  private int adaptedExonReferenceEnd   = -1;

  private String fragmentId = "";
  private String alignmentId = "";

  private int readAlignedLength = 0;
  private int numInsertions     = 0;
  private int numDeletions      = 0;
  
  private int readIndex = -1;

  private boolean pairedEnd = true;

  private int genomeAlignmentStart = 0;
  private int genomeAlignmentEnd   = 0;

  private int startDiff = -1;
  private int endDiff = -1;

  private int startDifference = 0;
  private int lengthDifference = 0;

  private boolean isFirstGenomicInterval = false;
  private boolean isLastGenomicInterval  = false;

  private boolean isTranscriptExonAlignment = false;

  private int maxDiff = 0;

  
  /***********************************************************************************
   *
   *                            Contructors
   *
   ***********************************************************************************/


  public BedRecord (String bedRecordString) throws IOException {

    parseBedRecordString (bedRecordString);
    
    int genomeExonStart = exon.getChrStart();
    int genomeExonEnd   = exon.getChrEnd();

    int referenceExonDiff = exonReferenceEnd - exonReferenceStart + 1;
    int genomeExonDiff    = genomeExonEnd - genomeExonStart + 1;

    if (debugLevel >= 2) {
      System.err.println ("exonReferenceStart: " + exonReferenceStart + " exon reference end - start: " + referenceExonDiff +
			  ", genome end - start: " + genomeExonDiff + ", is genomic record: " + isGenomicRecord());
    }


    /* Some 5' exons start before the annotated start of the transcript; since bedTools does not allow for negative
       coordinates this is cannot be captured in the BED file so we correct for it here */
    if (exonReferenceStart == 1 && referenceExonDiff < genomeExonDiff && ! isGenomicRecord()) {
      if (debugLevel >= 2) {
	System.err.println("Setting exon reference start");
      }
      exonReferenceStart = exonReferenceStart - (genomeExonDiff - referenceExonDiff);
    }

    if (debugLevel >= 2) {
      System.err.println ("exonReferenceStart: " + exonReferenceStart);
    }

    adaptedExonReferenceStart = exonReferenceStart;
    adaptedExonReferenceEnd   = exonReferenceEnd;

    /* Since the exons of some RefSeq transcripts are in the wrong order (the first exon is at the end of the
       transcript sequence etc), we need to change the strand of the exons of these transcripts. Otherwise the
       assumed order of the exons is incorrect. When creating the transcript exon file, we encode the inversion
       of the strand in the referenceSequenceId by appending the suffix invertStrandSuffix. */
    String invertStrandSuffix = UtilLib.getInvertStrandSuffix ();
    if (referenceSequenceId.endsWith (invertStrandSuffix)) {
      referenceSequenceId = referenceSequenceId.substring(0, referenceSequenceId.length() - invertStrandSuffix.length());
      if (debugLevel >= 1) {
	System.err.println ("Inverting strand for BED record: " + this);
      }
      exon.invertStrand ();
    }

    /* A junction is composed of the concatenation of the *plus* strand genome sequences of the (parts of)
       the exon intervals. So for a junction the exon strand is always "+" although the exon may be part
       of a transcript that aligns to the reverse strand of the genome. Of course, in this case the exon strand
       for the transcript is still "-". This plays an important role in the interpretation of the overlap intervals
       in computeGenomicInterval. */
    if (referenceSequenceId.indexOf("junction") != -1) {
      exon.setStrand("+");
    }
    
    /* Remove the version to match the chromosome in Exon */
    if (! exon.isGenomicExon ()) {
      Matcher versionMatcher = versionPattern.matcher(referenceSequenceId);
      referenceSequenceId = versionMatcher.replaceAll("");
    }

    Pattern pattern = null;
    Matcher matcher = null;

    int lastAIndex = bedRecordReadId.lastIndexOf("-A");
    int lastLIndex = bedRecordReadId.indexOf("-L", lastAIndex);

    if (debugLevel >= 2) {
      System.out.println ("bedRecordReadId: " + bedRecordReadId + ", last A index: " + lastAIndex + ", last L index: " + lastLIndex);
    }

    if (lastLIndex > 0 && bedRecordReadId.indexOf("-I", lastLIndex) > 0 && bedRecordReadId.indexOf("-D", lastLIndex) > 0) {
    
      int lastIIndex      = bedRecordReadId.indexOf("-I", lastLIndex);
      int lastDIndex      = bedRecordReadId.indexOf("-D", lastIIndex);
      int lastFIndex      = bedRecordReadId.indexOf("-F", lastDIndex);
      int lastHyphenIndex = bedRecordReadId.indexOf("-",  lastFIndex+1);
      int lastSlashIndex  = bedRecordReadId.indexOf("/",  lastHyphenIndex);

      fragmentId  = bedRecordReadId.substring(0, lastAIndex);
      alignmentId = bedRecordReadId.substring(0, lastLIndex);

      readAlignedLength = Integer.parseInt (bedRecordReadId.substring(lastLIndex + 2, lastIIndex));
      numInsertions     = Integer.parseInt (bedRecordReadId.substring(lastIIndex + 2, lastDIndex));

      if (debugLevel >= 2) {
	System.out.println ("Bed record last D index + 2: " + (lastDIndex + 2) + ", lastSlashIndex: " + (lastSlashIndex) + " : " +
			  bedRecordReadId.substring(lastHyphenIndex+1, lastSlashIndex));
      }
      
      numDeletions  = Integer.parseInt (bedRecordReadId.substring(lastDIndex + 2, lastFIndex));
      fieldIndex    = Integer.parseInt (bedRecordReadId.substring(lastFIndex + 2, lastHyphenIndex));
      maxFieldIndex = Integer.parseInt (bedRecordReadId.substring(lastHyphenIndex + 1, lastSlashIndex));

      if (fieldIndex == 0) {
	isFirstGenomicInterval = true;
      }

      if (fieldIndex == maxFieldIndex) {
	isLastGenomicInterval   = true;
      }

    } else {

      pattern = Pattern.compile("-A[0-9]*/[SP][12]$");
      matcher = pattern.matcher(bedRecordReadId);
    
      fragmentId = matcher.replaceAll("");

      pattern = Pattern.compile("/[SP][12]$");
      matcher = pattern.matcher(bedRecordReadId);

      alignmentId = matcher.replaceAll("");

      isFirstGenomicInterval = true;
      isLastGenomicInterval  = true;
      
    }

    if (bedRecordReadId.endsWith("1") || bedRecordReadId.endsWith("2")) {
      readIndex = UtilLib.toInt(bedRecordReadId.substring(bedRecordReadId.length()-1, bedRecordReadId.length()));
    }

    if (debugLevel >= 2) {
      System.out.println("Read index for read " + bedRecordReadId + ": " + readIndex);
    }

    String bedRecordReadIdEnd = bedRecordReadId.substring(bedRecordReadId.length()-2, bedRecordReadId.length());
    
    if (bedRecordReadIdEnd.equals("S1") || bedRecordReadIdEnd.equals("S2")) {
      readId = fragmentId + "/S" + readIndex;
      readAlignmentId = alignmentId + "/S" + readIndex;
      pairedEnd = false;
    } else {
      readId = fragmentId + "/P" + readIndex;
      readAlignmentId = alignmentId + "/P" + readIndex;
    }

    alignmentBaseId = readId;
    if (pairedEnd) {
      alignmentBaseId = fragmentId;
    }

  }

  /***********************************************************************************
   *
   *                       parseBedRecordString
   *
   ***********************************************************************************/

  private void parseBedRecordString (String bedRecordString) throws IOException {
    
    int i = 0;
    try {
      /* StringTokenizer st = new StringTokenizer (bedRecordString, "\t");
	 String token = "";
	 while (st.hasMoreTokens()) {
	   token = st.nextToken();
      */

      int pos = 0;
      while (pos < bedRecordString.length()) {

	int end = bedRecordString.indexOf ("\t", pos);
	if (end == -1) {
	  end = bedRecordString.length();
	}
	String token = bedRecordString.substring(pos, end);
	pos = end + 1;

	if (i == 0) {
	  referenceSequenceId = token;
	  i++;
	} else if (i == 1) {
	  readAlignStart = UtilLib.toInt (token) + 1;
	  i++;
	} else if (i == 2) {
	  readAlignEnd = UtilLib.toInt (token);
	  i++;
	} else if (i == 3) {
	  bedRecordReadId = token;
	  i++;
	} else if (i == 4) {
	  alignScore = UtilLib.toDouble(token);
	  i++;
	} else if (i == 5) {
	  alignStrand = token;
	  i++;
	} else if (i == 6) {
	  if (! referenceSequenceId.equals(token)) {
	    throw new IOException ("Reference sequence id: " + referenceSequenceId + " does not equal: " + token);
	  };
	  i++;
	} else if (i == 7) {
	  exonReferenceStart = UtilLib.toInt (token) + 1;
	  i++;
	} else if (i == 8) {
	  exonReferenceEnd = UtilLib.toInt (token);
	  i++;
	} else if (i == 9) {
	  exonId = token;
	  exon = new Exon (exonId, "/", referenceSequenceId);
	  isTranscriptExonAlignment = exon.isTranscriptExon ();
	  i++;
	} else if (i == 10) {
	  exonReferenceScore = UtilLib.toDouble(token);
	  i++;
	} else if (i == 11) {
	  exonReferenceStrand = token;
	  i++;
	} else if (i == 12) {
	  overlap = UtilLib.toInt (token);
	  i++;
	}

      }

      if (i < 13) {
	throw new IOException ("Too few BED entries: " + i);
      }

    } catch (Exception e) {
      throw new IOException ("ERROR: problem with bed record: " + bedRecordString + "  when reading field: " + i + " -> " + (e==null?"Null message":e.getMessage()));
    }

  }

  /***********************************************************************************
   *
   *                              compare
   *
   ***********************************************************************************/

  public int compareTo (BedRecord b) {

    if (! getAlignmentId ().equals(b.getAlignmentId ())) {
      return alignmentId.compareTo(b.getReferenceSequenceId ());
    }

    if (! referenceSequenceId.equals(b.getReferenceSequenceId ())) {
      return referenceSequenceId.compareTo(b.getReferenceSequenceId ());
    }

    if (readIndex != b.getReadIndex ()) {
      return UtilLib.compareInt (readIndex, b.getReadIndex ());
    }

    /* if (fieldIndex != b.getFieldIndex ()) {
      return UtilLib.compareInt (fieldIndex, b.getFieldIndex ());
      } */

    if (readAlignStart != b.getReadAlignStart ()) {
      return UtilLib.compareInt (readAlignStart, b.getReadAlignStart ());
    }

    if (readAlignEnd != getReadAlignEnd ()) {
      return UtilLib.compareInt (readAlignEnd, getReadAlignEnd ());
    }

    if (exonReferenceStart != b.getExonReferenceStart ()) {
      return UtilLib.compareInt (exonReferenceStart, b.getExonReferenceStart ());
    }

    if (exonReferenceEnd != b.getExonReferenceEnd ()) {
      return UtilLib.compareInt (exonReferenceEnd, b.getExonReferenceEnd ());
    }

    return exonReferenceStrand.compareTo(b.getExonReferenceStrand ());

  }

      
  /***********************************************************************************
   *
   *                Compute the genomic alignment of a read
   *
   *  Use the exon coordinates to infer the genomic coordinates of the read alignment.
   *  Differences in the length of the exon on the genome vs the transcript pose a
   *  problem.
   *
   *  index: current index in sequence of bedRecords that belong to the GenomicAlignment
   *     index == 0            => start of sequence
   *     index == maxIndex - 1 => end of sequence
   *
   ***********************************************************************************/

  public GenomicInterval computeGenomicInterval (int index, int maxIndex, GenomicAlignment genomicAlignment)
    throws IOException {

    int genomeExonStart = exon.getChrStart();
    int genomeExonEnd   = exon.getChrEnd();

    int diffTranscript = adaptedExonReferenceEnd - adaptedExonReferenceStart + 1;
    int diffGenome     = genomeExonEnd - genomeExonStart + 1;

    int lengthDifference = diffTranscript - diffGenome;
    int startDifference  = 0;

    /* Check for a length difference of the exon mapped to the transcript vs the length of the genomic mapping */
    if (lengthDifference != 0) {
	
      if (! isJunction ()) {
	/* Note that we should have that the number of aligned genome bases plus the genome insertions equal the number of aligned transcript bases
	   plus the transcript insertsions (= genome deletions) */
	if (diffTranscript + numDeletions > diffGenome + numInsertions  && (debugLevel >= 2 || Math.abs(diffTranscript - diffGenome) > 100)) {
	  String alignmentString = referenceSequenceId + "-" + exon + "-" +  diffTranscript + "-" +  diffGenome;
	  if (! incorrectStartCoordinates.contains (alignmentString)) {
	    System.err.println("WARNING: Incorrect start coordinate on transcript for read " + readAlignmentId + " on transcript: "
			       + referenceSequenceId + ", exon " + exon + " length on transcript: " + diffTranscript + " vs on genomeExon " + diffGenome + ".");
	    System.err.println("Alignment string: " + alignmentString);
	    incorrectStartCoordinates.add(alignmentString);
	  }
	  maxDiff = Math.max(maxDiff, Math.abs(diffTranscript - diffGenome));
	}
      }

      if (debugLevel >= 2) {
	System.out.println("Index: " + index + ", maxIndex: " + maxIndex + ".");
      }

      startDifference = 0;
      if (isJunction () && index == 0) {
	startDifference = diffTranscript - diffGenome;
	adaptedExonReferenceStart = adaptedExonReferenceEnd - diffGenome + 1;
      } else if (index == maxIndex - 1) {
	adaptedExonReferenceEnd = adaptedExonReferenceStart + diffGenome - 1;
      } else if (index == 0) {
	startDifference = diffTranscript - diffGenome;
	adaptedExonReferenceStart = adaptedExonReferenceEnd - diffGenome + 1;
      } else {
	startDifference = diffTranscript - diffGenome/2;
	adaptedExonReferenceStart = adaptedExonReferenceEnd   - diffGenome/2 + 1;
	adaptedExonReferenceEnd   = adaptedExonReferenceStart + (diffGenome - diffGenome/2) - 1;
      }
    }


    /********************************************************************************
     **
     ** The overlap on the transcripts starts with the maximum of the read start and 
     ** exon start and ends with the minimum of the read end and exon end
     **
     ********************************************************************************/

    int overlapStart = Math.max(readAlignStart, adaptedExonReferenceStart);
    int overlapEnd   = Math.min(readAlignEnd,   adaptedExonReferenceEnd);

    startDiff = Math.abs(overlapStart - adaptedExonReferenceStart);
    endDiff   = Math.abs(overlapEnd   - adaptedExonReferenceEnd);

    if (debugLevel >= 2) {
      System.out.println (startDiff + " = " + overlapStart + " - " + adaptedExonReferenceStart);
      System.out.println (endDiff   + " = " + overlapEnd   + " - " + adaptedExonReferenceEnd);
    }
    

    /********************************************************************************
     *
     *  Compute the genomic interval of the intersection of the alignment of the read
     *  and the exon.
     *
     *  The strand of the exon is the genomic orientation of the exon and, thus, of the
     *  alignment of the read. Note that the read coordinates are in the correct
     *  orientation if the read maps to the genome (start is before end) or the
     *  read maps to a transcript (start on transcript is before end on transcript)
     *  and the orientation of the exon on the genome is "+". If the orientation
     *  of the exon on the genome is "-", then the order of start and end of the
     *  genomic read alignment needs to get reversed in the process of mapping the
     *  read to the genome.
     *
     *  Since a junction is composed of the concatenation of the *plus* strand genome
     *  sequences, we already set the exon strand to "+" and, therefore, the read
     *  alignment strand is just the alignment strand.
     *
     ********************************************************************************/

    String exonStrand = "+";
    if (isGenomicRecord() || referenceSequenceId.indexOf("junction") != -1) {
      genomicIntervalStrand = getAlignStrand();
    } else {
      exonStrand = exon.getStrand();
      genomicIntervalStrand = exon.getStrand ().equals(getAlignStrand())?"+":"-";
    }
    
    int overlapDiff = overlap - 1;

    if (startDiff <= endDiff) {
      if (exonStrand.equals("+")) {
	genomeAlignmentStart = genomeExonStart + startDiff;
	genomeAlignmentEnd   = genomeExonStart + startDiff + overlapDiff;
      } else {
	genomeAlignmentStart = Math.max(1, genomeExonEnd - startDiff - overlapDiff);
	genomeAlignmentEnd   = genomeExonEnd - startDiff;
      }
    } else {
      if (exonStrand.equals("+")) {
	genomeAlignmentStart = Math.max(1, genomeExonEnd - endDiff - overlapDiff);
	genomeAlignmentEnd   = genomeExonEnd - endDiff;
      } else {
	genomeAlignmentStart = genomeExonStart + endDiff;
	genomeAlignmentEnd   = genomeExonStart + endDiff + overlapDiff;
      }
    }

    boolean conformingAlignment = isConformingAlignment ();

    if (debugLevel >= 2 || ComputeCounts.getSpecialReadIdSet ().contains(getFragmentId ())) {
      System.out.println(referenceSequenceId + "/" + readAlignStart  + "/" + readAlignEnd + " vs. " + exon + " conforming: " + conformingAlignment);
    }

    if (debugLevel >= 2) {
      System.out.println("Exon of BedRecord before creating GenomicInterval: " + exon);
    }

    /* Note that the strand of the exon defines the strand of the alignment of the transcript against the genome.
       So if the strand of the alignment of the read against the transcript equals the strand of alignment of the
       transcript against the genome (+ and + or - and -), then the strand of the genomic interval on the genome
       is plus, otherwise, it is minus. */
    GenomicInterval genomicInterval =
	new GenomicInterval(exon.getChromosome (), getGenomeAlignmentStart (), getGenomeAlignmentEnd (), genomicIntervalStrand, getOverlap(), exon,
			    lengthDifference, readAlignedLength, numInsertions, numDeletions, index, genomicAlignment, this, conformingAlignment);

    if (debugLevel >= 2) {
      System.out.println ("Creating a genomic interval " + genomicInterval + " for fieldIndex: " + index);
    }


    return genomicInterval;
  }


  /***********************************************************************************
   *
   *       Compute the genomic alignment of a read (for genomic BED records)
   *
   ***********************************************************************************/

  public GenomicInterval computeGenomicInterval (GenomicAlignment genomicAlignment) throws IOException {

    genomeAlignmentStart  = readAlignStart;
    genomeAlignmentEnd    = readAlignEnd;
    if (! isGenomicRecord()) {
      throw new IOException ("computeGenomicInterval called on none genome BED record: " + this);
    } 
    genomicIntervalStrand = getAlignStrand();
    
    if (debugLevel >= 2) {
      System.out.println ("Creating new genomic interval: " + exon.getChromosome () + "/" + getGenomeAlignmentStart () + "/" + getGenomeAlignmentEnd () +
			  " with overlap " + getOverlap () + " for exon: " + exon);
    }

    boolean conformingAlignment = isConformingAlignment ();

    if (debugLevel >= 2 || ComputeCounts.getSpecialReadIdSet ().contains(getFragmentId ())) {
      System.out.println(referenceSequenceId + "/" + readAlignStart  + "/" + readAlignEnd + " conforming: " + conformingAlignment);
    }
        
    return new GenomicInterval(exon.getChromosome (), getGenomeAlignmentStart (), getGenomeAlignmentEnd (), getAlignStrand (),
			       getOverlap(), exon, lengthDifference, readAlignedLength, numInsertions,
			       numDeletions, fieldIndex, genomicAlignment, this, conformingAlignment);
  }

  

  /***********************************************************************************
   *
   *                         getCountObjectIds
   *
   ***********************************************************************************/

  public HashSet<String> getCountObjectIds (int overlapThreshold, Hashtable<String, Integer> countObjectLengthTable) throws IOException {

    
    if (Exon.getCountObjectTable () == null || overlap < overlapThreshold) {
      return new HashSet<String> ();
    }


    HashSet<String> countObjectIds = Exon.getCountObjectTable ().get (exon);

    if (countObjectIds == null) {
      countObjectIds = new HashSet<String> ();
      countObjectIds.add (exon.toString());
    }
  
    return countObjectIds;
    
  }
	
    
  /***********************************************************************************
   *
   *                             alignmentEquals
   *
   ***********************************************************************************/

  public boolean alignmentEquals (BedRecord b) {

    if (b == null) {
      return false;
    }
    
    return referenceSequenceId.equals(b.getReferenceSequenceId ()) && readAlignStart == b.getReadAlignStart() &&
      readAlignEnd == b.getReadAlignEnd();
	    
  }


  /***********************************************************************************
   *
   *                             isConformingAlignment
   *
   *  A conforming alignment is an alignment of (a spliced part of) the read against
   *  the genomic interval of the exon where either:
   *
   *  - The end points of the (spliced part of the) read alignment and the exon
   *    coincide or
   *  - The read alignment starts at the beginning of the read and the start position
   *    is contained within the exon or
   *  - The read alignment ends at the end of the read and the end position
   *    is contained within the exon.
   *
   *
   *   Read:   |---------|      <-----...    ...----->
   *   Exon:   |---------|    |----------    -----------|
   *
   ***********************************************************************************/

  public boolean isConformingAlignment () throws IOException {

    /* For stranded libraries we want the genomic interval to have the same strand
       as the exon. */
    if (UtilLib.strandedMode () && ! exon.getOriginalStrand ().equals(genomicIntervalStrand)) {
      if (UtilLib.warningsOn ()) {
	System.err.println(this.toCompleteString ());
	System.err.println("exon strand: " + exon.getStrand () + " genomicIntervalStrand: " + genomicIntervalStrand);
      }
      return false;
    }

    /* Check if we aligned against a different reference sequence set, e.g. the transcripts or junctions; in this
       case the reads are correctly spliced since we map them back to the genome */
    if (! isGenomicRecord ()) {
      return true;
    }


    if (! UtilLib.spliceConformingCountMode() && ! UtilLib.containmentMode()) {
      return true;
    }
      
    if (! (exon.getChrStart() == getExonReferenceStart()) || exon.getChrEnd() != getExonReferenceEnd()) {
      throw new IOException ("Alignment against genome and exon coodinates disagree in BedRecord:\n" +
			     toString () + "\t" + getReferenceSequenceId () + "\t" + getExonReferenceStart() + "\t" +
			     getExonReferenceEnd());
    }

    if (debugLevel >= 2) {
      System.out.println ("isFirstGenomicInterval (): " + isFirstGenomicInterval () + ", " + exon.getChrStart() + " > " + getReadAlignStart () + ", " +
			  "Math.abs(exon.getChrStart() - getReadAlignStart ()): " + Math.abs(exon.getChrStart() - getReadAlignStart ()) +
			  ", isLastGenomicInterval (): " + isLastGenomicInterval () + ", Math.abs(exon.getChrEnd() - getReadAlignEnd ()): " +
			  Math.abs(exon.getChrEnd() - getReadAlignEnd ()));
    }

    
    if (getReadAlignedLength () - getNumInsertions () +  getNumDeletions () < getOverlap ()) {
      throw new IOException ("Computed overlap and intersectBed overlap disagree in BedRecord:\n" +
			     toString () + "\n" +
			     getReadAlignedLength () + " - " + getNumInsertions () + " + " + getNumDeletions () + " != " + getOverlap());
    }

    if (debugLevel >= 2) {
      System.err.println ("Checking conforming mode for [" + getReadAlignStart () + ", " + getReadAlignEnd () + "]");
      System.err.println ("spliceConformingCountMode: " + UtilLib.spliceConformingCountMode() + ", containmentMode: " + UtilLib.containmentMode());
    }
    
    if (UtilLib.spliceConformingCountMode()) {
      /* It may happen that an exon intersects only a part of the aligned genomic interval */
      if (getReadAlignedLength () - getNumInsertions () +  getNumDeletions () > getOverlap ()) {	
	if (debugLevel >= 2) {
	  System.out.println ("Length of aligned genomic interval is larger than the overlap with the exon in BedRecord:\n" +
			      toString () + "\n" +
			      getReadAlignedLength () + " - " + getNumInsertions () + " + " + getNumDeletions () + " != " + getOverlap());
	}
	return false;
      }
    } else { /* UtilLib.containmentMode() */
      if (debugLevel >= 2) {
	System.err.println ("exon.getChrStart: " + exon.getChrStart() + ", getReadAlignStart: " + getReadAlignStart () + "; getReadAlignEnd: " + getReadAlignEnd () +
			    ", exon.getChrEnd: " + exon.getChrEnd());
      }
      return exon.getChrStart() <= getReadAlignStart () && getReadAlignEnd () <= exon.getChrEnd();
    }


    if (debugLevel >= 2) {
      System.out.println("Is first genomic interval: " + isFirstGenomicInterval () + ", exon chr start: " + exon.getChrStart() + ", read align start: " + getReadAlignStart ());
      System.out.println("Is last genomic interval: " + isLastGenomicInterval () + ", exon chr end: " + exon.getChrEnd() + ", read align end: " + getReadAlignEnd ());
    }

    /* Note that isFirstGenomicInterval () == true just implies this BED record represents the first genomic interval of the genomic alignment */
    if (isFirstGenomicInterval () && exon.getChrStart() > getReadAlignStart ()) {
      return false;
    }

    if (! isFirstGenomicInterval () && Math.abs(exon.getChrStart() - getReadAlignStart ()) > UtilLib.getExonStartSlack ()) {
      return false;
    }

    if (isLastGenomicInterval () && exon.getChrEnd() < getReadAlignEnd ()) {
      return false;
    }

    if (! isLastGenomicInterval () && Math.abs(exon.getChrEnd() - getReadAlignEnd ()) > UtilLib.getExonEndSlack()) {
      return false;
    }
      
    
    return true;
    
  }

  
  /***********************************************************************************
   *
   *                              Get methods
   *
   ***********************************************************************************/

  public String getReferenceSequenceId () {
    return referenceSequenceId;
  }

  public boolean isJunction () {
    return referenceSequenceId.indexOf("junction") >= 0;
  }

  public int getReadAlignStart () {
    return readAlignStart;
  }

  public int getReadAlignEnd () {
    return readAlignEnd;
  }

  public String getReadAlignmentId() {
    return readAlignmentId;
  }

  public double getAlignScore () {
    return alignScore;
  }

  public String getAlignStrand () {
    return alignStrand;
  }

  public int getExonReferenceStart () {
    return exonReferenceStart ;
  }

  public int getExonReferenceEnd () {
    return exonReferenceEnd;
  }

  public int getAdaptedExonReferenceStart () {
    return adaptedExonReferenceStart ;
  }

  public int getAdaptedExonReferenceEnd () {
    return adaptedExonReferenceEnd;
  }

  public String getExonId () {
    return exonId;
  }

  public Exon getExon () {
    return exon;
  }

  public String getReferenceGenomicAlignmentOrientation () {
    if (isGenomicRecord ()) {
      return "+";
    }

    return getExon().getStrand ();
    
  }

  public double getExonReferenceScore () {
    return exonReferenceScore;
  }

  public String getExonReferenceStrand () {
    return exonReferenceStrand;
  }

  public int getOverlap () {
    return overlap;
  }

  public String getFragmentId () {
    return fragmentId;
  }

  public String getAlignmentId () {
    return alignmentId;
  }

  public int getReadIndex () {
    return readIndex;
  }

  public String getReadId () {
    return readId;
  }

  public String getAlignmentBaseId () {
    return alignmentBaseId;
  }

  public boolean isPairedEnd () {
    return pairedEnd;
  }

  public String getAlignment () {
    return referenceSequenceId + "/" + readAlignStart + "/" + readAlignEnd;
  }


  public int getGenomeAlignmentStart () {
    return genomeAlignmentStart;
  }

  public int getGenomeAlignmentEnd () {
    return genomeAlignmentEnd;
  }

 
  public boolean isTranscriptExonAlignment () {
    return isTranscriptExonAlignment;
  }

  public int getStartDiff () {
    return startDiff;
  }

  public int getEndDiff () {
    return endDiff;
  }

  public int getMaxDiff () {
    return maxDiff;
  }

  public int getReadAlignedLength () {
    return readAlignedLength;
  }

  public int getNumInsertions () {
    return numInsertions;
  }

  public int getNumDeletions () {
    return numDeletions;
  }

  public int getFieldIndex () {
    return fieldIndex;
  }

  public int getMaxFieldIndex () {
    return maxFieldIndex;
  }

  public void setFieldIndex (int value) {
    fieldIndex = value;
    
    if (debugLevel >= 2) {
      System.out.println("Setting fieldIndex to " + value);
    }

    if (fieldIndex == 0) {
      isFirstGenomicInterval = true;
    }

    if (fieldIndex == maxFieldIndex) {
      isLastGenomicInterval = true;
    }
    
  }

  public void setMaxFieldIndex (int value) {
    maxFieldIndex = value;

    if (debugLevel >= 2) {
      System.out.println("Setting maxFieldIndex to " + value);
    }

    if (fieldIndex == maxFieldIndex) {
      isLastGenomicInterval = true;
    }

  }

  public boolean isFirstGenomicInterval () {
    return isFirstGenomicInterval;
  }

  public boolean isLastGenomicInterval () {
    return isLastGenomicInterval;
  }

  public boolean isGenomicRecord () {
    return referenceSequenceId.equals(exon.getChromosome ());
  }

  public String toString() {
    return (referenceSequenceId + "\t" + readAlignStart + "\t" + readAlignEnd + "\t" + readAlignmentId + "\t" + exonId);
  }

  public String toCompleteString() {
    return (referenceSequenceId + "\t" + readAlignStart + "\t" + readAlignEnd + "\t" + bedRecordReadId + "\t" + alignScore + "\t" +
	    alignStrand + "\t" + referenceSequenceId + "\t" + exonReferenceStart + "\t" + exonReferenceEnd + "\t" + exonId + "\t" + 
	    exonReferenceScore + "\t" + exonReferenceStrand + "\t" + overlap);
  }

  

}

