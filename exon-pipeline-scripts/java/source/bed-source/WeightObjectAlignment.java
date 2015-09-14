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
 *                              Class WeightObjectAlignment
 *
 *  One weight object has one or more weight object alignments associated to it.
 *  A weight object alignment is identified by an alignment id which is given by the
 *  alignment id of the BedRecords it stores. It stores two sets of BedRecords (one
 *  for the first read and one for the second read).
 *
 *  When the file with the intersection BedRecords is read, each BedRecord is added
 *  to the weight object which in turn adds it to the current weightObjectAlignment
 *  if the alignment id of the BedRecord matches, otherwise a new weightObjectAligment
 *  is created.
 *
 *  The BedRecords are added to the first or second set depending on their readIndex.
 *
 *  One of he main function of a WeightObjectAlignment is to provide a representation of
 *  genomic alignment in terms of the genomicAlignmentString. The genomic alignment
 *  string has the format <genomic alignment1>[-<genomic alignment2>].
 *
 *  It also offers the function computeOverlaps which computes 
 *
 ***********************************************************************************/

class WeightObjectAlignment {

  private static int debugLevel = UtilLib.getDebugLevel ();
  
  /***********************************************************************************
   *
   *                Object variables and constructors
   *
   ***********************************************************************************/

  private String alignmentId = "";
  private String alignmentBaseId = "";
  private String genomicAlignmentString = null;
  private String fragmentId = "";
  private int    numReads = 0;
  private boolean pairedEnd = false;

  GenomicAlignment genomicAlignment1 = null;
  GenomicAlignment genomicAlignment2 = null;

  private TreeSet<BedRecord> bedRecords1 = null;
  private TreeSet<BedRecord> bedRecords2 = null;

  boolean fromGenomicBedRecords = false;

  private String referenceId = null;

  public WeightObjectAlignment () {
    alignmentId = "";
    numReads    = 0;

    bedRecords1 = new TreeSet<BedRecord> ();
    bedRecords2 = new TreeSet<BedRecord> ();
    
  }

  public WeightObjectAlignment (String alignmentId, String alignmentBaseId, String fragmentId) {

    debugLevel = UtilLib.getDebugLevel ();

    this.alignmentId     = alignmentId;
    this.alignmentBaseId = alignmentBaseId;
    numReads             = 0;
    this.fragmentId      = fragmentId;

    bedRecords1 = new TreeSet<BedRecord> ();
    bedRecords2 = new TreeSet<BedRecord> ();
   
  }

 
  /***********************************************************************************
   *
   *                              addBedRecord
   *
   ***********************************************************************************/

  public void addBedRecord (BedRecord bedRecord) throws Exception {

    if ((bedRecords1.size() > 0 || bedRecords2.size() > 0) && bedRecord.isPairedEnd () != pairedEnd) {
      throw new Exception ("BedRecord: " + bedRecord + " differs in paired-end status from previous records of the same alignment.");
    }

    if (referenceId == null) {
      referenceId = bedRecord.getReferenceSequenceId ();
    } else if (! referenceId.equals(bedRecord.getReferenceSequenceId ())) {
      throw new Exception ("BedRecord: " + bedRecord + " has a different reference sequence id than the weightObjectAlignment: " + referenceId);
    }
    
    pairedEnd = bedRecord.isPairedEnd ();

    int readIndex = bedRecord.getReadIndex ();
    if (readIndex == 1) {
      if (debugLevel >= 2) {
	System.out.println("Adding " + bedRecord.toString () + " to bedRecords1:");
	for (BedRecord curBedRecord: bedRecords1) {
	  System.out.println(curBedRecord);
	}
      }
      bedRecords1.add(bedRecord);
      if (debugLevel >= 2) {
	System.out.println("New bedRecords1:");
	for (BedRecord curBedRecord: bedRecords1) {
	  System.out.println(curBedRecord);
	}
      }

    } else if (readIndex == 2) {
      if (debugLevel >= 2) {
	System.out.println("Adding " + bedRecord.toString () + " to bedRecords2:");
	for (BedRecord curBedRecord: bedRecords2) {
	  System.out.println(curBedRecord);
	}
      }
      bedRecords2.add(bedRecord);      
      if (debugLevel >= 2) {
	System.out.println("New bedRecords2:");
	for (BedRecord curBedRecord: bedRecords2) {
	  System.out.println(curBedRecord);
	}
      }
    } else {
      throw new Exception ("Unknown read index for bed record: " + bedRecord.toString());
    }
    
  }


  /***********************************************************************************
   *
   *                 getGenomicAlignmentString
   *
   *  compute and/or return the genomic alignment string
   *
   ***********************************************************************************/

  public String getGenomicAlignmentString () {

    if (genomicAlignmentString != null) {
      return genomicAlignmentString;
    }

    return getAlignmentId () + ": No genomic alignment string";
    
  }

  
  /***********************************************************************************
   *
   *                 setGenomicAlignmentString
   *
   *  compute and return the genomic alignment string
   *
   ***********************************************************************************/

  public void setGenomicAlignmentString () throws IOException {

    if (genomicAlignment1 == null && genomicAlignment2 == null) {
      return;
    }

    String genomicAlignment1String = "";
    if (genomicAlignment1 != null) {
      genomicAlignment1String = genomicAlignment1.toGenomicAlignmentString ();
    }

    String genomicAlignment2String = "";
    if (genomicAlignment2 != null) {
      genomicAlignment2String = genomicAlignment2.toGenomicAlignmentString ();
    }

    if (genomicAlignment1String != "" && genomicAlignment2String != "") {
      if (genomicAlignment1.getStart() <= genomicAlignment2.getStart()) {
	genomicAlignmentString = genomicAlignment1String + "-" + genomicAlignment2String;
      } else {
	genomicAlignmentString = genomicAlignment2String + "-" + genomicAlignment1String;
      }
    } else {
      if (genomicAlignment2 == null) {
	genomicAlignmentString = "l-" + genomicAlignment1String;
      } else if (genomicAlignment1 == null) {
	genomicAlignmentString = "r-" + genomicAlignment2String;
      } else {
	genomicAlignmentString = genomicAlignment1String + genomicAlignment2String;
      }
    }

    genomicAlignmentString = alignmentBaseId + ":" + genomicAlignmentString;
    
  }

  
  /***********************************************************************************
   *
   *                 computeGenomicAlignmentString
   *
   *  compute and return the genomic alignment string
   *
   ***********************************************************************************/

  public void initializeGenomicAlignmentString () throws IOException {

    if (bedRecords1.size() > 0) {
      genomicAlignment1 = new GenomicAlignment (bedRecords1, this);
      fromGenomicBedRecords = genomicAlignment1.fromGenomicBedRecords ();
    }
      
    if (bedRecords2.size() > 0) {
      genomicAlignment2 = new GenomicAlignment (bedRecords2, this);
      fromGenomicBedRecords = genomicAlignment2.fromGenomicBedRecords ();
    }

    if (bedRecords1.size() > 0 && bedRecords2.size() > 0 && genomicAlignment1.fromGenomicBedRecords () != genomicAlignment2.fromGenomicBedRecords ()) {
      System.err.println ("Left and right alignments come from different BED record types for " + getAlignmentId ());
    }
 
    setGenomicAlignmentString ();
    
  }

  
  /***********************************************************************************
   *
   *                 adjustGenomicAlignmentString
   * Since the mapping of the transcript to the exons on the genome is sometimes
   * ambiguous (the exon can be too large or too small), there can be some ambiguity 
   * (or slack) in the genomic coordinates of the read. If the read is mapped to
   * a different exon, the alignment may look different though using the slack caused
   * caused by the ambiguity, it may be possible to unify the alignments.
   *
   * Called by adjustGenomicAlignments:
   * First getGenomicAlignmentString is called for each weightGenomicAlignment of a
   * weight object and then the genomic alignment strings are adjusted.
   *
   ***********************************************************************************/

  public void adjustGenomicAlignmentString (Vector<WeightObjectAlignment> weightObjectAlignments) throws IOException {

    if (genomicAlignment1 != null) {
      if (debugLevel >= 2) {
	System.out.println("Alignment 1 of " + alignmentId);
	System.out.println("Genomic intervals before adjustment:");
	System.out.println(genomicAlignment1.getGenomicIntervalString());
	System.out.print("Adjusting " + genomicAlignment1 + " to ");
      }
      genomicAlignment1.adjustGenomicAlignmentString (weightObjectAlignments, 1);
      if (debugLevel >= 2) {
	System.out.println (genomicAlignment1.toGenomicAlignmentString());
	System.out.println("Genomic intervals after adjustment:");
	System.out.println(genomicAlignment1.getGenomicIntervalString());
      }
    }

    if (genomicAlignment2 != null) {
      if (debugLevel >= 2) {
	System.out.println("Alignment 2 of " + alignmentId);
	System.out.println("Genomic intervals before adjustment:");
	System.out.println(genomicAlignment2.getGenomicIntervalString());
	System.out.print("Adjusting alignment 2 of " + alignmentId+ ": " + genomicAlignment2 + " to ");
      }
      genomicAlignment2.adjustGenomicAlignmentString (weightObjectAlignments, 2);
      if (debugLevel >= 2) {
	System.out.println (genomicAlignment2.toGenomicAlignmentString());
	System.out.println("Genomic intervals after adjustment:");
	System.out.println(genomicAlignment2.getGenomicIntervalString());
      }
     }

    setGenomicAlignmentString ();

  }


  /***********************************************************************************
   *
   *  combineExons:
   *
   *  combineExons is called when the genomic alignments are normalized according to their
   *  alignmentString. It merges the exon sets of the genomic intervals.
   *
   ***********************************************************************************/
  
  public void combineExons (WeightObjectAlignment weightObjectAlignment) throws IOException {

    if (genomicAlignment1 != null && genomicAlignment2 != null &&
	weightObjectAlignment.getGenomicAlignment1 () != null && weightObjectAlignment.getGenomicAlignment2 () != null) {
      if (debugLevel >= 2) {
	System.out.println("genomicAlignment1: " + genomicAlignment1 + ", genomicAlignment2: " + genomicAlignment2);
	System.out.println("weightObjectAlignment - genomicAlignment1: " + weightObjectAlignment.getGenomicAlignment1 () +
			   ", genomicAlignment2: " + weightObjectAlignment.getGenomicAlignment2 ());
      }
      if (genomicAlignment1.equals(weightObjectAlignment.getGenomicAlignment1 ()) && genomicAlignment2.equals(weightObjectAlignment.getGenomicAlignment2 ())) {
	genomicAlignment1.combineExons(weightObjectAlignment.getGenomicAlignment1 ());
	genomicAlignment2.combineExons(weightObjectAlignment.getGenomicAlignment2 ());
	return;
      } else if (genomicAlignment1.equals(weightObjectAlignment.getGenomicAlignment2 ()) && genomicAlignment2.equals(weightObjectAlignment.getGenomicAlignment1 ())) {
	if (debugLevel >= 2) {
	  System.out.println("Combining: " + genomicAlignment1 + " with " + weightObjectAlignment.getGenomicAlignment2 () + " and " +
			     genomicAlignment2 + " with " + weightObjectAlignment.getGenomicAlignment1 ());
	}
	if (debugLevel >= 2) {
	  System.out.println("Combining: " + genomicAlignment1.getGenomicIntervals() + " with " + weightObjectAlignment.getGenomicAlignment2 ().getGenomicIntervals() + " and " +
			     genomicAlignment2.getGenomicIntervals() + " with " + weightObjectAlignment.getGenomicAlignment1 ().getGenomicIntervals());
	}
 	genomicAlignment1.combineExons(weightObjectAlignment.getGenomicAlignment2 ());
	genomicAlignment2.combineExons(weightObjectAlignment.getGenomicAlignment1 ());
	return;
      }
    }

    /* We merge the exon sets even if only the first and second read have the same alignment since
       the genomic alignment in these two cases is the same and this is used for normalization */
    if (genomicAlignment1 != null && weightObjectAlignment.getGenomicAlignment1 () != null  &&
	genomicAlignment1.equals(weightObjectAlignment.getGenomicAlignment1 ())) {
      genomicAlignment1.combineExons(weightObjectAlignment.getGenomicAlignment1 ());
      return;
    }

    if (genomicAlignment2 != null && weightObjectAlignment.getGenomicAlignment2 () != null  &&
	genomicAlignment2.equals(weightObjectAlignment.getGenomicAlignment2 ())) {
      genomicAlignment2.combineExons(weightObjectAlignment.getGenomicAlignment2 ());
      return;
    }

    if (genomicAlignment1 != null && weightObjectAlignment.getGenomicAlignment2 () != null &&
	genomicAlignment1.equals(weightObjectAlignment.getGenomicAlignment2 ())) {
      genomicAlignment1.combineExons(weightObjectAlignment.getGenomicAlignment2 ());
      return;
    }
    
    if (genomicAlignment2 != null && weightObjectAlignment.getGenomicAlignment1 () != null &&
	genomicAlignment2.equals(weightObjectAlignment.getGenomicAlignment1 ())) {
      genomicAlignment2.combineExons(weightObjectAlignment.getGenomicAlignment1 ());
      return;
    }

    throw new IOException ("Inconsistent weightObjectAlignments when combining the exons of " + this.getAlignmentId() +
			   " with: " + weightObjectAlignment.getAlignmentId() + "\n" +
			   (bedRecords1==null?"bedRecords1 null":bedRecords1) + "\n" +
			   (bedRecords2==null?"bedRecords2 null":bedRecords2) + "\n" +
			   genomicAlignment1 + "\n" +			   
			   genomicAlignment2 + "\n" +			   
			   weightObjectAlignment.getGenomicAlignment1 ()  + "\n" +
			   weightObjectAlignment.getGenomicAlignment2 ());
  }


  
  /***********************************************************************************
   *
   *                     computeOverlaps
   *
   * computeOverlaps computes for each countObjectId which has an overlap with the
   * weight object alignment the overlapping genomic alignments and the amount of
   * overlap with each of the genomic alignments.
   *
   ***********************************************************************************/

  public void computeOverlaps (Hashtable<String, Hashtable<GenomicAlignment, Integer>> overlapTable1,
			       Hashtable<String, Hashtable<GenomicAlignment, Integer>> overlapTable2) {


    Hashtable<String, Hashtable<GenomicAlignment, Integer>> alignmentOverlapTable1 = new Hashtable<String, Hashtable<GenomicAlignment, Integer>> ();
    Hashtable<String, Hashtable<GenomicAlignment, Integer>> alignmentOverlapTable2 = new Hashtable<String, Hashtable<GenomicAlignment, Integer>> ();

    debugLevel = UtilLib.getDebugLevel();
    if (ComputeCounts.getSpecialReadIdSet ().contains(getFragmentId())) {
      debugLevel = 2;
    }
    
    if (genomicAlignment1 == null && genomicAlignment2 == null) {
      genomicAlignmentString = getGenomicAlignmentString ();

      
      if (debugLevel >= 2) {
	System.out.println ("Genomic alignment string: " + genomicAlignmentString + " computed.");
      }
    }

    if (genomicAlignment1 != null) {
      if (debugLevel >= 1) {
	System.out.println ("Computing genomic overlaps for: " + genomicAlignment1);
      }
      genomicAlignment1.computeOverlaps (overlapTable1);
      numReads++;
    }
    if (genomicAlignment2 != null) {
      if (debugLevel >= 1) {
	System.out.println ("Computing genomic overlaps for: " + genomicAlignment2);
      }
      genomicAlignment2.computeOverlaps (overlapTable2);
      numReads++;
    }

    if (debugLevel >= 2) {
      System.out.println ("Genomic overlaps computed.");
    }

    debugLevel = UtilLib.getDebugLevel();

    
  }

  
  /***********************************************************************************
   *
   *                     getBedEntries
   *
   ***********************************************************************************/
  
  public Vector<BedEntry> getBedEntries () {

    if (genomicAlignment1 == null && genomicAlignment2 == null) {
      genomicAlignmentString = getGenomicAlignmentString ();
    }

    Vector<BedEntry> bedEntries = null;
    if (genomicAlignment1 != null) {
      bedEntries = genomicAlignment1.getBedEntries ();
    }

    if (genomicAlignment2 != null) {
      if (bedEntries == null) {
	bedEntries = genomicAlignment2.getBedEntries ();
      } else {
	bedEntries.addAll(genomicAlignment2.getBedEntries ());
      }
    }

    return (bedEntries);
    
  }


  /***********************************************************************************
   *
   *                         getCountObjectIds
   *
   ***********************************************************************************/

  public HashSet<String> getCountObjectIds (int overlapThreshold, Hashtable<String, Integer> countObjectLengthTable) throws IOException {

    HashSet<String> countObjectIds = new HashSet<String> ();
    
    for (BedRecord bedRecord: bedRecords1) {
      countObjectIds.addAll (bedRecord.getCountObjectIds (overlapThreshold, countObjectLengthTable));
    }

    for (BedRecord bedRecord: bedRecords2) {
      countObjectIds.addAll (bedRecord.getCountObjectIds (overlapThreshold, countObjectLengthTable));
    }

    return countObjectIds;
    
  }

  
  /***********************************************************************************
   *
   *                         printExonSet
   *
   ***********************************************************************************/


  public void printExonSet () {

    System.out.println("  Exon set of " + this);

    if (genomicAlignment1 != null) {
      genomicAlignment1.printExonSet ();
    }

    if (genomicAlignment2 != null) {
      genomicAlignment2.printExonSet ();
    }
  }

  
  /***********************************************************************************
   *
   *  isTranscriptExonAlignment
   *
   **********************************************************************************/

  public boolean isTranscriptExonAlignment () throws IOException {

    boolean isTranscriptExonAlignment1 = true;
    if (genomicAlignment1 != null) {
      isTranscriptExonAlignment1 = genomicAlignment1.isTranscriptExonAlignment ();
    }

    boolean isTranscriptExonAlignment2 = true;
    if (genomicAlignment2 != null) {
      isTranscriptExonAlignment2 = genomicAlignment2.isTranscriptExonAlignment ();
    }

    return isTranscriptExonAlignment1 && isTranscriptExonAlignment2;
    
  }

  /***********************************************************************************
   *
   *                     Basic methods
   *
   ***********************************************************************************/

  public String getFragmentId () {
    return fragmentId;
  }

  public String getAlignmentId () {
    return alignmentId;
  }

  public String getAlignmentBaseId () {
    return alignmentBaseId;
  }
  
  public int getNumReads () {
    return numReads;
  }

  public String getReferenceId () {
    return referenceId;
  }

  public GenomicAlignment getGenomicAlignment1 () {
    return genomicAlignment1;
  }


  public GenomicAlignment getGenomicAlignment2 () {
    return genomicAlignment2;
  }

  public int getLengthDifference () {
    return genomicAlignment1.getLengthDifference () + genomicAlignment2.getLengthDifference ();
  }


  public boolean isPairedEnd () {
    return pairedEnd;
  }

  public TreeSet<BedRecord> getBedRecords1 () {
    return bedRecords1;
  }


  public TreeSet<BedRecord> getBedRecords2 () {
    return bedRecords2;
  }

  public boolean fromGenomicBedRecords () {
    return fromGenomicBedRecords;
  }


  public boolean equals (Object o) {
    
    WeightObjectAlignment w = (WeightObjectAlignment) o;

    if (genomicAlignmentString == null) {
      genomicAlignmentString = getGenomicAlignmentString ();
    }
    
    if (getGenomicAlignmentString ().equals(w.getGenomicAlignmentString ())) {
      return true;
    }

    return false;
    
  }

  public int hashCode () {
    
    return getGenomicAlignmentString ().hashCode ();
    
  }


  public String toString () {
    return getAlignmentId() + " > " + getGenomicAlignmentString ();
  }


}
