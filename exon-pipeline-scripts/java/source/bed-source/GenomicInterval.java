/**File: GenomicAligment.java 

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
 *                    Class WeightObjectAlignmentComparator1
 *
 ***********************************************************************************/

class WeightObjectAlignmentComparator1 implements Comparator<WeightObjectAlignment> {

  @Override
  public int compare (WeightObjectAlignment weightObjectAlignment1, WeightObjectAlignment weightObjectAlignment2) {

    GenomicAlignment genomicAlignment1 = weightObjectAlignment1.getGenomicAlignment1 ();
    GenomicAlignment genomicAlignment2 = weightObjectAlignment2.getGenomicAlignment1 ();
    
    if (! genomicAlignment1.getChromosome ().equals(genomicAlignment2.getChromosome ())) {
      return genomicAlignment1.getChromosome ().compareTo (genomicAlignment2.getChromosome ());
    }

    if (genomicAlignment1.getLengthDifference () != genomicAlignment2.getLengthDifference ()) {
      (new Integer (genomicAlignment1.getLengthDifference ())).compareTo (new Integer (genomicAlignment2.getLengthDifference ()));
    }

    return (new Integer (genomicAlignment1.getChrStart ())).compareTo (new Integer (genomicAlignment2.getChrStart ()));
  }
}

/***********************************************************************************
 *
 *                    Class WeightObjectAlignmentComparator2
 *
 ***********************************************************************************/

class WeightObjectAlignmentComparator2 implements Comparator<WeightObjectAlignment> {

  @Override
  public int compare (WeightObjectAlignment weightObjectAlignment1, WeightObjectAlignment weightObjectAlignment2) {

    GenomicAlignment genomicAlignment1 = weightObjectAlignment1.getGenomicAlignment2 ();
    GenomicAlignment genomicAlignment2 = weightObjectAlignment2.getGenomicAlignment2();
    
    if (! genomicAlignment1.getChromosome ().equals(genomicAlignment2.getChromosome ())) {
      return genomicAlignment1.getChromosome ().compareTo (genomicAlignment2.getChromosome ());
    }

    if (genomicAlignment1.getLengthDifference () != genomicAlignment2.getLengthDifference ()) {
      (new Integer (genomicAlignment1.getLengthDifference ())).compareTo (new Integer (genomicAlignment2.getLengthDifference ()));
    }

    return (new Integer (genomicAlignment1.getChrStart ())).compareTo (new Integer (genomicAlignment2.getChrStart ()));
  }
}


/***********************************************************************************
 *
 *                              Class GenomicInterval
 *
 ***********************************************************************************/

public class GenomicInterval extends Exon implements Comparable {

  private static int debugLevel = UtilLib.getDebugLevel ();
  
  public static void setDebugLevel () {
    debugLevel = UtilLib.getDebugLevel ();
  }


  /***********************************************************************************
   *
   *                    Object variables and methods
   *
   ***********************************************************************************/

  private int overlap = 0;
  private HashSet<Exon> exonSet = null;
  
  private GenomicAlignment genomicAlignment = null;
  private BedRecord bedRecord = null;
  private int lengthDifference = 0;

  private int readAlignedLength    = 0;
  private int numInsertions = 0;
  private int numDeletions  = 0;

  private boolean chrStartEndAdjusted = false;

  private boolean isConformingAlignment = false;
  private int fieldIndex = -1;

  private int originalChromosomeStart = -1;
  private int originalChromosomeEnd = -1;


  GenomicInterval (String chromosome, int start, int end, String strand, int overlap, Exon exon, int lengthDifference,
		   int readAlignedLength, int numInsertions, int numDeletions, int fieldIndex , GenomicAlignment genomicAlignment,
		   BedRecord bedRecord, boolean isConformingAlignment) {
    
    this.chromosome              = chromosome;
    this.chromosomeStart         = start;
    this.chromosomeEnd           = end;
    this.originalChromosomeStart = start;
    this.originalChromosomeEnd   = end;
    this.strand                  = strand;

    this.overlap           = overlap;
    this.lengthDifference  = lengthDifference;

    this.readAlignedLength = readAlignedLength;
    this.numInsertions     = numInsertions;
    this.numDeletions      = numDeletions;

    this.fieldIndex = fieldIndex;

    this.genomicAlignment = genomicAlignment;
    this.bedRecord        = bedRecord;
 
    exonSet = new HashSet<Exon> ();

    this.isConformingAlignment = isConformingAlignment && (overlap >= Math.min(UtilLib.getOverlapThreshold (), exon.getLength ()));
    if (debugLevel >= 2) {
      System.out.println("this.isConformingAlignment = " + isConformingAlignment + " && (" + overlap + " >= Math.min(" + UtilLib.getOverlapThreshold () +
			 ", " + exon.getLength () + "))");
    }
    
    if (this.isConformingAlignment) {
      if (debugLevel >= 2) {
	System.out.println ("Adding exon: " + exon + " to " + exonSet);
      }
      exonSet.add(exon);
    }
    
  }


  /***********************************************************************************
   *
   *                           addToTableEntry  
   *
   ***********************************************************************************/

  public void addToTableEntry (Hashtable<String,  Hashtable<GenomicAlignment, Integer>> entryTable, String key, int value) {

    if (debugLevel >= 1) {
      System.out.println("Inserting " + key + "=" + value + " into " + entryTable);
      System.out.println("Retrieving key: " + key + " from entry table: " + entryTable);
    }
    
    Hashtable<GenomicAlignment, Integer> genomicAlignmentTable = entryTable.get(key);
    if (genomicAlignmentTable == null) {
      genomicAlignmentTable = new Hashtable<GenomicAlignment, Integer> ();
      entryTable.put (key, genomicAlignmentTable);

      if (debugLevel >= 1) {
	System.out.println("Inserting 0 into " + genomicAlignmentTable);
      }

      genomicAlignmentTable.put(genomicAlignment, new Integer (0));

      if (debugLevel >= 1) {
	System.out.println("New genomicAlignmentTable:" + genomicAlignmentTable);
      }
    }
    
    if (debugLevel >= 1) {
      System.out.println("Retrieving genomic alignment: " + genomicAlignment + " from genomic alignment table: " + genomicAlignmentTable);
    }

    Integer tableValueInt = genomicAlignmentTable.get(genomicAlignment);
    int tableValue = 0;
    if (tableValueInt != null) {
      tableValue = tableValueInt.intValue();
    }
	    
    tableValue = tableValue + value;
    genomicAlignmentTable.put(genomicAlignment, new Integer(tableValue));

    if (debugLevel >= 1) {
      System.out.println(tableValue + " inserted into " + genomicAlignmentTable);
    }    
    
  }

  /***********************************************************************************
   *
   *                    computeOverlaps
   *
   ***********************************************************************************/


  public void computeOverlaps (Hashtable<String, Hashtable<GenomicAlignment, Integer>> overlapTable) {


    debugLevel = UtilLib.getDebugLevel();

    HashSet<String> countObjectSet =  new HashSet<String> ();
    for (Exon exon: exonSet) {

     if (Exon.getCountObjectTable ().getSet (exon) == null) {
	
	if (UtilLib.warningsOn ()) {
	  System.err.println ("WARNING: No count object for " + exon.toString () + " found.");
	}

	if (debugLevel >= 1) {
	  System.out.println ("Adding exon in computeOverlaps: " + exon.toString ());
	}
	Exon.getCountObjectTable ().putValue (exon, exon.toString());
	
      }

      countObjectSet.addAll(Exon.getCountObjectTable ().getSet (exon));

      if (debugLevel >= 1) {
	System.out.println ("Exon: " + exon.toString () + " adding count objects for: " + countObjectSet);
      }
    }
      

    for (String countObjectId: countObjectSet) {

      if (debugLevel >= 1) {	      
	System.out.println ("GenomicInterval: " + this + " with exon set: " + exonSet);
	System.out.println ("genomicAlignment: " + genomicAlignment + ", weightObjectAlignment: " +
			    genomicAlignment.getWeightObjectAlignment () + ", countObject: " + countObjectId);
      }
	
      addToTableEntry (overlapTable, countObjectId, getOverlap ());

    }
      

    if (debugLevel >= 1) {
      System.out.println ("Overlap table after adding " + this + ": " + overlapTable);
    }
    
  }


  /***********************************************************************************
   *
   *                    adjustChrStartEnd
   *
   *  We define a partial ordering on the GenomicIntervals of one WeightObject which map
   *  to the same chromosome and have the same length and are comparable in the sense
   *  that the difference in start coordinates is not larger than the maximum of the
   *  length differences of both GenomicIntervals (where the length difference is
   *  defined as difference between the genomic and transcript exon length).
   *  The ordering is given by the ordering of the length difference and if the length
   *  difference is the same then the ordering is given by the start on the chromosome.
   *  In this way we have one canoncial element for each connected component of the
   *  partial order graph which is then used as the representative. (This is not totally
   *  stringent as it is easy to imagine cases in which the order of execution leads
   *  to the elimination or introduction of comparable pairs but in practice it is
   *  good enough.)
   *
   ***********************************************************************************/


  public boolean adjustChrStartEnd (GenomicInterval genomicInterval) throws IOException {

    if (genomicInterval == null) {
      return false;
    }

    if (! chromosome.equals(genomicInterval.getChromosome ())){
      return false;
    }

    /* If both alignments are to the same reference (i.e. junction or transcript), the two intervals need not be adjusted according to each other */
    if (getReferenceId().equals(genomicInterval.getReferenceId ())) {
      return false;
    }

    if (chromosomeStart == genomicInterval.getChrStart() && chromosomeEnd == genomicInterval.getChrEnd()) {
      return false;
    }

    if (this.lengthDifference < genomicInterval.getLengthDifference() ||
	(this.lengthDifference == genomicInterval.getLengthDifference() && chromosomeStart < genomicInterval.getChrStart())) {
      if (debugLevel >= 3) {
	System.out.println ("Changing from adjusting " + this + " to adjusting " + genomicInterval);
      }
      
      if (genomicInterval.adjustChrStartEnd (this)) {
	genomicInterval.getGenomicAlignment().setGenomicAlignmentString ();
	genomicInterval.getGenomicAlignment().getWeightObjectAlignment().setGenomicAlignmentString ();
	return true;
      } else {
	return false;
      }
    }

    if (debugLevel >= 3) {
      System.out.println ("Checking genomic interval: " + this + " against: " + genomicInterval);
      System.out.println (chromosome + " == " + genomicInterval.getChromosome () + " && Math.abs(" + chromosomeEnd + " - " +
			  chromosomeStart + ") = " + Math.abs(chromosomeEnd - chromosomeStart) + " == Math.abs(" + genomicInterval.getChrEnd() + " - " +
			  genomicInterval.getChrStart() + ") = " + Math.abs(genomicInterval.getChrEnd() - genomicInterval.getChrStart()) + " && Math.abs(" +
			  genomicInterval.getChrStart() + " - " + chromosomeStart + ") == " + Math.abs(genomicInterval.getChrStart() - chromosomeStart) +
			  " <= " + lengthDifference);

      System.out.println ("And: (" + this.lengthDifference + " > " + genomicInterval.getLengthDifference() + ") || (" + this.lengthDifference + " == " +
			  genomicInterval.getLengthDifference() + " && " + chromosomeStart + " > " + genomicInterval.getChrStart() + ")");
    }

    if (chromosome.equals(genomicInterval.getChromosome ()) &&
	Math.abs(chromosomeEnd - chromosomeStart) == Math.abs(genomicInterval.getChrEnd() - genomicInterval.getChrStart()) &&
	Math.abs(genomicInterval.getOriginalChromosomeStart() - originalChromosomeStart) <= lengthDifference) {
      if (this.lengthDifference > genomicInterval.getLengthDifference() ||
	  (this.lengthDifference == genomicInterval.getLengthDifference() && chromosomeStart > genomicInterval.getChrStart())) {

	if (debugLevel >= 3) {
	  System.out.print ("Adjusting genomic interval: " + this + " to " + genomicInterval + " (differences: " + this.lengthDifference + " vs. " +
			    genomicInterval.getLengthDifference() + ")");
	}
		
	chromosomeStart  = genomicInterval.getChrStart();
	chromosomeEnd    = genomicInterval.getChrEnd();
	lengthDifference = genomicInterval.getLengthDifference();
	// lengthDifference = lengthDifference - Math.abs(chromosomeStart - genomicInterval.getChrStart());

	if (debugLevel >= 2) {
	  System.out.println (" to: " + this);
	}

	return true;
	
      }
     
    }

    return false;
    
  }


  /***********************************************************************************
   *
   * adjustChrStartEnd
   *
   ***********************************************************************************/
  
  public void adjustChrStartEnd (Vector<WeightObjectAlignment> weightObjectAlignments, int genomicAlignmentIndex) throws IOException {
    
    if (weightObjectAlignments == null || lengthDifference == 0 || chrStartEndAdjusted) {      
      return;
    }

    chrStartEndAdjusted = true;

    if (genomicAlignmentIndex == 1) {
      // Collections.sort(weightObjectAlignments, new WeightObjectAlignmentComparator1 ());
    }

    if (genomicAlignmentIndex == 2) {
      // Collections.sort(weightObjectAlignments, new WeightObjectAlignmentComparator2 ());
    }

    /* Note that we always have to go through all weightObjectAlignments since we have to
       make sure that we change the representative element for all weightObjectAlignments
       if necessary - it could be the last element considered, for instance. */
    for (WeightObjectAlignment weightObjectAlignment: weightObjectAlignments) {

      if (genomicAlignmentIndex == 1) {
	GenomicAlignment genomicAlignment1 = weightObjectAlignment.getGenomicAlignment1 ();
	if (genomicAlignment1 != null && ! genomicAlignment.equals(genomicAlignment1)) {
	  for (GenomicInterval genomicInterval: genomicAlignment1.getGenomicIntervals()) {
	    String oldGenomicIntervalString = genomicInterval.toString();
	    adjustChrStartEnd (genomicInterval);
	    if (debugLevel >= 2) {
	      if (! genomicInterval.toString ().equals(oldGenomicIntervalString)) {
		System.out.println ("Readjusted genomic interval " + oldGenomicIntervalString + " of alignment 1 of " + weightObjectAlignment.getAlignmentId() +
				    " to " + genomicInterval);
	      }
	    }
	  }
	}
      }
      
      if (genomicAlignmentIndex == 2) {
	GenomicAlignment genomicAlignment2 = weightObjectAlignment.getGenomicAlignment2 ();
	if (genomicAlignment2 != null && ! genomicAlignment.equals(genomicAlignment2)) {
	  for (GenomicInterval genomicInterval: genomicAlignment2.getGenomicIntervals()) {
	    String oldGenomicIntervalString = genomicInterval.toString();
	    adjustChrStartEnd (genomicInterval);
	    if (debugLevel >= 2) {
	      if (! genomicInterval.toString ().equals(oldGenomicIntervalString)) {
		System.out.println ("Readjusted genomic interval " + oldGenomicIntervalString + " of alignment 2 of " + weightObjectAlignment.getAlignmentId() +
				    " to " + genomicInterval);
	      }
	    }
	  }
	}
      }
    }
    
  }



  /***********************************************************************************
   *
   *                    Exon set methods
   *
   ***********************************************************************************/

  public void addExonSet (HashSet<Exon> exonSet) {

    if (debugLevel >= 2) {
      System.out.println ("Adding exon set: " + exonSet + " to this.exonSet: " + this.exonSet);
    }
    
    this.exonSet.addAll (exonSet);

    if (debugLevel >= 2) {
      System.out.println ("New exon set: " + this.exonSet);
    }
    
  }


  public HashSet<Exon> getExonSet () {
    return exonSet;
  }


  /***********************************************************************************
   *
   *                    getOverlap: compute the overlap with an exon
   *
   *  Note that GenomicInterval is created in BedRecord, however, the start and end
   *  of the genomic interval are determined from the genomic exon coordinates so
   *  that all of the genomic coordinates are one based (i.e. they are not BED
   *  coordinates).
   *
   ***********************************************************************************/
  
  public int getOverlap (Exon exon) {
    return Math.max(0, Math.min(exon.getChrEnd(), getChrEnd()) - Math.max(getChrStart (), getChrStart ()) + 1);
  }
  
  public int getOverlap () {
    return overlap;
  }

  /***********************************************************************************
   *
   *  removeExons
   *  
   *  Remove exons with less overlap than the maximum overlap of an exon for a
   *  given genomic alignment
   *
   ***********************************************************************************/
  
  public HashSet<Exon> updateExonSet () {

    int maxOverlap = 0;

    if (exonSet.size () > 1) {
      for (Exon exon: exonSet) {
	if (! exon.isGenomicExon ()) {
	  if (getOverlap (exon) > maxOverlap) {
	    maxOverlap = getOverlap (exon);
	  }
	}
      }

      if (maxOverlap > 0) {
	overlap = Math.max(overlap, maxOverlap);
	
	HashSet<Exon> exonDeletionSet = new HashSet <Exon> ();
	for (Exon exon: exonSet) {
	  if (! exon.isGenomicExon ()) {
	    if (getOverlap (exon) < maxOverlap) {
	      if (UtilLib.warningsOn ()) {
		System.err.println("Exon with too little overlap: " + exon + " with overlap " + getOverlap (exon) + " vs " + maxOverlap +
				   " for " + this + " of alignment " + genomicAlignment + " - removing.");
		System.err.println("BED records:");
		for (BedRecord bedRecord: genomicAlignment.getBedRecords ()) {
		  System.err.println (bedRecord.toCompleteString());
		}
		System.err.println("Other exons:");
		for (Exon newExon: exonSet) {
		  System.err.println (newExon);
		}
	      }
	      exonDeletionSet.add (exon);
	    }
	  }
	}
	
	exonSet.removeAll(exonDeletionSet);
	
      }
    }

    return exonSet;
    
  }
  
  
  /***********************************************************************************
   *
   *                         printExonSet
   *
   ***********************************************************************************/

  public void printExonSet () {

    System.out.println("      Exon set of " + this + ": " + exonSet);

  }
  
  /***********************************************************************************
   *
   *                    Basic methods
   *
   ***********************************************************************************/


  public int getOriginalChromosomeStart () {
    return originalChromosomeStart;
  }

  public int getOriginalChromosomeEnd () {
    return originalChromosomeEnd;
  }

  public int getReadAlignedLength () {
    return readAlignedLength;
  }

  public int getLengthDifference () {
    return lengthDifference;
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

  public String getReferenceId () {
    return genomicAlignment.getReferenceId();
  }


  public GenomicAlignment getGenomicAlignment () {
    return genomicAlignment;
  }

  public BedRecord getBedRecord () {
    return bedRecord;
  }

  public boolean isTranscriptExon() {
    return UtilLib.isTranscriptExon(chromosome);
  }

  
  /***********************************************************************************
   *
   *                    compareTo
   *
   ***********************************************************************************/
  
  public int compareTo (Object o) {

    GenomicInterval g = (GenomicInterval) o;

    if (UtilLib.isTranscriptExon(chromosome) || UtilLib.isTranscriptExon(g.getChromosome ())) {
      BedRecord bedRecord1 = getBedRecord ();
      BedRecord bedRecord2 = g.getBedRecord ();
      if (bedRecord1.getReferenceSequenceId ().equals(bedRecord2.getReferenceSequenceId ())) {
	Exon exon1 = bedRecord1.getExon ();
	Exon exon2 = bedRecord2.getExon ();
	if (exon1.getStrand().equals (exon2.strand)) {
	  int compValue = (exon1.getStrand().equals("+")?1:-1) * UtilLib.compareInt (bedRecord1.getExonReferenceStart (), bedRecord2.getExonReferenceStart ());
	  if (compValue == 0 && ! this.equals(g)) {
	    System.err.println("compareTo (" + this + ", " + g + ") = 0 but no equality");
	    System.err.println("BED record 1: " + bedRecord1 + ", " + bedRecord1.getExonReferenceStart () + ", exon: " + exon1);
	    System.err.println("BED record 2: " + bedRecord2 + ", " + bedRecord2.getExonReferenceStart () + ", exon: " + exon2);
	    System.exit (1);
	  }								    
	  return compValue;
	} else {
	  System.err.println ("WARNING: Comparison of " + this + " with " + g + " where the exons of the BED records have different strands.\n" +
			      bedRecord1 + "\n" + bedRecord2);
	}
      } else {
	System.err.println ("WARNING: Comparison of " + this + " with " + g + " where BED records have different reference sequence ids.\n" +
			    bedRecord1 + "\n" + bedRecord2);
      }
    }

    if (! chromosome.equals(g.getChromosome ())) {
      System.err.println ("WARNING: Comparison of " + this + " with " + g + " where the exons of the BED records have chromosomes.\n");
      return chromosome.compareTo(g.getChromosome ());
    }

    if (chromosomeStart != g.getChrStart ()) {
      if (chromosomeStart > g.getChrStart ()) {
	return 1;
      } else {
	return -1;
      }
    }

    if (chromosomeEnd != g.getChrEnd ()) {
      if (chromosomeEnd > g.getChrEnd ()) {
	return 1;
      } else {
	return -1;
      }
    }
    
    return getStrand ().compareTo(g.getStrand ());
    
  }
 

  public BedEntry toBedEntry () {
    
    return new BedEntry (chromosome, chromosomeStart, chromosomeEnd, "", overlap, "+");
    
  }


  /***********************************************************************************
   *
   * toGenomicBedEntryString ():
   *
   ***********************************************************************************/
  
  public String toGenomicBedEntryString (String bedEntryId) {

    String genomicBedEntryString = "";
    for (Exon exon: exonSet) {
      if (! genomicBedEntryString.equals("")) {
	genomicBedEntryString = genomicBedEntryString + "\n";
      }
      
      genomicBedEntryString = genomicBedEntryString + toBedString("\t") + "\t" + bedEntryId + "\t" + "255" + "\t" + strand + "\t" + 
	exon.toBedString("\t") + "\t" + exon.toString() + "/" + exon.getStrand () + "\t" + "0" + "\t" + "+" + "\t" + getOverlap(exon);
      
    }

    return genomicBedEntryString;

  }


  /***********************************************************************************
   *
   * toString:
   *
   * If the exon does not align perfectly against the transcript there is some
   * uncertainty in the genomic location of a read (due to the slack in the mapping
   * from the transcript to the genome.) In such a case we can try to identify other
   * genomic mappings to which our mapping can be converted (given the slack).
   *
   ***********************************************************************************/
  
  public String toString (Vector<WeightObjectAlignment> weightObjectAlignments, int genomicAlignmentIndex) throws IOException {

    adjustChrStartEnd (weightObjectAlignments, genomicAlignmentIndex);

    return toString ();
  }

  
  public String toCompleteString () {
    
    return super.toString();

  }


  public String toString () {
    
    if (UtilLib.isTranscriptExon(getChromosome ()) && ! UtilLib.getGeneCountMode ()) {
      return "";
    }

    if (debugLevel >= 3) {
      System.out.println ("toString " + super.toString() + " - exon set " + exonSet);
    }

    return super.toString();

  }

}


