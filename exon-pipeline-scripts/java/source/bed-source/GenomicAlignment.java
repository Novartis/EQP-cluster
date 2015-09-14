/**File: GenomicAlignment.java 

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
 *                              Class GenomicAlignment
 *
 ***********************************************************************************/

public class GenomicAlignment implements Comparable {

  private static int debugLevel = UtilLib.getDebugLevel ();


  /***********************************************************************************
   *
   *                     Object variables and methods
   *
   ***********************************************************************************/

  private TreeSet<GenomicInterval> genomicIntervals = null;
  private Hashtable<GenomicInterval, GenomicInterval> genomicIntervalTable = null;
  private int start = Integer.MAX_VALUE;
  private WeightObjectAlignment weightObjectAlignment = null;
  private String alignmentString = null;
  private String strand = null;
  private String referenceGenomicAlignmentOrientation = null;

  private int maxFieldIndex = -1;
  private int correctionDifference = 0;

  private int readAlignedLength = -1;
  private int readAlignStart  = -1;
  private boolean genomicBedRecords = false;

  private int numInsertions = -1;
  private int numDeletions  = -1;

  private String chromosome = "";
  private int lengthDifference = 0;

  private boolean isTranscriptExonAlignmentIsSet = false;
  private boolean isTranscriptExonAlignment      = false;
  
  GenomicAlignment () {
    debugLevel = UtilLib.getDebugLevel ();
    genomicIntervals     = new TreeSet<GenomicInterval> ();
    genomicIntervalTable = new Hashtable <GenomicInterval, GenomicInterval> ();
  }


  /***********************************************************************************
   *
   *                              computeGenomicAlignments
   *
   * Note that the BedRecords are traversed in sorted order. If there are overlapping
   * exons, then this could cause problems with counting consecutive exons, for instance,
   * along one junction. This is not a problem for spliced genomic records but for
   * unspliced genomic records this can be a problem if we want to compute exon-intron
   * junction counts. In order to alleviate this we give the same field index to
   * overlapping exons until we encounter an exon that does not overlap with the first
   * overlapping exon.
   *
   * If we allow non-conforming alignments (for the computation of intron-exon junctions),
   * then we treat each contiguous part of a spliced alignment separately since only
   * contiguous parts can overlap such a junctions meaningfully.
   *
   ***********************************************************************************/

  GenomicAlignment (TreeSet<BedRecord> bedRecords, WeightObjectAlignment weightObjectAlignment) throws IOException {

    this ();

    this.weightObjectAlignment = weightObjectAlignment;

    int sumExonOverlap = 0;
    int sumStartDiff   = 0;
    int sumEndDiff     = 0;
    int readCovered    = 0;

    BedRecord oldBedRecord = null;
    String    oldReadAlignmentId = "";
    Exon      oldExon = null;
  
    int oldExonTranscriptEnd = -1;
    int oldEndDiff           = -1;
    
    int sumReadAlignedLength = 0;
    int sumInsertions        = 0;
    int sumDeletions         = 0;
    int oldOverlap           = 0;

    GenomicInterval genomicInterval    = null;
    GenomicInterval oldGenomicInterval = null;
    
    int i = -1;
    int numBedRecords = bedRecords.size ();
      
    for (BedRecord bedRecord: bedRecords) {

      if (strand == null) {
	strand = bedRecord.getAlignStrand ();
      } else if (! strand.equals (bedRecord.getAlignStrand ())) {
	throw new IOException ("BED record: " + bedRecord + " has a different strand than the old BED record: " + oldBedRecord);
      }

      if (referenceGenomicAlignmentOrientation == null) {
	referenceGenomicAlignmentOrientation = bedRecord.getReferenceGenomicAlignmentOrientation ();
      } else if (! referenceGenomicAlignmentOrientation.equals (bedRecord.getReferenceGenomicAlignmentOrientation ())) {
	referenceGenomicAlignmentOrientation = "*";
	if (UtilLib.warningsOn ()) {
	  System.err.println ("BED record: " + bedRecord + " has a different referenceGenomicAlignmentOrientation than the old BED record: " + oldBedRecord);
	}
      }
      
      if (debugLevel >= 2) {
	System.err.println("Bed record: " + bedRecord);
      }

      if (! oldReadAlignmentId.equals("") && (genomicBedRecords != bedRecord.isGenomicRecord () ||
					      (oldBedRecord.getFieldIndex () == bedRecord.getFieldIndex () &&
					       (numInsertions != bedRecord.getNumInsertions () || numDeletions != bedRecord.getNumDeletions ())))) {
        throw new IOException ("Change of BED record type from " + oldBedRecord + " to " + bedRecord + "\n" +
			       genomicBedRecords + " != " + bedRecord.isGenomicRecord () + " or " + numInsertions + " != " + bedRecord.getNumInsertions () +
			       " or " + numDeletions + " != " + bedRecord.getNumDeletions ());
      }

      genomicBedRecords = bedRecord.isGenomicRecord ();
      numDeletions      = bedRecord.getNumDeletions ();
      numInsertions     = bedRecord.getNumInsertions ();

      String readAlignmentId = bedRecord.getReadAlignmentId();
      Exon   exon            = bedRecord.getExon();
      int    overlap         = bedRecord.getOverlap ();

      if (readAlignStart == -1 || bedRecord.getReadAlignStart () < readAlignStart) {
	readAlignStart = bedRecord.getReadAlignStart ();
      }
      	
      if (! oldReadAlignmentId.equals("") && ! readAlignmentId.equals(oldReadAlignmentId)) {
	throw new IOException ("Different read alignment ids: " + readAlignmentId + " vs " +  oldReadAlignmentId);
      }

      if (bedRecord.isGenomicRecord ()) {

	if (oldBedRecord != null && oldBedRecord.getMaxFieldIndex() != bedRecord.getMaxFieldIndex ()) {
	  throw new IOException ("Differing max field index: " + oldBedRecord.getMaxFieldIndex() + " vs. " +  bedRecord.getMaxFieldIndex () +
				 " for bed records: \n" + oldBedRecord + "\n" + bedRecord);	  
	}	

	maxFieldIndex = bedRecord.getMaxFieldIndex ();

	if (debugLevel >= 2) {
	  System.err.println( "Before computeGenomicInterval.");
	}


	genomicInterval = bedRecord.computeGenomicInterval (this);

	if (debugLevel >= 2) {
	  System.err.println( "Genomic interval computed.");
	  System.err.println (genomicIntervals);
	  System.err.println( "Before if");
	}


	if (! genomicIntervals.contains (genomicInterval)) {
	  if (debugLevel >= 2) {
	    System.err.println( "Adding " + genomicInterval + " to genomicIntervals: " + genomicIntervals);
	  }
	  genomicIntervals.add (genomicInterval);
	  genomicIntervalTable.put (genomicInterval, genomicInterval);
	} else {
	  if (debugLevel >= 2) {
	    System.err.println ("Else branch.");
	    System.err.println (genomicIntervalTable);
	  }
	  GenomicInterval genomicTableInterval = genomicIntervalTable.get (genomicInterval);
	  if (debugLevel >= 2) {
	    System.err.println ("genomicTableInterval: " + genomicTableInterval);
	  }
	  if (debugLevel >= 2) {
	    System.err.println( "Adding exon set " + genomicInterval.getExonSet () + " to genomicTableInterval: " + genomicTableInterval + " with exon set: " +
				genomicTableInterval.getExonSet ());
	  }
	  genomicTableInterval.addExonSet (genomicInterval.getExonSet ());
	}

	if (oldBedRecord != null && oldBedRecord.getFieldIndex () == bedRecord.getFieldIndex () && ! genomicInterval.genomicCoordinatesEquals(oldGenomicInterval)) {
	  throw new IOException ("Differing genomic intervals for " + oldBedRecord + " and " + bedRecord + "\n" +
				 "Genomic interval 1: " + genomicInterval + ", genomic interval 2: " + oldGenomicInterval);	  
	}

	if (bedRecord.getReadAlignedLength () > 0) {
	  sumReadAlignedLength = sumReadAlignedLength + bedRecord.getReadAlignedLength ();
	}

      } else {

	if (i == -1) {
	  maxFieldIndex = bedRecords.size () - 1;
	  bedRecord.setMaxFieldIndex (maxFieldIndex);
	}
	
	/* Non-genomic reference sequences are transcripts and junctions; these do not have overlapping exons. So we do not need
	   to take overlapping exons into account and can just increase the index for each genomic interval. The genomic intervals
	   are ordered along the transcript sequence - with one caveat: Note that we want to count *consecutive* exon aligments to
	   quantify junctions; for this special case we need to disregard BED records that do not represent a genomic alignment, that
	   is, alignments to tr-exons. For quantifying exons it does not influence the result and for genes we definitely need to
	   make sure that all alignments are captured as we want to ensure that the complete read aligns to the gene. */
	if (bedRecord.isTranscriptExonAlignment () && UtilLib.getCountMode ().equals("junction")) {
	  bedRecord.setFieldIndex (-1);
	} else {
	  i++;
	  bedRecord.setFieldIndex (i);
	}
	
	if (debugLevel >= 2) {
	  System.err.println ("Max field index: " + bedRecord.getMaxFieldIndex () + " for BED record: " + bedRecord);
	}


	if (sumReadAlignedLength == 0) {
	  sumReadAlignedLength = bedRecord.getReadAlignedLength ();
	} else if (sumReadAlignedLength != bedRecord.getReadAlignedLength ()) {
	  throw new IOException ("Aligned length: " + bedRecord.getReadAlignedLength () + " of BED record:\n" + bedRecord + "\n" + 
				 "differs from previous aligned length: " + sumReadAlignedLength);
	}

	if (sumInsertions == 0) {
	  sumInsertions = bedRecord.getNumInsertions ();
	}  else if (sumInsertions != bedRecord.getNumInsertions ()) {
	  throw new IOException ("Number of insertions: " + bedRecord.getNumInsertions () + " of BED record:\n" +
				 bedRecord.toCompleteString() + "\n" + 
				 "differs from previous insertions: " + sumInsertions);
	}

	if (sumInsertions < 0) {
	  throw new IOException ("ERROR: negative number of insertions for BED record: " + bedRecord.toCompleteString ());
	}


	if (sumDeletions == 0) {
	  sumDeletions = bedRecord.getNumDeletions ();
	} else if (sumDeletions != bedRecord.getNumDeletions ()) {
	  throw new IOException ("Number of deletions: " + bedRecord.getNumDeletions () + " of BED record:\n" +
				 bedRecord.toCompleteString() + "\n" + 
				 "differs from previous deletions: " + sumDeletions);
	}
	
	if (sumDeletions < 0) {
	  throw new IOException ("ERROR: negative number of deletions for BED record: " + bedRecord.toCompleteString ());
	}

	if (debugLevel >= 2) {
	  System.err.println ("old readCovered: " + readCovered + ", overlap: " + overlap);
	}

	readCovered = readCovered + overlap;
  
	/***************************************************************************
	 *
	 *  Assertion: the same readAlignmentId and a new bed entry implies the
	 *   the same alignment on the transcript and a new exon
	 *
	 ***************************************************************************/

	if (oldBedRecord != null) {
    
	  if (! bedRecord.alignmentEquals(oldBedRecord)) {
	    throw new IOException ("BedRecord: " + bedRecord + "\n" + "Read id: " + readAlignmentId + " -  changed alignment: " +  bedRecord.getAlignment () +
				   " vs " + oldBedRecord.getAlignment ());
	  }

	  if (exon.equals(oldExon)) {
	    throw new IOException ("old BedRecord: " + oldBedRecord + "\n" + "BedRecord: " + bedRecord + "\n" + "Unchanged exon: " + exon + " vs " + oldBedRecord.getExon());
	  }
	}

	/***************************************************************************
	 *
	 *  Some times there are overlaps between consecutive exons; we compute
	 *  the difference between the end of the old exon and the beginning of
	 *  the new exon (as we order according to readAlignStart and End as well as
	 *  exonTranscriptStart and End this is the correct orientation).
	 *
	 ***************************************************************************/
    
	int exonTranscriptStart = bedRecord.getExonReferenceStart ();      
	int exonOverlap = Math.max(0, oldExonTranscriptEnd - exonTranscriptStart + 1);
	if (oldExonTranscriptEnd != -1) {
	  sumExonOverlap  = sumExonOverlap + exonOverlap;
	}

	if ((debugLevel >= 2) && sumExonOverlap > 0) {
	  System.err.println("sum exon overlap: " + sumExonOverlap + ", bedRecord: " + bedRecord.toString ());
	  System.err.println("oldExonTranscriptEnd: " + oldExonTranscriptEnd + "exonTranscriptStart: " + exonTranscriptStart);
	}

	genomicInterval = bedRecord.computeGenomicInterval (i, numBedRecords, this);

	String genomicIntervalChromosome = genomicInterval.getChromosome ();
	if (! chromosome.substring(Math.max(0, chromosome.length() - genomicIntervalChromosome.length())).equals (genomicIntervalChromosome)) {
	  chromosome = chromosome + ":" + genomicIntervalChromosome;
	}

	lengthDifference = lengthDifference + genomicInterval.getLengthDifference ();
		
	if (debugLevel >= 2) {
	  System.err.println("Adding genomic interval " + genomicInterval + " to genomic alignment");
	  System.err.println(genomicIntervals + ".contains (" + genomicInterval + "): " + genomicIntervals.contains (genomicInterval));
	  if (genomicIntervals.contains (genomicInterval)) {
	    Iterator it = genomicIntervals.iterator ();
	    while (it.hasNext ()) {
	      GenomicInterval g = (GenomicInterval) it.next();
	      System.err.println(genomicInterval + ".equals (" + g + "): " + genomicIntervals.equals (g));
	    }
	  }
	}

	if (! genomicIntervals.contains (genomicInterval)) {
	  if (debugLevel >= 2) {
	    System.err.println("Genomic interval " + genomicInterval + " not contained in: " + genomicIntervals);
	  }
	  genomicIntervals.add (genomicInterval);
	  genomicIntervalTable.put (genomicInterval, genomicInterval);
	} else {
	  if (debugLevel >= 2) {
	    System.err.println("Genomic interval " + genomicInterval + " contained in: " + genomicIntervals);
	  }
	  GenomicInterval genomicTableInterval = genomicIntervalTable.get (genomicInterval);
	  genomicTableInterval.addExonSet (genomicInterval.getExonSet ());
	}
	 
	
	/********************************************************************************
	 **
	 ** Usually consecutive intersection intervals will have endDiff = 0 for i - 1 and
	 ** startDiff = 0 for i. However, if the start of the alignment falls into an
	 ** overlap region between two exons:
	 ** 
	 ** read:            |-------------------|
	 ** exon1: ..----------|
	 ** exon2:       |---------------------|
	 **
	 ** then startDiff may be larger than 0 for i (the overlap with exon 2). Similarly,
	 ** if the end falls into an overlap:
	 **
	 ** read:       |-------------------|
	 ** exon1:        ..--------------------|
	 ** exon2:                     |---------------------|
	 **
	 ** then endDiff may be larger than 0 for i - 1 (the overlap with exon 1).
	 **
	 ********************************************************************************/
	
	int endDiff = bedRecord.getEndDiff();

	if (debugLevel >= 2) {
	  System.err.println ("endDiff for " + bedRecord.toString () + ": " + endDiff);
	  System.err.flush();
	}
	
	if (numBedRecords > 1) {

	  int startDiff = bedRecord.getStartDiff();
	  if (i > 0) {
	    sumStartDiff  = sumStartDiff + startDiff;
	  }

	  if (i < numBedRecords - 1) {
	    sumEndDiff = sumEndDiff + endDiff;
	  }

	  if (! bedRecord.isGenomicRecord () && UtilLib.warningsOn ()) {
	    if ((i == 1 && startDiff > exonOverlap) || (i > 1 && startDiff > 0)) {
	      System.err.println("WARNING: at interval " + i + " of read " + readAlignmentId + " and startDiff: " + startDiff);
	    }
	    
	    if (((i == numBedRecords - 1) && oldEndDiff > exonOverlap) || (i < numBedRecords - 1 && oldEndDiff > 0)) {
	      System.err.println("WARNING: at interval " + i + " of read " + readAlignmentId + " and endDiff: " + endDiff);
	    }
	  }
	}

	/* Save some values for the next iteration */
	oldExon              = exon;
	oldExonTranscriptEnd = bedRecord.getExonReferenceEnd ();
	oldEndDiff           = endDiff;

      }

      if (genomicInterval.getChrStart () < start) {
	start = genomicInterval.getChrStart ();
      }

      /* Save some values for the next iteration */
      oldBedRecord         = bedRecord;
      oldReadAlignmentId   = readAlignmentId;
      oldGenomicInterval   = genomicInterval;
          
    }

    if (debugLevel >= 3) {
      System.err.println ("genomicIntervals: " + genomicIntervals);
      System.err.println ("genomicIntervalTable: " + genomicIntervalTable);
    }

    if (! genomicBedRecords) {
      if (UtilLib.warningsOn () && sumReadAlignedLength > 0) {
	if (readCovered + sumStartDiff + sumEndDiff < sumReadAlignedLength - sumInsertions + sumDeletions + sumExonOverlap) {
	  System.err.println("Read " + oldReadAlignmentId + " covers too little: " + readCovered + " + " + sumStartDiff + " + " + sumEndDiff +
			     " vs " + sumReadAlignedLength + " - " + sumInsertions + " + " + sumDeletions + " + " + sumExonOverlap);
	  for (BedRecord bedRecord: bedRecords) {
	    System.err.println (bedRecord.toCompleteString());
	  }
	}
      
	if (readCovered + sumStartDiff + sumEndDiff > sumReadAlignedLength - sumInsertions + sumDeletions + sumExonOverlap) {
	  System.err.println("Read " + oldReadAlignmentId + " covers too much: " + readCovered +  " + " + sumStartDiff + " + " + sumEndDiff +
			     " vs " + sumReadAlignedLength + " - " + sumInsertions + " + " + sumDeletions + " + " + sumExonOverlap);
	  for (BedRecord bedRecord: bedRecords) {
	    System.err.println (bedRecord.toCompleteString());
	  }
	}
      }

      if (readCovered < sumReadAlignedLength) {
	correctionDifference = sumReadAlignedLength - readCovered;
      }
    }   
  }


  /***********************************************************************************
   *
   *                              addToTableEntry
   *
   ***********************************************************************************/

  
  public void addToTableEntry (Hashtable<String, Hashtable<GenomicAlignment, Integer>> overlapTable, String countObjectId, int value) {

    Hashtable<GenomicAlignment, Integer> overlapTableGenomicAlignment = overlapTable.get (countObjectId);

    if (debugLevel >= 2) {
      System.err.println("overlapTableGenomicAlignment for " + countObjectId + ": " + overlapTableGenomicAlignment);
    }

    if (overlapTableGenomicAlignment == null) {
      return;
    }

    Integer overlap = overlapTableGenomicAlignment.get(this);

    if (debugLevel >= 2) {
      System.err.println("overlap for " + this + ": " + overlap);
    }


    if (overlap == null) {
      return;
    }

    overlapTableGenomicAlignment.put(this, new Integer (overlap.intValue () + value));
    
  }


  
  
  /***********************************************************************************
   *
   *                              computeOverlaps
   *
   ***********************************************************************************/


  public void computeOverlaps (Hashtable<String, Hashtable<GenomicAlignment, Integer>> overlapTable) {

    debugLevel = UtilLib.getDebugLevel();
    
    int oldDebugLevel = debugLevel;

    if (ComputeCounts.getSpecialReadIdSet().contains(weightObjectAlignment.getFragmentId())) {
      debugLevel = 1;
      UtilLib.setDebugLevel(1);
    }
    
    if (debugLevel >= 1) {
      System.err.println ("genomicIntervals: " + genomicIntervals);
    }


    for (GenomicInterval genomicInterval: genomicIntervals) {
      if (debugLevel >= 2) {
	System.err.println("Computing overlaps for: " + genomicInterval);
      }
      genomicInterval.computeOverlaps (overlapTable);
    }

    if (numInsertions - numDeletions != 0) {
      for (String countObjectId: overlapTable.keySet()) {
	if (debugLevel >= 2) {
	  System.err.println("Adding " + numInsertions + " - " + numDeletions + " for " + countObjectId + " and " + this + " to " + overlapTable);
	}
	addToTableEntry(overlapTable, countObjectId, numInsertions - numDeletions);

	if (debugLevel >= 2) {
	  System.err.println("New overlap table: " + overlapTable);
	}

      }
    }

    if (debugLevel >= 2) {
      System.err.println ("Done.");
    }

    UtilLib.setDebugLevel(oldDebugLevel);
    debugLevel = UtilLib.getDebugLevel();
    GenomicInterval.setDebugLevel ();
    
  }


  /***********************************************************************************
   *
   *                    getConformingCountObjects
   *
   *  Computes the countObjects that have conforming alignments with the GenomicAlignment.
   *
   *  In order to do this we first compute, for each countObject, the number of
   *  genomicIntervals of the GenomicAlignment which have a conforming alignment
   *  with the exons that are associated to the countObject.
   *
   *  A countObject has a conforming alignment with the GenomicAlignment if the
   *  number of genomicIntervals associated to the countObject is greater than
   *  countObjectThreshold. So, for instance, if we are interested in exons,
   *  then countObjectThreshold = 1; if we are interested in junctions, then
   *  countObjectThreshold = 2; if we are interested in genes, then all GenomicIntervals
   *  should have a conforming alignment with an exon associated to the countObject -
   *  - this is indicated by a negative countObjectThreshold (note that in this case
   *  the overlap threshold needs to be set to one).
   *
   ***********************************************************************************/

  public HashSet<String> getConformingCountObjects (int countObjectThreshold, HashSetTable<Exon, String> countObjectTable,
						    boolean countConsecutive, boolean excludeAmbiguousGenomicIntervals) throws IOException {


    /* A negative value of countObjectThreshold indicates that the countObject needs to have conforming alignment with all
       genomic intervals */
    if (UtilLib.getGeneCountMode () || countObjectThreshold < 0) {
      countObjectThreshold = maxFieldIndex + 1;
    }

    /* countObjectGenomicIntervalTable contains the array of indices of the genomic intervals for which there is a
       conforming alignment for a given count object */
    Hashtable<String, int []> countObjectGenomicIntervalTable = new Hashtable <String, int []> (500);

    /* conformingCountObjectSet contains the countObjects for whom the number of genomic intervals with a conforming
       alignment to an exon associated to the countObject is at least countObjectThreshold */
    HashSet<String> conformingCountObjectSet = new HashSet <String> (500);

    /* For each count object compute the number of genomic intervals for which there is a conforming alignment */
    for (GenomicInterval genomicInterval: genomicIntervals) {

      /* curCountObjects is the set of count objects associated to exons having a conforming alignment with the current genomicInterval.
	 This is used to determine consecutive genomicIntervals having conforming alignments with the exons associated to a countObject. */
      HashSet<String> curCountObjects = new HashSet <String> (500);

      HashSet<Exon> exonSet = genomicInterval.updateExonSet ();

      if (debugLevel >= 2) {
	System.err.println("Genomic interval: " + genomicInterval + " exon set: " + exonSet + ", field index: " + genomicInterval.getFieldIndex () +
			   ", excludeAmbiguousGenomicIntervals: " + excludeAmbiguousGenomicIntervals);
      }

      /* We exclude genomic intervals if more than two exons conform */
      if (excludeAmbiguousGenomicIntervals && exonSet.size() >= 2) {
	if (debugLevel >= 2) {
	  System.err.println("Excluding exon set: " + exonSet + " for genomic interval: " + genomicInterval);
	}
	exonSet.clear ();
      }

      for (Exon exon: exonSet) {

	HashSet<String> countObjectSet = null;
	if (countObjectTable != null) {
	  countObjectSet = countObjectTable.get (exon);
	  if (countObjectSet == null) {
	    if (debugLevel >= 1) {
	      System.err.println ("No count object associated to exon " + exon);
	    }
	  }
	}

	if (countObjectSet == null) {
	  countObjectSet = new HashSet<String> (10);
	  countObjectSet.add(exon.toString());
	  countObjectTable.putValue (exon, exon.toString());
	}

	if (debugLevel >= 2) {
	  System.err.println("count object set for exon " + exon + ": " + countObjectSet);
	}

	for (String countObject : countObjectSet) {

	  if (debugLevel >= 2) {
	    System.err.println("Count object: " + countObject);
	  }

	  if (! curCountObjects.contains (countObject)) {
	    curCountObjects.add(countObject);
	    int [] fieldIndexArray = countObjectGenomicIntervalTable.get(countObject);
	    if (fieldIndexArray == null) {
	      fieldIndexArray = new int [maxFieldIndex + 1];
	      for (int j = 0; j < fieldIndexArray.length; j++) {
		fieldIndexArray[j] = 0;
	      }
	      countObjectGenomicIntervalTable.put(countObject, fieldIndexArray);
	    }

	    if (debugLevel >= 2) {
	      System.err.println("Setting " + genomicInterval.getFieldIndex () + " of " + maxFieldIndex);
	    }

	    if (genomicInterval.getFieldIndex () >= 0) {
	      fieldIndexArray[genomicInterval.getFieldIndex ()] = 1;
	    }
	  }
	}
      }
    }

    /* Add the countObjects whose number of (consecutive) genomic intervals with a conforming alignment to an associated exon is at least
       countObjectThreshold to conformingCountObjectSet */
    for (String countObject: countObjectGenomicIntervalTable.keySet ()) {

      if (debugLevel >= 2) {
	System.err.println("Count object: " + countObject);
      }

      int [] fieldIndexArray = countObjectGenomicIntervalTable.get(countObject);

      if (debugLevel >= 2) {
	for (int j = 0; j < fieldIndexArray.length; j++) {
	  System.err.println("fieldIndexArray["+j+"]: " + fieldIndexArray[j]);
	}
      }


      
      if (fieldIndexArray != null) {
	int sum = 0;
	int maxConsecutiveSum = -1;
	for (int j = 0; j < fieldIndexArray.length; j++) {
	  if (countConsecutive && fieldIndexArray[j] == 0 && sum > 0) {
	    maxConsecutiveSum = Math.max(maxConsecutiveSum, sum);
	    sum = 0;
	  }
	  sum = sum + fieldIndexArray[j];
	}

	if (debugLevel >= 2) {
	  System.err.println("sum: " + sum + ", maxConsecutiveSum: " + maxConsecutiveSum);
	}

	if (Math.max(sum, maxConsecutiveSum) >= countObjectThreshold) {
	  
	  if (debugLevel >= 2) {
	    System.err.println("Adding count object: " + countObject + ", sum: " + sum + ", threshold: " + countObjectThreshold);
	  }
	  
	  conformingCountObjectSet.add(countObject);
	}
      } else {
	throw new IOException ("null countObjectGenomicIntervalNum in genomic alignment: " + this + "\n" +
			       "for count object: " + countObject + "\n" +
			       "countObjectGenomicIntervalTable: " + countObjectGenomicIntervalTable);	
      }
      
    }

    return conformingCountObjectSet;
    
  }


  /***********************************************************************************
   *
   *     combine the exon set that are associated to each genomic interval
   *
   ***********************************************************************************/

  public void combineExons (GenomicAlignment genomicAlignment) throws IOException {

    if (genomicAlignment == null) {
      throw new IOException ("ERROR: merging with null genomicAlignment: " + this);
    }

    if (debugLevel >= 2) {
      System.err.println ("Merging - genomicIntervals: " + genomicIntervals);
    }

    Iterator<GenomicInterval> it = genomicAlignment.getGenomicIntervals().iterator ();

    GenomicInterval genomicAlignmentInterval = null;
    int i = 0;
    for (GenomicInterval genomicInterval: genomicIntervals) {
      if (! genomicInterval.toString().equals("")) {
	/* Note that we walk through the intervals in genomicIntervals and in genomicAlignment.getGenomicIntervals() in lock step */
	if (it.hasNext()) {
	  genomicAlignmentInterval = it.next ();
	  while (it.hasNext() && genomicAlignmentInterval.toString().equals("")) {
	    genomicAlignmentInterval = it.next ();
	  }
	} else {
	  throw new IOException ("ERROR: too few genomic intervals in: " + genomicAlignment + " when merged with: " + this);
	}
	
	if (! genomicInterval.toString().equals(genomicAlignmentInterval.toString ())) {
	  if (UtilLib.warningsOn ()) {
	    System.err.println("Alignment: " + weightObjectAlignment.getAlignmentId () + " - " + toGenomicAlignmentString ());
	    System.err.println ("Genomic intervals:");
	    for (GenomicInterval gi: genomicIntervals) {
	      System.err.println (gi.toCompleteString());
	    }
	    
	    System.err.println ("Genomic intervals to be merged:");
	    System.err.println("From alignment: " + genomicAlignment.getWeightObjectAlignment ().getAlignmentId () + " - " + genomicAlignment.toGenomicAlignmentString ());
	    for (GenomicInterval gi: genomicAlignment.getGenomicIntervals()) {
	      System.err.println (gi.toCompleteString());
	    }
	  }
	  System.err.println ("Warning: " + i + "th genomic interval " + genomicAlignmentInterval + " of alignment " +
				 genomicAlignment + " of " + weightObjectAlignment.getAlignmentId () +
				 " is different from " + genomicInterval.toCompleteString() + " of " + this + " when merging.");
	}
	
	genomicInterval.addExonSet (genomicAlignmentInterval.getExonSet ());

	if (genomicInterval.getExonSet ().size () > 1 && debugLevel >= 2) {
	  System.err.println ("GA Weight object alignment: " + weightObjectAlignment);
	  System.err.println ("Genomic alignment: " + this);
	  System.err.println ("Genomic interval: " + genomicInterval);
	  System.err.println ("Exons:");
	  for (Exon exon: genomicInterval.getExonSet ()) {
	    System.err.println(exon);
	  }
	  // System.exit(0);
	}
           
      }
    }

    if (debugLevel >= 2) {
      System.err.println ("Done.");
    }
    
  }

  /***********************************************************************************
   *
   *  isTranscriptExonAlignment
   *
   **********************************************************************************/

  public boolean isTranscriptExonAlignment () throws IOException {

    if (! isTranscriptExonAlignmentIsSet) {

      isTranscriptExonAlignmentIsSet = true;
      for (GenomicInterval genomicInterval: genomicIntervals) {
	if (! genomicInterval.isTranscriptExon()) {
	  return false;
	}
      }

      isTranscriptExonAlignment = true;
    }
    
    return isTranscriptExonAlignment;

  }


  /***********************************************************************************
   *
   *                         printExonSet
   *
   ***********************************************************************************/

  public void printExonSet () {

    System.out.println("    Exon set of " + this);
    
    for (GenomicInterval genomicInterval: genomicIntervals) {
      genomicInterval.printExonSet ();
    }
    
  }

  
  /***********************************************************************************
   *
   *                            getGenomicBedEntryString
   *
   **********************************************************************************/

  public Vector<String> getGenomicBedEntryString (String bedEntryId) throws IOException {

    Vector<String> genomicBedEntries = new Vector<String> (10);
    String oldStrand = null;
    String oldChromosome = null;
    for (GenomicInterval genomicInterval: genomicIntervals) {
      if (oldChromosome != null && genomicInterval.getChromosome() != null) {
	if (! UtilLib.isTranscriptExon(oldChromosome) && !  UtilLib.isTranscriptExon(genomicInterval.getChromosome ())) {
	  if (! oldChromosome.equals(genomicInterval.getChromosome ())) {
	    throw new IOException ("Different chromosomes for genomic intervals of " + bedEntryId);
	  }
	  if (oldStrand != null && ! oldStrand.equals(genomicInterval.getStrand ())) {
	    throw new IOException ("Different strands for genomic intervals of " + bedEntryId);
	  }
	}
      }
      if (genomicInterval.getChromosome() != null) {
	String genomicIntervalBedEntryString = genomicInterval.toGenomicBedEntryString (bedEntryId);
	if (genomicIntervalBedEntryString != null && ! genomicIntervalBedEntryString.equals("")) {
	  genomicBedEntries.add(genomicIntervalBedEntryString);
	}
      }
      oldChromosome = genomicInterval.getChromosome ();
      oldStrand = genomicInterval.getStrand ();
    }

    return genomicBedEntries;
    
  }


  /***********************************************************************************
   *
   *                            adjustGenomicAlignmentString
   *
   ***********************************************************************************/

  public void adjustGenomicAlignmentString (Vector<WeightObjectAlignment> weightObjectAlignments, int genomicAlignmentIndex) throws IOException {

    alignmentString = "";
    for (GenomicInterval genomicInterval: genomicIntervals) {
      if (genomicInterval.getChromosome() == null) {
	throw new IOException ("Null chromosome in genomice interval");
      }
      if (alignmentString == "") {
	alignmentString = genomicInterval.toString (weightObjectAlignments, genomicAlignmentIndex);
      } else {
	alignmentString = alignmentString + ":" + genomicInterval.toString (weightObjectAlignments, genomicAlignmentIndex);
      }
    }
    
  }
  

  /***********************************************************************************
   *
   *                            setGenomicAlignmentString
   *
   ***********************************************************************************/

  public void setGenomicAlignmentString () throws IOException {

    alignmentString = null;
    for (GenomicInterval genomicInterval: genomicIntervals) {
      if (genomicInterval.getChromosome() == null) {
	throw new IOException ("Null chromosome in genomice interval");
      }

      if (alignmentString == null) {
	alignmentString = genomicInterval.toString ();
      } else {
	alignmentString = alignmentString + ":" + genomicInterval.toString ();
      }

    }
    
  }

  
  /***********************************************************************************
   *
   *                            toGenomicAlignmentString
   *
   ***********************************************************************************/

    public String toGenomicAlignmentString () throws IOException {

    if (alignmentString == null) {
      setGenomicAlignmentString ();
    }

    return alignmentString;
    
  }


  /***********************************************************************************
   *
   *                            toString
   *
   ***********************************************************************************/
  
  public String toString () {

    return alignmentString;
    
  }


  /* public boolean equalsNew (Object o) {

    if (debugLevel >= 2) {
      System.out.println ("Equals start");
      System.out.println (this + " vs " + o);
    }

    boolean value = super.equals(o);

    if (debugLevel >= 2) {
      System.out.println ("Equals end");
    }

    return value;
  }

  public int hashCodeNew () {

    if (debugLevel >= 2) {
      System.out.println ("Hashcode start");
      System.out.println ("Hashcode for " + this);
    }

    int value = super.hashCode ();

    if (debugLevel >= 2) {
      System.out.println ("Hashcode end");
    }

    return value;
    
    } */

  
  /***********************************************************************************
   *
   *                              Basic methods
   *
   ***********************************************************************************/
  
  public int getStart () {
    return start;
  }

  public int getChrStart () {
    return start;
  }


  public String getChromosome () {
    return chromosome;
  }

  public String getStrand () {
    return strand;
  }

  public String getReferenceGenomicAlignmentOrientation () {
    return referenceGenomicAlignmentOrientation;
  }

  public int getLengthDifference () {
    return lengthDifference;
  }

  public int getCorrectionDifference () {
    return correctionDifference;
  }

  public int getReadAlignedLength () {

    if (readAlignedLength >= 0) {
      return readAlignedLength;
    }

    readAlignedLength = 0;   
    for (GenomicInterval genomicInterval: genomicIntervals) {
      readAlignedLength = readAlignedLength + genomicInterval.getReadAlignedLength ();
    }

    return readAlignedLength;

  }
  
  public WeightObjectAlignment getWeightObjectAlignment () {
    return weightObjectAlignment;
  }

  public String getReferenceId () {
    return weightObjectAlignment.getReferenceId();
  }

  public int getReadAlignStart () {
    return readAlignStart;
  }

  public boolean isFirst () {

    if (weightObjectAlignment.getGenomicAlignment1 () == null) {
      return false;
    }
    
    return this.equals(weightObjectAlignment.getGenomicAlignment1 ());
  }
  

  public boolean isPairedEnd () {
    return weightObjectAlignment.isPairedEnd ();
  }

  
  public TreeSet<GenomicInterval> getGenomicIntervals () {
    return genomicIntervals;
  }


  public GenomicAlignment getPairedGenomicAlignment () {

    if (isFirst ()) {
      return weightObjectAlignment.getGenomicAlignment2 ();
    }

    return weightObjectAlignment.getGenomicAlignment1 ();
    
  }

  public boolean fromGenomicBedRecords () {
    return genomicBedRecords;
  }

  public String getGenomicIntervalString () {
    String genomicIntervalString = null;
    for (GenomicInterval genomicInterval: getGenomicIntervals ()) {
      if (genomicIntervalString == null) {
	genomicIntervalString = genomicInterval.toString();
      } else {
	genomicIntervalString = genomicIntervalString + "\n" + genomicInterval.toString();
      }
    }

    return genomicIntervalString;
  }

  public int compareTo (Object o) {

    GenomicAlignment genomicAlignment = (GenomicAlignment) o;

    Iterator<GenomicInterval> it   = getGenomicIntervals ().iterator ();
    Iterator<GenomicInterval> gaIt = genomicAlignment.getGenomicIntervals ().iterator ();

    while (it.hasNext () && gaIt.hasNext ()) {

      GenomicInterval gi   = it.next();
      GenomicInterval gaGi = gaIt.next();

      if (gi.compareTo(gaGi) != 0) {
	return gi.compareTo(gaGi);
      }
      
    }

    if (! it.hasNext () && ! gaIt.hasNext ()) {
      return 0;
    }

    if (! it.hasNext ()) {
      return -1;
    }

    if (! gaIt.hasNext ()) {
      return 1;
    }

    return 0;
    
  }


  public boolean equals (Object o) {

    GenomicAlignment genomicAlignment = (GenomicAlignment) o;

    try {
      return toGenomicAlignmentString ().equals (genomicAlignment.toGenomicAlignmentString ());
    } catch (IOException e) {
      System.err.println (e==null?"Null message":e.getMessage());
      System.exit (1);
    }

    return false;
    
  }


  public int hashCode () {

    try {
      return toGenomicAlignmentString ().hashCode();
    } catch (IOException e) {
      System.err.println (e==null?"Null message":e.getMessage());
      System.exit (1);
    }

    return -1;
    
  }


  public Vector<BedEntry> getBedEntries () {
    
    Vector<BedEntry> bedEntries = new Vector<BedEntry> (genomicIntervals.size());
    for (GenomicInterval genomicInterval: genomicIntervals) {
      bedEntries.add(genomicInterval.toBedEntry());
    }
    return bedEntries;
    
  }



  public TreeSet<BedRecord> getBedRecords () {

    if (weightObjectAlignment.getGenomicAlignment1 () != null && this.equals(weightObjectAlignment.getGenomicAlignment1 ())) {
      return weightObjectAlignment.getBedRecords1 ();
    }

    if (weightObjectAlignment.getGenomicAlignment2 () != null && this.equals(weightObjectAlignment.getGenomicAlignment2 ())) {
      return weightObjectAlignment.getBedRecords2 ();
    }

    System.err.println("weightObjectAlignment: " + weightObjectAlignment + " does not have genomic alignment for " + this);

    return new TreeSet<BedRecord> ();
    
  }

  
}
