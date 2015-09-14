/**File: ComputeCounts.java 

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
 * 
 *                           Class ComputeCounts
 *
 * This is the main method to compute exon, gene, and junction counts. It works as
 * follows:
 *
 * countTable = Hashtable(<String>, <Double>)
 *
 * for each bedRecord with overlap > 0 do
 *   if weightObject is null:
 *     weightObject = WeightObject (bedRecord)
 *   elsif oldFragmentId == bedRecord.getFragmentId:
 *     add bedRecord to weightObject
 *   else:
 *     get readWeight of oldFragmentId
 *     numGenomicAlignments = weightObject.getNumGenomicAlignments();
 *     objectWeight = min(readWeight, 1 / numGenomicAlignments)
 *     if objectWeight > readWeightThreshold:
 *       processWeightObject (weightObject, countTable)
 *     weightObject = WeightObject (bedRecord)
 *     oldFragmentId = bedRecord.getFragmentId
 *
 * process last weightObject:
 *   processWeightObject (weightObject, countTable)
 *
 * for countObjectId in countTable.keys union Exon.getCountObjectIds do
 *   print countObjectId, countTable[countObjectId]
 *
 ***********************************************************************************/

public class ComputeCounts {

  private static int debugLevel = 2;

  private static Hashtable<String, Hashtable<GenomicAlignment, Integer>> overlapTable1 = null;
  private static Hashtable<String, Hashtable<GenomicAlignment, Integer>> overlapTable2 = null;
  
  /* private static Hashtable<String, Integer> countObjectLengthTable = null; */
  private static HashSetTable<String, Exon> countObjectExonTable   = null;

  private static Hashtable<String, Double> readWeightTable = new  Hashtable<String, Double> (50 * 1000 * 1000);

  private static final String [] specialCountObjectIdSetValues = new String [] {};
  
  private static HashSet<String> specialCountObjectIdSet = new HashSet<String> (Arrays.asList(specialCountObjectIdSetValues));

  private static double specialCount = 0;
  private static String specialReadIds = "";
  
  private static HashSet<String> specialReadIdSet = new HashSet<String> ();

  public static String getSpecialReadId () {
    return specialReadIds;
  }

  public static HashSet<String> getSpecialReadIdSet () {
    return specialReadIdSet;
  }
  
  private static HashSet<WeightObject> countObjectIdWeightObjects    = new HashSet<WeightObject> ();
  private static HashSet<WeightObject> countObjectIdWeightObjectsNew = new HashSet<WeightObject> ();

  private static PriorityQueue<String> weightFragmentIdQueue = new PriorityQueue<String> (1000);
  private static Hashtable<String, Double> weightFragmentIdTable = new Hashtable<String, Double> (1000);
  

  /***********************************************************************************
   *
   *                           addToTableEntry  
   *
   ***********************************************************************************/

  public static void addToTableEntry (Hashtable<String, Double> entryTable, String key, double value) {

    Double tableValueDouble = entryTable.get(key);
    double tableValue = 0;

    if (tableValueDouble != null) {
      tableValue = tableValueDouble.doubleValue();
    }
	    
    tableValue = tableValue + value;
    entryTable.put(key, new Double(tableValue));
  }


  /***********************************************************************************
   * 
   *                           getReadWeight
   *
   *  Retrieve the weight for fragmentId from weightReader: 
   *    Read entries from weightReader until fragmentId is found
   *    if an entry with a larger id is found, raise an exception
   *
   ***********************************************************************************/

  private static double getReadWeight (String fragmentId, BufferedReader weightReader) throws IOException {

    if (debugLevel >= 1) {
      System.err.println ("Calling getReadWeight for: " +  fragmentId);
    }

    if (debugLevel >= 2 || debugLevel >= 1) {
      System.err.println ("Looking for fragmentId: " + fragmentId);
    }

    if (weightReader == null) {
      return 1.0;
    }

    String line = weightReader.readLine ();

    if (line == null) {
      return 1.0;
    }

    StringTokenizer st = new StringTokenizer (line, "\t");

    String weightFragmentId = "";
    if (st.hasMoreTokens()) {
      weightFragmentId = st.nextToken ();	
    }

    Double alignmentNum = null;
    if (st.hasMoreTokens()) {
      alignmentNum = new Double(st.nextToken ());	
    } else {
      throw new IOException ("No weight for read " + weightFragmentId + " found.");
    }

    if (debugLevel >= 2 || debugLevel >= 1) {
      System.err.println ("Weight fragment id: " + weightFragmentId + ", alignmentNum: " + alignmentNum);
    }

    if (weightFragmentIdQueue.size () > 1000) {
      String minFragmentId = weightFragmentIdQueue.poll ();
      if (debugLevel >= 2 || debugLevel >= 1) {
	System.err.println ("Min fragment id: " + minFragmentId + ", weightFragmentIdQueue.size (): " + weightFragmentIdQueue.size ());
      }

      if (minFragmentId != null) {
	weightFragmentIdTable.remove (minFragmentId);
      }
    }

    weightFragmentIdQueue.add (weightFragmentId);
    weightFragmentIdTable.put (weightFragmentId, alignmentNum);

    while (! weightFragmentId.equals (fragmentId)) {

      int compValue = weightFragmentId.compareTo (fragmentId);

      if (compValue == 1) {
	alignmentNum = weightFragmentIdTable.get (fragmentId);
	if (alignmentNum == null) {
	  String minFragmentId = weightFragmentIdQueue.poll ();
	  throw new IOException ("Fragment id: " + fragmentId + " not found in read weight file - current weight fragment id: " + weightFragmentId +
				 ", minimum stored weight fragment id: " + minFragmentId);
	} else {
	  return 1.0 / alignmentNum.doubleValue ();
	}
      }

      line = weightReader.readLine ();

      if (line == null) {
	throw new IOException ("Fragment id: " + fragmentId + " not found in read weight file.");
	/* return 1.0; */
      }

      st = new StringTokenizer (line, "\t");
      
      if (st.hasMoreTokens()) {
	weightFragmentId = st.nextToken ();	
      }

      if (st.hasMoreTokens()) {
	alignmentNum = new Double(st.nextToken ());	
      } else {
	throw new IOException ("No weight for read " + weightFragmentId + " found.");
      }

      if (debugLevel >= 2 || debugLevel >= 1) {
	System.err.println ("Weight fragment id: " + weightFragmentId + ", alignmentNum: " + alignmentNum);
      }

      if (weightFragmentIdQueue.size () > 1000) {
	String minFragmentId = weightFragmentIdQueue.poll ();
	if (debugLevel >= 2 || debugLevel >= 1) {
	  System.err.println ("min fragment id: " + minFragmentId);
	}

	if (minFragmentId != null) {
	  weightFragmentIdTable.remove (minFragmentId);
	}
      }
      
      weightFragmentIdQueue.add (weightFragmentId);
      weightFragmentIdTable.put (weightFragmentId, alignmentNum);
      

    }

    if (debugLevel >= 1) {
      System.err.println ("Weight: " + alignmentNum + " from " + line + " for " + fragmentId);
    }

    return 1.0 / alignmentNum.doubleValue ();

  }
  
  
  /***********************************************************************************
   * 
   *                          processWeightObject
   *
   *  processWeightObject takes a weightObject which is essentially a collection
   *  of bedRecords for one fragment and computes the conforming count objects, that is,
   *  exons, genes, or junctions which overlap with the fragment in way that respects
   *  the splicing pattern of the read and that respects the structure of the count object.
   *
   *  Preprocess genomic alignments:
   *    weightObject.adjustGenomicAlignments: Compute genomic alignments and adjust them if necessary
   *    weightObject.mergeWeightAlignments: Combine the exons of weight object alignments with the same genomic alignment
   *
   *  
   *
   ***********************************************************************************/

  private static double processWeightObject (WeightObject weightObject, double objectWeight, int minExonNum, boolean countConsecutive,
					     boolean pairedEndOnlyMode, Hashtable<String, Double> countTable, boolean excludeAmbiguousReads,
					     boolean genomicIntervalMode, boolean checkWeights) throws IOException {

    if (debugLevel >= 1) {
      System.err.println ();
      System.err.println("processWeightObject new start.");
    }
    
    if (debugLevel >= 2 || debugLevel >= 1) {
      System.err.println("weightObject: " + weightObject + ", objectWeight: " + objectWeight);
    }

    HashSetTable<Exon, String> countObjectTable = Exon.getCountObjectTable ();

    if (getSpecialReadIdSet().contains(weightObject.getFragmentId ())) {
      debugLevel = 2;
    }


    /* Compute genomic alignment strings and adjust them if necessary */
    weightObject.adjustGenomicAlignments ();
        
    /* Combine the exons of weight object alignments with the same genomic alignment string */
    weightObject.mergeWeightAlignments ();

    /* We store the contribution of the weightObject to a countObject in the weightObjectCountTable. We do this in order to account for
       weightObjects that have more than one (conforming) alignment to the countObject */
    Hashtable<String, Double> weightObjectCountTable  = new Hashtable<String, Double> (100);
    Hashtable<String, Double> weightObjectCountTable1 = new Hashtable<String, Double> (100);
    Hashtable<String, Double> weightObjectCountTable2 = new Hashtable<String, Double> (100);
    boolean excludeAmbiguousGenomicIntervals = excludeAmbiguousReads && genomicIntervalMode;
    
    if (debugLevel >= 2 || debugLevel >= 1) {
      System.err.println("processWeightObject - excludeAmbiguousGenomicIntervals: " + excludeAmbiguousGenomicIntervals + ", excludeAmbiguousReads: " + 
			 excludeAmbiguousReads + ", genomicIntervalMode: " +  genomicIntervalMode);
    }

    for (WeightObjectAlignment weightObjectAlignment: weightObject.getWeightObjectAlignmentSet()) {

      if (debugLevel >= 2 || debugLevel >= 1) {
	System.err.println("processWeightObject - WeightObjectAlignment: " + weightObjectAlignment);
      }

      GenomicAlignment genomicAlignment1 = weightObjectAlignment.getGenomicAlignment1();
      HashSet<String> conformingCountObjectIds1 = new HashSet<String> ();
      if (genomicAlignment1 != null) {
	conformingCountObjectIds1 =
	  genomicAlignment1.getConformingCountObjects (minExonNum, countObjectTable, countConsecutive, excludeAmbiguousGenomicIntervals);
      }

      if (debugLevel >= 2 || debugLevel >= 1) {
	System.err.println("conformingCountObjectIds1: " + conformingCountObjectIds1);
      }

      
      GenomicAlignment genomicAlignment2 = weightObjectAlignment.getGenomicAlignment2();
      HashSet<String> conformingCountObjectIds2 = new HashSet<String> ();
      if (genomicAlignment2 != null) {
	conformingCountObjectIds2 =
	  genomicAlignment2.getConformingCountObjects (minExonNum, countObjectTable, countConsecutive, excludeAmbiguousGenomicIntervals);
      }

      if (debugLevel >= 2 || debugLevel >= 1) {
	System.err.println("conformingCountObjectIds2: " + conformingCountObjectIds2);
      }

      /* Note that retainAll and addAll also work if genomicAlignment1 == null or genomicAlignment2 == null */
      if (weightObjectAlignment.isPairedEnd ()) {
	if (pairedEndOnlyMode) {
	  conformingCountObjectIds1.retainAll (conformingCountObjectIds2);
	} else {
	  conformingCountObjectIds1.addAll (conformingCountObjectIds2);
	}
	conformingCountObjectIds2.clear ();
	if (conformingCountObjectIds2.size() > 0) {
	  throw new IOException ("Clearing conformingCountObjectIds2 did not work ... exiting");
	}
      }

      for (String countObjectId: conformingCountObjectIds1) {
	if (specialCountObjectIdSet.contains(countObjectId)) {
	  System.err.println("Putting " + countObjectId + " into the weightObjectCountTable for " + weightObject + " with weight " + objectWeight +
			     " based on alignment 1 of " + weightObjectAlignment);

	}
	addToTableEntry (weightObjectCountTable1, countObjectId, objectWeight);
      }
	
      for (String countObjectId: conformingCountObjectIds2) {
	if (specialCountObjectIdSet.contains(countObjectId)) {
	  System.err.println("Putting " + countObjectId + " into the weightObjectCountTable for " + weightObject + " with weight " + objectWeight +
			     " based on alignment 1 of  " + weightObjectAlignment);

	}
	addToTableEntry (weightObjectCountTable2, countObjectId, objectWeight);
      }

      if (debugLevel >= 2 || debugLevel >= 1) {
	System.err.println("Resulting conformingCountObjectIds1: " + conformingCountObjectIds1);
      }

    }

    if (debugLevel >= 2 || debugLevel >= 1) {
      if (weightObjectCountTable.keySet ().size () == 0 && weightObjectCountTable1.keySet ().size () == 0 && weightObjectCountTable2.keySet ().size () == 0) {
	System.err.println ("No count object ids found for: " + weightObject.getFragmentId ());
      } else if (weightObjectCountTable.keySet ().size () > 1 / objectWeight) {
	System.err.println ("More than one count object id for " + weightObject.getFragmentId () + ": " + weightObjectCountTable.keySet () +
			    ", objectWeight: " + objectWeight);
      }
    }


    if (excludeAmbiguousReads && ! genomicIntervalMode) {
      /* Form the union of the count object ids for the single read and paired-end alignments */
      HashSet<String> countObjectIds = new HashSet<String> ();
      countObjectIds.addAll (weightObjectCountTable.keySet ());
      countObjectIds.addAll (weightObjectCountTable1.keySet ());
      countObjectIds.addAll (weightObjectCountTable2.keySet ());
      
      /* Check whether the weight object contributes to more than one count object */
      if (countObjectIds.size() >= 2) {
	return 0;
      }
    }

    /* Note that sumWeightAdded depends on the number of identifiers that are associated with a genomic locus; so even for
       a uniquely mapping read there may be a number of count objects that it contributes to, for instance, the exons it overlaps
       of which there can be several, for instance, overlapping exons, or anti-sense exons. In addition, as in the case of EMBL exons,
       several exon ids may be associated with one gene/genomic interval/strand tuple. */
    double sumWeightAdded = 0;
    for (String countObjectId: weightObjectCountTable.keySet ()) {
      if (debugLevel >= 2 || debugLevel >= 1) {
	System.err.println(weightObject.getFragmentId () + ": adding " + weightObjectCountTable.get(countObjectId).doubleValue () + " to " + countObjectId);
      }

      /* Disregarding a margin for rounding errors a weight object should contribute at most 1 to a count object */
      if (weightObjectCountTable.get(countObjectId).doubleValue () > 1.001 && checkWeights) {
	throw new IOException ("Weight " + weightObjectCountTable.get(countObjectId).doubleValue () + " of weight object " + weightObject + " for count object " +
			       countObjectId + " is larger than 1.");
      }
      addToTableEntry (countTable, countObjectId, Math.min (1, weightObjectCountTable.get(countObjectId).doubleValue ()));
      if (specialCountObjectIdSet.contains(countObjectId)) {
	System.err.println("Adding pe weight object: " + weightObject + " with weight " + weightObjectCountTable.get(countObjectId).doubleValue () + " to " + countObjectId +
			   " giving a total weight of " + countTable.get(countObjectId));
      }
      sumWeightAdded = sumWeightAdded + Math.min (1, weightObjectCountTable.get(countObjectId).doubleValue ());
    }

    HashSet<String> singleReadCountObjectIds = new HashSet<String> ();
    singleReadCountObjectIds.addAll (weightObjectCountTable1.keySet ());
    singleReadCountObjectIds.addAll (weightObjectCountTable2.keySet ());

    for (String countObjectId: singleReadCountObjectIds) {

      double weight1 = 0;
      if (weightObjectCountTable1.keySet ().contains(countObjectId)) {
	weight1 = weightObjectCountTable1.get(countObjectId).doubleValue ();
      }
      
      double weight2 = 0;
      if (weightObjectCountTable2.keySet ().contains(countObjectId)) {
	weight2 = weightObjectCountTable2.get(countObjectId).doubleValue ();
      }

      if (debugLevel >= 2 || debugLevel >= 1) {
	System.err.println(weightObject.getFragmentId () + ": adding max(" + weight1 + ", " + weight2 + ") to " + countObjectId);
      }
      	  
      double weight = Math.max (weight1, weight2);
      addToTableEntry (countTable, countObjectId, Math.min (1, weight));
      if (specialCountObjectIdSet.contains(countObjectId)) {
	System.err.println("Adding sr weight object: " + weightObject + " with weight max (" + weight1 + ", " + weight2 + ") to " + countObjectId +
			   " yielding a count of " + countTable.get(countObjectId));
      }
      sumWeightAdded = sumWeightAdded + Math.min (1, weight);
    }

    if (getSpecialReadIdSet().contains(weightObject.getFragmentId ())) {
      debugLevel = UtilLib.getDebugLevel ();
    }

    if (debugLevel >= 1) {
      System.err.println("processWeightObject new done with " + sumWeightAdded + " added.");
    }

    return sumWeightAdded;
       
  }

 

  /***********************************************************************************/

   private static void printHelp () {
    System.out.println("ComputeCounts\n" +                                              
    "USAGE: ComputeCounts [-W <read weight thresh.>] [-O <overlap thresh.>] [-u] \n" +
    "   [-g|-e|-j] [-p] [-N] [-m <count object map file>] [-w <read weight file>|none] \n" +
    "    [-M <count object file>] -b <intersect. bed file> -o <outputFile>\n" +
    "\n" +
    " -U: Count only unambiguous reads, i.e. reads that map to at most one gene,\n" +
    "     exon, or junction.\n" +
    " -u: Do not use read weights in counting. Read weights can still be used for filtering.\n" +
    " -b STRING: intersect. bed file - the bed file containing the intersection of\n" +
    "     the exons transcript intervals and the reads mapped to the transcripts.\n" +
    "     (- for STDIN) [-]\n" +
    " -o STRING: output file - the file to which the output is written (- for STDOUT)\n" + 
    "    [-]\n" +
    " -g: Compute gene counts: set paired-end mode (option -p), set the minimal overlap\n" +
    "    equal to the <read length> (option -o), and consider non-genomic exons\n" +
    "    as well (\"tr-\" exons). Note that one of -g, -e, or -j must set.\n" +
    " -e: Compute exon counts: ensure that a read is only counted if it overlaps\n" +
    "     an exon in a splicing compatible manner .\n" +
    " -j: Compute junction counts: ensure that a read is only counted if it overlaps\n" +
    "     both exons with the minimal overlap (option -o) and the exons occur\n" +
    "     consecutively in the alignment of the read.\n" +
    " -i: Compute intron counts. All reads are counted for a genomic interval. Spliced" +
    "     and unspliced reads." + 
    " -N: count all genomic alignments - even if they do not respect\n" +
    "     splicing patterns (non-splice conforming).\n" +
    " -C: count all genomic alignments - even if they do not respect\n" +
    "     splicing patterns.\n" +
    " -O INT: Minimal overlap of a read with a count object (for paired-end\n" +
    "    alignments twice the overlap is required.) [8]\n" +
    " -W FLOAT: Minimal weight of a read; reads with a lower weight are disregarded \n" +
    "    [0.01].\n" +
    " -w read weight file: a file containing the read weights with two columns:\n" +
    "      <fragment id>, <read weight>. Reads are weighted by the minimum of\n" +
    "      the read weight contained in this file,if present, and the read weight\n" +
    "      computed from the number of genomic alignments in the bed file. If the\n" +
    "      keyword \"none\" is selected, then the reads are not weighted.\n" +
    " -m count object map file: contains the exon to count object mapping. It\n" +
    "    consists of two columns, the exon id (<gene id>/<chr>/<start>/\n" +
    "    <end>/<strand>) and the id of the count object.\n" +
    " -M count object file: Contains the sequence of count objects which is used\n" +
    "    when outputting the results.\n" +
    " -s: set strand specific mode - only count reads on the same strand as the.\n" +
    "        gene, exon, or junction.\n" +
    " -S STRING: special count object id [<empty string>].\n" +
    " -n: output only non-zero counts (otherwise output all counts).\n" +
    "\n");
  }
                                     
  /***********************************************************************************/

  public static void main (String [] args) {

    int numDeletions  = 0;
    int numInsertions = 0;

    String separator = "#";
    
    String intersectionFilename = "-";
    String weightFilename = "none";
    String countObjectMapFilename = "";
    String countObjectFilename = "";
    String outputFilename = "-";

    double readWeightThreshold = 0.01;
    
    int  overlapThreshold = 1;

    String oldAlignmentId = "";
    String oldFragmentId = "";

    WeightObjectAlignment weightObjectAlignment = null;
    WeightObject          weightObject = null;

    boolean printLines             = false;
    boolean useReadWeights         = true;
    boolean excludeAmbiguousReads  = false;

    /* if pairedEndOnlyMode is set to true, then only the paired-end alignments are considered *if* a read has at least one;
       otherwise the single-read alignments are considered as well. */
    boolean pairedEndOnlyMode = false;
    boolean togglePairedEndOnlyMode = false;

    /* minExonNum is the minimum number of exons that a readWeightObject must overlap */
    int minExonNum = 1;
    /* If countConsecutive is set to true, minExonNum applies to consecutive exons */
    boolean countConsecutive = false;
    /* outputZeroes controls whether only non-zero count objects are output or all count objects */
    boolean outputZeroes     = true;

    int countUnit = 2 * 1250 * 1000;
    boolean checkWeights = false;
    boolean genomicIntervalMode = false;
    
    Getopt g = new Getopt("ComputeCounts.java", args, "ab:cCd:D:egGijm:M:nNo:O:pr:sS:uUw:W:zh");
    
    int c;
    String arg = "";
    
    c = g.getopt();
    
    while (c  != -1) {
      switch(c) {
      case 'a':
	UtilLib.setWarningsOn (true);
	break;	
      case 'b':
	intersectionFilename = g.getOptarg();
	break;
      case 'c':
	checkWeights = true;
	break;
      case 'C':
	UtilLib.setContainmentMode (true);
	break;
      case 'd':
	UtilLib.setDebugLevel (Integer.parseInt(g.getOptarg()));
	break;
      case 'D':
	specialReadIds = g.getOptarg();
	break;
      case 'e':
	UtilLib.setCountMode("exon");
	genomicIntervalMode = true;
	break;
      case 'g':
	UtilLib.setCountMode("gene");
	pairedEndOnlyMode = true;
	break;
      case 'i':
	/* Set intron count mode */
	UtilLib.setCountMode("junction");
	UtilLib.setSpliceConformingCountMode (false);
	break;
      case 'j':
	UtilLib.setCountMode("junction");
	minExonNum = 2;
	countConsecutive = true;
	break;
      case 'm':
	countObjectMapFilename = g.getOptarg();
	break;
      case 'M':
	countObjectFilename = g.getOptarg();
	System.err.println("Option -M used ... exiting.");
	System.exit (1);
	break;
      case 'n':
	outputZeroes = false;
	break;
      case 'N':
	UtilLib.setSpliceConformingCountMode (false);
	break;	
      case 'o':
	outputFilename = g.getOptarg();
	break;
      case 'O':
	overlapThreshold = Integer.parseInt(g.getOptarg());
	UtilLib.setOverlapThreshold (overlapThreshold);
	break;
      case 'p':
	togglePairedEndOnlyMode = true;
	break;
      case 's':
	UtilLib.setStrandedMode ();
	break;
      case 'S':
	specialCountObjectIdSet.add(g.getOptarg());
	break;
      case 'u':
	useReadWeights = false;
	break;
      case 'U':
	excludeAmbiguousReads = true;
	break;
      case 'w':
	weightFilename = g.getOptarg();
	break;	
      case 'W':
	readWeightThreshold = Double.parseDouble(g.getOptarg());
	break;
      case 'h':
	printHelp();
	System.exit(0);
	break;
      default:
	System.err.print("Error: getopt() returned unknown option: " + c + "\n");
      }
      c = g.getopt();
    }

    if (! UtilLib.getCountMode ().equals("gene") && ! UtilLib.getCountMode ().equals("exon")  && ! UtilLib.getCountMode ().equals("junction") ) {
      System.out.println ("Please specify the count mode: Options -g, -e, or -j.");
      printHelp();
      System.exit(0);
    }

    debugLevel = UtilLib.getDebugLevel ();
    UtilLib.setOverlapThreshold(overlapThreshold);
    
    StringTokenizer st = new StringTokenizer (specialReadIds, ":");
    while (st.hasMoreTokens ()) {
      specialReadIdSet.add(st.nextToken ());
    }

    if (togglePairedEndOnlyMode) {
      pairedEndOnlyMode = ! pairedEndOnlyMode;
    }

    if (excludeAmbiguousReads && genomicIntervalMode) {
      System.err.println("Computing unambiguous exon counts.");
    }
    
    String line = "";
    try {
      
      BufferedReader reader = UtilLib.getBufferedReader (intersectionFilename);
      
      System.err.println("Writing to " + (outputFilename.equals("-")?"stdout":outputFilename));
      System.err.flush();
      PrintWriter    outputWriter = UtilLib.getPrintWriter (outputFilename);

      System.err.println("Loading exon count object map file " + countObjectMapFilename);
      System.err.flush();
      BufferedReader countObjectMapReader = UtilLib.getBufferedReader (countObjectMapFilename);
      Exon.loadCountObjectFile (countObjectMapReader);
      if (debugLevel >= 3) {
	System.err.println("countObjectTable:");
	HashSetTable<Exon, String> countObjectTable = Exon.getCountObjectTable ();
	int i = 0;
	for (Exon exon: countObjectTable.keySet()) {
	  System.err.println(exon + ": " + countObjectTable.get(exon));
	  i = i + 1;
	  if (i > 10000) {
	    break;
	  }
	}
      }

      if (debugLevel >= 1) {
	System.err.println("useReadWeights: " + useReadWeights + ", readWeightThreshold: " + readWeightThreshold + ", weightFilename: " + weightFilename);
      }

      if (readWeightThreshold < 0 || weightFilename.equals("none")) {
	System.err.println("Computing unweighted counts.");
	useReadWeights = false;
      }

      if (debugLevel >= 1) {
	System.err.println("useReadWeights: " + useReadWeights);
      }

      BufferedReader weightReader = null;
      if (weightFilename != "" && ! weightFilename.equals("none")) {
	weightReader = UtilLib.getBufferedReader (weightFilename);
      }

      Hashtable<String, Double> countTable = new Hashtable<String, Double> (4 * 1000 * 1000);
      overlapTable1 = new Hashtable<String, Hashtable<GenomicAlignment, Integer>> (400);
      overlapTable2 = new Hashtable<String, Hashtable<GenomicAlignment, Integer>> (400);

      if (debugLevel >= 1) {
	System.err.println("Reading a line from file " + intersectionFilename);
      }
      line = reader.readLine();

      int lineNumber = 0;
      double objectWeight = 1.0;

      int numWeightObjects = 0;
      int numWeightObjectsIncluded = 0;
      int numWeightObjectsExcluded = 0;
      double sumObjectWeight = 0.0;
      double weightAdded = 0.0;

      System.err.println("Reading bed file: " + intersectionFilename + " (. = " + countUnit + " lines)");
      while (line != null) {
	
	BedRecord bedRecord = new BedRecord (line);

	if (bedRecord.getOverlap () > 0) {
	  
	  /* We consider fragments as weight objects. If there are two single-read alignments of one fragment
	     against a count object, then we count this as one. */
	  String  fragmentId = bedRecord.getFragmentId ();
	  
	  if (debugLevel >= 2 || debugLevel >= 1) {
	    System.err.println ("BedRecord Weight object id: " + fragmentId);
	  }
	  
	  if (printLines) {
	    System.err.println (bedRecord.toString());
	  }
	  
	  if (weightObject == null) {

	    weightObject = new WeightObject (bedRecord);
	    numWeightObjects++;
	    
	  } else if (! fragmentId.equals(oldFragmentId)) {

	    int numGenomicAlignments = weightObject.getNumGenomicAlignments();
	    if (weightReader != null) {
	      objectWeight = getReadWeight (oldFragmentId, weightReader);	      
	      if (debugLevel >= 1 || (! specialReadIds.equals("") && oldFragmentId.indexOf(specialReadIds) >= 0)) {
		System.err.println ("Old weight object id: " + oldFragmentId + " with weight: " + objectWeight);
	      }	      
	      if (debugLevel >= 1|| (! specialReadIds.equals("") && oldFragmentId.indexOf(specialReadIds) >= 0)) {
		System.err.println ("weightObject " + weightObject + " - numGenomicAlignments: " + numGenomicAlignments + " vs " + objectWeight);
	      }	      
	      if (objectWeight > 1.0 / numGenomicAlignments) {
		objectWeight = 1.0 / numGenomicAlignments;
	      }
	    }
	    
	    /* Finish old weight object */
	    if (objectWeight >= readWeightThreshold) {
	      double objectAddWeight = useReadWeights?objectWeight:1.0;
	      weightAdded = processWeightObject (weightObject, objectAddWeight, minExonNum, countConsecutive, pairedEndOnlyMode, countTable, excludeAmbiguousReads,
						 genomicIntervalMode, checkWeights);
	      sumObjectWeight = sumObjectWeight + weightAdded;
	      if (weightAdded > 0) {
		numWeightObjectsIncluded++;
	      } else {
		numWeightObjectsExcluded++;
	      }
	      	      
	    } else {
	      numWeightObjectsExcluded++;
	    }
	    
	    /* Create new weight object */
	    numWeightObjects++;
	    weightObject = new WeightObject (bedRecord);
	    
	  } else {
	    
	    if (debugLevel >= 2 || debugLevel >= 1) {
	      System.err.println ("adding bed record: " + bedRecord);
	    }
	    weightObject.addBedRecord (bedRecord);
	    
	  }
	  
	  if (! oldFragmentId.equals(fragmentId)) {
	    debugLevel = UtilLib.getDebugLevel ();
	  }
	  
	  oldFragmentId = fragmentId;
	}
	  
	lineNumber++;
	if (lineNumber % countUnit == 0) {
	  System.err.print(".");
	  System.err.flush();
	}
	  
	line = reader.readLine();

	if (debugLevel >= 1) {
	  System.err.println ("Line: " + line);
	}
	if (debugLevel >= 2 || debugLevel >= 1) {
	  System.err.println ("count table: " + countTable);
	}
	
      }

      if (debugLevel >= 1) {
	System.err.println ("Processing last fragment");
      }

      /* process last fragment */
      if (weightObject != null) {

	int numGenomicAlignments = weightObject.getNumGenomicAlignments();	  	  
	if (weightReader != null) {
	  objectWeight = getReadWeight (oldFragmentId, weightReader);	
	  if (debugLevel >= 1) {
	    System.err.println ("Old weight object id: " + oldFragmentId + " with weight: " + objectWeight);
	  }	  
	  if (debugLevel >= 1) {
	    System.err.println ("numGenomicAlignments: " + numGenomicAlignments);
	  }	  
	  if (objectWeight > 1.0 / numGenomicAlignments) {
	    objectWeight = 1.0 / numGenomicAlignments;
	  }
	}

	if (debugLevel >= 1) {
	  System.err.println ("Old weight object id: " + oldFragmentId + " with weight: " + objectWeight);
	}

	if (objectWeight >= readWeightThreshold) {
	  double objectAddWeight = useReadWeights?objectWeight:1.0;
	  weightAdded = processWeightObject (weightObject, objectAddWeight, minExonNum, countConsecutive, pairedEndOnlyMode, countTable, excludeAmbiguousReads,
					     genomicIntervalMode, checkWeights);
	  sumObjectWeight = sumObjectWeight + weightAdded;
	  if (weightAdded > 0) {
	    numWeightObjectsIncluded++;
	  } else {
	    numWeightObjectsExcluded++;
	  }
	} else {
	  numWeightObjectsExcluded++;
	}
      }

      if (lineNumber >= countUnit) {
	System.err.println();
      }

      System.err.println(lineNumber + " lines read.");
      
      reader.close();

      if (debugLevel >= 1) {
	if (debugLevel >= 3) {
	  System.err.println ("count table final: " + countTable);
	}
	System.err.println ("Number of columns for file " + countObjectMapFilename + ": " + Exon.getCountObjectMapColumn ());
      }

      if (countObjectFilename.equals ("")) {

	if (Exon.getCountObjectMapColumn () == 2) {
	  TreeSet<String> countObjectIds = null;
	  if (countObjectMapFilename != ""  && Exon.getCountObjectIds() != null) {
	    if (Exon.getCountObjectMapColumn () == 2)
	      countObjectIds = new TreeSet<String> (Exon.getCountObjectIds().keySet ());	
	  } else {
	    countObjectIds = new TreeSet<String> (countTable.keySet ());
	  }
	
	  for (String countObjectId: countObjectIds) {
	    if (specialCountObjectIdSet.contains(countObjectId)) {
	      System.err.println ("Count for " + countObjectId + ": " + countTable.get(countObjectId));
	    }
	    if (countTable.get(countObjectId) != null) {
	      if (outputZeroes || countTable.get(countObjectId).doubleValue() != 0) {
		outputWriter.println(countObjectId + "\t" + countTable.get(countObjectId).doubleValue());
	      }
	    } else if (outputZeroes) {
	      outputWriter.println(countObjectId + "\t" + 0);
	    }
	  }
	} else {
	  /* The count object ids of the third column of countObjectMapFilename correspond to
	     external count objects ids in the second column which are going to be used for output */
	  countObjectMapReader = UtilLib.getBufferedReader (countObjectMapFilename);
	  
	  countUnit = 500 * 1000;
	  lineNumber = 0;
	  line = countObjectMapReader.readLine();
	  
	  int numCountIds = 0;
	  int numCountIdsQuantified = 0;
	  while (line != null) {
	    st = new StringTokenizer (line, "\t");
	    if (! st.hasMoreTokens()) {
	      throw new IOException ("No exon id found.");
	    }
	    String exonId = st.nextToken ();
	
	    if (! exonId.equals("Exon Id")) {
	      if (! st.hasMoreTokens()) {
		throw new IOException ("No second column found for exon " + exonId + ".");
	      }
	      String outputCountObjectId = st.nextToken ();

	      if (! st.hasMoreTokens()) {
		throw new IOException ("No third column found for exon " + exonId + ".");
	      }
	      String countObjectId = st.nextToken ();
	      
	      numCountIds++;
	      if (countTable.get(countObjectId) != null) {
		if (outputZeroes || countTable.get(countObjectId).doubleValue() != 0) {
		  outputWriter.println(outputCountObjectId + "\t" + countTable.get(countObjectId).doubleValue());
		}
		numCountIdsQuantified++;
	      } else if (outputZeroes) {
		outputWriter.println(outputCountObjectId + "\t" + 0.0);
	      }
	    }

	    if (lineNumber % countUnit == 0) {
	      System.err.print (".");
	    }
	    	  
	    line = countObjectMapReader.readLine();
	    lineNumber++;
	  }
	  if (lineNumber > countUnit) {
	    System.err.println (".");
	  }
	  
	  countObjectMapReader.close ();
	}
      } else {

	System.err.println ("Reading count ids from file " + countObjectFilename);
	BufferedReader countObjectReader = UtilLib.getBufferedReader (countObjectFilename);
	int numCountIds = 0;
	int numCountIdsQuantified = 0;
	line = countObjectReader.readLine();
	while (line != null) {
	  int tabPos = line.indexOf ("\t");
	  String countObjectId = "";
	  if (tabPos == -1) {
	    countObjectId = line;
	  } else {
	    countObjectId = line.substring (0, tabPos);
	  }
	  numCountIds++;
	  
	  if (countTable.get(countObjectId) != null) {
	    if (outputZeroes || countTable.get(countObjectId).doubleValue() != 0) {
	      outputWriter.println(countObjectId + "\t" + countTable.get(countObjectId).doubleValue());
	    }
	    numCountIdsQuantified++;
	  } else if (outputZeroes) {
	    outputWriter.println(countObjectId + "\t" + 0.0);
	  }
	  line = countObjectReader.readLine();
	}
	countObjectReader.close ();
	
	System.err.println (numCountIds + " count ids written to file " + outputFilename);
	System.err.println ("of which " + numCountIdsQuantified + " were quantified.");
      }

      outputWriter.close ();

      System.err.format ("Num weight objects in BED file: " + numWeightObjects + ", num weight objects included: " + numWeightObjectsIncluded +
			 " (%.2f%%), num weight objects excluded: " + numWeightObjectsExcluded + " (%.2f%%).%n", numWeightObjectsIncluded * 100.0 / numWeightObjects,
			 numWeightObjectsExcluded * 100.0 / numWeightObjects);
      System.err.format ("Sum of weights contributing to count objects: %.2f%n", sumObjectWeight);

      if (UtilLib.getCountMode().equals ("gene")) {
	System.err.println("NUMBER_EXPRESSED_READS=" + numWeightObjectsIncluded);
      }

    }
    catch (Exception e) {
      System.err.println ("Problem in line: " + line + ": " + e==null?"No error message":e.getMessage());
      System.exit (1);
    }
  }
}

