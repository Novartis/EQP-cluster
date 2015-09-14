/**File: ComputeExonCountsBed.java 

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
 *                           Class ComputeExonCountsBed
 *
 *  
 ***********************************************************************************/


public class ComputeExonCountsBed {

  private static int debugLevel = 0;

  private static Hashtable<String, Hashtable<GenomicAlignment, Integer>> overlapTable1 = null;
  private static Hashtable<String, Hashtable<GenomicAlignment, Integer>> overlapTable2 = null;
  
  private static Hashtable<String, Double>  countTable = null;
  private static Hashtable<String, Integer> countObjectLengthTable = null;
  private static HashSetTable<String, Exon> countObjectExonTable   = null;

  private static Hashtable<String, Double> readWeightTable = new  Hashtable<String, Double> (50 * 1000 * 1000);

  private static String previousFragmentId = "";
  

  /***********************************************************************************
   *
   *                       Compute count object lengths  
   *
   ***********************************************************************************/

  public static void computeCountObjectLengths () {

    HashSetTable<Exon, String> exonCountObjectTable = Exon.getCountObjectTable ();

    countObjectExonTable = new HashSetTable<String, Exon> (2 * exonCountObjectTable.size());
    for (Exon exon: exonCountObjectTable.keySet ()) {
      for (String countObjectId: exonCountObjectTable.getSet (exon)) {
	countObjectExonTable.putValue (countObjectId, exon);
      }
    }

    countObjectLengthTable = new Hashtable<String, Integer> (countObjectExonTable.size ());

    for (String countObjectId: countObjectExonTable.keySet ()) {
      int countObjectLength = Exon.computeExonSetLength (countObjectExonTable.getSet (countObjectId));
      countObjectLengthTable.put (countObjectId, new Integer (countObjectLength));
    }
    
  }


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
   ***********************************************************************************/

  private static double getReadWeight (String fragmentId, BufferedReader weightReader, boolean useReadWeights) throws IOException {

    if (weightReader == null) {
      return 1.0;
    }

    if (! useReadWeights) {
      return 1.0;
    }

    /* Store the last used fragment id and delete it if it changes; note that getReadWeight can also be
       called if the read id changes */
    if (! previousFragmentId.equals(fragmentId) && previousFragmentId != "") {
      readWeightTable.remove (previousFragmentId);
    }
    previousFragmentId = fragmentId;

    Double readWeight = readWeightTable.get (fragmentId);
    if (readWeight != null) {
      return 1.0 / readWeight.doubleValue ();
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
    
    if (st.hasMoreTokens()) {
      readWeight = new Double(st.nextToken ());	
    } else {
      throw new IOException ("No weight for read " + weightFragmentId + " found.");
    }

    while (! weightFragmentId.equals (fragmentId)) {

      readWeightTable.put (weightFragmentId, readWeight);

      line = weightReader.readLine ();

      if (line == null) {
	return 1.0;
      }

      st = new StringTokenizer (line, "\t");
      
      if (st.hasMoreTokens()) {
	weightFragmentId = st.nextToken ();	
      }

      if (st.hasMoreTokens()) {
	readWeight = new Double(st.nextToken ());	
      } else {
	throw new IOException ("No weight for read " + weightFragmentId + " found.");
      }

    }

    readWeightTable.put (weightFragmentId, readWeight);

    if (debugLevel >= 1) {
      System.out.println ("Weight: " + readWeight + " from " + line + " for " + fragmentId);
    }

    return 1.0 / readWeight;

  }



  /***********************************************************************************
   * 
   *                          processWeightObject
   *
   ***********************************************************************************/

  private static void processWeightObject (WeightObject weightObject, double objectWeight, int globalOverlapThreshold) throws IOException {


    HashSet<String> countObjectIdSet = weightObject.getCountObjectIds (globalOverlapThreshold, countObjectLengthTable);
    
    if (debugLevel >= 1) {
      System.out.println("Processing weight object: " + weightObject + " (weight: " + objectWeight + ") with count object set: " + countObjectIdSet);
    }
    
    
    for (String countObjectId: countObjectIdSet) {

      addToTableEntry (countTable, countObjectId, objectWeight);
	      
    }
    
  }

 

  /***********************************************************************************/

   private static void printHelp () {
    System.out.println("ComputeExonCountsBed\n" +                                              
    "USAGE: ComputeExonCountsBed [-p] [-O <overlap thresh.>] [-R <read weight thresh.>]\n" +
    "   [-g] [-C <count object map file>] [-W <read weight file>|none] \n" +
    "     -b <intersect. bed file> -r <read length> -o <outputFile>\n" +
    "\n" +
    " -b STRING: intersect. bed file - the bed file containing the intersection of\n" +
    "     the exons transcript intervals and the reads mapped to the transcripts.\n" +
    "     (- for STDIN) [-]\n" +
    " -r INT: read length - the length of the reads that were aligned\n" +
    " -o STRING: output file - the file to which the output is written (- for STDOUT)\n" + 
    "    [-]\n" +
    " -g: Compute gene counts: set paired-end mode (option -p) and the minimal\n" +
    "    overlap to the <read length> (option -o) and consider non-genomic exons\n" +
    "    as well (\"tr-\" exons)\n" +
    " -O INT: Minimal overlap of a read with a count object (for paired-end\n" +
    "    alignments twice the overlap is required.) [8]\n" +
    " -R FLOAT: Minimal weight of a read; reads with a lower weight are disregarded \n" +
    "    [0.01]. If it is negative, then the reads are not weighted.\n" +
    " -p: paired end mode - if both ends of a fragment align to exons which belong\n" +
    "    to a collection of exons C, then *both* alignments must have a length of \n" +
    "    at least the overlap threshold.\n" +
    " -W read weight file: a file containing the read weights with two columns:\n" +
    "      <fragment id>, <read weight>. Reads are weighted by the minimum of\n" +
    "      the read weight contained in this file,if present, and the read weight\n" +
    "      computed from the number of genomic alignments in the bed file. If the\n" +
    "      keyword \"none\" is selected, then the reads are not weighted.\n" +
    " -m count object map file: contains the exon to count object mapping. It\n" +
    "    consists of two columns, the exon id (<gene id>/<chr>/<start>/\n" +
    "    <end>/<strand>) and the id of the count object.\n" +
    "\n" +
    "Example: ComputeExonCountsBed -s sma-sample1-refseq-mm9.bed -r 51\n" + 
    "  -o sma-sample1-refseq-mm9.cnt\n");


  }
                                     

  /***********************************************************************************/

  public static void main (String [] args) {

    int numDeletions  = 0;
    int numInsertions = 0;

    String separator = "#";
    
    boolean geneCountMode = false;
    
    String intersectionFilename = "-";
    String weightFilename = "none";
    String countObjectMapFilename = "";
    String outputFilename = "-";

    double readWeightThreshold = 0.01;
    
    int  readLength = -1;
    int  overlapThreshold = 3;

    String oldAlignmentId = "";
    String oldFragmentId = "";
    String oldReadId = "";

    boolean oldIsPairedEnd = false;

    WeightObjectAlignment weightObjectAlignment = null;
    WeightObject          weightObject = null;

    boolean printLines    = false;
    boolean pairedEndMode = false;
    boolean useReadWeights = true;
    boolean strandSpecificProcessing = false;

    final int countUnit = 2 * 1250 * 1000;
    
    Getopt g = new Getopt("ComputeExonCountsBed", args, "b:d:gL:m:o:O:psw:W:h");
    
    int c;
    String arg = "";
    
    c = g.getopt();
    
    while (c  != -1) {
      switch(c) {
      case 'b':
	intersectionFilename = g.getOptarg();
	break;
      case 'd':
	UtilLib.setDebugLevel (Integer.parseInt(g.getOptarg()));
	break;
      case 'g':
	geneCountMode = true;
	break;
      case 'L':
	readLength = Integer.parseInt(g.getOptarg());
	break;
      case 'm':
	countObjectMapFilename = g.getOptarg();
	break;
      case 'o':
	outputFilename = g.getOptarg();
	break;
      case 'O':
	overlapThreshold = Integer.parseInt(g.getOptarg());
	break;
      case 'p':
	pairedEndMode = true;
	break;
      case 's':
	strandSpecificProcessing = true;
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
	System.out.print("Error: getopt() returned " + c + "\n");
      }
      c = g.getopt();
    }

    debugLevel = UtilLib.getDebugLevel ();

    if (readLength <= 0) {
      printHelp();
      System.exit(0);	
    }

    String line = "";
    try {
      
      BufferedReader reader       = UtilLib.getBufferedReader (intersectionFilename);
      PrintWriter    outputWriter = UtilLib.getPrintWriter    (outputFilename);

      BufferedReader countObjectReader = UtilLib.getBufferedReader (countObjectMapFilename);
      System.err.println("Loading exon count object file " + countObjectMapFilename);
      Exon.loadCountObjectFile (countObjectReader);
      computeCountObjectLengths ();
      System.err.println("Done.");

      if (readWeightThreshold < 0 || weightFilename.equals("none")) {
	useReadWeights = false;
      }
     
      BufferedReader weightReader = null;
      if (weightFilename != "" && useReadWeights) {
	weightReader = UtilLib.getBufferedReader (weightFilename);
      }

      countTable    = new Hashtable<String, Double> (400 * 1000);

      line = reader.readLine();
      int lineNumber = 1;
      Double objectWeightDouble = null;
      double objectWeight = 1.0;

      System.err.println("Reading bed file: " + intersectionFilename + " (. = " + countUnit + " lines)");

      while (line != null) {
	
	BedRecord bedRecord = new BedRecord (line);

	if (bedRecord.getOverlap () > 0) {
	
	  /* We consider fragments as weight objects. If there are two single-read alignments of one fragment
	     against a count object, then we count this as one. */
	  String  fragmentId  = bedRecord.getFragmentId ();
	  String  readId      = bedRecord.getReadId ();
	  boolean isPairedEnd = bedRecord.isPairedEnd ();

	  if (debugLevel >= 2) {
	    System.out.println ("BedRecord Weight object id: " + fragmentId);
	  }

	  if (printLines) {
	    System.out.println (bedRecord.toString());
	  }

	  if (weightObject == null) {
	    
	    weightObject = new WeightObject (bedRecord);
	  
	  } else if ((! fragmentId.equals(oldFragmentId) && oldIsPairedEnd) || (! oldIsPairedEnd && ! readId.equals(oldReadId))){

	    objectWeight = getReadWeight (oldFragmentId, weightReader, useReadWeights);
	  
	    if (debugLevel >= 1) {
	      System.out.println ("Weight for " + oldFragmentId + ": " + objectWeight);
	    }

	    if (objectWeight >= readWeightThreshold) {
	      processWeightObject (weightObject, objectWeight, overlapThreshold);
	    }

	    /* Create new weight object */
	    weightObject = new WeightObject (bedRecord);
	  
	  } else {
	  
	    weightObject.addBedRecord (bedRecord);
	  
	  }

	  oldFragmentId  = fragmentId;
	  oldReadId      = readId;
	  oldIsPairedEnd = isPairedEnd;
	}
	  
	if (lineNumber % countUnit == 0) {
	  System.err.print(".");
	}
	lineNumber++;
	  
	
	line = reader.readLine();

	if (debugLevel >= 2) {
	  System.out.println ("Line: " + line);
	  System.out.println ("count table: " + countTable);
	}
	
      }

      if (debugLevel >= 1) {
	System.out.println ("Processing last fragment");
      }

      /* process last fragment */
      if (weightObject != null) {
	objectWeight = getReadWeight (oldFragmentId, weightReader, useReadWeights);
	
	if (objectWeight >= readWeightThreshold) {
	  processWeightObject (weightObject, objectWeight, overlapThreshold);
	}
      }

      if (lineNumber >= countUnit) {
	System.err.println();
      }
      
      reader.close();

      if (debugLevel >= 1) {
	System.out.println ("count table final: " + countTable);
      }

      TreeSet<String> countObjectIds = new TreeSet<String> (countTable.keySet ());
      if (countObjectMapFilename != "" && Exon.getCountObjectIds() != null) {
	countObjectIds = new TreeSet<String> (Exon.getCountObjectIds().keySet ());	
      }

      for (String countObjectId: countObjectIds) {
	if (debugLevel >= 1) {
	  System.out.println ("Count object id: " + countObjectId);
	}
	if (countTable.get(countObjectId) != null) {
	  outputWriter.println(countObjectId + "\t" + countTable.get(countObjectId).doubleValue());
	} else {
	  outputWriter.println(countObjectId + "\t" + 0);
	}
      }

      outputWriter.close ();

    }
    catch (Exception e) {
      System.out.println ("Problem in line: " + line + ": " + e==null?"No error message":e.getMessage());
    }
  }
}

