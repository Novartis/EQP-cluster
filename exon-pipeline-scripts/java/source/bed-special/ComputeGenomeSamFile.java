/**File: ComputeGenomeSamFile.java 

Original Author: Sven Schuierer
Date: 08/10/2013

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
 *                           Class ComputeGenomeSamFile
 *
 *  
 ***********************************************************************************/

public class ComputeGenomeSamFile {

  private static int debugLevel = 0;

  private static String bedLine          = "";
  private static String oldBedFragmentId = "";
  private static int bedLineNumber       = 1;
  
  private static WeightObject classWeightObject = null;
  
  private static String combinedSamLine = "";
  private static String oldCombinedSamFragmentName = "";
  private static int    combinedSamLineNumber = 1;

  private static Vector<SamRecord> classCombinedSamRecords = null;
  private static Vector<SamRecord> classGenomeSamRecords   = null;
  
  private static final int countUnit = 1 * 1000 * 1000;



  /***********************************************************************************/

  private static WeightObject readNextBedRecord (BufferedReader bedReader, HashSet<String> chromosomeIds) throws Exception {

    if (bedLine == null) {
      return null;
    }

    while (bedLine != null) {
      BedRecord bedRecord = new BedRecord (bedLine);
	
      if (bedRecord.getOverlap () > 0 && ! chromosomeIds.contains(bedRecord.getReferenceSequenceId ())) {
	  
	String  fragmentId = bedRecord.getFragmentId ();
	  
	if (classWeightObject == null) {	    
	  classWeightObject = new WeightObject (bedRecord);	    
	} else if (! fragmentId.equals(oldBedFragmentId)) {	    
	  WeightObject oldWeightObject = classWeightObject;
	  classWeightObject = null;
	  oldBedFragmentId = fragmentId;
	  return oldWeightObject;	  
	} else {
	  classWeightObject.addBedRecord (bedRecord);
	}
	  
	oldBedFragmentId = fragmentId;
      }
	
      if (bedLineNumber % countUnit == 0) {
	/* System.err.print("."); */
      }
      bedLineNumber++;
	
      bedLine = bedReader.readLine();
	
      if (debugLevel >= 1) {
	System.err.println ("BedLine: " + bedLine);
      }
	
    }
      
    return classWeightObject;
    
    
  }
  

  /***********************************************************************************/

  private static Vector<SamRecord> readNextCombinedSamRecord (BufferedReader samReader, String samFilename) throws IOException {

    if (combinedSamLine == null) {
      return null;
    }
    
    while (combinedSamLine != null) {

      if (! combinedSamLine.startsWith("@") && ! combinedSamLine.equals("")) {
	
	SamRecord samRecord    = new SamRecord (combinedSamLine, combinedSamLineNumber, samFilename);
	String    queryName    = samRecord.getQueryName();
	String    fragmentName = samRecord.getFragmentName();

	if (classCombinedSamRecords == null) {
	  classCombinedSamRecords = new Vector<SamRecord> (500);
	  classCombinedSamRecords.add(samRecord);
	} else if (! fragmentName.equals(oldCombinedSamFragmentName)) {
	  Vector<SamRecord> oldSamRecords = classCombinedSamRecords;
	  classCombinedSamRecords = null;
	  oldCombinedSamFragmentName = fragmentName;
	  return oldSamRecords;
	} else {
	  classCombinedSamRecords.add(samRecord);
	}
		
	oldCombinedSamFragmentName = fragmentName;
	  
	if (combinedSamLineNumber % countUnit == 0) {
	  /* System.err.print("."); */
	}	  
	combinedSamLineNumber++;

      }
	
      combinedSamLine = samReader.readLine();
	
    }

    return classCombinedSamRecords;
    
  }
  

  /***********************************************************************************/

  private static String readSamHeader (BufferedReader headerReader) throws IOException {

    String samHeader = "";
    String headerLine = headerReader.readLine();
    while (headerLine != null) {

      if (headerLine.startsWith("@")) {
	if (headerLine.startsWith("@PG")) {
	  headerLine = "@PG\tID:EQP\tVN:1.0";
	}
	if (samHeader.equals("")) {
	  samHeader = headerLine;
	} else {
	  samHeader = samHeader + "\n" + headerLine;
	}
      }

      headerLine = headerReader.readLine();
	
    }

    return samHeader;
    
  }

  
  /***********************************************************************************/

  private static HashSet<String> readChromosomeIds (BufferedReader chromosomeIdReader) throws IOException {

    HashSet<String> chromosomeIds = new HashSet<String> (200);
    String chromosomeIdLine = chromosomeIdReader.readLine();
    while (chromosomeIdLine != null) {

      StringTokenizer st = new StringTokenizer (chromosomeIdLine, "\t");
      if (st.hasMoreTokens()) {
	chromosomeIds.add(st.nextToken());
      }
      chromosomeIdLine = chromosomeIdReader.readLine();
	
    }

    UtilLib.setChromosomeIds (chromosomeIds);

    return chromosomeIds;
    
  }


  /***********************************************************************************/

  private static Vector<SamRecordPair>[] sortSamRecordPairsByCategory (Vector<SamRecordPair> samRecordPairVector, boolean pairedEndAlignments) throws IOException {
    
    /* We divide a vector of non-paired end sets of SAM records into contiguous stretches of alignments for
       the first read and for the second read */

    @SuppressWarnings("unchecked")
    Vector <SamRecordPair> [] samRecordPairVectorArray = (Vector <SamRecordPair>[]) new Vector [3];
    
    int allocationSize = 1;
        
    Vector <SamRecordPair> samRecordPairVectorPaired = new Vector <SamRecordPair> (allocationSize);
    Vector <SamRecordPair> samRecordPairVector1 = new Vector <SamRecordPair> (allocationSize);
    Vector <SamRecordPair> samRecordPairVector2 = new Vector <SamRecordPair> (allocationSize);

    for (SamRecordPair samRecordPair: samRecordPairVector) {
      if (debugLevel >= 2) {
	System.err.println("Sorting SAM record pair: " + samRecordPair);
      }
      if (samRecordPair.isPairedEndAlignment () || ! pairedEndAlignments) {
	if (debugLevel >= 2) {
	  System.err.println("Adding SAM record pair to paired end.");
	}
	samRecordPairVectorPaired.add (samRecordPair);
      } else if (samRecordPair.hasFirstRead ()) {
	if (debugLevel >= 2) {
	  System.err.println("Adding SAM record pair to first read set.");
	}
	samRecordPairVector1.add(samRecordPair);
      } else if (samRecordPair.hasSecondRead ()) {
	if (debugLevel >= 2) {
	  System.err.println("Adding SAM record pair to second read set.");
	}
	samRecordPairVector2.add(samRecordPair);
      }
    }

    samRecordPairVectorArray[0] = samRecordPairVectorPaired;
    samRecordPairVectorArray[1] = samRecordPairVector1;
    samRecordPairVectorArray[2] = samRecordPairVector2;
    
    return samRecordPairVectorArray;    
  }


  /***********************************************************************************/

  private static SamRecordPair createDefaultSamRecordPair (Vector<SamRecord> samRecords, boolean pairedEndAlignments) throws IOException {

    SamRecord firstRead  = null;
    SamRecord secondRead = null;

    SamRecordPair defaultSamPair = new SamRecordPair ();
    for (SamRecord samRecord: samRecords) {
      if (firstRead == null && samRecord.getReadIndex () == 1) {
	firstRead = samRecord;
      }
      
      if (secondRead == null && samRecord.getReadIndex () == 2)  {
	secondRead = samRecord;
      }

      if (firstRead != null && (! pairedEndAlignments || secondRead != null)) {
	break;
      }

    }

    if ((pairedEndAlignments && (firstRead == null || secondRead == null)) || (! pairedEndAlignments && firstRead == null)) {
      throw new IOException ("ERROR: first or second read missing from SAM records: " + samRecords);
    }

    defaultSamPair.setFirstRead  (firstRead);
    defaultSamPair.setSecondRead (secondRead);

    return defaultSamPair;

  }


  
  /***********************************************************************************/

  private static Vector<SamRecordPair>[] createSamRecordPairs (Vector<SamRecord> samRecords, boolean pairedEndAlignments, HashSet<String> chromosomeIds) throws IOException {

    if (debugLevel >= 2) {
      System.err.println("Creating SAM record pairs for " + samRecords.size() + " SAM records: " + samRecords);
    }

    HashSet<SamRecordPair> samRecordPairSet    = new HashSet<SamRecordPair> (samRecords.size());
    Vector <SamRecordPair> samRecordPairVector = new Vector<SamRecordPair> (samRecords.size());

    HashSet<Integer> processedIndices = new HashSet<Integer> (2 * samRecords.size() + 1);
    for (int i = 0; i < samRecords.size(); i++) {

      Integer iInt = Integer.valueOf (i);
      if (! processedIndices.contains(iInt)) {

	processedIndices.add(iInt);

	SamRecord samRecord = samRecords.get(i);
	if (debugLevel >= 2) {
	  System.err.println("Processing alignment to: " + samRecord.getReferenceName () + "; is a chromosome: " + chromosomeIds.contains(samRecord.getReferenceName ()));
	}
	if (chromosomeIds.contains(samRecord.getReferenceName ())) {
	  SamRecord mateSamRecord = null;
	  if (samRecord.isMapped () && samRecord.mateIsMapped ()) {
	    mateSamRecord = samRecord.findMate(samRecords, i, processedIndices);

	    if (mateSamRecord != null && ! chromosomeIds.contains(mateSamRecord.getReferenceName ())) {
	      mateSamRecord = null;
	    }

	    if (debugLevel >= 2) {
	      System.err.println("Mate for: " + samRecord + " is " + mateSamRecord);
	    }
	  } 

	  SamRecordPair samRecordPair = null;
	  if (mateSamRecord == null) {
	    if (samRecord.isMapped ()) {
	      samRecordPair = new SamRecordPair (samRecord);
	    }
	  } else {
	    samRecordPair = new SamRecordPair (samRecord, mateSamRecord);
	  }

	  if (samRecordPair != null) {
	    if (! samRecordPairSet.contains(samRecordPair)) {
	      samRecordPairSet.add(samRecordPair);
	      samRecordPairVector.add(samRecordPair);
	    }
	  }
	}
      }
    }

    if (debugLevel >= 1) {
      System.err.println(samRecordPairVector.size() + " SAM record pairs created.");
    }

    return sortSamRecordPairsByCategory (samRecordPairVector, pairedEndAlignments);
    
  }

  
  /***********************************************************************************/

  private static Vector<SamRecordPair> removeOverlapping (Vector<SamRecordPair> testSamRecordPairs, Vector <SamRecordPair> controlSamRecordPairs,
							  HashSet<SamRecordPair> excludedControlSamRecordPairs) throws IOException {

    if (debugLevel >= 1) {
      System.err.println ("Removing overlapping SAM records against " + controlSamRecordPairs);
    }

    Vector<SamRecordPair> retainedSamRecordPairs = new Vector<SamRecordPair> (testSamRecordPairs.size ());
    for (SamRecordPair testSamRecordPair: testSamRecordPairs) {

      if (debugLevel >= 2) {
	System.err.println ("Test SAM record pair: " + testSamRecordPair);
      }
      
      boolean addPair = true;
      if (testSamRecordPair.isUnmapped ()) {
	addPair = false;
      } else {
	for (SamRecordPair controlSamRecordPair: controlSamRecordPairs) {
	  /* Note that unmapped reads pairs overlap with any other read pair */
	  if (! excludedControlSamRecordPairs.contains(testSamRecordPair) && controlSamRecordPair.overlapsWithAndIsCloseTo (testSamRecordPair, 0.3, 0.93)) {
	    addPair = false;
	    break;
	  }
	}
      }

      if (addPair) {
	if (debugLevel >= 2) {
	  System.err.println ("Test SAM record pair added.");
	}
	retainedSamRecordPairs.add (testSamRecordPair);
      } else {
	if (debugLevel >= 2) {
	  System.err.println ("Test SAM record pair removed.");
	}

      }
    }

    testSamRecordPairs.clear ();
    testSamRecordPairs.addAll(retainedSamRecordPairs);

    if (debugLevel >= 1) {
      System.err.println ("Overlapping SAM records removed.");
    }


    return testSamRecordPairs;
  }
  
  /***********************************************************************************/

  private static void printGenomeSamRecordsByCategory (Vector <SamRecordPair> computedGenomeSamPairs, Vector<SamRecordPair> genomeSamRecordPairs, int level,
						       boolean [] alignmentsFound, boolean [] alignmentsPrinted, PrintWriter outputWriter) throws IOException {

    for (SamRecordPair genomeSamPair: computedGenomeSamPairs) {
      if (debugLevel >= 2) {
	System.err.println ("Printing genome SAM pair of computedGenomeSamPairs:");
	System.err.println ("First read: " + genomeSamPair.getFirstRead ());
	System.err.println ("Second read: " + genomeSamPair.getSecondRead ());
      }
      outputWriter.println(genomeSamPair.toString (alignmentsFound, alignmentsPrinted, level));
    }
    
    for (SamRecordPair genomeSamPair: genomeSamRecordPairs) {
      if (debugLevel >= 2) {
	System.err.println ("Printing genome sam pair:");
	System.err.println ("First read: " + genomeSamPair.getFirstRead ());
	System.err.println ("Second read: " + genomeSamPair.getSecondRead ());
	System.err.println ("alignmentsFound: " + Arrays.toString(alignmentsFound) + ", alignmentsPrinted: " + Arrays.toString(alignmentsPrinted) + ", level: " + level);
      }
      outputWriter.println(genomeSamPair.toString (alignmentsFound, alignmentsPrinted, level));
    }
    
  }

    
  /*****************************************************************************
   *
   *                           printGenomeSamRecords
   *
   *****************************************************************************/

  private static void printGenomeSamRecords (Vector <SamRecordPair> computedGenomeSamPairs, Vector <SamRecord> originalGenomeSamRecords, boolean pairedEndAlignments,
					      HashSet<String> chromosomeIds, PrintWriter outputWriter) throws IOException {

    Vector<SamRecordPair>[] originalGenomeSamRecordPairs = createSamRecordPairs(originalGenomeSamRecords, pairedEndAlignments, chromosomeIds);
    if (debugLevel >= 2) {
      System.err.println (originalGenomeSamRecordPairs[0].size() + " + " + originalGenomeSamRecordPairs[1].size() + " + " + originalGenomeSamRecordPairs[2].size() +
			  " genome SAM pairs created.");
    }
    
    HashSet<SamRecord>     genomeSamRecordSet     = new HashSet<SamRecord>     (originalGenomeSamRecords);
    HashSet<SamRecordPair> genomeSamRecordPairSet = new HashSet<SamRecordPair> (originalGenomeSamRecordPairs[0]);
    genomeSamRecordPairSet.addAll (originalGenomeSamRecordPairs[1]);
    genomeSamRecordPairSet.addAll (originalGenomeSamRecordPairs[2]);
    
    Vector <SamRecordPair>  genomeSamPairsRetained   = new Vector <SamRecordPair> (computedGenomeSamPairs.size());
    HashSet<SamRecordPair>  genomeSamPairSetExcluded = new HashSet<SamRecordPair> (computedGenomeSamPairs.size());

    /* Merge the computed genome SAM Pairs with the original genome SAM pairs:
       As a first step remove all genome SAM pairs that are already contained in the original genome alignments */
    for (SamRecordPair samRecordPair: computedGenomeSamPairs) {
      if (samRecordPair.isPairedEndAlignment () && genomeSamRecordSet.contains(samRecordPair.getFirstRead()) && genomeSamRecordSet.contains(samRecordPair.getSecondRead())) {
	if (debugLevel >= 2) {
	  System.err.println ("Excluding paired-end genome SAM pair: " + samRecordPair);
	}
	genomeSamPairSetExcluded.add(samRecordPair);
      } else if (! samRecordPair.isPairedEndAlignment () &&
		 ((samRecordPair.hasFirstRead () && genomeSamRecordSet.contains(samRecordPair.getFirstRead())) ||
		  (samRecordPair.hasSecondRead () && genomeSamRecordSet.contains(samRecordPair.getSecondRead())))) {
	if (debugLevel >= 2) {
	  System.err.println ("Excluding genome SAM pair: " + samRecordPair);
	}
	genomeSamPairSetExcluded.add(samRecordPair);
      } else {
	genomeSamPairsRetained.add(samRecordPair);
      }

    }

    Vector<SamRecordPair>[] genomeSamPairsArray = sortSamRecordPairsByCategory (genomeSamPairsRetained, pairedEndAlignments);
    boolean[] alignmentsFound  = {false, false};
    boolean pairedEndFound     = false;
    boolean pairedEndPreferred = true;
    
    /* Print genome SAM pairs and non-overlapping original genome alignments */
    for (int i = 0; i < 3; i++) {

      if (debugLevel >= 1) {
	System.err.println ("Before removeOverlapping level " + i + " number of genome SAM pairs: " + genomeSamPairsArray[i].size () + " and original genome pairs: " +
			    originalGenomeSamRecordPairs[i].size());
      }
      
      /* Compute first and second read and remove overlapping original genome alignments */
      if (i == 0) {
	originalGenomeSamRecordPairs[i] = removeOverlapping (originalGenomeSamRecordPairs[i], genomeSamPairsArray[i], genomeSamPairSetExcluded);
      } else {
	originalGenomeSamRecordPairs[i] = removeOverlapping (originalGenomeSamRecordPairs[i], genomeSamPairsRetained, genomeSamPairSetExcluded);
      }

      if (debugLevel >= 1) {
	System.err.println ("Level " + i + " number of genome SAM pairs: " + genomeSamPairsArray[i].size () + " and original genome pairs: " +
			    originalGenomeSamRecordPairs[i].size());
      }

      if (genomeSamPairsArray[i].size () + originalGenomeSamRecordPairs[i].size() > 0) {	
	if (i == 0) {
	  alignmentsFound[0] = true;
	  alignmentsFound[1] = true;
	  pairedEndFound     = true;
	} else {
	  alignmentsFound[i - 1] = true;
	}
      }
    }

    /* Count the number of alignments and add NH field */
    int numberOfAlignments = 0;
    for (int i = 0; i < 3; i++) {
      numberOfAlignments += genomeSamPairsArray[i].size() + originalGenomeSamRecordPairs[i].size();
      if (pairedEndFound && pairedEndPreferred) {
	break;
      }
    }

    for (int i = 0; i < 3; i++) {
      for (SamRecordPair genomeSamPair: genomeSamPairsArray[i]) {
	genomeSamPair.addOptionalField("NH", "i", String.valueOf (numberOfAlignments));
      }
      for (SamRecordPair genomeSamPair: originalGenomeSamRecordPairs[i]) {
	genomeSamPair.addOptionalField("NH", "i", String.valueOf (numberOfAlignments));
      }
      if (pairedEndFound && pairedEndPreferred) {
	break;
      }
    }
    
    /* alignmentsPrinted is needed to keep track of primary and secondary alignments */
    boolean[] alignmentsPrinted = {false, false};
    for (int i = 0; i < 3; i++) {
      if (debugLevel >= 1) {
	System.err.println ("Printing genome SAM reads for level " + i);
      }
      printGenomeSamRecordsByCategory (genomeSamPairsArray[i], originalGenomeSamRecordPairs[i], i, alignmentsFound, alignmentsPrinted, outputWriter);
      if (pairedEndFound && pairedEndPreferred) {
	break;
      }
    }

    
    if (debugLevel >= 2) {
      System.err.println ("alignmentsPrinted: (" + alignmentsPrinted[0] + ", " + alignmentsPrinted[1] + "), alignmentsFound: (" +
			  alignmentsFound[0] +  ", " + alignmentsPrinted[1] + ")");
    }

    if (! alignmentsPrinted[0] || (! alignmentsPrinted[1] && pairedEndAlignments)) {

      if (debugLevel >= 2) {
	System.err.println ("Creating default SAM records from originalGenomeSamRecords: " + originalGenomeSamRecords);
      }

      SamRecordPair defaultSamRecordPair = createDefaultSamRecordPair (originalGenomeSamRecords, pairedEndAlignments);
      SamRecord firstRead  = defaultSamRecordPair.getFirstRead ();
      SamRecord secondRead = defaultSamRecordPair.getSecondRead ();
    
      if (debugLevel >= 2) {
	System.err.println ("firstRead: " + firstRead + ", secondRead: " + secondRead);
      }

      if (firstRead == null) {
	throw new IOException ("No first read found for: " + secondRead);
      }
      
      if (pairedEndAlignments && secondRead == null) {
	throw new IOException ("No second read found for: " + firstRead);
      }
        
      if (! alignmentsPrinted[0]) {
	outputWriter.println(firstRead.toUnmappedString (alignmentsPrinted[1])); 
      }

      if (pairedEndAlignments && ! alignmentsPrinted[1]) {
	outputWriter.println(secondRead.toUnmappedString (alignmentsPrinted[0])); 
      }
    }

    if (debugLevel >= 2) {
      System.err.println("printGenomeSamRecords done.");
    }
  }
  
  /***********************************************************************************/

  private static void printHelp () {
    System.err.println("ComputeGenomeSamFile\n" +                                              
    "USAGE: ComputeGenomeSamFile -b <intersect. bed file> -s <Combined SAM file>\n" +
    "   -c <chromosome id file> -H <SAM header file> -o <outputFile>\n" +
    "\n" +
    " -b STRING: intersect. bed file - the bed file containing the intersection of\n" +
    "     the exons transcript intervals and the reads mapped to the transcripts.\n" +
    "     (- for STDIN) [-]\n" +
    " -s STRING: file containing the best SAM entries for each read chosen between\n" +
    "     the alignments agains the transcripts, junctions, and genome\n" +
    " -c STRING: file containing all possible chromosome identifiers in the first\n" +
    "     column (can be a fai file).\n" +
    " -H STRING: file with genome SAM header (used to create a header for the output\n" +
    "     file).\n" +
    " -o STRING: output file - the file to which the output is written (- for STDOUT)\n" + 
    "    [-]\n" +
    " -S: alignments are single read alignments\n" +
    "\n");
  }


  /***********************************************************************************/

  public static void main (String [] args) {

    String intersectionBedFilename = "-";
    String combinedSamFilename     = "";
    String genomeSamFilename       = "";
    String chromosomeIdFilename    = "";
    String headerFilename          = "";
    String outputFilename          = "-";
    boolean quiet = false;
    boolean singleReadAlignments = false;
    
    Getopt g = new Getopt("ComputeGenomeSamFile.java", args, "b:c:d:g:H:o:qs:SWh");
    
    int c;
    String arg = "";
    
    c = g.getopt();
    
    while (c  != -1) {
      switch(c) {
      case 'b':
	intersectionBedFilename = g.getOptarg();
	break;
      case 'c':
	chromosomeIdFilename = g.getOptarg();
	break;
      case 'd':
	debugLevel = Integer.parseInt(g.getOptarg());
	UtilLib.setDebugLevel (debugLevel);
	break;
      case 'g':
	genomeSamFilename = g.getOptarg();
	break;
      case 'H':
	headerFilename = g.getOptarg();
	break;
      case 'o':
	outputFilename = g.getOptarg();
	break;
      case 'q':
	quiet = true;
	break;
      case 's':
	combinedSamFilename = g.getOptarg();
	break;
      case 'S':
	singleReadAlignments = false;
      case 'W':
	UtilLib.setWarningsOn (true);
	break;	
      case 'h':
	printHelp();
	System.exit(0);
	break;
      default:
	System.err.print("Error: getopt() returned " + c + "\n");
      }
      c = g.getopt();
    }

    /* We consider all alignment, splice conforming or not */
    UtilLib.setSpliceConformingCountMode (false);

    /* We are interested in complete reads, not only the parts of the reads that map to exons */
    UtilLib.setCountMode ("gene");

    debugLevel = UtilLib.getDebugLevel ();
    int lineNumber = 0;
    try {
      
      BufferedReader chromosomeIdReader = UtilLib.getBufferedReader (chromosomeIdFilename, true, "Chromosome Id file");
      BufferedReader bedReader          = UtilLib.getBufferedReader (intersectionBedFilename, true, "Intersection BED file");
      BufferedReader combinedSamReader  = UtilLib.getBufferedReader (combinedSamFilename, true, "Combined SAM file");
      PrintWriter    outputWriter       = UtilLib.getPrintWriter    (outputFilename);

      String samHeader = "";
      if (! headerFilename.equals("")) {
	BufferedReader headerReader = UtilLib.getBufferedReader (headerFilename);
	System.err.println("Reading SAM header from " + (headerFilename.equals("-")?"stdin":"file " + headerFilename));
	samHeader = readSamHeader (headerReader);
	headerReader.close();
      }
	
      System.err.println("Reading chromosome ids: " + chromosomeIdFilename);
      HashSet<String> chromosomeIds = readChromosomeIds (chromosomeIdReader);
      chromosomeIdReader.close ();
      
      System.err.println("Reading BED file: " + intersectionBedFilename);
      bedLine = bedReader.readLine();

      System.err.println("Reading combined SAM file: " + combinedSamFilename);
      combinedSamLine = combinedSamReader.readLine();
      
      System.err.println("Writing to genome SAM file: " + outputFilename + (quiet?"":" (. = " + countUnit + " reads)"));
      outputWriter.println(samHeader);
      outputWriter.flush();

      WeightObject weightObject = readNextBedRecord (bedReader, chromosomeIds);
      int numBedRecords = 1;
      int numFragments  = 0;
      
      Vector<SamRecord> combinedSamRecords = readNextCombinedSamRecord (combinedSamReader, combinedSamFilename);
      boolean pairedEndAlignments = combinedSamRecords.get(0).isPairedInSequencing ();
      if (debugLevel >= 1) {
	System.err.println("pairedEndAlignments: " + pairedEndAlignments);
      }
      
      while (weightObject != null || combinedSamRecords != null) {
	/* Note that bedLine and combinedSamLine can be null at this point */

	if (debugLevel >= 1) {
	  System.err.println ("Number of combined SAM records: " + combinedSamRecords.size ());
	  System.err.flush();
	}
	
	String combinedSamFragmentName = combinedSamRecords.get(0).getFragmentName();
	numFragments++;
	
	if (debugLevel >= 1) {
	  System.err.println ("SAM fragment: " + combinedSamFragmentName + ", BED fragment: " + weightObject.getFragmentId ());
	  System.err.flush();
	}

	Vector <SamRecordPair> computedGenomeSamPairs = new Vector <SamRecordPair> ();;	
	if (weightObject != null && weightObject.getFragmentId ().compareTo(combinedSamFragmentName) == 0) {
	  
	  /* Compute the genome SAM read for the current read that is contained in the BED file */
	  computedGenomeSamPairs = weightObject.computeGenomeSamPairs (combinedSamRecords, chromosomeIds);
	  
	  if (debugLevel >= 1) {
	    System.err.println (computedGenomeSamPairs.size () + " SAM record pairs computed from the BED entries for fragment: " + combinedSamFragmentName);
	    for (SamRecordPair samRecordPair: computedGenomeSamPairs) {
	      System.err.println ("First read: " + samRecordPair.getFirstRead ());
	      System.err.println ("Second read: " + samRecordPair.getSecondRead ());
	    }
	  }

	  weightObject = readNextBedRecord (bedReader, chromosomeIds);
	  numBedRecords++;

	} else if (weightObject != null && weightObject.getFragmentId ().compareTo(combinedSamFragmentName) < 0) {
	  throw new IOException ("BED fragment: " + weightObject.getFragmentId () + " without a SAM Record.\n" +
				 "Next SAM Record: " + combinedSamFragmentName);
	}

	if (debugLevel >= 1) {
	  System.err.println ("Print at most " + combinedSamRecords.size () + " genome SAM records for fragment: " + combinedSamFragmentName);
	}
	
	printGenomeSamRecords (computedGenomeSamPairs, combinedSamRecords, pairedEndAlignments, chromosomeIds, outputWriter);
	
	/* Read next SAM records */
	combinedSamRecords = readNextCombinedSamRecord (combinedSamReader, combinedSamFilename);
       
	lineNumber++;
	if (lineNumber % countUnit == 0 && ! quiet) {
	  System.err.print(".");
	}
      }

      if (bedLine != null && combinedSamLine == null) {
	throw new IOException ("BED fragment: " + weightObject.getFragmentId () + " without a SAM Record (at end of file).\n");
      }

      if (debugLevel >= 2) {
	System.err.println("Weight object: " + weightObject);
	System.err.println("SAM records: " + combinedSamRecords);
      }


      /* Print the genome alignments for the fragments at the end of the SAM that are not contained in the BED file */

      if (debugLevel >= 2 || false) {
	System.err.println ("Sam records before while: " + combinedSamRecords);
	System.err.flush();
      }

      while (combinedSamRecords != null) {
	for (SamRecord samRecord: combinedSamRecords) {
	  if (chromosomeIds.contains(samRecord.getReferenceName ()) || samRecord.getReferenceName ().equals("*")) {
	    outputWriter.println (samRecord);
	  }
	}
	combinedSamRecords = readNextCombinedSamRecord (combinedSamReader, combinedSamFilename);
      }

      if (debugLevel >= 2 || false) {
	System.err.println ("Sam records after while: " + combinedSamRecords);
	System.err.flush();
      }

      if (debugLevel >= 1) {
	System.err.println ("Processing last fragment");
      }

      combinedSamReader.close();
      bedReader.close();
      outputWriter.close ();

      if (lineNumber > countUnit  && ! quiet) {
	System.err.println();
      }
      
      System.err.println (numBedRecords + " BED records read for " + numFragments + " fragments.");

    }
    catch (Exception e) {
      System.err.println ("ERROR: Problem with BED record: " + bedLine + "\n" +
			  "Combined SAM record: " + combinedSamLine + "\n" +
			  "Error message: " + (e==null?"No error message":e.getMessage()));
      System.exit(1);
    }

    System.err.println("\nGenome SAM file " + outputFilename + " created.");
  }
}
