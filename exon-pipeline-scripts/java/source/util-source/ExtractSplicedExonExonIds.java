/**File:ExtractSplicedExonExonIds.java 

Original Author: Sven Schuierer
Date: 29/09/2014

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
 *                              Class GtfEntry
 *
 *   Stores the values of a GTF entry
 *
 ***********************************************************************************/

class GtfEntry {

  private static int debugLevel = 0;

  
  /***********************************************************************************
   *
   *                             Object variables
   *
   ***********************************************************************************/

  private String referenceName = "";
  private String source        = "";
  private String feature       = "";
  private int    start         = 0;
  private int    end           = 0;
  private String score         = "";
  private String strand        = "";
  private String frame         = "";
  
  private Hashtable<String, String> attributeValues = new Hashtable<String, String> (8);
  
  private String geneId         = "";
  private String geneName       = "";
  private String transcriptId   = "";
  private String transcriptName = "";
  private String exonId         = "";
  private int    exonNumber     = 0;

  
  /***********************************************************************************
   *
   *                             Constructors
   *
   ***********************************************************************************/

  GtfEntry () {
    exonId = "unknown";
  }



  GtfEntry (String line) throws IOException {

    StringTokenizer st = new StringTokenizer (line, "\t");
    int i = 0;
    String attributes = "";
    while (st.hasMoreTokens ()) {
      String token = st.nextToken ();
      switch (i) {
      case 0:
	referenceName = token;
	break;
      case 1:
	source = token;
	break;
      case 2:
	feature = token;
	break;
      case 3:
	start = Integer.parseInt(token);
	break;
      case 4:
	end = Integer.parseInt(token);
	break;
      case 5:
	score = token;
	break;
      case 6:
	strand = token;
	break;
      case 7:
	frame = token;
	break;
      case 8:
	attributes = token;
	break;
      default:
	System.err.println ("Too many fields for GFT entry: " + line);
      }

      i += 1;
      
    }

    String idField = "";
    StringTokenizer attSt = new StringTokenizer (attributes, ";");
    while (attSt.hasMoreTokens ()) {
      String token = attSt.nextToken ();
      
      StringTokenizer tokenSt = new StringTokenizer (token, "\"");
      if (tokenSt.hasMoreTokens ()) {
	idField = tokenSt.nextToken ().trim();
	if (tokenSt.hasMoreTokens ()) {
	  attributeValues.put(idField, tokenSt.nextToken ().trim ());
	  if (debugLevel >= 2) {
	    System.out.println (idField + ": " + attributeValues.get(idField));
	  }
	  
	  if (tokenSt.hasMoreTokens ()) {
	    throw new IOException ("Attributes of GTF entry " + line + " wrongly formatted for field " + token);
	  }
	} else {
	  throw new IOException ("Attributes of GTF entry " + line + " wrongly formatted for field " + token);
	}
      } else {
	throw new IOException ("Attributes of GTF entry " + line + " wrongly formatted for field " + token);
      }
    }

    geneId         = attributeValues.get("gene_id");
    geneName       = attributeValues.get("gene_name");
    transcriptId   = attributeValues.get("transcript_id");
    transcriptName = attributeValues.get("transcript_name");
    exonId         = attributeValues.get("exon_id");

    if (debugLevel >= 2) {
      System.out.println ("exon id: " + exonId);
    }
    
    String exonNumberString = attributeValues.get("exonNumber");
    if (exonNumberString != null) {
      exonNumber = Integer.parseInt (exonNumberString);
    }
  }

  
  /***********************************************************************************
   *
   *                             Get functions
   *
   ***********************************************************************************/

  public String getReferenceName () {
   return referenceName;
  }

  public String getSource () {
   return source;
  }

  public String getFeature () {
   return feature;
  }

  public int getStart () {
   return start;
  }

  public int getEnd () {
   return end;
  }

  public String getScore () {
   return score;
  }

  public String getStrand () {
   return strand;
  }

  public String getFrame () {
   return frame;
  }

  public Hashtable<String, String> getAttributes () {
    return attributeValues;
  }

  public String getGeneId  () {
   return geneId;
  }

  public String getGeneName () {
   return geneName;
  }

  public String getTranscriptId () {
   return transcriptId;
  }

  public String getTranscriptName () {
   return transcriptName;
  }

  public String getExonId () {
    if (exonId != null) {
      return exonId;
    }
    return referenceName + "/" + start  + "/" + end  + "/" + strand;
  }
 
  public int getExonNumber () {
   return exonNumber;
  }

  
  /***********************************************************************************
   *
   *                          hashCode and equals
   *
   ***********************************************************************************/

  public int hashCode () {
    return getExonId().hashCode ();
  }

  public boolean equals (Object o) {
    GtfEntry g = (GtfEntry) o;
    return getExonId().equals(g.getExonId ());
  }

}



/***********************************************************************************
 *
 *                              Class ExtractSplicedExonExonIds
 *
 *   Reads a SAM file line by line and saves the genomic intervals spanned by spliced
 *   reads.
 *
 ***********************************************************************************/

public class ExtractSplicedExonExonIds {

  private static int debugLevel   = UtilLib.getDebugLevel ();
  private static final int countUnit    = 2500 * 1000;
  private static final int gtfCountUnit =  100 * 1000;


  /***********************************************************************************
   *
   *  Read GTF file
   *
   ***********************************************************************************/
  
  private static ArrayList<HashSetTable <String, GtfEntry>> readGtfFile (String gtfFilename) throws IOException {
      
    System.err.println("Reading GTF file " + (gtfFilename.equals("-")?"stdin":gtfFilename));
    System.err.println("(. = " + gtfCountUnit + " entries.)");
    System.err.flush();
    BufferedReader gtfReader = UtilLib.getBufferedReader (gtfFilename);
    String line = gtfReader.readLine();
    int lineNumber = 0;

    HashSetTable <String, GtfEntry> leftBoundarySet  = new HashSetTable <String, GtfEntry> (100000);
    HashSetTable <String, GtfEntry> rightBoundarySet = new HashSetTable <String, GtfEntry> (100000);
    while (line != null) {
	
      GtfEntry gtfEntry = new GtfEntry (line);
      leftBoundarySet.putValue(gtfEntry.getReferenceName () + "/" + gtfEntry.getStart (), gtfEntry);
      rightBoundarySet.putValue(gtfEntry.getReferenceName () + "/" +  gtfEntry.getEnd (), gtfEntry);
	
      lineNumber++;
      if (lineNumber % gtfCountUnit == 0) {
	System.err.print(".");
      }
	  
      line = gtfReader.readLine();
    }
      
    gtfReader.close ();
    if (lineNumber > gtfCountUnit) {
      System.err.print("\n");
    }

    ArrayList<HashSetTable <String, GtfEntry>> boundarySetArrayList = new ArrayList<HashSetTable <String, GtfEntry>> ();
    boundarySetArrayList.add(leftBoundarySet);
    boundarySetArrayList.add(rightBoundarySet);
   
    return boundarySetArrayList;
    
  }
  
  /***********************************************************************************
   *
   *  Read SAM file (and weight file if specified)
   *
   ***********************************************************************************/

  private static TreeSet<Interval> readSamFile (String samFilename, String weightFilename, boolean quantify, boolean quiet) throws IOException {

    Counter fragmentCounter = new Counter (9);
    Pattern pattern = Pattern.compile("^F[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]$");
    
    System.err.println("Reading SAM file " + (samFilename.equals("-")?"stdin":samFilename));
    System.err.println("(. = " + countUnit + " entries.)");
    System.err.flush();

    int lineNumber = 0;
    TreeSet<Interval> intronIntervals = new TreeSet <Interval> ();

    /* numSplicedSamRecords is the number of SAM records that are spliced */
    int numSplicedSamRecords  = 0;
      
    /* numSamJunctions is the number of SAM record/junction pairs or the sum SAM record splice multiplicities */
    int numSamJunctions = 0;

    int duplicateCoverage = 0;

    /* Note that one fragment can cover one junction twice, once for each read of the fragment; a fragment can cover
       multiple junction either if one SAM entry has a splice multiplicity > 1 or if the two reads cover different
       junctions. */
    Hashtable<String, Integer> splicedFragmentTable = new Hashtable<String, Integer> ();
    boolean counterChecked = false;
    boolean createNewFragmentIds = false;
    String oldMappedFragmentName = "";
    String fragmentName = "";
    FragmentEntry fragmentEntry = null;
    int numMismatches = 0;
    int numUniquelyMappingFragments = 0;
    int numMultiMappers = 0;
    int softClippingLengthTotal = 0;

    
    /***********************************************************************************
     *
     * Loop over lines in SAM file
     *
     ***********************************************************************************/

    String line = "";
    String weightLine = "";
    try {
      BufferedReader samReader = UtilLib.getBufferedReader (samFilename);
      line = samReader.readLine();

      BufferedReader weightReader = null;
      if (weightFilename != "") {
	weightReader = UtilLib.getBufferedReader (weightFilename);
	weightLine = weightReader.readLine();
      }
      
      while (line != null) {

	if (debugLevel >= 3) {
	  System.err.println ("Reading line: " + line);
	  System.err.flush();
	}

	if (! line.startsWith("@") && ! line.equals("")) {
	
	  SamRecord samRecord = new SamRecord (line, lineNumber, samFilename);
	  String originalFragmentName = samRecord.getOriginalFragmentName ();
	  
	  if (! counterChecked) {
	    counterChecked = true;
	    Matcher match = pattern.matcher(samRecord.getFragmentName());
	    if (! match.find()) {
	      System.err.println("Changing fragment id format from original id format to: Fnnnnnnnnn.");
	      SamRecord.setCounter (fragmentCounter);
	      createNewFragmentIds = true;
	    }
	  }

	  if (! samRecord.isMapped()) {
	    line = samReader.readLine();
	    continue;
	  }
	
	  if (! originalFragmentName.equals(oldMappedFragmentName)) {
	    if (oldMappedFragmentName != "") {
	      fragmentCounter.inc();
	      if (weightFilename != "" && numMismatches != fragmentEntry.getSumEditDistance () && ! quiet) {
		throw new IOException ("Number of SAM mismatches for " + oldMappedFragmentName + ": " + numMismatches +
				       " does not match weight file (" + fragmentEntry.getFragmentName() + "): " + fragmentEntry.getSumEditDistance ());
	      }
	    }

	    fragmentName = samRecord.getFragmentName ();
	    if (weightFilename != "") {
	      fragmentEntry = new FragmentEntry (weightLine);
	      if (debugLevel >= 2) {
		System.err.println ("SAM fragment name: " + fragmentName + ", new fragment entry name: " + fragmentEntry.getFragmentName ());
	      }
	      if (! createNewFragmentIds) {
		while (fragmentEntry.getFragmentName ().compareTo (fragmentName) < 0) {
		  weightLine = weightReader.readLine ();
		  fragmentEntry = new FragmentEntry (weightLine);
		  if (debugLevel >= 2) {
		    System.err.println ("New fragment entry name: " + fragmentEntry.getFragmentName ());
		  }
		}
	      }
	      weightLine = weightReader.readLine ();
	      if (! fragmentEntry.getFragmentName ().equals(fragmentName)) {
		throw new IOException ("No weight entry for " + fragmentName + " found - weight fragment: " + fragmentEntry.getFragmentName ());
	      }
	    }
	    oldMappedFragmentName = originalFragmentName;

	    if (debugLevel >= 2) {
	      System.err.println("OldMappedFragmentName: " +  oldMappedFragmentName);
	    }
	    numMismatches = samRecord.getNumMismatches ();
	  } else {
	    numMismatches += samRecord.getNumMismatches ();
	  }
	
	
	  fragmentName = samRecord.getFragmentName ();
	  CiagrString ciagrString = new CiagrString (samRecord, true);
	  if (debugLevel >= 2) {
	    System.err.println ("ciagrString: " + ciagrString);
	  }
	  softClippingLengthTotal += ciagrString.getSoftClippingLength ();
	
	  if (debugLevel >= 2) {
	    System.err.println ("samRecord.isSpliced (): " + samRecord.isSpliced ());
	  }
	  
	  if (samRecord.isSpliced ()) {
	    if (debugLevel >= 2) {
	      System.err.println ("SAM record " + samRecord + " is spliced.");
	    }
	    int position = samRecord.getPosition ();
	    Vector<int []> localIntronCoordinates = ciagrString.getIntronCoordinates ();

	    if (debugLevel >= 2) {
	      System.err.println ("localIntronCoordinates: " + localIntronCoordinates);
	    }
	    if (localIntronCoordinates.size() > 0) {
	      numSplicedSamRecords += 1;

	      if (quantify) {
		Integer numFragmentJunctions = splicedFragmentTable.get(fragmentName);
		if (numFragmentJunctions == null) {
		  if (weightFilename != "") {
		    if (fragmentEntry.getNumAlignments () == 1) {
		      numUniquelyMappingFragments += 1;
		    } else {
		      numMultiMappers += 1;
		    }
		  }
		  splicedFragmentTable.put(fragmentName, new Integer (1));
		} else {
		  splicedFragmentTable.put(fragmentName, new Integer (numFragmentJunctions.intValue () + 1));
		}
	      }

	      /* A relative intron interval [I0, I1] corresponds to the genomic interval [position - 1 + I0, position - 1 + I1] since the relative
		 intron intervals start with position 1, that is position corresponds to 1; for instance, [1,2] would correspond to [position, position + 1] */
	      for (int [] localIntronInterval: localIntronCoordinates) {
		localIntronInterval[0] += position - 1;
		localIntronInterval[1] += position - 1;
		Interval interval = null;
		if (quantify) {
		  interval = new Interval(samRecord.getReferenceName(), localIntronInterval, fragmentName);
		} else {
		  interval = new Interval(samRecord.getReferenceName(), localIntronInterval);
		}
		if (! intronIntervals.contains (interval)) {
		  intronIntervals.add (interval);
		} else {
		  duplicateCoverage++;
		}
		numSamJunctions++;

		if (debugLevel >= 2) {
		  int numIntronIntervalFragments = 0;
		  for (Interval curInterval: intronIntervals) {
		    numIntronIntervalFragments += curInterval.getFragmentNameSet().size();
		  }
		
		  int numSplicedSamFragments = 0;
		  if (quantify) {
		    for (String curFragmentName: splicedFragmentTable.keySet()) {
		      numSplicedSamFragments += splicedFragmentTable.get(curFragmentName).intValue();
		    }
		  }
		  System.out.println((fragmentName.length()>36?fragmentName.substring(36):fragmentName) + "\t" + samRecord.getCiagrString() + "\t" +
				     localIntronInterval[0] + "-" + localIntronInterval[1]);
		  /* System.out.println("Num. spliced SAM: " + numSplicedSamRecords + ", num SAM record/junction pairs: " + numSamJunctions +
		     ", num. spliced frag: " + splicedFragmentTable.keySet().size() + ", num. intron intervals: " + intronIntervals.size() +
		     ", num intron int frag: " + numIntronIntervalFragments + ", sum frag splice mult: " + numSplicedSamFragments); */
		}
	      }
	    }
	  }
	}
	
	lineNumber++;
	if (lineNumber % countUnit == 0) {
	  System.err.print(".");
	}
	  
	line = samReader.readLine();
      }

      if (weightFilename != "" && weightLine != null) {
	System.err.println ("Surplus weight entry: " + weightLine);
	weightLine = weightReader.readLine ();
	if (weightLine != null) {
	  int numAdditionalWeightEntries = 0;
	  while (weightLine != null) {
	    weightLine = weightReader.readLine ();
	    numAdditionalWeightEntries += 1;
	  }
	  System.err.println (numAdditionalWeightEntries + " additional entries in the weight file.");
	}
      }

      samReader.close();
      if (weightFilename != "") {
	weightReader.close ();
      }
      if (lineNumber > countUnit) {
	System.err.print("\n");
      }
    }
    catch (Exception e) {
      throw new IOException ("Problem with reading SAM file " + samFilename + " in line:\n" + line + "\n" + " and weight line of file " + weightFilename +
			     ": " + weightLine + "\nMessage: " + e);
    }

    System.err.println("Number of spliced SAM records: " + numSplicedSamRecords);
    System.err.println("Number of read/junction pairs: " + numSamJunctions);
    System.err.println("Number of intron intervals: " + intronIntervals.size());
    System.err.println("Number of SAM records that cover the same junction as the other read of the fragment: " + duplicateCoverage);
    if (weightFilename != "") {
      System.err.println("Number of uniquely mapping fragments: " + numUniquelyMappingFragments);
      System.err.println("Number of multi-mappers: " + numMultiMappers);
    }

    if (quantify) {
      System.err.println("Number of spliced fragments: " + splicedFragmentTable.keySet().size());
      int numIntronIntervalFragments = 0;
      for (Interval interval: intronIntervals) {
	if (interval.getFragmentNameSet() == null) {
	  throw new IOException ("Intron interval " + interval + " has a null fragment name set.");
	}
	numIntronIntervalFragments += interval.getFragmentNameSet().size();
      }
      System.err.println("Number of fragment/junction pairs: " + numIntronIntervalFragments);
      int numSplicedSamFragments = 0;
      for (String curFragmentName: splicedFragmentTable.keySet()) {
	numSplicedSamFragments += splicedFragmentTable.get(curFragmentName).intValue();
      }
      System.err.println("Sum of spliced SAM records per fragment: " + numSplicedSamFragments);
    }

    System.err.println("Soft clipped bases: " + softClippingLengthTotal);
    System.err.println();
    if (quantify) {
      System.err.println("NUM_SPLICED_FRAGMENTS=" + splicedFragmentTable.keySet().size());
    }
    System.err.println("NUM_SPLICED_SAM_RECORDS=" + numSplicedSamRecords);
    System.err.println("NUM_INTRON_INTERVALS=" + intronIntervals.size());
    if (weightFilename != "") {
      System.err.println("NUM_UNIQUE_MAPPERS=" + numUniquelyMappingFragments);
      System.err.println("NUM_MULTI_MAPPERS=" + numMultiMappers);
    }    

    return intronIntervals;
    
  }
     

  /***********************************************************************************/

   private static void printHelp () {
    System.out.print("ExtractSplicedExonExonIds\n" +                                              
    "USAGE: ExtractSplicedExonExonIds [-w <weight file name>] [-j <junction covering\n" +
    "          threshold>] [-Q] -s <SAM file> -g <GTF file> -o <exon exon pair file>\n");
  }
                                     
  /***********************************************************************************/

  public static void main (String [] args) {

    String samFilename = "-";
    String gtfFilename = "none";
    String outputFilename = "-";
    String weightFilename = "";

    boolean quantify = false;
    boolean quiet = false;
    boolean printNonMatchingFragments = false;

    int junctionCoveringThreshold = 5;

    Getopt g = new Getopt("ExtractSplicedExonExonIds.java", args, "d:g:j:o:qQs:w:h");
    
    int c;
    String arg = "";
    
    c = g.getopt();
    
    while (c  != -1) {
      switch(c) {
      case 'd':
	debugLevel = Integer.parseInt(g.getOptarg());
	break;	
      case 'g':
	gtfFilename = g.getOptarg();
	break;	
      case 'j':
	junctionCoveringThreshold = Integer.parseInt(g.getOptarg());
	break;	
      case 'o':
	outputFilename = g.getOptarg();
	break;
      case 'p':
	printNonMatchingFragments = true;
	break;	
      case 'q':
	quiet = true;
	break;	
      case 'Q':
	quantify = true;
	Interval.setQuantify();
	break;	
      case 's':
	samFilename = g.getOptarg();
	break;
      case 'w':
	weightFilename = g.getOptarg();
	break;
      case 'h':
	printHelp();
	System.exit(0);
	break;
      default:
	System.out.print("Error: getopt() returned unknown option: " + c + "\n");
      }
      c = g.getopt();
    }

    String line = "";
    int lineNumber = 0;
    try {

      ArrayList<HashSetTable <String, GtfEntry>> boundarySetArrayList = readGtfFile (gtfFilename);
      
      HashSetTable <String, GtfEntry> leftBoundarySet  = boundarySetArrayList.get(0);
      HashSetTable <String, GtfEntry> rightBoundarySet = boundarySetArrayList.get(1);
      
      TreeSet<Interval> intronIntervals = readSamFile (samFilename, weightFilename, quantify, quiet);
      
      /***********************************************************************************
       *
       *  Write exon - exon file
       *
       ***********************************************************************************/
      
      System.err.println("Writing exon - exon pairs to file:\n" + "  " + (outputFilename.equals("-")?"stdout":outputFilename));
      System.err.flush();
      PrintWriter outputWriter = UtilLib.getPrintWriter (outputFilename);

      int numNoMatchingBoundary     = 0;
      int numLeftMatchingBoundary   = 0;
      int numRightMatchingBoundary  = 0;
      int numBothMatchingBoundaries = 0;
      
      int[] highCoverageJunctions   = new int [junctionCoveringThreshold + 1];
      for (int i = 1; i <= junctionCoveringThreshold; i++) {
	highCoverageJunctions[i] = 0;
      }

      for (Interval interval: intronIntervals) {
	String intronIntervalString = interval.getChromosome()  + ":" + interval.getStart() + "-" + interval.getEnd();
	
	// System.err.println(interval);
	Set<GtfEntry> leftGtfEntrySet  = rightBoundarySet.get(interval.getChromosome() + "/" + (interval.getStart() - 1));	
	Set<GtfEntry> rightGtfEntrySet = leftBoundarySet.get(interval.getChromosome()  + "/" + (interval.getEnd() + 1));

	String additionalFieldsString = "";
	if (quantify) {
	  if (leftGtfEntrySet == null && rightGtfEntrySet == null) {
	    numNoMatchingBoundary += interval.getFragmentNameSet().size();
	    if (printNonMatchingFragments) {
	      for (String curFragmentName: interval.getFragmentNameSet()) {
		System.err.println ("Non matching: " + curFragmentName);
	      }
	    }
	    
	    leftGtfEntrySet = new HashSet<GtfEntry> ();
	    leftGtfEntrySet.add(new GtfEntry());
	    rightGtfEntrySet = new HashSet<GtfEntry> ();
	    rightGtfEntrySet.add(new GtfEntry());
	    
	  } else if (leftGtfEntrySet == null) {
	    numLeftMatchingBoundary += interval.getFragmentNameSet().size();
	    if (printNonMatchingFragments) {
	      for (String curFragmentName: interval.getFragmentNameSet()) {
		System.err.println ("Left matching: " + curFragmentName);
	      }
	    }
	    leftGtfEntrySet = new HashSet<GtfEntry> ();
	    leftGtfEntrySet.add(new GtfEntry());
	  } else if (rightGtfEntrySet == null) {
	    numRightMatchingBoundary += interval.getFragmentNameSet().size();
	    if (printNonMatchingFragments) {
	      for (String curFragmentName: interval.getFragmentNameSet()) {
		System.err.println ("Right matching: " + curFragmentName);
	      }
	    }
	    rightGtfEntrySet = new HashSet<GtfEntry> ();
	    rightGtfEntrySet.add(new GtfEntry());
	  } else {
	    numBothMatchingBoundaries += interval.getFragmentNameSet().size();
	    int upperThreshold = Math.min (interval.getFragmentNameSet().size(), junctionCoveringThreshold);
	    for (int i = 2; i <= upperThreshold; i++) {
	      highCoverageJunctions[i] += 1;
	    }
	  }
	  
	  additionalFieldsString = "\t" + interval.getChromosome()  + ":" + interval.getStart() + "-" + interval.getEnd() + "\t" + interval.getFragmentNameSet().size();
	  interval.clearFragmentNameSet ();
	} else if (leftGtfEntrySet == null || rightGtfEntrySet == null) {
	  continue;
	}

	for (GtfEntry leftGtfEntry: leftGtfEntrySet) {
	  for (GtfEntry rightGtfEntry: rightGtfEntrySet) {
	    outputWriter.println(leftGtfEntry.getExonId ()  + "\t" + rightGtfEntry.getExonId () + additionalFieldsString);
	    if (! rightGtfEntry.getExonId ().equals(leftGtfEntry.getExonId ())) {
	      outputWriter.println(rightGtfEntry.getExonId () + "\t" + leftGtfEntry.getExonId ()  + additionalFieldsString);
	    }
	  }
	}	
      }
      
      outputWriter.close ();

      if (quantify) {
	System.err.println("Number of spliced fragment/junction pairs with no matching boundary: "     + numNoMatchingBoundary);
	System.err.println("Number of spliced fragment/junction pairs with left matching boundary: "   + numLeftMatchingBoundary);
	System.err.println("Number of spliced fragment/junction pairs with right matching boundary: "  + numRightMatchingBoundary);
	System.err.println("Number of spliced fragment/junction pairs with both matching boundaries: " + numBothMatchingBoundaries);
	System.err.println("Number of junctions with at least " + junctionCoveringThreshold + " reads coverage: " + highCoverageJunctions);
	System.err.println();
	System.err.println("NON_BOUNDARY_MATCHING_FRAG_JUNC="   + numNoMatchingBoundary);
	System.err.println("LEFT_BOUNDARY_MATCHING_FRAG_JUNC="  + numLeftMatchingBoundary);
	System.err.println("RIGHT_BOUNDARY_MATCHING_FRAG_JUNC=" + numRightMatchingBoundary);
	System.err.println("BOTH_BOUNDARY_MATCHING_FRAG_JUNC="  + numBothMatchingBoundaries);
	for (int i = 2; i <= junctionCoveringThreshold; i++) {
	  System.err.println("NUM_HIGH_COVERAGE_JUNCTIONS_" + i + "=" + highCoverageJunctions[i]);
	}
      }

      Runtime runtime = Runtime.getRuntime();
      System.err.println("Memory used: " + (runtime.totalMemory() - runtime.freeMemory()));

      
    }
    catch (IOException e) {
      System.err.println ("IO ERROR: Problem in line: " + line + ", message: " + (e==null?"No error message":e.getMessage()));
    }
    catch (Exception e) {
      System.err.println ("ERROR: Problem in line: " + line + ", message: " + (e==null?"No error message":e.getMessage()));
    }
  }
}
