/**File: ComputeGeneCountsSam.java 

Original Author: Sven Schuierer
Date: 06/01/2012

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
 *                           Class ComputeGeneCountsSam
 *
 *  
 ***********************************************************************************/


public class ComputeGeneCountsSam {

  private static int debugLevel;
  private static HashSetTable<String, String> transcriptGeneMapTable     = new HashSetTable<String, String>  (7000000);
  private static Hashtable<String, Integer>   transcriptPositionMapTable = null;
  private static Hashtable<String, Integer>   transcriptLengthTable      = new Hashtable<String, Integer> (100 * 1000);;
  private static int minEffectiveLength = 100;

  
  /***********************************************************************************
   * 
   *                           loadGeneTranscriptFile
   *
   ***********************************************************************************/

  private static void loadGeneTranscriptFile (String transcriptGeneMapFilename) throws IOException {

    if (transcriptGeneMapFilename == "") {
      return;
    }
    
    BufferedReader reader = UtilLib.getBufferedReader (transcriptGeneMapFilename);

    System.err.println("Loading gene transcript file " + transcriptGeneMapFilename);

    String line = reader.readLine();
    int lineNumber = 0;
    int numTranscripts = 0;
    while (line != null) {
      StringTokenizer st = new StringTokenizer (line, "\t");

      String transcriptId = "";
      if (st.hasMoreTokens()) {
	transcriptId = st.nextToken ();	
      }

      String geneId = "";
      if (st.hasMoreTokens()) {
	geneId = st.nextToken ();
      } else {
	throw new IOException ("No gene for transcript " + transcriptId + " found.");
      }

      transcriptGeneMapTable.putValue (UtilLib.getSimplifiedReferenceId(transcriptId), geneId);

      if (st.hasMoreTokens()) {
	
	if (transcriptPositionMapTable == null) {
	  transcriptPositionMapTable = new  Hashtable<String, Integer> (7000000);
	}
	
	String positionToken = null;
	try {
	  positionToken = st.nextToken ();
	  transcriptPositionMapTable.put(transcriptId, new Integer(Integer.parseInt(positionToken)));
	  numTranscripts++;
	}
	catch (NumberFormatException e) {
	  throw new IOException ("Position field for transcript" + transcriptId  + ": " + positionToken + " is not an integer" + "\n" +
				 "Error message: " + (e==null?"Null message":e.getMessage()));
	}
      }
      
      line = reader.readLine();
      lineNumber++;
      
    }

    System.err.println(lineNumber + " lines read and " + numTranscripts + " transcripts added to the position table.");
    reader.close ();

    if (debugLevel >= 2) {
      if (transcriptPositionMapTable != null) {
	System.err.println("Size of transcriptPositionMapTable: " + transcriptPositionMapTable.size());
      }
    }
  }

  
  /***********************************************************************************
   * 
   *                           loadTranscriptLengthFile
   *
   ***********************************************************************************/

  private static Hashtable<String, BitSet> loadTranscriptLengthFile (BufferedReader transcriptLengthReader, String transcriptLengthReaderType) throws IOException {

    Hashtable<String, BitSet> transcriptStartPositionTable = new Hashtable<String, BitSet> (100 * 1000);
    boolean isSamFile = transcriptLengthReaderType.startsWith("SAM file");

    String line = transcriptLengthReader.readLine();
    int lineNumber = 0;
    int numTranscripts = 0;
    boolean sequenceHeaderStarted = false;
    while (line != null) {

      if (isSamFile && (! line.startsWith ("@") || (sequenceHeaderStarted && ! line.startsWith ("@SQ")))) {
	break;
      }
      
      StringTokenizer st = new StringTokenizer (line, "\t ");

      String transcriptId = "";
      if (st.hasMoreTokens()) {
	String idToken = st.nextToken ();
	if (isSamFile) {
	  if (! idToken.equals("@SQ")) {
	    continue;
	  } else {
	    sequenceHeaderStarted = true;
	  }
	  
	  if (st.hasMoreTokens()) {
	    idToken = st.nextToken ();
	    if (! idToken.startsWith("SN:")) {
	      throw new IOException ("SN field missing for line: " + line);
	    }
	    idToken = idToken.substring(3);
	  } else {
	    throw new IOException ("SN field missing for line: " + line);
	  }
	}
	transcriptId = idToken;
      } else {
	throw new IOException ("Transcript id field missing in line: " + line);
      }

      String lengthToken = null;
      if (st.hasMoreTokens()) {
	try {
	  lengthToken = st.nextToken ();
	  if (isSamFile) {
	    if (! lengthToken.startsWith("LN:")) {
	      throw new IOException ("LN field missing in line: " + line);
	    }
	    lengthToken = lengthToken.substring(3);
	  }
	  int transcriptLength = Integer.parseInt(lengthToken);
	  BitSet bitSet = new BitSet (transcriptLength);
	  transcriptStartPositionTable.put(transcriptId, bitSet);

	  /* Note that the size of bitSet is not transcriptLength but the next multiple of 128
	     that is at least as large as transcriptLength. Hence, we need to store transcriptLength
	     separately. */
	  transcriptLengthTable.put(transcriptId, new Integer (transcriptLength));

	  numTranscripts++;
	}
	catch (NumberFormatException e) {
	  throw new IOException ("Length field for transcript" + transcriptId  + ": " + lengthToken + " is not an integer" + "\n" +
				 "Error message: " + (e==null?"Null message":e.getMessage()));
	}
      } else {
	throw new IOException ("Length field missing in line: " + line);
      }
      
      line = transcriptLengthReader.readLine();
      lineNumber++;
      
    }

    System.err.println(lineNumber + " lines read and " + numTranscripts + " transcripts added to the start position table.");    
    if (transcriptLengthReaderType.equals("SAM file stdin") && ! line.startsWith ("@")) {
      throw new IOException ("ERROR: Header of SAM file ends directly after sequence length part and SAM file is read from stdin.");
    }

    if (debugLevel >= 2) {
      if (transcriptStartPositionTable != null) {
	System.err.println("Size of transcriptStartPositionTable: " + transcriptStartPositionTable.size());
      }
    }

    return transcriptStartPositionTable;
    
  }
  
  
  /***********************************************************************************
   * 
   *                           estimateEffectiveLength
   *
   * The number of expected distinct values d when drawing n times from a pool 
   * of m values (with replacement) is E[d] = m (1 - (1 - 1/m)^n). We observe d, in our
   * case the number of different starting positions and n, the number of reads,
   * and would like to estimate m. In order to do this we try to find a value for
   * m that approximately satisfies the above equation (using binary search starting
   * with m_max = length of the transcript; we can use binary search since the 
   * function m (1 - (1 - 1/m)^n) is monotone in m).
   *
   ***********************************************************************************/

  private static int estimateEffectiveLength (int numStartPositions, int numReads, int transcriptLength) {

    if (numStartPositions == 0) {
      return 0;
    }
    
    if (numStartPositions >= numReads) {
      return numStartPositions;
    }
    
    double Ed;
    int m;

    if (debugLevel >= 2) {
      System.err.println ("estimateEffectiveLength - numStartPositions: " + numStartPositions + ", numReads: " + numReads + ", transcriptLength: " + transcriptLength);
    }

    int l = 1;
    int r = transcriptLength;
    while (l <= r) {
       m = l + (r - l) / 2;

       Ed = m * (1 - Math.pow(1 - 1.0/m, numReads));
       
       if (numStartPositions < Ed) {
	 r = m - 1;
       } else if (numStartPositions > Ed) {
	 l = m + 1;
       } else {
	 return m;
      }
    }

    m = l + (r - l) / 2;
    
    return Math.min(m, transcriptLength);
    
  }


  /***********************************************************************************/

   private static void printHelp () {
    System.out.println("ComputeGeneCountsSam.java\n" +                                              
      "   -- Script to count fragments aligning to genes\n" +
      "\n" +
    "USAGE: ComputeGeneCountsSam [-A] [-W <read weight thresh.>] [-a]\n" +
    "      [-m <transcript gene map file>] [-w <read weight file>] [-s <sam file>]\n" +
    "      [-o <output file>] [-O <overlap>]\n" +
    "\n" +
    "sam file: file with the sam alignments (- for STDIN) [-]\n" +
    "transcript gene map file: file storing the transcript - gene mappings and\n" +
    "  optionally a position required for the overlap computation (see option -O)\n" +
    "read weight file: file storing the read weights\n" +
    "output file: gene count file (- for STDOUT) [-]\n" +
    "-A: use all genes (even those not contained in <transcript gene file>\n" +
    "-r: count read alignments instead of fragments (single reads are always\n" +
    "    counted separately - even if both ends align)\n" +
    "-W FLOAT: Minimal weight of a read; reads with a lower weight are disregarded \n" +
    "    [0.01]. If it is negative, then the reads are not weighted.\n" +
    "-a: turn warnings on\n" +
    "-O INT: required overlap of a read to the left and right of transcript position\n" +
    "   provided in the transcript gene map file (see option -m)\n" +
    "-P: consider only primary alignment (as specified in the flag field of the SAM\n" +
    "    entries)\n" +
    "-p STRING: number of start positions per transcript file\n" +
    "-l STRING: transcript length file (if STRING = SAM, then the SAM header is\n" +
    "      used to determine the transcript lengths)\n" +
    " -n: output only non-zero counts (otherwise output all counts).\n" +
    "\n" +
    "Reads a SAM file from STDIN or from <sam file> and output a file with gene counts.\n");
  }
                                     

  /***********************************************************************************/

  public static void main (String [] args) {

    String  transcriptGeneMapFilename   = "";
    String  geneCountFilename           = "-";
    String  readWeightFilename          = "none";
    String  startPositionNumberFilename = "none";
    String  transcriptLengthFilename    = "none";
    String  samFilename                 = "-";

    double readWeightThreshold = 0.01;
    
    boolean useAllGenes         = false;
    boolean countReadAlignments = false;
    boolean saveGeneWeights     = false;

    String geneWeightFilename = null;

    String sampleName = "none";

    int overlap = -1;

    boolean primaryAlignmentsOnly = false;
    boolean outputZeroes = true;

    boolean printLines = false;
    boolean saveProcessedReads = false;
    boolean unambiguous = false;
    
    Getopt g = new Getopt("ComputeGeneCountsSam.java", args, "aAd:g:G:l:m:no:O:p:rR:s:S:uw:W:h");
    
    int c;
    String arg = "";
    
    c = g.getopt();
    
    while (c  != -1) {
      switch(c) {
      case 'a':
	UtilLib.setWarningsOn (true);
	break;	
      case 'A':
	useAllGenes = true;
	break;	
      case 'd':
	arg = g.getOptarg();
	UtilLib.setDebugLevel(Integer.parseInt(arg));
	break;
      case 'g':
	geneWeightFilename = g.getOptarg();
	saveGeneWeights    = true;
	break;
      case 'G':
	SamProcessorCount.setSpecialGeneId(g.getOptarg());
	break;
      case 'l':
	transcriptLengthFilename = g.getOptarg();
	break;
      case 'm':
	transcriptGeneMapFilename = g.getOptarg();
	break;
      case 'n':
	outputZeroes = false;
	break;
      case 'o':
	geneCountFilename = g.getOptarg();
	break;
      case 'O':
	overlap = Integer.parseInt(g.getOptarg());
	break;
      case 'p':
	startPositionNumberFilename = g.getOptarg();
	break;
      case 'P':
	primaryAlignmentsOnly = true;
	break;
      case 'r':
	countReadAlignments = true;
	break;
      case 'R':
	/* legacy option */
	readWeightThreshold = Double.parseDouble(g.getOptarg());
	break;
      case 's':
	samFilename = g.getOptarg();
	break;
      case 'S':
	sampleName = g.getOptarg();
	break;
      case 'u':
	unambiguous = true;
	break;
      case 'w':
	readWeightFilename = g.getOptarg();
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
    
    SamRecord.init ();
    if (debugLevel >= 2) {
      System.out.println("Options read.");
    }

    int countUnit = 5 * 1000 * 1000;

    String line = "";
    try {

      UtilLib.setDebugLevel(0);
      if (geneCountFilename == "") {
	throw new IOException ("Gene count file not specified.");
      }

      FragmentEntry.setWarningsOn (UtilLib.warningsOn ());

      BufferedReader readWeightReader = null;
      if (readWeightFilename != "" && ! readWeightFilename.equals("none")) {
	readWeightReader = new BufferedReader (new FileReader (readWeightFilename));
      }

      /* Fill transcriptGeneMapTable and transcriptPositionMapTable */
      if (transcriptGeneMapFilename != "") {
	loadGeneTranscriptFile (transcriptGeneMapFilename);
	
	if (transcriptPositionMapTable == null) {
	  System.err.println("No position information provided in the transcript gene map file " + transcriptGeneMapFilename +
			     ". All mapped reads will be counted.");
	} else if (overlap == -1) {
	  System.err.println("No overlap threshold given (see option -O). All mapped reads will be counted.");
	}
      }
      
      Hashtable<String, Double> geneCountTable = new Hashtable<String, Double> (10000000);
      if (transcriptGeneMapTable != null) {
	for (String transcriptId: transcriptGeneMapTable.keySet ()) {
	  for (String geneId: transcriptGeneMapTable.getSet(transcriptId)) {
	    geneCountTable.put (geneId, new Double(0));
	  }
	}
      }

      if (debugLevel >= 2) {
	System.err.println ("Transcript gene count table with " + geneCountTable.size() + " different gene ids loaded.");
      }

      UtilLib.setDebugLevel(debugLevel);

      /* Read the transcripts lengths and initialize the transcriptStartPositionTable */
      BufferedReader transcriptLengthReader = null;
      Hashtable<String, BitSet> transcriptStartPositionTable = null;
      if (startPositionNumberFilename != "" && ! startPositionNumberFilename.equals("none")) {

	String transcriptLengthReaderType = "tab-delimited file";
	if (transcriptLengthFilename.toLowerCase().equals("sam")) {
	  transcriptLengthFilename   = samFilename;
	  transcriptLengthReaderType = "SAM file";
	  if (samFilename.equals("-")) {
	    transcriptLengthReaderType = transcriptLengthReaderType + " stdin";
	  }
	} else if (transcriptLengthFilename.toLowerCase().endsWith("sam")) {
	  transcriptLengthReaderType = "SAM file";
	}

	if (! transcriptLengthFilename.equals("-")) {
	  System.err.println ("Reading transcript lengths from file " + transcriptLengthFilename);
	} else {
	  System.err.println ("Reading transcript lengths from stdin");
	}
	transcriptLengthReader = UtilLib.getBufferedReader (transcriptLengthFilename);
	transcriptStartPositionTable = loadTranscriptLengthFile (transcriptLengthReader, transcriptLengthReaderType);
	if (! transcriptLengthFilename.equals(samFilename) || ! transcriptLengthFilename.equals("-")) {
	  transcriptLengthReader.close ();
	}

	if (debugLevel >= 2) {
	  for (String transcriptId: transcriptLengthTable.keySet ()) {
	    System.err.println(transcriptId + ": " + transcriptLengthTable.get(transcriptId));
	  }
	}
      }

      /* Intialitze the SAM file reader */
      BufferedReader samFileReader = null;
      if (transcriptLengthFilename.equals(samFilename) && samFilename.equals("-")) {
	samFileReader = transcriptLengthReader;
      } else {
	samFileReader = UtilLib.getBufferedReader (samFilename);
      }

      /* Intialitze the SAM reader and the SAM processor classes */
      SamReader samReader = new SamReader (samFileReader, saveProcessedReads, samFilename.equals("-")?"std in":samFilename, countReadAlignments);
      SamProcessorCount samProcessorCount =
	new SamProcessorCount (transcriptGeneMapTable, transcriptPositionMapTable, readWeightReader, readWeightThreshold, overlap, geneCountTable,
			       transcriptStartPositionTable, useAllGenes, countReadAlignments, primaryAlignmentsOnly, unambiguous);

      if (saveGeneWeights) {
	samProcessorCount.setGeneWeightsWriter(new PrintWriter (geneWeightFilename));
      }

      /* Read the input SAM file */
      samReader.readSamFile (samProcessorCount);

      
      /* Write the gene counts */
      PrintWriter geneCountWriter = UtilLib.getPrintWriter (geneCountFilename);
      TreeSet<String> geneIds = null;
      if (transcriptGeneMapFilename != "") {
	geneIds = new TreeSet<String> ();
	for (String transcriptId: transcriptGeneMapTable.keySet()) {	  
	  geneIds.addAll (transcriptGeneMapTable.getSet (transcriptId));
	}		  
      } else {
	geneIds = new TreeSet<String> (geneCountTable.keySet());
      }

      for (String geneId: geneIds) {
	if (geneCountTable.get(geneId) != null) {
	  if (outputZeroes || geneCountTable.get(geneId).doubleValue() != 0.0) {
	    geneCountWriter.println(geneId + "\t" + geneCountTable.get(geneId).doubleValue());
	  }
	} else if (outputZeroes) {
	  geneCountWriter.println(geneId + "\t" + 0);
	}
      }

      geneCountWriter.close();
      System.err.println (samProcessorCount.getFragmentsCounted () + " fragments counted." );
      
      
      /* Write the number of start position */
      if (startPositionNumberFilename != "" && ! startPositionNumberFilename.equals("none")) {
	PrintWriter startPositionCountWriter = UtilLib.getPrintWriter (startPositionNumberFilename);
	
	Hashtable<String, Integer> transcriptStartPositionCountTable = new Hashtable<String, Integer> (transcriptStartPositionTable.size ());
	Hashtable<String, Integer> transcriptMinStartPositionTable   = new Hashtable<String, Integer> (transcriptStartPositionTable.size ());
	
	/* Compute the cardinality of each transcript start position bit set and aggregate on
	   the gene level (to the maximum) if required */
	for (String transcriptId: transcriptStartPositionTable.keySet()) {

	  /* Get the id for the start positions: either the transcript id or the gene id */
	  String startPositionId  = transcriptId;

	  /* Aggregate to maximum - only needed if the transcript to gene id mapping is provided */
	  BitSet startPositionSet = transcriptStartPositionTable.get (transcriptId);
	  if (transcriptStartPositionCountTable.keySet().contains (startPositionId)) {
	    int currentCardinality = transcriptStartPositionCountTable.get (startPositionId).intValue ();
	    if (currentCardinality < startPositionSet.cardinality()) {
	      transcriptStartPositionCountTable.put (startPositionId, new Integer(Math.max (currentCardinality, startPositionSet.cardinality())));
	      transcriptMinStartPositionTable.put (startPositionId, new Integer(startPositionSet.nextSetBit(0)));
	      transcriptLengthTable.put (startPositionId, transcriptLengthTable.get(transcriptId));
	    }												       
	  } else {
	    transcriptStartPositionCountTable.put (startPositionId, new Integer(startPositionSet.cardinality()));
	    transcriptMinStartPositionTable.put (startPositionId, new Integer(startPositionSet.nextSetBit(0)));
	    transcriptLengthTable.put (startPositionId, transcriptLengthTable.get(transcriptId));
	  }
	}

	String prefix = "Id";
	if (! sampleName.equals("none")) {
	  prefix = "Sample name" + "\t" + prefix;
	}
	startPositionCountWriter.println(prefix + "\t" + "Number of start positions" + "\t" + "Min start position" + "\t" + "Length" + "\t" + "Count" + "\t" +
					 "Est. eff. length" + "\t" + "Est. min start pos.");
	for (String startPositionId: transcriptStartPositionCountTable.keySet()) {
	  prefix = startPositionId;
	  if (! sampleName.equals("none")) {
	    prefix = sampleName + "\t" + prefix;
	  }

	  if (geneCountTable.get(startPositionId) == null) {
	    geneCountTable.put(startPositionId, new Double (0));
	  }

	  int transcriptStartPositionCount = transcriptStartPositionCountTable.get(startPositionId).intValue();
	  int transcriptLength             = transcriptLengthTable.get(startPositionId).intValue();
	  double geneCount                 = geneCountTable.get(startPositionId).doubleValue();
	
	  int estimatedEffectiveLength = estimateEffectiveLength (transcriptStartPositionCount, (int) geneCount, transcriptLength);

	  int minTranscriptStartPosition = transcriptMinStartPositionTable.get(startPositionId).intValue() + 1;

	  /* According to the "German tank problem" (which is similar but not identical to the problem we consider here
	     since we allow a position to be occupied more than once) the minimum variance unbiased estimator for the
	     maximum of a uniform distribution if the observed maximum is m and k value have been seen is
	     m + m/k - 1. */
	  long estimatedMinTranscriptStartPosition = Math.max(1, minTranscriptStartPosition -
							      Math.max(0, Math.round((transcriptLength - minTranscriptStartPosition + 1) * 1.0 / geneCount) - 1));

	  startPositionCountWriter.println(prefix + "\t" + transcriptStartPositionCount + "\t" + minTranscriptStartPosition + "\t" + transcriptLength + "\t" +
					   geneCount + "\t" + estimatedEffectiveLength + "\t" + estimatedMinTranscriptStartPosition);
	  
	}

	startPositionCountWriter.close();
      }
    }
    catch (Exception e) {
      System.err.println ("Problem in line: " + line + ": " + e==null?"No error message":e.getMessage());
      System.exit(1);
    }
  }
}

