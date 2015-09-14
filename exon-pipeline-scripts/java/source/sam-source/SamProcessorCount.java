/**File: SamProcessorCount.java 

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
 *                              Class SamProcessorCount
 *
 ***********************************************************************************/

class SamProcessorCount implements SamProcessor {

  private static int debugLevel = UtilLib.getDebugLevel ();

  private static String specialGeneId = "";
  public static void setSpecialGeneId (String s) {
    specialGeneId = s;
  }


  /***********************************************************************************
   *
   *                     Object variables and methods
   *
   ***********************************************************************************/

  private HashSetTable<String, String>     transcriptGeneMapTable       = null;
  private Hashtable<String, Integer>       transcriptPositionTable      = null;
  private Hashtable<String, FragmentEntry> readWeightTable              = null;
  private Hashtable<String, Double>        geneCountTable               = null;
  private Hashtable<String, BitSet>        transcriptStartPositionTable = null;

  private BufferedReader readWeightReader = null;
  
  private HashSet<String> geneIds = new HashSet<String>  (2000);

  private String readId               = "";
  private String fragmentName         = "";
  private boolean useAllGenes         = false;
  private boolean countReadAlignments = false;
  private double readWeightThreshold  = 0.01;
  
  private int overlap          = 0;
  private int fragmentsCounted = 0;

  private String readWeightFragmentName    = "";
  private String oldReadWeightFragmentName = "";
  private String readWeightLine            = null;

  private boolean primaryAlignmentsOnly = false;
  private boolean unambiguous           = false;

  
  /***********************************************************************************
   *
   *                     Constructors
   *
   ***********************************************************************************/
  
  SamProcessorCount () {
  }

  SamProcessorCount (boolean useAllGenes) {
    this.useAllGenes = useAllGenes;
  }

  SamProcessorCount (HashSetTable<String, String> transcriptGeneMapTable, Hashtable<String, Integer> transcriptPositionTable, BufferedReader readWeightReader,
		     double readWeightThreshold, int overlap, Hashtable<String, Double> geneCountTable, Hashtable<String, BitSet> transcriptStartPositionTable,
		     boolean useAllGenes, boolean countReadAlignments, boolean primaryAlignmentsOnly, boolean unambiguous) {
    
    this.transcriptGeneMapTable       = transcriptGeneMapTable;
    this.transcriptPositionTable      = transcriptPositionTable;
    this.readWeightReader             = readWeightReader;
    this.geneCountTable               = geneCountTable;
    this.transcriptStartPositionTable = transcriptStartPositionTable;
    
    this.overlap                 = overlap;
    
    this.useAllGenes             = useAllGenes;
    this.countReadAlignments     = countReadAlignments;
    this.primaryAlignmentsOnly   = primaryAlignmentsOnly;
    this.readWeightThreshold     = readWeightThreshold;
    this.unambiguous             = unambiguous;
  }

  public void init (String readId) {
        
  }


  public void init (SamRecord samRecord) {

    init (samRecord.getQueryName ());
    
  }

  
  /***********************************************************************************
   *
   *                           Using gene weights
   *
   * There is the possibility to use the number of genes that a read maps to as 
   * read weights and to write out these weights to a file.
   *
   ***********************************************************************************/

  private boolean saveGeneWeights = false;
  private PrintWriter geneWeightsWriter = null;
  
  public  void setGeneWeightsWriter (PrintWriter w) {
    saveGeneWeights = true;
    geneWeightsWriter = w;
  }

  /***********************************************************************************
   *
   *                           retrieveFragmentEntry
   *
   ***********************************************************************************/
  
  public FragmentEntry retrieveFragmentEntry (String fragmentName, BufferedReader readWeightReader) throws IOException {

    if (fragmentName == null) {
      return null;
    }

    /* Ensure non-null readWeightLine */
    if (readWeightLine == null) {
      if (readWeightReader != null) {
	readWeightLine = readWeightReader.readLine ();
      }
      if (readWeightLine == null) {
	return null;
      }
    }

    FragmentEntry fragmentEntry = new FragmentEntry (readWeightLine);
    readWeightFragmentName = fragmentEntry.getFragmentName ();

    while (readWeightLine != null && fragmentName.compareTo (readWeightFragmentName) > 0) {
      readWeightLine = readWeightReader.readLine ();
      if (readWeightLine != null) {
	fragmentEntry = new FragmentEntry (readWeightLine);
	readWeightFragmentName = fragmentEntry.getFragmentName ();
	
	if (oldReadWeightFragmentName.compareTo(readWeightFragmentName) > 0) {
	  throw new IOException ("Entries in read weight file not sorted: Fragment " + oldReadWeightFragmentName + " occurs before " + readWeightFragmentName);
	}
	
	oldReadWeightFragmentName = readWeightFragmentName;
      }
    }

    if (readWeightLine == null) {
      return null;
    }

    if (fragmentName.equals(readWeightFragmentName)) {
      return fragmentEntry;
    } else {
      return null;
    }
    
  }


  /***********************************************************************************
   *
   *                           addReadWeight
   *
   ***********************************************************************************/

  public void addReadWeight (String fragmentName, HashSet<String> geneIds, BufferedReader readWeightReader, Hashtable<String, Double> geneCountTable)
    throws IOException {
    
    if (debugLevel >= 1) {
      System.out.println("Adding read counts for read id: " + readId);
    }
	    
    double readWeight = 1;
    if (readWeightReader != null) {
      FragmentEntry fragmentEntry = retrieveFragmentEntry (fragmentName, readWeightReader);
      if (fragmentEntry != null && fragmentEntry.getNumAlignments() > 0) {
	readWeight = 1.0 / fragmentEntry.getNumAlignments();
      } else {
	if (UtilLib.warningsOn ()) {
	  System.err.println("WARNING: No weight for fragmentName: " + fragmentName + " found.");
	}
      }
      
      if (debugLevel >= 2) {
	System.out.println("Read weight for fragment " + fragmentName + ": " + readWeight);
      }
      
    } else if (saveGeneWeights) {
      readWeight = 1.0 / geneIds.size();
      geneWeightsWriter.println(fragmentName + "\t" + readWeight);
    }

    if (readWeight >= readWeightThreshold) {
      for (String geneId: geneIds) {      
	Double geneWeightDouble = geneCountTable.get(geneId);
	
	if (geneWeightDouble != null) {

	  if (debugLevel >= 2 || geneId.equals(specialGeneId)) {
	    System.out.println("Gene: " + geneId + " adding " + fragmentName + " with weight " + readWeight + " to " + geneWeightDouble.doubleValue());
	  }
	  
	  geneCountTable.put(geneId, new Double(geneWeightDouble.doubleValue() + readWeight));
	} else {
	  if (UtilLib.warningsOn ()) {
	    System.err.println("WARNING: geneId " + geneId + " not in geneCountTable.");
	  }
	  geneCountTable.put(geneId, new Double(readWeight));
	}
      }

      fragmentsCounted++;
    }

    geneIds.clear ();
  }

  
  /***********************************************************************************
   *
   *                           getFragmentsCounted
   *
   ***********************************************************************************/

  public int getFragmentsCounted () {
    return fragmentsCounted;
  }


  /***********************************************************************************
   *
   *                           getGeneId
   *
   ***********************************************************************************/

  public HashSet<String> getGeneIdSet (SamRecord samRecord, HashSetTable<String, String> transcriptGeneMapTable, boolean useAllGenes) {
    
    String transcriptId = UtilLib.getSimplifiedReferenceId(samRecord.getReferenceName());
    if (debugLevel >= 2) {
      System.out.println("Transcript: " + transcriptId);
    }

    HashSet<String> geneIdSet = transcriptGeneMapTable.get(transcriptId);
    if (debugLevel >= 2) {
      System.out.println("Gene id set: " + geneIdSet);
    }
    
    if (geneIdSet == null) {
      if (useAllGenes) {
	transcriptGeneMapTable.putValue (transcriptId, transcriptId);
	geneCountTable.put(transcriptId, new Double(0));
	
	geneIdSet = new HashSet<String> ();
	geneIdSet.add (transcriptId);
      } else if (UtilLib.warningsOn()) {
	System.out.println("No gene id for transcript: " + transcriptId + " found.");
      }
    }

    return geneIdSet;
     
  }


  /***********************************************************************************
   *
   *                           processSamRecords   
   *
   ***********************************************************************************/

  public void processSamRecords (Vector<SamRecord> samRecords) throws IOException {
      
    HashSet<Integer> processedIndices = new HashSet<Integer> (2 * samRecords.size());

    if (transcriptGeneMapTable == null) {
      transcriptGeneMapTable = new HashSetTable<String, String> (500 * 1000);
    }

    if (samRecords.size() > 0) {
      SamRecord samRecord = samRecords.get(0);
      fragmentName = samRecord.getFragmentName();
      readId       = samRecord.getQueryName ();      
    } else {
      throw new IOException ("Empty samRecords set.");
    }

    if (debugLevel >= 2) {
      System.out.println("Processing " + samRecords.size() + " sam records for fragment " + fragmentName + ".");
    }

    SamRecord mateSamRecord = null;
    for (int i = 0; i < samRecords.size(); i++) {

      Integer iInt = Integer.valueOf (i);
      if (! processedIndices.contains(iInt)) {

	processedIndices.add(iInt);

	SamRecord samRecord = samRecords.get(i);
	
	if (debugLevel >= 1) {
	  System.out.println("Processing read: " + samRecord.getQueryName ());
	}

	if (! primaryAlignmentsOnly || samRecord.isPrimary ()) {
	  if (debugLevel >= 1) {
	    System.out.println("Processing read: " + samRecord.getQueryName ());
	  }

	  if (countReadAlignments && ! readId.equals(samRecord.getQueryName ())) {
	    throw new IOException ("Two different reads ids in samRecords collection: " + readId + " and " + samRecord.getQueryName ());
	  }

	  if (! fragmentName.equals(samRecord.getFragmentName ())) {
	    throw new IOException ("Two different fragment ids in samRecords collection: " + fragmentName + " and " + samRecord.getFragmentName ());
	  }

	  mateSamRecord = null;
	  if (samRecord.hasMate() && ! countReadAlignments) {	  
	    mateSamRecord = samRecord.findMate(samRecords, i, processedIndices);
	    if (primaryAlignmentsOnly && ! mateSamRecord.isPrimary ()) {
	      mateSamRecord = null;
	    }
	  }

	  /* Note that the transcriptGeneMapTable uses the simplified transcriptIds/SAM reference ids
	     of UtilLib.getSimplifiedReferenceId for the juncion ids whereas the transcriptPositionTable
	     does not as each position depends on the actual sequence and not just the exons defining
	     the junction. */
	  HashSet<String> geneIdSet = getGeneIdSet (samRecord, transcriptGeneMapTable, useAllGenes);
	  Integer position = null;
	  if (transcriptPositionTable != null) {
	    position = transcriptPositionTable.get(samRecord.getReferenceName());
	  }

	  int leftPosition  = -1;
	  int rightPosition = -1;
	  if (debugLevel >= 2) {
	    System.err.println ("overlap: " + overlap);
	    System.err.println ("Position: " + position);
	  }

	  if (position != null && overlap > 0) {
	    leftPosition  = Math.max(position.intValue () - overlap + 1, 1);
	    rightPosition = position.intValue () + overlap;
	  }

	  if (debugLevel >= 2) {
	    System.err.println ("Sam record: " + samRecord);
	    System.err.println ("containsInterval (" + leftPosition + ", " + rightPosition + "): " + samRecord.containsInterval (leftPosition, rightPosition));
	  }
	
	  if (samRecord.containsInterval (leftPosition, rightPosition) ||
	      (mateSamRecord != null && mateSamRecord.containsInterval (leftPosition, rightPosition))) {

	    if (debugLevel >= 2) {
	      System.out.println("Adding gene id set: " + geneIdSet);
	    }

	    if (geneIdSet != null) {
	      geneIds.addAll(geneIdSet);
	    }
	  
	    if (debugLevel >= 2) {
	      System.out.println("Genes: " + geneIds);
	    }
	  } else if (UtilLib.warningsOn ()) {
	    System.err.println ("WARNING: [" + leftPosition + "," + rightPosition + "] is not contained in " + samRecord +
				(mateSamRecord!=null?" nor in the mate " + mateSamRecord:""));
	  }

	  if (transcriptStartPositionTable != null) {
	    String transcriptId = samRecord.getReferenceName();
	    int minStartPos     = samRecord.getPosition ();
	    if (mateSamRecord != null && mateSamRecord.getPosition () != -1) {
	      if (minStartPos >= 1) {
		minStartPos = Math.min(minStartPos, mateSamRecord.getPosition ());
	      } else {
		minStartPos = mateSamRecord.getPosition ();
	      }
	    }

	    if (minStartPos >= 1) {
	      BitSet bitSet = transcriptStartPositionTable.get(transcriptId);
	      if (bitSet != null) {
		bitSet.set (minStartPos - 1);
	      } else {
		throw new IOException ("Transcript id " + transcriptId + " not found in transcript length file.");
	      }
	    }
	  }
	}
      }
    }

    if (geneIds.size () > 0 && (! unambiguous || geneIds.size () == 1)) {
      addReadWeight (fragmentName, geneIds, readWeightReader, geneCountTable);
    }

    if (debugLevel >= 3) {
      System.out.println("Gene count table: " + geneCountTable);
    }

  }
  
}

