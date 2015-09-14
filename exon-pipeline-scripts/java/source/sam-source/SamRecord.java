/**File: SamRecord.java 

Original Author: Sven Schuierer
Date: 19/12/2011

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
 *                              Class SamRecord
 *
 ***********************************************************************************/

public class SamRecord implements Comparable<SamRecord> {

  private static int debugLevel = 0;
  
  private static boolean noSplice = false;

  private static boolean warningsOn = false;
  private static boolean printLines = false;
  private static Counter fragmentCounter = null;
  private static int bedEntryNum = 0;

  public static void setCounter (Counter counter) {
    fragmentCounter = counter;
  }
  
  private Pattern asPattern = Pattern.compile("AS:i:([+-]?[0-9]+)");
  private Pattern nmPattern = Pattern.compile("[nN]M:i:([0-9]+)");

  public static void setWarningsOn (boolean value) {
    warningsOn = value;
  }

  public static void setPrintLines (boolean value) {
    printLines = value;
  }

  private static Hashtable<String, Integer> processedReads = null;
  
  public static void setNoSplice (boolean value) {
    noSplice = value;
  }

  public static void init () {
  }

  public static int processedReads () {
    return(processedReads.size());
  }
  

  /***********************************************************************************
   *
   *                     Object variables and methods
   *
   ***********************************************************************************/


  private String samRecordString = null;
  private String originalQueryName;
  private String originalFragmentName;
  private int    flag;
  private String referenceName;
  private int    position;
  private int    mapQuality;
  private String ciagrString;
  private String mateReferenceName;
  private int    matePosition;
  private String insertSize;
  private String sequence;
  private String qualityString;
  private String optionalFields = "";

  private int editDistance   = -2;
  private int alignmentScore = Integer.MAX_VALUE;

  private int readLength    = 0;
  private int numInsertions = 0;   // Insertions and deletions are with respect to the reference
  private int numDeletions  = 0;   // see above (1I -> a "*" in the reference sequence)

  private int hitIndex = -1;
  private int numMismatches = -1;
  
  private int readIndex = -1;
  private String strand = "";
  private int flagReadIndex = -1;
  private boolean hasMate = false;
  private int endPosition = -1;

  private boolean useOriginalFragmentName = false;
  
  Vector<BedEntry> bedEntries = new Vector<BedEntry> (5);
  

  /***********************************************************************************
   *
   *                     Constructor
   *
   ***********************************************************************************/
  
  public SamRecord (String line, int lineNumber, String filename) throws IOException {

    debugLevel = UtilLib.getDebugLevel ();

    samRecordString = line;
    
    if (debugLevel >= 2) {
      System.out.println ("Creating SAM record for: " + samRecordString);
    }

    int i = 0;
    String token = "";
    
    try {

      if (samRecordString == null) {
	throw new Exception ("SamRecord constructor called on null string");
      }
      
      /* StringTokenizer st = new StringTokenizer (samRecordString, "\t");
	 while (st.hasMoreTokens()) {
	   token = st.nextToken(); */
      
      int pos = 0;
      while (pos < samRecordString.length()) {

	int end = samRecordString.indexOf ("\t", pos);
	if (end == -1) {
	  end = samRecordString.length();
	}
	token = samRecordString.substring(pos, end);
	pos = end + 1;

	if (i == 0) {
	  originalQueryName = token;

	  int spaceIndex = originalQueryName.indexOf (' ');
	  if (spaceIndex > 0) {
	    originalQueryName = originalQueryName.substring (0, spaceIndex);
	    if (token.endsWith ("/1")) {
	      originalQueryName = originalQueryName + "/1";
	    }
	    if (token.endsWith ("/2")) {
	      originalQueryName = originalQueryName + "/2";
	    }
	  }

	  originalFragmentName = originalQueryName;
	  if (originalFragmentName.endsWith ("/1")) {
	    readIndex = 1;
	    originalFragmentName = originalFragmentName.substring (0, originalFragmentName.length () - 2);
	  } else if (originalFragmentName.endsWith ("/2")) {
	    readIndex = 2;
	    originalFragmentName = originalFragmentName.substring (0, originalFragmentName.length () - 2);
	  }
	  i++;
	} else if (i == 1) {
	  flag = Integer.parseInt(token);

	  /* Determine read index from SAM flag */
	  flagReadIndex = getFlagReadIndex ();
	  if (readIndex == -1) {
	    readIndex = flagReadIndex;
	  } else if (readIndex != flagReadIndex) {
	    throw new IOException ("Read index in flag: " + flagReadIndex + " differs from name read index:" + readIndex);
	  }

	  /* Determine strand from SAM flag */
	  if (flag % 8 < 4) {
	    strand = "+";
	    if (flag % 32 >= 16) {
	      strand = "-";
	    }
	  }
	  
	  i++;
	} else if (i == 2) {
	  referenceName = token;
	  i++;
	} else if (i == 3) {
	  position = token.equals("*")?-1:Integer.parseInt(token);
	  i++;
	} else if (i == 4) {
	  mapQuality = Integer.parseInt(token);
	  i++;
	} else if (i == 5) {
	  ciagrString = token;
	  i++;	} else if (i == 6) {
	  mateReferenceName = token;
	  i++;
	} else if (i == 7) {
	  matePosition = token.equals("*")?-1:Integer.parseInt(token);
	  i++;
	} else if (i == 8) {
	  insertSize = token;
	  i++;
	} else if (i == 9) {
	  sequence = token;
	  i++;
	} else if (i == 10) {
	  qualityString = token;
	  i++;
	} else if (i >= 11) {
	  if (optionalFields != "") {
	    optionalFields = optionalFields + "\t" + token;
	  } else {
	    optionalFields = token;
	  }
	  i++;
	  if (token.indexOf(":") >= 0) {
	    i--;
	  }
	  if (token.startsWith("HI:i:")) {
	    hitIndex = new Integer(token.substring(5)).intValue ();
	  } else if (token.toLowerCase().startsWith("nm:i:")) {
	    numMismatches = new Integer(token.substring(5)).intValue ();
	    if (debugLevel >= 2) {
	      System.err.println("Setting number of mismatch for " + getFragmentName () + " to " + getNumMismatches() + " based on " + token);
	    }
	  }
	}
      }

      readLength = sequence.length();

      if (! qualityString.equals("*") && sequence.length() != qualityString.length()) {
	throw new IOException ("sequence and qualities have different lengths: " + sequence + " (" + sequence.length() + "), " + qualityString + " (" +
			       qualityString.length() + ")");
      }

      if (i < 11) {
	throw new IOException ("Too few fields in SAM record (" + i + ")");
      }

      if (i > 12) {
	throw new IOException ("Too many fields in SAM record (" + i + ")");
      }

      if (debugLevel >= 2) {
	System.err.println("Read index for " + this + ": " + readIndex);
      }

    } catch (Exception e) {
      throw new IOException ("ERROR: Problem parsing sam file in sam record:\n" + samRecordString + "\n" +
			     "at position: " + i + " with token: " + token + "\n" +
			     ((lineNumber >= 0)?("at line number " + lineNumber):"") +
			     ((filename != null)?("of file: " + filename):"") +
			     ((lineNumber >= 0 || filename != null)?"\n":"") +
			     (e==null?"Null message":e.getMessage ()));
    }
  }


  public SamRecord (String samRecordString) throws IOException {

    this (samRecordString, -1, null);
    
  }


  public SamRecord (String queryName, String fragmentName, int flag, String referenceName, int position, int mapQuality, String ciagrString,
		    String mateReferenceName, int matePosition, String insertSize, String sequence, String qualityString, String optionalFields) {
    
      this.originalQueryName    = queryName;
      this.originalFragmentName = fragmentName;
      
      this.flag              = flag;
      this.referenceName     = referenceName;
      this.position          = position;
      this.mapQuality        = mapQuality;
      this.ciagrString       = ciagrString;
      this.mateReferenceName = mateReferenceName;
      this.matePosition      = matePosition;
      this.insertSize        = insertSize;
      this.sequence          = sequence;
      this.qualityString     = qualityString;
      this.optionalFields    = optionalFields;

      this.readLength = sequence.length ();

  }


  /***********************************************************************************
   *
   *             Processing the SAM record which may contain
   *             spliced reads - note that the fragmentId can be different
   *             from the fragment name (if new fragment ids are created).
   *
   *  Note that in order to intersect the BED entries in a strand-specific manner
   *  the alignment of the second read of a pair to the reverse complement
   *  needs to be inverted.
   *
   ***********************************************************************************/

  public Vector<BedEntry> computeBedEntries (String fragmentId, String strandSpecificDirection) throws IOException {

    if (debugLevel >= 2) {
      System.out.println ("Computing bed entries for fragmentId: " + fragmentId);
    }

    if (bedEntries.size () > 0) {
      return bedEntries;
    }

    int flagReadIndex = getFlagReadIndex ();
    int bedStart = Math.max(0, position - 1);

    Vector<String> ciagrFields   = new Vector<String> (20);
    Vector<Integer> intronFields = new Vector<Integer> (20);

    if (hasMate()) {
      fragmentId = fragmentId + "/P" + flagReadIndex;
    } else {
      fragmentId = fragmentId + "/S" + flagReadIndex;
    }

    if (debugLevel >= 2) {
      System.out.println ("fragmentId: " + fragmentId);
    }

    if (noSplice) {
      ciagrFields.add(ciagrString);
    } else {

      Pattern pattern = Pattern.compile("[0-9]+N");
      Matcher match = pattern.matcher(ciagrString);
      
      int ciagrPos = 0;
      while(match.find()) {
	ciagrFields.add(ciagrString.substring(ciagrPos, match.start()));
	intronFields.add(new Integer(ciagrString.substring(match.start(), match.end()-1)));
	ciagrPos = match.end();
      }

      if (ciagrPos == 0) {
	ciagrFields.add(ciagrString);
      } else {
	ciagrFields.add(ciagrString.substring(ciagrPos, ciagrString.length()));
      }
    }

    /* We need as many intron fields as ciagr fields */
    intronFields.add (new Integer (0));

    if (intronFields.size() > 0 && false) {
      fragmentId = "Splice-" + fragmentId;
    }

    int ciagrFieldNum  = ciagrFields.size();
    
    int curCiagrFieldNum     = 0;
    int sumRefAlignedLength  = 0;
    int sumReadAlignedLength = 0;
    int sumIntronLength      = 0;

    for (String ciagrField: ciagrFields) {

      if (debugLevel >= 1) {
	System.out.println("ciagrField: " + ciagrField + ", bedStart: " + bedStart + ", sumRefAlignedLength: " + sumRefAlignedLength +
			   ", sumIntronLength: " + sumIntronLength);
      }
      
      int start = bedStart + sumRefAlignedLength + sumIntronLength;

      if (debugLevel >= 1) {
	System.out.println("start: " + start + ", sumIntronLength  = "   + sumIntronLength   + " + " + (new Integer(intronFields.get(curCiagrFieldNum))).intValue());
	sumIntronLength = sumIntronLength + (new Integer(intronFields.get(curCiagrFieldNum))).intValue();
      }

      if ((ciagrFieldNum > 1 && debugLevel >= 1)) {
	System.err.println("CIAGR field for " + originalQueryName + ": " + ciagrField);
      }

      
      Pattern pattern2 = Pattern.compile("[0-9]+[A-Z=]"); // 
      Matcher match    = pattern2.matcher(ciagrField);
      
      numInsertions = 0;
      numDeletions  = 0;
      int readPartLength = 0;
      int softClippingLength = 0;
      int hardClippingLength = 0;
      while(match.find()) {
	int    fieldLength    = Integer.parseInt (ciagrField.substring(match.start(), match.end()-1));
	String fieldOperation = ciagrField.substring(match.end()-1, match.end());
	if (fieldOperation.equals("M") || fieldOperation.equals("X") || fieldOperation.equals("=")) {
	  readPartLength = readPartLength + fieldLength;
	} else if (fieldOperation.equals("I")) {
	  numInsertions  = numInsertions + fieldLength;
	} else if (fieldOperation.equals("D")) {
	  numDeletions = numDeletions + fieldLength;
	} else if (fieldOperation.equals("S")) {
	  softClippingLength = softClippingLength + fieldLength;
	} else if (fieldOperation.equals("H")) {
	  hardClippingLength = hardClippingLength + fieldLength;
	} else if (fieldOperation.equals("N")) {
	  throw new IOException ("ERROR: N discovered in unspliced region of " + fragmentId);
	} else if (! fieldOperation.equals("P")) {
	  throw new IOException ("ERROR: Unknown operation in CIAGR string for fragment " + fragmentId + ": " + fieldLength + fieldOperation);
	}

	if (debugLevel >= 2) {
	  System.out.println("fieldLength: " + fieldLength + ", fieldOperation: " + fieldOperation);
	}
      }

      int refAlignedLength  = readPartLength + numDeletions;
      int readAlignedLength = readPartLength + numInsertions + softClippingLength;

      String strand = getStrand();
      if (strandSpecificDirection.equals("forward")) {
	strand = getFragmentStrand ();
      } else if (strandSpecificDirection.equals("backward")) {
	strand = getReverseFragmentStrand ();
      }

      if (readAlignedLength > 0) {
	bedEntries.add(new BedEntry(referenceName, start, start + refAlignedLength, fragmentId, bedEntryNum, strand, readAlignedLength,
				    numInsertions + softClippingLength, numDeletions, curCiagrFieldNum, ciagrFieldNum));
	bedEntryNum += 1;
      }

      if (debugLevel >= 2 ) {
	System.out.println("sumRefAlignedLength = "   + sumRefAlignedLength  + " + " + refAlignedLength +
			   ", sumReadAlignedLength: " + sumReadAlignedLength + " + " + readAlignedLength);
      }

      sumRefAlignedLength  = sumRefAlignedLength  + refAlignedLength;
      sumReadAlignedLength = sumReadAlignedLength + readAlignedLength;
      sumIntronLength      = sumIntronLength + intronFields.get(curCiagrFieldNum).intValue ();

      curCiagrFieldNum++;	

    }
    

    if (debugLevel >= 2 ) {
      System.err.println("Done: " + sumReadAlignedLength + " vs " + readLength);
      System.err.flush();
    }
    
    if (sumReadAlignedLength != readLength) {
      throw new IOException ("Fragment " + fragmentId + " mapped to " + referenceName + ": read length " + readLength +
			     " differs from sum of matches, insertions, and soft clipped bases " + sumReadAlignedLength);
    }

    if (debugLevel >= 2) {
      System.err.println (bedEntries.size() + " BED entries computed.");
      System.err.flush();
    }
    
    return bedEntries;

  }


  /***********************************************************************************
   *
   *                     printBedRecords (possibly with mate)
   *
   ***********************************************************************************/

  public void printBedRecords (SamRecord mateSamRecord, String fragmentId, Counter alignmentCounter,
			       HashSetTable<String, String> transcriptVersionIds, int alignmentScoreThreshold,
			       int editDistanceThreshold, String strandSpecificDirection, PrintWriter out) throws IOException {

    if (mateSamRecord != null) {
      if (mateSamRecord.getEditDistance () > editDistanceThreshold || mateSamRecord.getAlignmentScore () < alignmentScoreThreshold) {
	mateSamRecord = null;
      }
    }
	
    if (getEditDistance () > editDistanceThreshold || getAlignmentScore () < alignmentScoreThreshold) {
      if (mateSamRecord != null) {
	mateSamRecord.printBedRecords(null, fragmentId, alignmentCounter, transcriptVersionIds, alignmentScoreThreshold, editDistanceThreshold, strandSpecificDirection, out);
      }
      return;
    }
    
    
    if (debugLevel >= 2) {
      System.out.println("Fragment id: " + fragmentId + " printing bed records for " + getQueryName() + " and " +
			 (mateSamRecord==null?"null":mateSamRecord.getQueryName()));
    }

    Vector<BedEntry> bedEntries = computeBedEntries (fragmentId, strandSpecificDirection);

    if (debugLevel >= 2) {
      System.out.println("Number of bed entries: " + bedEntries.size());
    }
    
    Vector<BedEntry> mateBedEntries = new Vector<BedEntry> ();
    if (mateSamRecord != null) {
      mateBedEntries = mateSamRecord.computeBedEntries(fragmentId, strandSpecificDirection);
    }

    HashSet<String> transcriptVersionIdSet = null;
    if (transcriptVersionIds != null) {
      transcriptVersionIdSet = transcriptVersionIds.getSet(referenceName);
    }

    if (transcriptVersionIdSet == null) {
      transcriptVersionIdSet = new HashSet<String> (1);
      transcriptVersionIdSet.add(referenceName);
    }


    if (debugLevel >= 2) {
      System.out.println("Number of transcriptVersionIds: " + transcriptVersionIdSet.size() + ", number of bedEntries: " + bedEntries.size());
    }

    if (debugLevel >= 2) {
      System.out.println(referenceName + ", " + mateSamRecord + ", " + mateReferenceName);
      if (mateSamRecord != null) {
	System.out.println(referenceName + ", " + mateSamRecord.getMateReferenceName ()  + ", " + mateReferenceName  + ", " + mateSamRecord.getReferenceName ());
      }
    }
    

    if (mateSamRecord == null || (referenceName.equals (mateSamRecord.getMateReferenceName ()) &&
				  getMateReferenceName().equals (mateSamRecord.getReferenceName ()))) {
      for (String transcriptVersionId: transcriptVersionIdSet) {

	String alignmentId = "-A" + alignmentCounter.getZeroFilledCount();

	if (debugLevel >= 2) {
	  System.out.println("Printing bed records for transcriptVersionId: " + transcriptVersionId + " and alignmentId: " + alignmentId);
	}

	for (BedEntry bedEntry: bedEntries) {
	  out.println(bedEntry.toString(alignmentId, transcriptVersionId));
	}

	for (BedEntry bedEntry: mateBedEntries) {
	  out.println(bedEntry.toString(alignmentId, transcriptVersionId));
	}

	out.flush();

	alignmentCounter.inc();
      
      }

      alignmentCounter.dec();
      
    } else {
      
      throw new IOException ("Differing referenceName mate found for: " + this + "\n" + "mate: " + mateSamRecord);
      
    }
    
  }


  /***********************************************************************************
   *
   *                          findMate
   *
   ***********************************************************************************/
  
  public SamRecord findMate (Vector<SamRecord> samRecords, int i, HashSet<Integer> processedIndices) {

    boolean mateFound = false;
    SamRecord mateSamRecord = null;
    Integer jInt = null;
    for (int j = i + 1; j < samRecords.size() && ! mateFound; j++) {

      jInt = Integer.valueOf (j);
      if (! processedIndices.contains(jInt)) {
	mateSamRecord = samRecords.get(j);
	mateFound     = isMate(mateSamRecord);
      }	    
    }

    if (! mateFound) {
      if (warningsOn) {
	System.err.println("WARNING: No matching entry for " + getQueryName() + " at position " +
			   getReferenceName() + "/" + getPosition() + " found.");
      }
      return null;
    }
    
    processedIndices.add(jInt);

    return mateSamRecord;
    
  }

  
  /***********************************************************************************
   *
   *                     Getting and setting readIndex
   *
   ***********************************************************************************/

  public int getReadIndex () {
    return readIndex;
  }

  public int getFlagReadIndex () {

    int flagReadIndex = 1;
    if ((getFlag () % 2) == 1) {
      if (flag % 128 >= 64) {
	flagReadIndex = 1;
      } else if (flag % 256 >= 128) {
	flagReadIndex = 2;
      }
    } else if (originalQueryName.endsWith("/1")) {
      flagReadIndex = 1;
    } else if (originalQueryName.endsWith("/2")) {
      flagReadIndex = 2;
    } else if (getReadIndex () != -1) {
      flagReadIndex = getReadIndex ();
    }

    return (flagReadIndex);
    
  }


  /***********************************************************************************
   *
   *                     Checking for a mate
   *
   ***********************************************************************************/


  public boolean isMate (SamRecord samRecord) {

    if (debugLevel >= 2) {
      System.out.println("Testing: " + getMateReferenceName() + " = " + samRecord.getReferenceName()+ ", " + getMatePosition() + " = " +
			 samRecord.getPosition() + ", " + getReferenceName() + " = " + samRecord.getMateReferenceName()+ ", " +
			 getPosition() + " = " + samRecord.getMatePosition());
    }

    if (hitIndex >= 0) {
      return (samRecord.getHitIndex() == hitIndex);
    }

    return (samRecord.getFlagReadIndex () != getFlagReadIndex() && originalQueryName.equals(samRecord.getOriginalQueryName()) &&
	    getMateReferenceName().equals(samRecord.getReferenceName()) && getMatePosition() == samRecord.getPosition() &&
	    getReferenceName().equals(samRecord.getMateReferenceName()) && getPosition() == samRecord.getMatePosition());
    
  }


  /***********************************************************************************
   *
   *                          getEndPosition
   *
   ***********************************************************************************/
  
  public int getEndPosition () throws IOException {

    if (endPosition >= 0) {
      return endPosition;
    }

    String strandSpecificDirection = "none";
    Vector<BedEntry> bedEntries = computeBedEntries (getFragmentName(), strandSpecificDirection);

    if (debugLevel >= 2) {
      System.err.println (bedEntries.size() + " BED entries computed for getEndPosition.");
    }

    endPosition = bedEntries.get(bedEntries.size() - 1).getEnd ();

    return endPosition;
    
  }
  
  
  /***********************************************************************************
   *
   *                          containsInterval
   *
   ***********************************************************************************/
  
  public boolean containsInterval (int leftPosition, int rightPosition) throws IOException {

    if (rightPosition <= 0) {
      return true;
    }

    if (getPosition () > leftPosition) {
      return false;
    }
    
    if (rightPosition <= getEndPosition ()) {
      return true;
    }

    return false;
    
  }


  /***********************************************************************************
   *
   *                          overlap
   *
   * Compute the overlap with a genomic SAM record s (s should not be spliced)
   *
   ***********************************************************************************/
  
  public double overlap (SamRecord s) throws IOException {

    if (s == null) {
      return 0;
    }

    if (! getReferenceName ().equals(s.getReferenceName ())) {
      return 0;
    }

    int samStart = s.getPosition ();
    int samEnd   = s.getEndPosition ();
    
    if (getEndPosition () < samStart || samEnd < getPosition()) {
      return 0;
    }

    String strandSpecificDirection = "none";
    Vector<BedEntry> bedEntries = computeBedEntries (getFragmentName (), strandSpecificDirection);

    if (debugLevel >= 2) {
      System.err.println (bedEntries.size() + " BED entries computed for overlap.");
    }

    int overlap = 0;
    for (BedEntry bedEntry: bedEntries) {
      overlap += Math.max(0, Math.min(samEnd, bedEntry.getEnd()) - Math.max(samStart, bedEntry.getStart()) + 1);
    }

    if (debugLevel >= 2) {
      System.err.println ("overlap: " + overlap);
    }
    
    return overlap;
   
  }


  /***********************************************************************************
   *
   *                          getMappingDistance
   *
   ***********************************************************************************/
  
  public int getAlignmentDistance (SamRecord s) throws IOException {

    if (s == null) {
      return -1;
    }
    
    if (! getReferenceName ().equals(s.getReferenceName ())) {
      return -1;
    }

    return Math.max(Math.abs(getPosition () - s.getPosition ()), Math.abs(getEndPosition () - s.getEndPosition ()));
    
  }

  /***********************************************************************************
   *
   *                              Get methods
   *
   ***********************************************************************************/

  public Vector<BedEntry> getBedEntries () {
    return (bedEntries);
  }

  
  public String getReadId () {
    if (originalQueryName.length() == 10 && originalQueryName.startsWith("F")) {
      return originalQueryName + "/" + (isFirstRead()?"1":"2");
    } else {
      return getQueryName ();
    }
  }

  public String getQueryName () {
    if (fragmentCounter == null || useOriginalFragmentName) {
      return originalQueryName;
    }
    
    return "F" + fragmentCounter.getZeroFilledCount () + "/" + (isFirstRead()?"1":"2");
  }

  public String getOriginalQueryName () {
    return originalQueryName;
  }


  public String getFragmentName () {
    if (fragmentCounter == null || useOriginalFragmentName) {
      return originalFragmentName;
    }
    
    return "F" + fragmentCounter.getZeroFilledCount ();
  }

  public void fixFragmentName () {
    if (fragmentCounter != null) {
      originalFragmentName = "F" + fragmentCounter.getZeroFilledCount ();
      originalQueryName = "F" + fragmentCounter.getZeroFilledCount () + "/" + (isFirstRead()?"1":"2");
      useOriginalFragmentName = true;
    }
  }


  public String getOriginalFragmentName () {
    return originalFragmentName;
  }

  public int getFlag () {
    return flag;
  }

  public void setBitFlag (int i, boolean b) {
    flag = UtilLib.setBit (flag, i, b);
  }

  public String getReferenceName () {
    return referenceName;
  }

  public int getPosition () {
    return position;
  }

  public String getStrand () {
    return strand;
  }

  public String getFragmentStrand () {
    if (readIndex == 1) {
      return strand;
    } else {
      return (strand.equals("+")?"-":"+");
    }
  }

  public String getReverseFragmentStrand () {
    if (readIndex == 2) {
      return strand;
    } else {
      return (strand.equals("+")?"-":"+");
    }
  }

  public int getMapQuality () {
    return mapQuality;
  }

  public String getCiagrString () {
    return ciagrString;
  }

  public boolean isSpliced () {
    return ciagrString.indexOf("N") > 0;
  }


  public String getMateReferenceName () {
    if (mateReferenceName.equals("=")) {
      return referenceName;
    } else {
      return mateReferenceName;
    }
  }

  public String getMateReferenceNameForOutput () {
    if (mateReferenceName.equals("=")) {
      return "=";
    }

    if (mateReferenceName.equals(referenceName) && ! referenceName.equals("*")) {
      return "=";
    }
    
    return mateReferenceName;
  }
  
  public void setMateReferenceName (String referenceName) {
    mateReferenceName = referenceName;
  }


  public int getMatePosition () {
    return matePosition;
  }

  public void setMatePosition (int i) {
    matePosition = i;
  }

  public String getInsertSize () {
    return insertSize;
  }

  public void setInsertSize (String s) {
    insertSize = s;
  }
  
  public String getSequence () {
    return sequence;
  }

  public int getLength () {
    return sequence.length();
  }
  
  public String getQualityString () {
    return qualityString;
  }

  public String getOptionalParameters () {
    return optionalFields;
  }

  public int getHitIndex () {
    return hitIndex;
  }

  public int getNumMismatches () {
    return numMismatches;
  }

  
  /***********************************************************************************
   *
   *                              Flag values
   *
   ***********************************************************************************/

  public boolean isPairedInSequencing () {
    return (getFlag () % 2) == 1;
  }

  public boolean hasMate () {
    return ((getFlag () % 4) >= 2) && ((getFlag () % 2) == 1);
  }

  public void setHasMate (boolean b) {
    flag = UtilLib.setBit(flag, 1, b);
  }

  public boolean isUnmapped () {
    return ((getFlag () % 8) >= 4);
  }

  public void setIsUnmapped (boolean b) {
    flag = UtilLib.setBit(flag, 2, b);
  }

  public boolean isMapped () {
    return (getPosition() >= 0 && (getFlag () % 8) < 4);
  }

  public boolean mateIsMapped () {
    return (getFlag () % 16) < 8;
  }

  public void setMateIsUnmapped (boolean b) {
    flag = UtilLib.setBit(flag, 3, b);
  }

  public boolean isOnReverseStrand () {
    return (getFlag () % 32) >= 16;
  }

  public void setIsOnReverseStrand (boolean b) {
    flag = UtilLib.setBit(flag, 4, b);
  }
  
  public boolean mateIsOnReverseStrand () {
    return (getFlag () % 64) >= 32;
  }

  public void setMateIsOnReverseStrand (boolean b) {
    flag = UtilLib.setBit(flag, 5, b);
  }  

  public boolean isFirstRead () {
    return (getFlag () % 128) >= 64;
  }
  
  public boolean isSecondRead () {
    return (getFlag () % 256) >= 128;
  }

  public boolean isPrimary () {
    return ((getFlag () % 512) < 256);
  }

  public void setIsSecondaryAlignment (boolean b) {
    flag = UtilLib.setBit(flag, 8, b);
  }
  
  
  /***********************************************************************************
   * 
   *                          addOptionalField
   *
   ***********************************************************************************/

  public void addOptionalField (String fieldName, String fieldType, String fieldValue) {

    String separator = "\t";
    if (optionalFields.indexOf(fieldName) > 0) {
      Pattern fieldPattern = Pattern.compile(fieldName + ":[^" + separator + "]*");
      Matcher fieldMatch = fieldPattern.matcher(optionalFields);
      if (fieldMatch.find()) {
	optionalFields = optionalFields.substring(0, Math.max(0, fieldMatch.start(1) - 1)) + optionalFields.substring(fieldMatch.end(1));
      }
    }

    optionalFields = optionalFields + separator + fieldName + ":" + fieldType + ":" + fieldValue;
    
  }

  /***********************************************************************************
   *
   *                Edit distance and alignment score
   *
   ***********************************************************************************/
  
  public int getEditDistance () {

    if (editDistance != -2) {
      return editDistance;
    }

    Matcher nmMatch = nmPattern.matcher(optionalFields);
    if (nmMatch.find()) {
      editDistance = new Integer(optionalFields.substring (nmMatch.start(1), nmMatch.end(1)));
    }

    if (editDistance == -2) {
      editDistance = -1;
    }

    return editDistance;
  }


  public boolean isDefinedEditDistance () {
    return getEditDistance () >= 0;
  }


  public int getAlignmentScore () {

    if (alignmentScore < Integer.MAX_VALUE) {
      return alignmentScore;
    }

    Matcher asMatch = asPattern.matcher(optionalFields);
    if (asMatch.find()) {
      alignmentScore = new Integer(optionalFields.substring (asMatch.start(1), asMatch.end(1)));
    }

    if (alignmentScore == Integer.MAX_VALUE) {
      alignmentScore = Integer.MAX_VALUE - 1;
    }

    return alignmentScore;
    
  }

  public boolean isDefinedAlignmentScore () {
    return getAlignmentScore () < Integer.MAX_VALUE - 1;
  }

  
  /***********************************************************************************
   *
   *                toUnmappedString
   *
   ***********************************************************************************/
  
  public String toUnmappedString (boolean mateIsMapped) throws IOException {

    String queryName = getQueryName ();

    int flag = (isPairedInSequencing ()?1:0) + 4 + (hasMate()?(mateIsMapped?0:8):8) + (isFirstRead ()?64:0) +  + (isSecondRead ()?128:0);

    queryName = UtilLib.modifyFragmentId (queryName);
    
    return (queryName + "\t" + flag + "\t" + "*" + "\t" + 0 + "\t" + 0 + "\t" + "*" + "\t" +
	    "*" + "\t" + 0 + "\t" + 0 + "\t" + getSequence () + "\t" + getQualityString ());
  }


  /***********************************************************************************
   *
   *              compareTo
   *
   ***********************************************************************************/
  
  public int compareTo (SamRecord s) {
    return getReadId ().compareTo(s.getReadId ());
  }


  /***********************************************************************************
   *
   *                toString with mate information
   *
   ***********************************************************************************/
  
  public String toString (boolean mateIsMapped, boolean isSecondary, int level) throws IOException {

    String queryName = getQueryName ();

    queryName = UtilLib.modifyFragmentId (queryName);

    setIsSecondaryAlignment(isSecondary);

    String mateReferenceName = "*";
    int    matePosition      = 0;
    int    flag              = getFlag ();
    if (mateIsMapped) {
      if (mateIsMapped ()) {
	mateReferenceName = getMateReferenceNameForOutput ();
	matePosition      = getMatePosition ();
      } else if (level == 0) {
	throw new IOException ("Mate is mapped in toString though not in original alignment for " + this);
      }
    } else if (isPairedInSequencing() && mateIsMapped ()) {
      flag = flag + 8;
    }

    return ( queryName        + "\t" + flag              + "\t" + getReferenceName () + "\t" + getPosition () + "\t" +
	     getMapQuality () + "\t" + getCiagrString () + "\t" + mateReferenceName   + "\t" + matePosition   + "\t" +
	     getInsertSize () + "\t" + getSequence ()    + "\t" + getQualityString () + "\t" + getOptionalParameters ());

  }



  /***********************************************************************************
   *
   *                toString with a suffix to the query name
   *
   ***********************************************************************************/
  
  public String toString (String suffix) throws IOException {

    String queryName = getQueryName ();

    queryName = UtilLib.modifyFragmentId (queryName);

    if (suffix != "" && ! queryName.endsWith (suffix)) {
      queryName = queryName + suffix;
    }
    
    return ( queryName        + "\t" + getFlag ()        + "\t" + getReferenceName ()     + "\t" + getPosition ()     + "\t" +
	     getMapQuality () + "\t" + getCiagrString () + "\t" + getMateReferenceName () + "\t" + getMatePosition () + "\t" +
	     getInsertSize () + "\t" + getSequence ()    + "\t" + getQualityString ()     + "\t" + getOptionalParameters ());
  }

  

  /***********************************************************************************
   *
   *                toString
   *
   ***********************************************************************************/
  
  public String toString() {

    if (samRecordString != null) {
      return samRecordString;
    }
    
    try {
      return (toString(""));
    }
    catch (IOException e) {
      System.err.println ("ERROR: " + (e==null?"Null message":e.getMessage()));
      System.exit (1);
    }

    return null;
  }
  

  /***********************************************************************************
   *
   *                hashCode
   *
   ***********************************************************************************/

  public int hashCode () {
    
    return (originalQueryName + " - " + referenceName + " - " + position  + " - " + ciagrString).hashCode ();

  }


  /***********************************************************************************
   *
   *                Equals
   *
   ***********************************************************************************/
  
  public boolean equals (Object o) {

    if (o == null) {
      return false;
    }

    SamRecord s = (SamRecord) o;

    if (originalQueryName.equals(s.getOriginalQueryName ()) && referenceName.equals (s.getReferenceName ()) && getPosition () == s.getPosition () &&
	ciagrString.equals (s.getCiagrString ())) {
      return true;
    }

    return false;
  }


  /***********************************************************************************/

  public static void main (String [] args) {
    
    String samRecordString =
      "NRCHBS-WDL30299:95:D0902ACXX:1:1101:1096:154794 137     chr1    39951394       50       100M    *       0       0       NCCTTCCAAAATCCCAACCATGTCTAAGAAGACCACCACTGCCTCCCCCAGGACTCCAGGTCCCAAGCGATAACACTGTCTAAGCACCCCCAAGCCACTA    #1=BDDDFGFBFD<CBG@CGGIJIGEIIHJJGIDIIHHJJGICGGGDHIICFHHIGHHGHIIFH>;ADDBC>?;AC>;>;ABCCDDA9,899?AB??CCC    AS:i:255        XN:i:0  XM:i:1  XO:i:0  XG:i:0  NM:i:1  MD:Z:0A99      YT:Z:UU  XS:A:+  NH:i:1\n" +
      "F000000022	163	chr6	118024609	255	2S7M1I8P2D257N35M1D5=	chr6	118024833	274	GTGGGGCATGGAGACATACAGATAAAACTTCTGAAGAAAACAGGTAGTTG	HHHHD44555;@A??DDD>>>>+:=4445554455HHHHGBHFHE#####	XA:i:1 MD:Z:49T0 NM:i:1\n" + 
      "F000000022	83	chr6	118024833	255	50M	chr6	118024609	-274	TTTCTTCCCTGGACCATCAGATTGACTGAGATTGTGTAAGTAATTAAAAG	#@A<A44525555.555&4555555544.5@=/@>555555'55444535	XA:i:2 MD:Z:12C0A36 NM:i:2\n" + 
      "F000000029	99	chrM	1	255	20M100N30M	chrM	225	274	CCGCAACCTCAACACCACCTTCTTCGACCCCGCCGGAGGAGGAGACCCCA	HHHHHHHHGHHHHHHHHHHHEFFFFEHHHHFAF>7:C93<GGHGHHHHHH	XA:i:0 MD:Z:50 NM:i:0\n" + 
      "F000000029	147	chrM	225	255	50M	chrM	1	-274	TGTGAGCACACNATATATTTACAGTAGGAATAGACGTAGACACACGAGCG	HHHHHA@>?@G#GFGHHHHHIHHHHHHHHHHHHHHHHHHHHHHHHHHHHH	XA:i:2 MD:Z:11C37A0 NM:i:2\n" + 
      "F000000043	99	chr7	107398703	255	50M	chr7	107398936	283	ATCAACCATCGCCATATTAGAGCTGGAAAACCTGTTACCCGTGCTTCACT	5,55555555BFEF846145A:A=@HHHHHDD7A8+4-54+1445DDDD7	XA:i:1 MD:Z:17G32 NM:i:1\n" + 
      "F000000043	147	chr7	107398936	255	50M	chr7	107398703	-283	TCCAGCAGAATTGTCCATGGCTCCACCTCCACCTCGATCGGTCAGTCAGG	FFHHHHHHAHHHH8HHDHHH4444(@A?AA;HHHHHHHHHHHHFHHHHHH	XA:i:0 MD:Z:50 NM:i:0\n" + 
      "F000000046	163	chr15	101455195	255	50M	chr15	101455405	260	CTGGGGAGGGAGCTGTTGGCCATTTCTGTGTTTCCCTTNAAACCAGATCC	D:ADAFFF<BBHH7G###################################	XA:i:1 MD:Z:38T11 NM:i:1\n" + 
      "F000000046	83	chr15	101455405	255	50M	chr15	101455195	-260	GTGGACCTGCTTACAGAGCAAGCCACGCCTCTTTCCGAGGTGAAGGTGGG	########################################DD82CHHHHH	XA:i:0 MD:Z:50 NM:i:0\n" + 
      "F000000050	163	chr18	29128265	255	50M	chr18	29128481	266	AATATTTCCCACAACAATTTTCATAATTTTCAAAGACTAATTTCTTGACT	HHHHHHBHHHHHHHHHHHHHFGFGEHHHHHHHHHIGHHHHHHIHHHFFHH	XA:i:0 MD:Z:50 NM:i:0\n" + 
      "F000000050	83	chr18	29128481	255	50M	chr18	29128265	-266	TCCATCATCTAGAATTGTTTACTTAGTAATTGTTGTTTCTTTTATTATTA	HHHH@HHHHHHHHHBHHHHHHHHHHHHHHHHHHHBEEHEEHHHHFDADAD	XA:i:0 MD:Z:50 NM:i:0\n" +
      "NRCHBS-WDL30299:95:D0902ACXX:1:1206:20154:6524	137	chr11	123485490	50	54M3509N12D336N41M1I4M	*	0	0	CAACGGTTTCCACCTGCAGAGCGTGTCCAAGCTGCTGCTGGTTATCAGCTGTGTTCTGGTGCTGCTGGTCATCCTTAACATGATGCTCTTCTACAGATCG	BC@FFFFDHHHGHGIJIEGDIIJC@CGGIIIJIIGGECGGIGHIICII@@GGHIGCDAEHIJIIJJEEHFFFFE<CEEEEDCEEDDDCDD@CCDDDDDDB	AS:i:197        XN:i:0  XM:i:2  XO:i:2  XG:i:13 NM:i:15 MD:Z:54^GATCTGTTTCAG42A1T0      YT:Z:UU XS:A:+  NH:i:1\n" +
      "F000000001	137	chr11	123485490	50	54M3509N12M336N29M1I4M	*	0	0	CAACGGTTTCCACCTGCAGAGCGTGTCCAAGCTGCTGCTGGTTATCAGCTGTGTTCTGGTGCTGCTGGTCATCCTTAACATGATGCTCTTCTACAGATCG	BC@FFFFDHHHGHGIJIEGDIIJC@CGGIIIJIIGGECGGIGHIICII@@GGHIGCDAEHIJIIJJEEHFFFFE<CEEEEDCEEDDDCDD@CCDDDDDDB	AS:i:197        XN:i:0  XM:i:2  XO:i:2  XG:i:13 NM:i:15 MD:Z:54^GATCTGTTTCAG42A1T0      YT:Z:UU XS:A:+  NH:i:1";

    samRecordString = "3:1101:1307:2161	345	84134-junction-01-04::	6	32	90M2I3M1I4M	=	6	0	CAGTTCTGCCATCACTCAAGATGGCTGCCCCCATCAAGATGACCGGGGTGTGCCGGGGGGAAAGGGGCAGCATGATGGTCTGAGATGGTGTAGCGTCGGA	*	AS:i:-29";

    int leftPosition = -1;
    int rightPosition = -1;

    Getopt g = new Getopt("SamRecord", args, "p:P:h");
    
    int c;
    String arg = "";
    
    c = g.getopt();
    
    while (c  != -1) {
      switch(c) {
      case 'p':
	arg = g.getOptarg();
	leftPosition = Integer.parseInt(arg);
	break;
      case 'P':
	arg = g.getOptarg();
	rightPosition = Integer.parseInt(arg);
	break;
      case 'h':
	System.out.println("No help available");
	System.exit(0);
	break;
      default:
	System.out.print("Error: getopt() returned " + c + "\n");
      }
      c = g.getopt();
    }
    
    try {

      PrintWriter writer = new PrintWriter (System.out);
      StringTokenizer st = new StringTokenizer (samRecordString, "\n");

      Counter counter = new Counter (3);

      String[] chromosomes = {"chr6", "chrM", "chr7", "chr15", "chr18" };
      HashSetTable<String, String> transcriptVersionIds = new HashSetTable<String, String> ();

      for (int i = 0; i < chromosomes.length; i++) {
	transcriptVersionIds.putValue (chromosomes[i], chromosomes[i]);	
      }

      while (st.hasMoreTokens ()) {

	String token = st.nextToken ();

	SamRecord samRecord = new SamRecord (token);

	if (samRecord.containsInterval(leftPosition, rightPosition)) {
	  System.out.println("SamRecord contains the interval [" + leftPosition + "," + rightPosition + "]");
	} else {
	  System.out.println("SamRecord does not contain the interval [" + leftPosition + "," + rightPosition + "]");
	}

	// Vector<BedEntry> bedEntryVector = samRecord.computeBedEntries (samRecord.getFragmentName ());

	// System.out.println (token);

	
	
	samRecord.printBedRecords (null, samRecord.getFragmentName (), counter, null, Integer.MIN_VALUE, Integer.MAX_VALUE, "none", writer);
      
      }

      writer.close();
    }
    catch (IOException e) {
      System.err.println (e);
    }
  }
}
