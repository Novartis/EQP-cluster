/**File: UtilLib.java 

Original Author: Sven Schuierer
Date: 04/01/2012

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
import java.util.zip.*;

/***********************************************************************************
 *
 * 
 *                           Class UtilLib
 *
 *  
 ***********************************************************************************/


public class UtilLib {

  private static int debugLevel = 0;

  private static String idCutOffString = "";
  private static String commonPrefix = "";
  private static int    commonPrefixLen = -1;
  private static boolean warningsOn;
  
  public static void setDebugLevel (int value) {
    debugLevel = value;
  }

  public static int getDebugLevel () {
    return debugLevel;
  }

  public static void setWarningsOn (boolean value) {
    warningsOn = value;
  }

  public static boolean warningsOn () {
    return warningsOn;
  }
  

  /* The commonPrefix and the idCutOffString are only set in MergePairedSamFiles, ConvertSamBed, and
     ComputeReadWeightSam */
  public static void setCommonPrefix (String value) {
    commonPrefix = value;
  }

  public static void setIdCutOffString (String value) {
    idCutOffString = value;
  }

  /***********************************************************************************
   * 
   *  Exons whose chromosome starts with tr- are artificial exons used to cover transcripts
   *
   ***********************************************************************************/
  
  private static final String transcriptExonPrefix = "tr-";

  public static String getTranscriptExonPrefix () {
    return transcriptExonPrefix;
  }

  public static boolean isTranscriptExon (String chromosome) {
    return chromosome.startsWith(transcriptExonPrefix);
  }


  private static String invertStrandSuffix = ".inv-strand";
  public static String getInvertStrandSuffix () {
    return invertStrandSuffix;
  }


  /***********************************************************************************
   * 
   * Set the slack for the difference of an alignment against the genome with exon
   * start and end coordinates; important for determining if an alignment of (a spliced
   * part of) a read conforms with an exon or not.
   *
   ***********************************************************************************/
  
  private static int exonStartSlack = 0;
  private static int exonEndSlack   = 0;
  
  public static void setExonSlack (int value) {
    exonStartSlack = value;
    exonEndSlack = value;
  }
  
  public static void setExonStartSlack (int value) {
    exonStartSlack = value;
  }

  public static int getExonStartSlack () {
    return exonStartSlack;
  }

  public static void setExonEndSlack (int value) {
    exonEndSlack = value;
  }

  public static int getExonEndSlack () {
    return exonEndSlack;
  }


  /***********************************************************************************
   * 
   * Set the threshold for the overlap that is required for to count a read for an exon
   *
   ***********************************************************************************/
  
  private static int overlapThreshold = 1;
  
  public static void setOverlapThreshold (int value) {
    overlapThreshold = value;
  }

  public static int getOverlapThreshold () {
    return overlapThreshold;
  }


  /***********************************************************************************
   * 
   * Set strandedMode
   *
   ***********************************************************************************/
  
  private static Boolean strandedMode = false;
  public static void setStrandedMode () {
    strandedMode = true;
  }

  public static boolean strandedMode () {
    return strandedMode;
  }

  
  /***********************************************************************************
   * 
   * Set geneCountMode
   *
   ***********************************************************************************/
  
  private static String countMode = "none";

  public static void setCountMode (String mode) {
    countMode = mode;
  }

  public static String getCountMode () {
    return countMode;
  }
  
  public static boolean getGeneCountMode () {
    return countMode.equals("gene");
  }

  
  /***********************************************************************************
   * 
   * Set spliceConformingCountMode
   *
   ***********************************************************************************/
  
  private static boolean spliceConformingCountMode = true;
  public static void setSpliceConformingCountMode (boolean value) {
    spliceConformingCountMode = value;
  }

  public static boolean spliceConformingCountMode () {
    return spliceConformingCountMode;
  }

  
  /***********************************************************************************
   * 
   * Set containmentMode
   *
   ***********************************************************************************/
  
  private static boolean containmentMode = false;
  public static void setContainmentMode (boolean value) {
    containmentMode = value;
  }

  public static boolean containmentMode () {
    return containmentMode;
  }

  /***********************************************************************************
   * 
   * Set excludeAmbiguousReads
   *
   ***********************************************************************************/
  
  private static boolean excludeAmbiguousReads = false;

  public static void setExcludeAmbiguousReads (boolean value) {
    excludeAmbiguousReads = value;
  }

  public static boolean getExcludeAmbiguousReads () {
    return excludeAmbiguousReads;
  }


  /***********************************************************************************
   * 
   *                           compareInt
   *
   ***********************************************************************************/
  
  public static int compareInt (int a, int b) {
    
    if (a > b) {
      return 1;
    }

    if (b > a) {
      return -1;
    }
    
    return 0;
    
  }


  /***********************************************************************************
   * 
   *                    convert strings to numeric values
   *
   ***********************************************************************************/
  
  public static int toInt (String s) throws IOException {
    
    try {
      return Integer.parseInt(s);
    } catch(NumberFormatException nfe) {
      throw new IOException ("String " + s + " cannot be converted to an integer.");
    }
    
  }

  public static double toDouble (String s) throws IOException {
    
    try {
      return Double.parseDouble(s);
    } catch(NumberFormatException nfe) {
      throw new IOException ("String " + s + " cannot be converted to a numeric value.");
    }
    
  }


  /***********************************************************************************
   * 
   *                           set and unset a bit in an integer
   *
   ***********************************************************************************/
  
  public static int setBit (int f, int i) {

   int k = (1 << i);
   if (f % (2*k) < k) {
      return f + k;
    }

    return f;
    
  }

  
  public static int unsetBit (int f, int i) {

   int k = (1 << i);
    if (f % (2*k) >= k) {
      return f - k;
    }

    return f;
    
  }


  public static int setBit (int f, int i, boolean b) {

    if (b) {
      return setBit (f, i);
    }

    return unsetBit (f, i);
    
  }

  /***********************************************************************************
   * 
   *                           getPrintWriter
   *
   ***********************************************************************************/

  public static PrintWriter getPrintWriter (String filename) throws IOException {

    if (filename == "") {
      return null;
    }

    if (filename.equals("-")) {
      return new PrintWriter (System.out);
    }
    
    File  outputFile = new File (filename);            
    if (outputFile.exists() && ! outputFile.canWrite()) {
      throw new IOException ("Warning: file " + filename + " cannot be written to.");
    }
    
    if (filename.endsWith(".gz")) {
      return new PrintWriter (new OutputStreamWriter (new GZIPOutputStream (new FileOutputStream (outputFile))));
    }

    return (new PrintWriter (outputFile));

    
  }


  /***********************************************************************************
   * 
   *                           getBufferedReader
   *
   ***********************************************************************************/

  public static BufferedReader getBufferedReader (String filename) throws IOException {

    if (filename == "") {
      return null;
    }

    if (filename.equals("-")) {
      System.err.println ("Processing input from std in.");
      return new BufferedReader (new InputStreamReader (System.in));
    }
    
    File  inputFile = new File (filename);      
    if (! inputFile.exists()) {
      throw new IOException ("Warning: file " + filename + " not found.");
    }
      
    if (! inputFile.canRead()) {
      throw new IOException ("Warning: file " + filename + " cannot be read.");
    }

    if (filename.endsWith(".gz")) {
      if (debugLevel >= 2) {
	System.err.println("Return reader for gzipped file.");
      }
      return new BufferedReader (new InputStreamReader (new GZIPInputStream (new FileInputStream (inputFile))));
    }

    return new BufferedReader (new FileReader (inputFile));

    
  }


  /***********************************************************************************
   * 
   *                           getBufferedReader
   *
   ***********************************************************************************/

  public static BufferedReader getBufferedReader (String filename, boolean check, String fileDescription) throws IOException {

    BufferedReader bufferedReader = getBufferedReader (filename);

    if (bufferedReader == null) {
      throw new IOException ("ERROR: " + fileDescription + ": \"" + filename + "\" not found.");
    }

    return bufferedReader;
    
  }


  /***********************************************************************************
   * 
   *                           getReader
   *
   ***********************************************************************************/

  public static Reader getReader (String filename) throws IOException {

    if (filename == "") {
      return null;
    }

    if (filename.equals("-")) {
      return new BufferedReader (new InputStreamReader (System.in));
    }
    
    File  inputFile = new File (filename);      
    if (! inputFile.exists()) {
      throw new IOException ("Warning: file " + filename + " not found.");
    }
      
    if (! inputFile.canRead()) {
      throw new IOException ("Warning: file " + filename + " cannot be read.");
    }

    if (filename.endsWith(".gz")) {
      return new InputStreamReader (new GZIPInputStream (new FileInputStream (inputFile)));
    }

    return new FileReader (inputFile);

    
  }

  /***********************************************************************************
   * 
   *                     reverseComplement
   *
   ***********************************************************************************/
  
  public static String reverseComplement (String sequence) {
    
    char [] sequenceArray = sequence.toUpperCase().toCharArray ();
    char [] reverseArray  = new char[sequence.length()];

    for (int i = 0; i < reverseArray.length; i++) {
      switch (sequenceArray[sequenceArray.length - 1 - i]) {
      case 'A':
	reverseArray[i] = 'T';
	break;
      case 'C':
	reverseArray[i] = 'G';
	break;
      case 'G':
	reverseArray[i] = 'C';
	break;
      case 'N':
	reverseArray[i] = 'N';
	break;
      case 'T':
	reverseArray[i] = 'A';
	break;
      }
      
    }

    return new String(reverseArray);

  }


  /***********************************************************************************
   *
   *                           checkFastqFilename 
   *
   ***********************************************************************************/
  
  public static boolean checkFastqFilename (String filename) throws Exception {

    if (filename.equals("-")) {
      return true;
    }

    try {
      if (filename.equals("")) {
	throw new Exception ("No reads file given.");
      }
      
      if (! filename.substring(filename.length() - 6).equals(".fastq") &&
	  ! filename.substring(filename.length() - 9).equals(".fastq.gz")&&
	  ! filename.substring(filename.length() - 3).equals(".fq") &&
	  ! filename.substring(filename.length() - 6).equals(".fq.gz")) {
	throw new Exception ("File is not a fastq file:" + filename);
      }

      return true;
    } catch (Exception e) {
      System.err.println(e==null?"null":e.getMessage());
    }

    return false;
    
  }


  /***********************************************************************************
   *
   *                           getFastqFilenameBase
   *
   ***********************************************************************************/
  
  public static String getFastqFilenameBase (String filename) throws Exception {

    if (filename.equals("-")) {
      return "-";
    }

    try {
      if (filename == null || filename.equals("")) {
	throw new Exception ("No reads file given.");
      }
      
      if (! filename.substring(filename.length() - 6).equals(".fastq") &&
	  ! filename.substring(filename.length() - 9).equals(".fastq.gz")&&
	  ! filename.substring(filename.length() - 3).equals(".fq") &&
	  ! filename.substring(filename.length() - 6).equals(".fq.gz")) {
	throw new Exception ("File is not a fastq file with ending (fastq|fq)[.gz]:" + filename);
      }

      if (filename.substring(filename.length() - 6).equals(".fastq") || filename.substring(filename.length() - 6).equals(".fq.gz")) {
	return filename.substring(0, filename.length() - 6);
      }
      
      if (filename.substring(filename.length() - 9).equals(".fastq.gz")) {
	return filename.substring(0, filename.length() - 9);
      }
      
      if(filename.substring(filename.length() - 3).equals(".fq")) {
	return filename.substring(0, filename.length() - 3);
      }

    } catch (Exception e) {
      System.err.println(e==null?"null":e.getMessage());
    }

    return null;
    
  }


  /***********************************************************************************
   *
   *                           addToTableEntry  
   *
   ***********************************************************************************/

  public static void addToTableEntry (Hashtable<String, Integer> entryTable, String key, int value) {

    Integer tableValueInt = entryTable.get(key);
    int tableValue = 0;
    if (tableValueInt != null) {
      tableValue = tableValueInt.intValue();
    }
	    
    tableValue = tableValue + value;
    entryTable.put(key, new Integer(tableValue));
  }


  /***********************************************************************************
   *
   *                           modifyFragmentId
   * modifyFragmentId is called in SamRecord/toString (suffix) (to write SamRecords),
   * SamProcessorFragmentEntry/outputFragmentEntry () (to write weight records), and
   * SamProcessorBed/init () (to set the fragmentId which is used in 
   * SamRecord/printBedRecords being called by SamProcessorBed/processSamRecords)
   *
   ***********************************************************************************/

  public static String modifyFragmentId (String fragmentId) throws IOException {

    if (! fragmentIdModificationIsSet ()) {
      return fragmentId;
    }

    if (commonPrefixLen < 0) {
      commonPrefixLen = commonPrefix.length();
    }

    if (idCutOffString != "" && ! fragmentId.startsWith (idCutOffString) && fragmentId.indexOf (idCutOffString) >= 0) {
      fragmentId = fragmentId.substring(fragmentId.indexOf (idCutOffString));
    }
    
    if (commonPrefixLen > 0 && fragmentId.startsWith (commonPrefix)) {
      fragmentId = fragmentId.substring(commonPrefixLen, fragmentId.length ());
    }
	  
    return fragmentId;
  }


  /***********************************************************************************
   *
   *                      fragmentIdModificationIsSet 
   *
   ***********************************************************************************/

  public static boolean fragmentIdModificationIsSet () {
    
    if (commonPrefix != "" || idCutOffString != "") {
      return true;
    }

    return false;

  }
  

  /***********************************************************************************
   *
   *                      getSimplifiedReferenceId
   *
   ***********************************************************************************/

  public static String getSimplifiedReferenceId (String referenceId) {

    int junctionIndex = referenceId.indexOf ("-junction");
    if (junctionIndex > 0) {
      int colonIndex = referenceId.indexOf (":");
      if (colonIndex > 0) {
	return referenceId.substring (0, colonIndex);
      }
    }
    return referenceId;
    
  }
  
  
  /***********************************************************************************
   *
   *                      getCountObjectCombinationId
   *
   ***********************************************************************************/

  public static String getCountObjectCombinationId (String countObjectId) {

    int exonIndex = countObjectId.indexOf ("-exon");
    if (exonIndex > 0) {
      return countObjectId.substring (0, exonIndex);
    }

    int junctionIndex = countObjectId.indexOf ("-junction");
    if (junctionIndex > 0) {
      return countObjectId.substring (0, junctionIndex);
    }

    return countObjectId;
    
  }


  /***********************************************************************************
   *
   *                      getChromosomeIds
   *
   ***********************************************************************************/

  private static HashSet<String> chromosomeIds = null;
  
  public static HashSet<String> getChromosomeIds () {
    return chromosomeIds;
  }

  
  /***********************************************************************************
   *
   *                      setChromosomeIds
   *
   ***********************************************************************************/

  public static void setChromosomeIds (HashSet<String> chromosomeIdSet) {

    chromosomeIds = chromosomeIdSet;
  }


}
