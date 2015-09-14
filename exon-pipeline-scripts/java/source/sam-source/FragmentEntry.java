/**File: FragmentEntry.java 

Original Author: Sven Schuierer
Date: 03/01/2012

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


/***********************************************************************************
 *
 *                              Class FragmentEntry
 *
 *  A fragment entry contains
 *  - the fragment name
 *  - the read index (1 for the first read, 2 for the second read, and 3 if both
 *       are aligned but not as a pair)
 *  - the multiplicity (1 if only one read is aligned, 2 if both are aligned, as
 *       a pair or separately)
 *  - whether an alignment mate exists (hasMate)
 *  - the sum of the edit distances of the different alignments of a fragment and
 *  - the number of different alignments
 *
 *  The number of alignments is calculated as:
 *
 *      numAlignmentsPe + max(numAlignmentsSr1, numAlignmentsSr2)
 *
 *  Note that we assume that a read has either paired-end or single-read alignments.
 *
 *  and the edit distance:
 *        sum edit distance / number of alignments
 *
 *  where sum edit distance is the sum of the edit distance of all alignments
 *
 *  The latter two are counted separately for paired-end alignments, alignments
 *  of the first, and alignments of the second read. Furthermore, a fragment entry
 *  can store two sets of edit distances and number of alignments values which
 *  allows to compare between the two sets. The class also contains methods
 *  to compare the multiplicity and edit distance of two fragment entries.
 *
 *  The class is used in two different ways:
 *
 *  1. To construct a fragment entry from a collection of SamRecords for one
 *     fragment. See constructor "FragmentEntry (SamRecord s)" and method add;
 *     used by SamProcessorFragmentEntry.
 *  2. To construct a fragment entry from an entry in a weight file. The entry
 *     needs to consist of for tab-separated fields:
 *     - fragment name, number of alignments, sum edit distances, type of alignment
 *     used by FilterUnalignedFastqEntries, CombineSamFiles, SamProcessorCount,
 *       CompareReadWeightFiles, CombineReadWeightFiles
 *
 ***********************************************************************************/

class FragmentEntry implements Comparable {

  private static boolean warningsOn = false;

  private static int countUnit = 5 * 1000 * 1000;

  private static FragmentEntry fragmentEntry = null;
  private static String fragmentId = null;
  
  private static boolean outputNumSplicedFragments = false;
  public  static void setOutputNumSplicedFragments () {
    outputNumSplicedFragments = true;
  }

  private static boolean outputNumSoftMaskedBases = false;
  public  static void setOutputNumSoftMaskedBases () {
    outputNumSoftMaskedBases = true;
  }

  
  public static void setWarningsOn (boolean value) {
    warningsOn = value;
  }

  /***********************************************************************************
   *
   *                          Object variables and methods
   *
   ***********************************************************************************/

  private int debugLevel = 0;

  private boolean firstValuesSet = false;
  private int sumEditDistancePe;
  private int sumEditDistanceSr1;
  private int sumEditDistanceSr2;
  private int numAlignmentsPe;
  private int numAlignmentsSr1;
  private int numAlignmentsSr2;
  private int readIndex = 1;
  private int multiplicity = 0;
  private boolean hasMate = false;
  
  private int numSplicedFragments = 0;
  private int numSoftMaskedBases  = 0;

  private boolean secondValuesSet = false;
  private int sumEditDistance2Pe;
  private int sumEditDistance2Sr1;
  private int sumEditDistance2Sr2;
  private int numAlignments2Pe;
  private int numAlignments2Sr1;
  private int numAlignments2Sr2;
  private boolean hasMate2 =false;

  private String fragmentName = null;

  public FragmentEntry () {

    int debugLevel = UtilLib.getDebugLevel ();
    
    sumEditDistancePe  = 0;
    sumEditDistanceSr1 = 0;
    sumEditDistanceSr2 = 0;
    
    numAlignmentsPe    = 0;
    numAlignmentsSr1   = 0;
    numAlignmentsSr2   = 0;
    
  }

  public FragmentEntry (int numAlignments, int sumEditDistance, boolean hasMate, int readIndex) {

    this ();
    setValues (numAlignments, sumEditDistance, hasMate, readIndex);
    
  }


  /***********************************************************************************
   * 
   *                          Constructor
   *
   ***********************************************************************************/

  public FragmentEntry (String line) throws IOException {

    this ();
    
    if (line == null) {
      throw new IOException ("Fragment entry called on null line.");
    }

    if (debugLevel >= 1) {
      System.out.println("Line: " + line);
    }

    parseFragmentEntryLine (line);
    
  }


  /***********************************************************************************
   * 
   *          Parse a line containing the information for a FragmentEntry
   *
   ***********************************************************************************/

  public void parseFragmentEntryLine (String line) throws IOException {
      
    int pos = 0;
    int end = line.indexOf ("\t", pos);
    if (end == -1) {
      end = line.length();
    }
    String readId = line.substring(pos, end);;	
    pos = end + 1;

    if (readId.endsWith ("/1/1")) {
      readId = readId.substring(0, readId.length () - 2);
    } else if (readId.endsWith ("/2/2")) {
      readId = readId.substring(0, readId.length () - 2);
    } else if (readId.endsWith ("/1/2") || readId.endsWith ("/2/1")) {
      throw new IOException ("Wrong readId for read weight line: " + line);
    }

    int readIndex = 1;
    if (readId.endsWith ("/2")) {
      readIndex = 2;
    }

    fragmentName = readId;
    if (fragmentName.endsWith ("/1") || fragmentName.endsWith ("/2")) {
      fragmentName = fragmentName.substring (0, fragmentName.length () - 2);
    }

    String token = "";
    int numAlignments = 0;
    int sumEditDistance = 0;
    try {

      if (pos < line.length ()) {
	end = line.indexOf ("\t", pos);
	if (end == -1) {
	  end = line.length();
	}
	token = line.substring(pos, end);
	pos = end + 1;
    
	numAlignments = UtilLib.toInt (token);	
      } else {
	throw new IOException ("No weight for read " + readId + " found.");
      }

      if (pos < line.length ()) {
	end = line.indexOf ("\t", pos);
	if (end == -1) {
	  end = line.length();
	}
	token = line.substring(pos, end);
	pos = end + 1;
    
	sumEditDistance = UtilLib.toInt (token);	
      } else {
	if (warningsOn) {
	  throw new IOException ("No edit distance for read " + readId + " found.");
	}
      }
    }
    catch (NumberFormatException e) {
      throw new IOException ("Wrong number format for " + token + " in line: " + line);
    }

    boolean hasMate = false;
    if (pos < line.length ()) {
      end = line.indexOf ("\t", pos);
      if (end == -1) {
	end = line.length();
      }
      token = line.substring(pos, end);
      pos = end + 1;
    
      hasMate = token.equals("paired-end");
      if (token.indexOf ("both") >= 0) {
	readIndex = 3;
      }
    } else {
      if (warningsOn) {
	throw new IOException ("No paired end tag for read " + readId + " found.");
      }
    }

    setValues (numAlignments, sumEditDistance, hasMate, readIndex);
    
  }

  /***********************************************************************************
   * 
   *  Parse a line containing the information for a FragmentEntry (StringTokenizer
   *  version)
   *
   ***********************************************************************************/
  
  public void parseFragmentEntryLineStringTokenizer (String line) throws IOException {

    StringTokenizer st = new StringTokenizer (line, "\t");
      
    String readId = "";
    if (st.hasMoreTokens()) {
      readId = st.nextToken ();	
    }

    if (readId.endsWith ("/1/1")) {
      readId = readId.substring(0, readId.length () - 2);
    } else if (readId.endsWith ("/2/2")) {
      readId = readId.substring(0, readId.length () - 2);
    } else if (readId.endsWith ("/1/2") || readId.endsWith ("/2/1")) {
      throw new IOException ("Wrong readId for read weight line: " + line);
    }

    int readIndex = 1;
    if (readId.endsWith ("/2")) {
      readIndex = 2;
    }

    fragmentName = readId;
    if (fragmentName.endsWith ("/1") || fragmentName.endsWith ("/2")) {
      fragmentName = fragmentName.substring (0, fragmentName.length () - 2);
    }

    String token = "";
    int numAlignments = 0;
    int sumEditDistance = 0;
    try {
   
      if (st.hasMoreTokens()) {
	token = st.nextToken ();
	numAlignments = UtilLib.toInt (token);	
      } else {
	throw new IOException ("No weight for read " + readId + " found.");
      }

      if (st.hasMoreTokens()) {
	token = st.nextToken ();
	sumEditDistance = UtilLib.toInt (token);	
      } else {
	if (warningsOn) {
	  throw new IOException ("No edit distance for read " + readId + " found.");
	}
      }
    }
    catch (NumberFormatException e) {
      throw new IOException ("Wrong number format for " + token + " in line: " + line);
    }

    boolean hasMate = false;
    if (st.hasMoreTokens()) {
      token = st.nextToken ();
      hasMate = token.equals("paired-end");
      if (token.indexOf ("both") >= 0) {
	readIndex = 3;
      }
    } else {
      if (warningsOn) {
	throw new IOException ("No paired end tag for read " + readId + " found.");
      }
    }

    setValues (numAlignments, sumEditDistance, hasMate, readIndex);

  }
  
  /***********************************************************************************
   * 
   *                          setValues
   *
   ***********************************************************************************/

  public void setValues (int numAlignments, int sumEditDistance, boolean hasMate, int readIndex) {

    this.debugLevel = UtilLib.getDebugLevel ();

    if (hasMate) {
      this.numAlignmentsPe  = numAlignments;
      this.sumEditDistancePe = sumEditDistance;
    } else {

      if (readIndex == 1) {
	this.numAlignmentsSr1   = numAlignments;
	this.sumEditDistanceSr1 = sumEditDistance;
      } else {
	this.numAlignmentsSr2   = numAlignments;
	this.sumEditDistanceSr2 = sumEditDistance;
      }
    }
    
    this.hasMate = hasMate;
    this.readIndex   = readIndex;
    firstValuesSet   = true;

  }

  /***********************************************************************************
   *
   *                          add
   *
   ***********************************************************************************/


  public void add (int numAlignments, int sumEditDistance, boolean hasMate, int readIndex) {

    if (hasMate) {
      this.sumEditDistancePe = this.sumEditDistancePe + sumEditDistance;
      this.numAlignmentsPe   = this.numAlignmentsPe + numAlignments;
    } else {
      
      if (readIndex == 1) {
	this.sumEditDistanceSr1= this.sumEditDistanceSr1 + sumEditDistance;
	this.numAlignmentsSr1  = this.numAlignmentsSr1   +   numAlignments;
      } else {
	this.sumEditDistanceSr2= this.sumEditDistanceSr2 + sumEditDistance;
	this.numAlignmentsSr2  = this.numAlignmentsSr2   +   numAlignments;
      }

    }
    
    this.hasMate    = hasMate;
    this.readIndex  = readIndex;
    firstValuesSet  = true;

  }

  public void add (int numAlignments, int sumEditDistance, boolean hasMate) {
    add(numAlignments, sumEditDistance, hasMate, readIndex);
  }
  

  public void add (FragmentEntry f) throws IOException {

    fragmentName = f.getFragmentName ();

    if (debugLevel >= 1) {
      System.out.println ("Adding fragment entry: " + f + " to " + this);
    }

    if (f.hasMate ()) {
      this.sumEditDistancePe = this.sumEditDistancePe + f.getSumEditDistancePe ();
      this.numAlignmentsPe   = this.numAlignmentsPe   + f.getNumAlignments();
    } else {
      if (f.getReadIndex () != readIndex) {
	readIndex = 3;
      }
      this.sumEditDistanceSr1 = this.sumEditDistanceSr1 + f.getSumEditDistanceSr1 ();
      this.sumEditDistanceSr2 = this.sumEditDistanceSr2 + f.getSumEditDistanceSr2 ();

      this.numAlignmentsSr1 = this.numAlignmentsSr1 +  f.getNumAlignmentsSr1 ();
      this.numAlignmentsSr2 = this.numAlignmentsSr2 +  f.getNumAlignmentsSr2 ();

    }

  }

  public void add (SamRecord s, boolean readsMode) throws IOException {

    String curFragmentName;
    if (! readsMode) {
      curFragmentName = s.getFragmentName ();
    } else {
      curFragmentName = s.getReadId ();
    }

    if (fragmentName == null) {
      fragmentName = curFragmentName;
    } else if (! fragmentName.equals(curFragmentName)) {
      throw new IOException ("Adding SAM record for fragment " + curFragmentName + " to FragmentEntry of fragment " + fragmentName);
    }
    
    firstValuesSet = true;

    if (readsMode) {
      hasMate = s.hasMate();
    } else {
      hasMate = false;
    }

    if (s.getReadIndex () != readIndex) {
      readIndex = 3;
    }

    if (s.getReadIndex () == 1) {
      numAlignmentsSr1++;
      sumEditDistanceSr1 = sumEditDistanceSr1 + s.getEditDistance ();
    } else {
      numAlignmentsSr2++;
      sumEditDistanceSr2 = sumEditDistanceSr2 + s.getEditDistance ();
    }

    if (s.isSpliced ()) {
      numSplicedFragments += 1;
    }

    if (outputNumSoftMaskedBases) {
      CiagrString ciagrString = new CiagrString (s, true);
      numSoftMaskedBases += ciagrString.getSoftClippingLength ();;
    }

    if (debugLevel >= 2) {
      System.err.println("Adding SAM record: " + s + "\nNum alignments sr1: " + numAlignmentsSr1 + ", num alignments sr2: " + numAlignmentsSr2);
    }
  }

  public void add (SamRecord s1, SamRecord s2) throws IOException {

    if (fragmentName == null) {
      fragmentName = s1.getFragmentName ();
    } else if (! fragmentName.equals(s1.getFragmentName())) {
      throw new IOException ("Adding SAM record for fragment " + s1.getFragmentName() + " to FragmentEntry of fragment " + fragmentName);
    }

    if (! fragmentName.equals(s2.getFragmentName())) {
      throw new IOException ("Adding SAM record for fragment " + s2.getFragmentName() + " to FragmentEntry of fragment " + fragmentName);
    }

    firstValuesSet  = true;
    sumEditDistancePe = sumEditDistancePe + s1.getEditDistance () + s2.getEditDistance ();
    numAlignmentsPe++;
    hasMate = true;

    if (s1.isSpliced () || s2.isSpliced ()) {
      numSplicedFragments += 1;
    }

    if (debugLevel >= 2) {
      System.err.println("Adding SAM records: " + s1 + "\n" + s2 + "\nNum alignments pe: " + numAlignmentsPe);
    }
  }


  public void backup () {

    secondValuesSet     = true;
    sumEditDistance2Pe  = sumEditDistancePe;
    sumEditDistance2Sr1 = sumEditDistanceSr1;
    sumEditDistance2Sr2 = sumEditDistanceSr2;
    numAlignments2Pe    = numAlignmentsPe;
    numAlignments2Sr1   = numAlignmentsSr1;
    numAlignments2Sr2   = numAlignmentsSr2;
    hasMate2        = hasMate;
    
  }

  public void reset () {

    firstValuesSet     = false;
    sumEditDistancePe  = 0;
    sumEditDistanceSr1 = 0;
    sumEditDistanceSr2 = 0;
    numAlignmentsPe    = 0;
    numAlignmentsSr1   = 0;
    numAlignmentsSr2   = 0;
    hasMate        = false;
    
  }

  /***********************************************************************************
   *
   *                          compare value sets
   *
   ***********************************************************************************/

  public int getCorrectValueSetIndex (double distanceSlack) throws IOException {

    if (firstValuesSet) {	    
      if (! secondValuesSet) {
	return 1;
      } else {
	if (hasMate1 ()) {
	  if (! hasMate2 ()) {
	    return 1;
	  } else if (getEditDistance1 () <= getEditDistance2 () + 2 * distanceSlack) {
	    return 1;
	  } else {
	    return 2;
	  }
	} else if (! hasMate2 () && getEditDistance1 () <= getEditDistance2 () + distanceSlack) {
	  return 1;
	} else {
	  return 2;
	}
      }
    } else {
      return 2;
    }
    
  }
  

  /***********************************************************************************
   *
   *                          get methods
   *
   ***********************************************************************************/

  public String getFragmentName () throws IOException {
    return UtilLib.modifyFragmentId(fragmentName);
  }
  
  public int getSumEditDistance () {
    return sumEditDistancePe + sumEditDistanceSr1 + sumEditDistanceSr2;
  }

  public int getSumEditDistancePe () {
    return sumEditDistancePe;
  }

  public int getSumEditDistanceSr1 () {
    return sumEditDistanceSr1;
  }

  public int getSumEditDistanceSr2 () {
    return sumEditDistanceSr2;
  }
  
  public double getEditDistance () throws IOException {
    return getSumEditDistance () * 1.0 / getNumAlignments ();
  }

  public int getNumAlignments () throws IOException {
    if (getNumAlignmentsPe () > 0 &&  getNumAlignmentsSr () > 0) {
      System.err.println ("Warning for " + getFragmentName () + ": the number of PE alignments: " + getNumAlignmentsPe () + " and the number of single-read alignments: " +
			  getNumAlignmentsSr () + " are both greater than 0.");
    }
    return getNumAlignmentsPe () + getNumAlignmentsSr ();
  }

  public int getNumAlignmentsPe () {
    return numAlignmentsPe;
  }

  public int getNumAlignmentsSr () {
    return Math.max (numAlignmentsSr1,  numAlignmentsSr2);
  }

  
  public int getNumAlignmentsSr1 () {
    return numAlignmentsSr1;
  }

  public int getNumAlignmentsSr2 () {
    return numAlignmentsSr2;
  }

  public void setNumAlignments (int value) {
    if (hasMate()) {
      numAlignmentsPe = value;
    } else if (readIndex == 1) {
      numAlignmentsSr1 = value;
    } else {
      numAlignmentsSr2 = value;
    }
  }

  public void capNumAlignments (int value) {
    
    if (value < 0) {
      return;
    }
    
    if (hasMate()) {
      numAlignmentsPe = Math.min(value, numAlignmentsPe);
    } else if (readIndex == 1) {
      numAlignmentsSr1 = Math.min(value, numAlignmentsSr1);
    } else {
      numAlignmentsSr2 = Math.min(value, numAlignmentsSr2);
    }
  }

  public int getReadIndex () {
    return readIndex;
  }


  public int getSumEditDistance1 () {
    return getSumEditDistance ();
  }
  
  public int getSumEditDistance2 () {    
    return sumEditDistance2Pe + sumEditDistance2Sr1 + sumEditDistance2Sr2;
  }


  public double getEditDistance1 () throws IOException {
    return getEditDistance ();
  }

  public double getEditDistance2 () {
    return getSumEditDistance2 () * 1.0 / getNumAlignments2 ();
  }

  
  public int getNumAlignments1 () throws IOException {
    return getNumAlignments ();
  }

  public int getNumAlignments2 () {
    return numAlignments2Pe + Math.max(numAlignments2Sr1, numAlignments2Sr2);
  }


  public int getMultiplicity () {
    return ((hasMate() || readIndex > 2)?2:1);
  }
  
  public int getMultiplicity1 () {
    return getMultiplicity ();
  }
  
  public int getMultiplicity2 () {
    return (hasMate2()?2:1);
  }

  
  /***********************************************************************************
   *
   *                          test methods
   *
   ***********************************************************************************/

  public boolean hasMate () {
    return hasMate;
  }

  public boolean hasMate1 () {
    return hasMate;
  }

  public boolean hasMate2 () {
    return hasMate2;
  }
  
  public boolean firstValuesSet () {
    return firstValuesSet;
  }
  
  public boolean secondValuesSet () {
    return secondValuesSet;
  }

  public boolean identicalValues () throws IOException {
    return getSumEditDistance1 () == getSumEditDistance2 () && getNumAlignments1 () == getNumAlignments2 () &&
      hasMate1 () == hasMate2 ();
  }

  public boolean identicalNumAlignments () throws IOException {
    return getNumAlignments1 () == getNumAlignments2 ();
  }
  
  
  /***********************************************************************************
   *
   *                          toString
   *
   ***********************************************************************************/

  public String toString () {

    String alignmentType = hasMate?"paired-end":"single-read";
    if (! hasMate && getMultiplicity () == 2) {
      alignmentType = alignmentType + " (both)";
    }
    try {
      String returnString = fragmentName + "\t" + getNumAlignments () + "\t" + getSumEditDistance () + "\t" + alignmentType;
      if (outputNumSplicedFragments) {
	returnString = returnString + "\t" + numSplicedFragments;
      }
      if (outputNumSoftMaskedBases) {
	returnString = returnString + "\t" + numSoftMaskedBases;
      }
      return returnString;
    }
    catch (Exception e) {
      System.err.println ("ERROR: " + ((e==null)?"Null message":e.getMessage()));
      System.exit(1);
    }

    return fragmentName + "\t" + 0 + "\t" + getSumEditDistance () + "\t" + alignmentType;
  }
  
  public String toPrintString () throws IOException {

    String alignmentType = hasMate?"paired-end":"single-read";
    if (! hasMate && getMultiplicity () == 2) {
      alignmentType = alignmentType + " (both)";
    }
    String returnString = getFragmentName() + "\t" + getNumAlignments () + "\t" + getSumEditDistance () + "\t" + alignmentType;
    if (outputNumSplicedFragments) {
      returnString = returnString + "\t" + numSplicedFragments;
    }
    if (outputNumSoftMaskedBases) {
      returnString = returnString + "\t" + numSoftMaskedBases;
    }
    return returnString;
  }

  public String toString1 () throws IOException {
    return toPrintString ();
  }

  public String toString2 () {
    return getNumAlignments2 () + "\t" + getSumEditDistance2 () + "\t" + (hasMate2?"paired-end":"single-read");
  }

  
  /***********************************************************************************
   *
   *                         compareMultiplicity
   *
   ***********************************************************************************/

  public int compareMultiplicity (FragmentEntry f) {
    
    if (hasMate()) {
      if (! f.hasMate ()) {
	return 1;
      }
    } else if (f.hasMate ()) {
      return -1;
    }

    if (getMultiplicity () > 1) {
      if (f.getMultiplicity () <= 1) {
	return 1;
      }
    } else if (f.getMultiplicity () > 1) {
      return -1;
    }

    return 0;
    
  }

  /***********************************************************************************
   *
   *                         compareEditDistance
   *
   ***********************************************************************************/
  
  public int compareEditDistance (FragmentEntry f, double distanceSlack) throws IOException {
    
    int returnValue = compareMultiplicity (f);

    if (returnValue != 0) {
      return returnValue;
    }

    if (hasMate ()) {
      distanceSlack = 2 * distanceSlack;
    }

    if (f.getEditDistance () <= getEditDistance () && f.getNumAlignments() > 1.1 * getNumAlignments ()) {
      return -1;
    }
    
    if (f.getEditDistance () < getEditDistance () - distanceSlack) {
      return -1;
    }
  
    return 1;
    
  }


  /***********************************************************************************
   *
   *                         compareTo
   *
   ***********************************************************************************/

  public int compareTo (Object o) {

    if (o == null) {
      return 1;
    }

    FragmentEntry f = (FragmentEntry) o;

    try {
      return getFragmentName().compareTo(f.getFragmentName ());
    } catch (IOException e) {
      System.err.println ("Comparison of " + this + " to " + f + " failed.");
      System.exit (1);
    }

    return 0;
    
  }
  

}
