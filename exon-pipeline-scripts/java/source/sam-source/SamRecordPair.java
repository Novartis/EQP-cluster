/**File SamRecordPair.java 

Original Author: Sven Schuierer
Date: 28/04/2014

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
 *                           Class SamRecordPair
 *
 * Needed for the computation of genomic alignments (computeGenomicSamRecords) in order
 * to output only one alignment per genomic position.
 *
 ***********************************************************************************/

class SamRecordPair implements Comparable {

  private int debugLevel = 0;

  private SamRecord firstRead = null;
  private SamRecord secondRead = null;

  private int startPosition = -1;
  
  SamRecordPair () {
    debugLevel = UtilLib.getDebugLevel ();
  }

  SamRecordPair (SamRecord s1) throws IOException {
    this ();
    addSamRecord (s1);
  }  

  SamRecordPair (SamRecord s1, SamRecord s2) throws IOException {
    this (s1);
    addSamRecord (s2);
  }
  

  /***********************************************************************************
   * 
   *                        adding SAM records
   *
   ***********************************************************************************/

  public void addSamRecord (SamRecord s) throws IOException {
    
    if (s == null) {
      return;
    }
    
    if (s.isFirstRead () || ! s.isPairedInSequencing()) {
      if (! hasFirstRead ()) {
	firstRead = s;
      } else {
	throw new IOException ("Two first reads assigned to one SAM pair: " + firstRead + " and " + s);
      }
    } else {
      if (! hasSecondRead ()) {
	secondRead = s;
      } else {
	throw new IOException ("Two second reads assigned to one SAM pair: " + secondRead + " and " + s);
      }
    }
      
    if (startPosition != -1) {
      startPosition = Math.min(startPosition, s.getPosition ());
    } else {
      startPosition = s.getPosition ();
    }

    if (firstRead != null && secondRead != null) {

      /* Adjust the mate positions: Sometimes the start position of the mates needs to be adjusted due to overlapping
	 exons on a transcript. */
      firstRead.setMatePosition (secondRead.getPosition ());
      secondRead.setMatePosition (firstRead.getPosition ());
      
      if (! firstRead.getMateReferenceName ().equals(secondRead.getReferenceName ()) || ! secondRead.getMateReferenceName ().equals(firstRead.getReferenceName ())) {
	throw new IOException ("Inconsistent SAM pair:\n" + "First read: " + firstRead + "\n" + "Second read: " + secondRead);
      }
    }
    
  }

  
  /***********************************************************************************
   * 
   *                           get methods
   *
   ***********************************************************************************/
  
  public SamRecord getFirstRead () {
    return firstRead;
  }

  public SamRecord getSecondRead () {
    return secondRead;
  }
  
  public void setFirstRead (SamRecord s) {
    firstRead = s;
  }

  public void setSecondRead (SamRecord s) {
    secondRead = s;
  }

  public int getStartPosition () {
    return startPosition;
  }
  

  /***********************************************************************************
   * 
   *                         property methods
   *
   ***********************************************************************************/
  
  public boolean hasFirstRead () {
    return firstRead != null;
  }
  
  public boolean hasSecondRead () {
    return secondRead != null;
  }

  public boolean isPairedEndAlignment () {
    return firstRead != null && secondRead != null;
  }

  
  /***********************************************************************************
   * 
   *                        overlapsWithAndIsCloseTo
   *
   ***********************************************************************************/

  public boolean overlapsWithAndIsCloseTo (SamRecordPair samRecordPair, double lengthFraction, double overlapFraction) throws IOException {

    if (samRecordPair.isPairedEndAlignment ()) {
      if (! hasFirstRead () || ! hasSecondRead ()) {
	return false;
      }

      double firstReadOverlapFraction  = 1.0 * firstRead.overlap (samRecordPair.getFirstRead ()) / firstRead.getLength ();
      double secondReadOverlapFraction = 1.0 * secondRead.overlap (samRecordPair.getSecondRead ()) / secondRead.getLength ();

      if (firstReadOverlapFraction >= overlapFraction && secondReadOverlapFraction >= overlapFraction) {
	return true;
      }

      if (firstReadOverlapFraction >= lengthFraction && secondReadOverlapFraction >= lengthFraction &&
	  firstRead.getAlignmentDistance(samRecordPair.getFirstRead ()) <= lengthFraction * samRecordPair.getFirstRead ().getLength () &&
	  secondRead.getAlignmentDistance(samRecordPair.getSecondRead ()) <= lengthFraction * samRecordPair.getSecondRead ().getLength ()) {
	return true;
      }
    } else if (samRecordPair.hasFirstRead ()) {
      if (! hasFirstRead ()) {
	return false;
      }

      if (! firstRead.isMapped ()) {
	return false;
      }

      double firstReadOverlapFraction  = 1.0 * firstRead.overlap (samRecordPair.getFirstRead ()) / firstRead.getLength ();
      
      if (firstReadOverlapFraction >= overlapFraction) {
	return true;
      }

      if (firstReadOverlapFraction >= lengthFraction && 
	  firstRead.getAlignmentDistance(samRecordPair.getFirstRead ()) <= lengthFraction * firstRead.getLength ()) {
	return true;
      }
    } else if (samRecordPair.hasSecondRead ()) {
      if (! hasSecondRead ()) {
	return false;
      }

      double secondReadOverlapFraction  = 1.0 * secondRead.overlap (samRecordPair.getSecondRead ()) / secondRead.getLength ();
      
      if (secondReadOverlapFraction >= overlapFraction) {
	return true;
      }
      
      if (secondReadOverlapFraction >= lengthFraction && 
	  secondRead.getAlignmentDistance(samRecordPair.getSecondRead ()) <= lengthFraction * secondRead.getLength ()) {
	return true;
      }
    }
    
    return false;
  }


  /***********************************************************************************
   * 
   *                          isUnmapped
   *
   ***********************************************************************************/

   public boolean isUnmapped () throws IOException {

     boolean returnValue = false;
     if (firstRead != null && firstRead.isUnmapped ()) {
       if (secondRead == null) {
	 returnValue = true;
       } else {
	 throw new IOException ("Unmapped read paired with a second read in: " + this);
       }
     } else if (secondRead != null && secondRead.isUnmapped ()) {
       if (firstRead == null) {
	 returnValue = true;
       } else {
	 throw new IOException ("Unmapped read paired with a second read in: " + this);
       }
     } else if (firstRead == null && secondRead == null) {
       throw new IOException ("SAM pair entries with two null reads.");
     }

     if (debugLevel >= 2) {
       System.err.println ("isUnmapped returnValue: " + returnValue);
     }

     return returnValue;
    
   }

  
  /***********************************************************************************
   * 
   *                          addOptionalField
   *
   ***********************************************************************************/

  public void addOptionalField (String fieldName, String fieldType, String fieldValue) {
    
    if (firstRead != null) {
      firstRead.addOptionalField (fieldName, fieldType, fieldValue);
    }
    
    if (secondRead != null) {
      secondRead.addOptionalField (fieldName, fieldType, fieldValue);
    }
    
  }

  /***********************************************************************************
   * 
   *                          compareTo
   *
   ***********************************************************************************/

   public int compareTo (Object o) {

    SamRecordPair s = (SamRecordPair) o;

    if (startPosition < s.getStartPosition ()) {
      return -1;
    } else if (startPosition > s.getStartPosition ()) {
      return 1;
    }

    return 0;
    
   }

  /***********************************************************************************
   * 
   *                           toString with mate information
   *
   ***********************************************************************************/

  public String toString (boolean  [] alignmentsFound, boolean  [] alignmentsPrinted, int level) throws IOException {

    String returnString = "";
    if (firstRead != null) {
      returnString = firstRead.toString (alignmentsFound[1], alignmentsPrinted[0], level);
      alignmentsPrinted[0] = true;
    }
    
    if (secondRead != null) {
      if (returnString != "") {
	returnString = returnString + "\n" + secondRead.toString (alignmentsFound[0], alignmentsPrinted[1], level);
      } else {
	returnString = secondRead.toString (alignmentsFound[0], alignmentsPrinted[1], level);
      }
      alignmentsPrinted[1] = true;
    }

    return returnString;

  }

  /***********************************************************************************
   * 
   *                           toString
   *
   ***********************************************************************************/

  public String toString () {

    String returnString = "";
    if (firstRead != null) {
      returnString = firstRead.toString ();
    }
    
    if (secondRead != null) {
      if (returnString != "") {
	returnString = returnString + "\n" + secondRead.toString ();
      } else {
	returnString = secondRead.toString ();
      }
    }

    return returnString;

  }

  
  /***********************************************************************************
   * 
   *                           toFirstRead
   *
   ***********************************************************************************/

  public String toFirstRead () {

    String returnString = null;
    if (firstRead != null) {
      returnString = firstRead.toString ();
    }

    return returnString;

  }

  /***********************************************************************************
   * 
   *                           toSecondRead
   *
   ***********************************************************************************/

  public String toSecondRead () {

    String returnString = null;
    if (secondRead != null) {
      returnString = secondRead.toString ();
    }

    return returnString;

  }

  
  /***********************************************************************************
   * 
   *                           hashCode
   *
   ***********************************************************************************/

  public int hashCode () {

    long hashCode1 = 0;
    if (firstRead != null) {
      hashCode1 = firstRead.hashCode ();
    }
    
    long hashCode2 = 0;
    if (secondRead != null) {
      hashCode2 = secondRead.hashCode ();
    }

    return (int) ((hashCode1 + hashCode2) % Integer.MAX_VALUE);

  }


  /***********************************************************************************
   * 
   *                           equals
   *
   ***********************************************************************************/

  public boolean equals (Object o) {

    if (o == null) {
      return false;
    }

    SamRecordPair s = (SamRecordPair) o;

    if (firstRead == null) {
      if (secondRead == null) {
	return s.getFirstRead() == null && s.getSecondRead () == null;
      } else {
	return s.getFirstRead() == null && secondRead.equals (s.getSecondRead ());
      }
    } else {
      if (secondRead == null) {
	return firstRead.equals (s.getFirstRead()) && s.getSecondRead () == null;
      } else {
	return firstRead.equals (s.getFirstRead()) && secondRead.equals (s.getSecondRead ());
      }
    }

  }

}
