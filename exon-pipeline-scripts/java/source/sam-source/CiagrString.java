/**File CiagrString.java 

Original Author: Sven Schuierer
Date: 11/10/2013

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
 *                              Class CiagrString
 *
 *
 ***********************************************************************************/

class CiagrString {

  private static int debugLevel = UtilLib.getDebugLevel ();

  /***********************************************************************************
   *
   *                         Object variables
   *
   ***********************************************************************************/

  private Vector<Integer> ciagrOpLengthVector = new Vector<Integer> (10);
  private Vector<String>  ciagrOperations     = new Vector<String> (10);

  private int [] ciagrOpLengthArray   = null;
  private int [] ciagrRefLengthArray  = null;
  private int [] ciagrReadLengthArray = null;

  private int numInsertions = 0;
  private int numDeletions  = 0;
  private int softClippingLength = 0;
  private int hardClippingLength = 0;
  private int numIntrons = 0;
  private Vector<int []> relativeIntronCoordinates = null;
  
  private int matchingLength = 0;

  private int numCiagrFields = 0;
  private int readAlignedLength = -1;
  private int refAlignedLength = -1;

  private String ciagrString = null;
  private String sequence = null;
  private String qualityString = null;

  private boolean storeIntronCoordinates = true;
  
  
  /***********************************************************************************
   *
   *                         Constructors
   *
   ***********************************************************************************/

  public CiagrString () {
  }

  
  public CiagrString (String ciagrString) throws IOException {

    debugLevel = UtilLib.getDebugLevel ();
    this.ciagrString = ciagrString;
    analyseCiagrString (ciagrString);

  }


  public CiagrString (String ciagrString, String sequence, String qualityString) throws IOException {

    this (ciagrString);
    this.sequence      = sequence;
    this.qualityString = qualityString;

    if (debugLevel >= 2 || false) {
      System.out.println ("Sequence: " + sequence + " and quality string: " + qualityString + " added.");
    }
    
  }

  public CiagrString (SamRecord samRecord, boolean storeIntronCoordinates) throws IOException {

    this.ciagrString     = samRecord.getCiagrString();
    this.sequence        = samRecord.getSequence ();
    this.qualityString   = samRecord.getQualityString ();
    this.storeIntronCoordinates = storeIntronCoordinates;

    analyseCiagrString (ciagrString);
    
    if (debugLevel >= 2 || false) {
      System.out.println ("Sequence: " + sequence + " and quality string: " + qualityString + " added.");
    }
    
  }

  
  /***********************************************************************************
   *
   *          addField              
   *
   ***********************************************************************************/

  private void analyseCiagrString (String ciagrString) throws IOException {
    
    Pattern pattern = Pattern.compile("[0-9]+[MIDSHNPX=N]");
    Matcher match   = pattern.matcher(ciagrString);
    int sumIntronLength = 0;
    int relativeGenomePosition = 0;
    if (storeIntronCoordinates) {
      relativeIntronCoordinates = new Vector<int []> (3);
    }
      
    while (match.find()) {

      int fieldLength = Integer.parseInt (ciagrString.substring(match.start(), match.end()-1));
      String fieldOperation = ciagrString.substring(match.end()-1, match.end());
      addField(fieldLength, fieldOperation);
 
      if (fieldOperation.equals("M") || fieldOperation.equals("X") || fieldOperation.equals("=")) {
	matchingLength = matchingLength + fieldLength;
	relativeGenomePosition = relativeGenomePosition + fieldLength;
      } else if (fieldOperation.equals("I")) {
	numInsertions  = numInsertions + fieldLength;
      } else if (fieldOperation.equals("D")) {
	numDeletions = numDeletions + fieldLength;
	relativeGenomePosition = relativeGenomePosition + fieldLength;
      } else if (fieldOperation.equals("S")) {
	softClippingLength = softClippingLength + fieldLength;
      } else if (fieldOperation.equals("H")) {
	hardClippingLength = hardClippingLength + fieldLength;
      } else if (fieldOperation.equals("N")) {
	if (storeIntronCoordinates) {
	  sumIntronLength = sumIntronLength + fieldLength;
	  numIntrons += 1;
	  int [] intronInt = {relativeGenomePosition + 1, relativeGenomePosition + fieldLength};
	  relativeIntronCoordinates.add(intronInt);
	} else {
	  throw new IOException ("ERROR: N discovered in ciagrString: " + ciagrString);
	}
	relativeGenomePosition = relativeGenomePosition + fieldLength;
      } else if (! fieldOperation.equals("P")) {
	throw new IOException ("ERROR: Unknown operation in CIAGR string: " + ciagrString);
      }
    }

    refAlignedLength  = matchingLength + numDeletions;
    readAlignedLength = matchingLength + numInsertions + softClippingLength + hardClippingLength;
    relativeGenomePosition = matchingLength + numDeletions + sumIntronLength;

    numCiagrFields = ciagrOperations.size ();

    convertCiagrOpLengthToArray ();
  }


  /***********************************************************************************
   *
   *          addField              
   *
   ***********************************************************************************/

  public void addField (int length, String operation) {

    if (debugLevel >= 2) {
      System.out.println ("Adding CIAGR field " + length + operation + " to " + this);
      System.out.flush();
    }

    if (operation == null || operation.equals("") || length <= 0) {
      return;
    }
    
    if (numCiagrFields > 0) {
      String lastCiagrOperation = ciagrOperations.get(numCiagrFields - 1);
      if (lastCiagrOperation.equals("D") && operation.equals("N")) {
	ciagrOperations.set(numCiagrFields - 1, "N");
	lastCiagrOperation = "N";
      }
      if (operation.equals(lastCiagrOperation) || (lastCiagrOperation.equals("N") && operation.equals("D"))) {
	int opLength = getOpLength (numCiagrFields - 1);
	ciagrOpLengthVector.set(numCiagrFields - 1, new Integer (opLength + length));
	return;
      }
    }

    ciagrOpLengthVector.add(new Integer (length));
    ciagrOperations.add(operation);
    numCiagrFields++;
    
  }

  
  /***********************************************************************************
   *
   *          addSequence         
   *
   ***********************************************************************************/

  public void addSequence (String sequencePart, String qualityStringPart) {
    
    if (sequence != null) {
      sequence = sequence + sequencePart;
    } else {
      sequence = sequencePart;
    }
    
    if (qualityString != null) {
      qualityString = qualityString + qualityStringPart;
    } else {
      qualityString = qualityStringPart;
    }

  }


  /***********************************************************************************
   *
   *          transferInterval
   *
   * This function takes an interval [leftPosition, rightPosition] of relative coordinates
   * of the alignment of the read against the reference sequence/genome as input and
   * transfers the CIAGR entries, sequence and qualityString parts that fall into this
   * interval to ciagString.
   *
   ***********************************************************************************/

  public void transferInterval (int leftPosition, int rightPosition, CiagrString ciagrString, boolean isTranscriptExon, boolean isAtAlignmentMargin) throws IOException {

    if (leftPosition > rightPosition) {
      throw new IOException ("Right position " + rightPosition + " is smaller than leftPosition " + leftPosition + "  in transferInterval.");
    }
    
    if (ciagrString == null) {
      throw new IOException ("New ciagr string is null in transferInterval.");
    }
    
    if (ciagrOpLengthArray == null) {
      convertCiagrOpLengthToArray ();
    }

    if (debugLevel >= 2 || false) {
      System.out.println("Left position: " + leftPosition + ", right position: " + rightPosition);
    }	
    

    int leftIndex = 0;
    int leftSumOpRefLengths = 0;
    int leftSumOpReadLengths = 0;
    /* Find the first index i in ciagrRefLengthArray such that sum(ciagrRefLengthArray[0:i]) >= leftPosition */
     while (leftIndex < ciagrRefLengthArray.length - 1 && leftSumOpRefLengths + ciagrRefLengthArray[leftIndex] < leftPosition) {
      leftSumOpRefLengths  = leftSumOpRefLengths  + ciagrRefLengthArray[leftIndex];
      leftSumOpReadLengths = leftSumOpReadLengths + ciagrReadLengthArray[leftIndex];
      leftIndex = leftIndex + 1;
    }

    if (leftSumOpRefLengths + ciagrRefLengthArray[leftIndex] < leftPosition) {
      throw new IOException ("Left position " + leftPosition + " is not contained in CIAGR string " + this);
    }

    if (debugLevel >= 2 || false) {
      System.out.println("Left index: " + leftIndex + " for " + leftPosition + " in " + this + " computed.");
    }

    int rightIndex = leftIndex;
    int rightSumOpRefLengths  = leftSumOpRefLengths;
    int rightSumOpReadLengths = leftSumOpReadLengths;
    /* Find the last index i in ciagrRefLengthArray such that sum(ciagrRefLengthArray[0:i]) >= rightPosition and for the last index i'
       with sum(ciagrRefLengthArray[0:i']) < sum(ciagrRefLengthArray[0:i]) we also have sum(ciagrRefLengthArray[0:i']) < rightPosition.
       We do this in order to capture any trailing S or I fields. Note that this will not lead to duplicate assignment of S or I fields
       as long as  new leftPosition == old righPosition + 1. */
    while (rightIndex < ciagrRefLengthArray.length - 1 && rightSumOpRefLengths + ciagrRefLengthArray[rightIndex] <= rightPosition) {
      rightSumOpRefLengths  = rightSumOpRefLengths  + ciagrRefLengthArray[rightIndex];
      rightSumOpReadLengths = rightSumOpReadLengths + ciagrReadLengthArray[rightIndex];
      rightIndex = rightIndex + 1;
    }

    if (rightSumOpRefLengths + ciagrRefLengthArray[rightIndex] < rightPosition) {
      throw new IOException ("Right position " + rightPosition + " is not contained in CIAGR string " + this);
    }

    /* Go back one for rightIndex if we overshot */
    if (rightSumOpRefLengths + ciagrRefLengthArray[rightIndex] > rightPosition && rightSumOpRefLengths == rightPosition) {
      rightSumOpRefLengths  = rightSumOpRefLengths  - ciagrRefLengthArray[rightIndex - 1];
      rightSumOpReadLengths = rightSumOpReadLengths - ciagrReadLengthArray[rightIndex - 1];
      rightIndex = rightIndex - 1;
    }

    if (debugLevel >= 2 || false) {
      System.out.println("Right index: " + rightIndex + " for " + rightPosition + " computed (" + rightSumOpRefLengths + ", " +
			 ciagrRefLengthArray[rightIndex] + ")");
    }

    /* Compute the read sequence position that leftPosition corresponds to. Recall that leftPosition is a position on the reference sequence.
       Note that if ciagrReadLengthArray[leftIndex] > 0, then ciagrReadLengthArray[leftIndex] == ciagrOpLengthArray[leftIndex] and
       by the condition that ended the while loop for the leftIndex we also have ciagrReadLengthArray[leftIndex] == ciagrRefLengthArray[leftIndex] */
    int leftSequencePos = leftSumOpReadLengths + 1;
    if (ciagrReadLengthArray[leftIndex] > 0) {
      leftSequencePos = leftSequencePos - 1 + leftPosition - leftSumOpRefLengths;
    }
    
    /* The right sequence position is computed iteratively starting at the left sequence position */
    int rightSequencePos = leftSequencePos - 1;

    
    String defaultCiagrOperation = null;
    if (isTranscriptExon) {
      defaultCiagrOperation = "I";
      if (isAtAlignmentMargin) {
	defaultCiagrOperation = "S";
      }
    }
    
    if (leftIndex == rightIndex) {
      if (ciagrOpLengthArray[leftIndex] < rightPosition - leftPosition + 1) {
	throw new IOException ("Interval [" + leftPosition + "," + rightPosition + " is not contained in CIAGR operation " + leftIndex + ": " +
			       ciagrOpLengthArray[leftIndex] + getOperation(leftIndex, defaultCiagrOperation));
      }

      if (debugLevel >= 2 || false) {
	System.out.println("Adding field: " + rightPosition + " - " + leftPosition + " + 1 " + getOperation(leftIndex, defaultCiagrOperation));
      }	
      ciagrString.addField (rightPosition - leftPosition + 1, getOperation(leftIndex, defaultCiagrOperation));
      
      if (ciagrReadLengthArray[leftIndex] > 0) {
	rightSequencePos = rightSequencePos + rightPosition - leftPosition + 1;
      }
      
    } else {
      
      ciagrString.addField (leftSumOpRefLengths + ciagrOpLengthArray[leftIndex] - leftPosition + 1, getOperation(leftIndex, defaultCiagrOperation));
      if (ciagrReadLengthArray[leftIndex] > 0) {
	rightSequencePos = rightSequencePos + leftSumOpRefLengths + ciagrOpLengthArray[leftIndex] - leftPosition + 1;
      }

      if (debugLevel >= 2) {
	System.out.println ("1. rightSequencePos: " + rightSequencePos);
      }

      for (int i = leftIndex + 1; i < rightIndex; i++) {
	ciagrString.addField (ciagrOpLengthArray[i], getOperation(i, defaultCiagrOperation));
	rightSequencePos = rightSequencePos + ciagrReadLengthArray[i];

	if (debugLevel >= 2) {
	  System.out.println ("2. rightSequencePos: " + rightSequencePos + ", ciagrReadLengthArray["+i+"]: " + ciagrReadLengthArray[i]);
	}
      }

      if (ciagrRefLengthArray[rightIndex] > 0) {
	ciagrString.addField (rightPosition - rightSumOpRefLengths, getOperation(rightIndex, defaultCiagrOperation));
	if (ciagrReadLengthArray[rightIndex] > 0) {
	  /* In this case ciagrReadLengthArray[rightIndex] == ciagrRefLengthArray[rightIndex] == ciagrOpLengthArray[rightIndex] */
	  rightSequencePos = rightSequencePos + rightPosition - rightSumOpRefLengths;
	  if (debugLevel >= 2) {
	    System.out.println ("3. rightSequencePos: " + rightSequencePos);
	  }

	}
      } else {
	/* In this case ciagrReadLengthArray[rightIndex] == ciagrOpLengthArray[rightIndex] > 0 */
	ciagrString.addField (ciagrOpLengthArray[rightIndex], getOperation(rightIndex, defaultCiagrOperation));
	rightSequencePos = rightSequencePos + ciagrOpLengthArray[rightIndex];
	if (debugLevel >= 2) {
	  System.out.println ("4. rightSequencePos: " + rightSequencePos);
	}
      }

    }

    
    if (debugLevel >= 2 || false) {
      System.out.println("Adding substring: " + leftSequencePos + ", " + rightSequencePos + ": " + sequence.substring (leftSequencePos - 1, rightSequencePos));
    }
  
    ciagrString.addSequence (sequence.substring (leftSequencePos - 1, rightSequencePos), qualityString.substring (leftSequencePos - 1, rightSequencePos));

    
    if (debugLevel >= 2 || false) {
      System.out.println("Done.");
    }

  }
  

  /***********************************************************************************
   *
   *          convertCiagrOpLengthToArray
   *
   ***********************************************************************************/

  private void convertCiagrOpLengthToArray () {

    if (ciagrOpLengthArray != null) {
      return;
    }

    if (debugLevel >= 2 || false) {
      System.out.println ("Converting vectors of length " + ciagrOpLengthVector.size() + " to arrays.");
    }
    
    ciagrOpLengthArray = new int [ciagrOpLengthVector.size()];

    for (int i = 0; i < ciagrOpLengthArray.length; i++) {
      ciagrOpLengthArray[i] = ciagrOpLengthVector.get(i).intValue();
    }

    ciagrRefLengthArray  = new int [ciagrOpLengthVector.size()];
    ciagrReadLengthArray = new int [ciagrOpLengthVector.size()];

    refAlignedLength  = 0;
    readAlignedLength = 0;
    for (int i = 0; i < ciagrRefLengthArray.length; i++) {
      char ciagrOperation = ciagrOperations.get(i).charAt(0);
      if ("MX=".indexOf(ciagrOperation) >= 0){
	ciagrRefLengthArray [i] = ciagrOpLengthVector.get(i).intValue();
	ciagrReadLengthArray[i] = ciagrOpLengthVector.get(i).intValue();	
      } else if (ciagrOperation == 'D' || ciagrOperation == 'N') {
	ciagrRefLengthArray [i] = ciagrOpLengthVector.get(i).intValue();
	ciagrReadLengthArray[i] = 0;
      } else if (ciagrOperation == 'I' || ciagrOperation == 'S') {
	ciagrRefLengthArray [i] = 0;
	ciagrReadLengthArray[i] = ciagrOpLengthVector.get(i).intValue();
      } else {
	ciagrRefLengthArray [i] = 0;
	ciagrReadLengthArray[i] = 0;
      }

      refAlignedLength  = refAlignedLength  + ciagrRefLengthArray[i];	
      readAlignedLength = readAlignedLength + ciagrReadLengthArray[i];

    }
    
  }

  /***********************************************************************************
   *
   *           invertSequence         
   *
   ***********************************************************************************/

  public void invertSequence () {

    if (sequence != null) {
      sequence = UtilLib.reverseComplement (sequence);
    }

    if (qualityString != null) {
      qualityString = (new StringBuffer (qualityString)).reverse ().toString ();
    }
    
  }
    
    
  /***********************************************************************************
   *
   *           invert              
   *
   ***********************************************************************************/

  public void invert () {

    invertSequence ();

    if (numCiagrFields == 1) {
      return;
    }

    /* Invert ciagrOpLengthVector */
    Vector<Integer> ciagrOpLengthVectorTemp = ciagrOpLengthVector;
    ciagrOpLengthVector = new Vector<Integer> (numCiagrFields);
    for (int i = 0; i < numCiagrFields; i++) {
      ciagrOpLengthVector.add(ciagrOpLengthVectorTemp.get(numCiagrFields - 1 - i));
    }

    /* Invert ciagrOperations */
    Vector<String>  ciagrOperationsTemp = ciagrOperations;
    ciagrOperations = new Vector<String> (numCiagrFields);
    for (int i = 0; i < numCiagrFields; i++) {
      ciagrOperations.add(ciagrOperationsTemp.get(numCiagrFields - 1 - i));
    }
    
    /* Invert ciagrOpLengthArray */
    if (ciagrOpLengthArray != null) {
      int [] ciagrOpLengthArrayTemp = ciagrOpLengthArray;
      ciagrOpLengthArray = new int [numCiagrFields];  
      for (int i = 0; i < numCiagrFields; i++) {
      ciagrOpLengthArray[i] = ciagrOpLengthArrayTemp[numCiagrFields - 1 - i];
      }
    }
   
    /* Invert ciagRefLengthArray */
    if (ciagrRefLengthArray != null) {
      int [] ciagrRefLengthArrayTemp = ciagrRefLengthArray;
      ciagrRefLengthArray = new int [numCiagrFields];
      for (int i = 0; i < numCiagrFields; i++) {
	ciagrRefLengthArray[i] = ciagrRefLengthArrayTemp[numCiagrFields - 1 - i];
      }
    }

    /* Invert ciagReadLengthArray */
    if (ciagrReadLengthArray != null) {
      int [] ciagrReadLengthArrayTemp = ciagrReadLengthArray;
      ciagrReadLengthArray = new int [numCiagrFields];
      for (int i = 0; i < numCiagrFields; i++) {
	ciagrReadLengthArray[i] = ciagrReadLengthArrayTemp[numCiagrFields - 1 - i];
      }
    }    
  }

  
  
  /***********************************************************************************
   *
   *           getMethods              
   *
   ***********************************************************************************/

  public int getOpLength (int i) {
    if (i < numCiagrFields) {
      return ciagrOpLengthVector.get(i).intValue();
    } else {
      return -1;
    }
  }

  public String getOperation (int i, String s) {

    String op = null;
    if (i < numCiagrFields) {
      op = ciagrOperations.get(i);
    }

    /* Operations N and D can be replaced by an "S" operation as this will lead to inconsistencies with the sequence length */
    if (op == null || (s != null && ! op.equals("D") && ! op.equals("N"))) {
      return s;
    }

    return op;
    
  }
  
  public String getOperation (int i) {
    if (i < numCiagrFields) {
      return ciagrOperations.get(i);
    } else {
      return null;
    }
  }

  
  public int getLength () {
    return numCiagrFields;
  }


  public String getSequence () {
    return sequence;
  }
  
  
  public String getQualityString () {
    return qualityString;
  }
  

  public int getReadAlignedLength () {

    if (readAlignedLength < 0) {
      convertCiagrOpLengthToArray ();
    }

    return readAlignedLength;
  }

  
  public int getNumInsertions () {
    return numInsertions;
  }

  public int getNumDeletions () {
    return numDeletions;
  }
  
  public int getNumIntrons () {
    return numIntrons;
  }

  public int getSoftClippingLength () {
    return softClippingLength;
  }

  public int getHardClippingLength () {
    return hardClippingLength;
  }
  
  public Vector<int []> getIntronCoordinates () {
    if (relativeIntronCoordinates == null) {
      return new Vector<int []> ();
    }
    return relativeIntronCoordinates;
  }

  /***********************************************************************************
   *
   *           toString              
   *
   ***********************************************************************************/
  
  public String toString () {
    
    if (ciagrString != null) {
      return ciagrString;
    }

    String curCiagrString = "";
    for (int i = 0; i < numCiagrFields; i++) {
      curCiagrString = curCiagrString + ciagrOpLengthVector.get(i) + getOperation (i);
    }

    return curCiagrString;
  }

}
