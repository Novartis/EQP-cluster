/**File: FastqEntry

Original Author: Sven Schuierer
Date: 13/01/2012

Classes :   FastqEntry

$Author: Sven Schuierer $

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

/**********************************************************************************
 *
 *  Class FastqEntry
 *
 **********************************************************************************/

public class FastqEntry {

  private static final int debugLevel = 0;

  private String fullReadId;
  private String readId;
  private String fragmentId;
  private String sequence;
  private String qualities;

  private boolean isFastaEntry = false;

  private int readIndex;
  private int readMultiplicity = 1;

  FastqEntry (BufferedReader reader, boolean extractReadMultiplicities) throws IOException {

    fullReadId = reader.readLine ();
    if (fullReadId == null) {
      throw new IOException ("Fastq input complete.");
    }

    if (! fullReadId.startsWith ("@")) {
      if (fullReadId.startsWith (">")) {
	isFastaEntry = true;
      } else {
	throw new IOException ("Read id not found in ." + reader.toString () + " instead: " + fullReadId);
      }
    }

    int spaceIndex = fullReadId.indexOf(" ");
    if (spaceIndex != -1) {
      readId = fullReadId.substring(0, spaceIndex);
      if (extractReadMultiplicities) {
	try {
	  readMultiplicity = Integer.parseInt(fullReadId.substring(spaceIndex + 1));
	} catch (NumberFormatException e) {
	  throw new IOException ("ERROR: Fragment " + fullReadId + " does not contain a number after the first space.");
	}
      }
    } else {
      if (extractReadMultiplicities) {
	throw new IOException ("ERROR: Fragment " + fullReadId + " does not contain a space separated field with read multiplicities.");
      }
      readId = fullReadId;
    }

    readIndex = 1;
    if (fullReadId.endsWith("/2")) {
      readIndex = 2;
    }

    fragmentId = readId;
    if (fragmentId.endsWith ("/1") || fragmentId.endsWith ("/2")) {
      fragmentId = fragmentId.substring(0, fragmentId.length() - 2);
    }
	
    if (debugLevel >= 2) {
      System.out.println ("Reading: " + reader.toString () + " - " + readId);
    }

    sequence = reader.readLine ().toUpperCase ();
    if (sequence == null) {
      throw new IOException ("Sequence not found for the fastq entry of " + fullReadId + " in " + reader.toString ());
    }

    if (! isFastaEntry) {
      String line = reader.readLine ();
      if (line == null || ! line.startsWith ("+")) {
	throw new IOException ("Second read id not found for the fastq entry of " + fullReadId + " in " + reader.toString ());
      }
      
      qualities = reader.readLine ();
      if (qualities == null) {
	throw new IOException ("Qualities not found for the fastq entry of " + fullReadId + " in " + reader.toString ());
      }
    }
  }


  FastqEntry (BufferedReader reader) throws IOException {
    this (reader, false);
  }


  public String getReadId () {
    return readId;
  }

  public String getFragmentId () {
    return fragmentId;
  }

  public String getFullReadId () {
    return fullReadId;
  }

  public String getSequence () {
    return sequence;
  }

  public String getSequenceWithoutChar (String c) {
    return sequence.replace(c, "");
  }

  public String getQualities () {
    return qualities;
  }

  public String toMultiplicityString (Counter fragmentCounter) throws IOException {

    fragmentId = "F" + fragmentCounter.getZeroFilledCount ();
    return fragmentId + "\t" + readMultiplicity;
    
  }
  
  public boolean equals (Object o) {

    if (o == null) {
      return false;
    }

    FastqEntry f = (FastqEntry) o;

    return f.getReadId().equals(getReadId ());
    
  }

  public int hashCode () {

    return getReadId().hashCode();
    
  }

  public String toString () {
    return
      getFullReadId () + "\n" +
      getSequence ()   + "\n" +
      "+"              + "\n" +
      getQualities ();
  }


  public String toString (int length) {
    return
      getFullReadId () + "\n" +
      getSequence ().substring(0, length) + "\n" +
      "+"              + "\n" +
      getQualities ().substring(0, length) ;
  }

  public String toString (int length, int inc) {

    String sequence  = getSequence ();
    String qualities = getQualities ();
    int    seqLength = sequence.length();

    String returnString = "";
    for (int start = 0; start + length < seqLength; start = start + inc) {

      if (start > 0) {
	returnString = returnString + "\n";
      }
      
      returnString = returnString + 
	getFragmentId () + "-" + start + "/" + readIndex + "\n" +
	sequence.substring (start, start + length)   + "\n" +
	"+"                                          + "\n" +
	qualities.substring (start, start + length);
      
    }

    return returnString;
    
  }


  public String toString (Counter fragmentCounter, int readIndex) throws IOException {

    fragmentId = "@F" + fragmentCounter.getZeroFilledCount ();
    readId     = fragmentId + "/" + readIndex;
    fragmentCounter.inc();

    return
      readId + "\n" + 
      getSequence ()   + "\n" +
      "+"              + "\n" +
      getQualities ();
    
  }

  public String toString (String fragmentId, int readIndex) throws IOException {

    this.fragmentId = fragmentId;
    readId = fragmentId + "/" + readIndex;
    return
      readId + "\n" + 
      getSequence ()   + "\n" +
      "+"              + "\n" +
      getQualities ();
    
  }
  
}
