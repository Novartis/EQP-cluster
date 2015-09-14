/**File: BedEntry.java 

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
 *                              Class BedEntry
 *
 *  BedEntry represents a BED6 entry and is used to convert genomic alignments
 *  to bed entries. Used in SamRecord to construct the BedEntries a SamRecord
 *  consists of.
 *
 ***********************************************************************************/

class BedEntry {

  private static int debugLevel = 0;

  private String chromosome;
  private int start;
  private int end;
  private String name;
  private int score;
  private String strand;
  private String additionalParameter = "";

  private int numCiagrFields = 0;
  private int ciagrFieldIndex = 0;

  public BedEntry (String chromosome, int start, int end, String name, int score, String strand) {
    
    this.chromosome = chromosome;
    this.start      = start;
    this.end        = end;
    this.name       = name;
    this.score      = score;
    this.strand     = strand;

  }
  

  /***********************************************************************************
   *
   *  score is the mapQuality when called by SamRecord
   *
   ***********************************************************************************/

  public BedEntry (String chromosome, int start, int end, String name, String score, String strand) {

    this (chromosome, start, end, name, score.equals("*")?-1:Integer.parseInt (score), strand);
 
  }


  public BedEntry (String chromosome, int start, int end, String name, int score, String strand, int readAlignedLength, int numInsertions,
		   int numDeletions, int ciagrFieldIndex, int numCiagrFields) throws IOException {

    this (chromosome, start, end, name, score, strand);

    if (readAlignedLength <= 0 || numInsertions < 0 || numDeletions < 0) {
      throw new IOException ("ERROR: Wrong parameters for name " + name + " mapping to " + chromosome + " - read length: " + readAlignedLength +
			     ", num insertions: " + numInsertions + ", numDeletions: " + numDeletions);
    }

    this.additionalParameter = "-L" + readAlignedLength + "-I" + numInsertions + "-D" + numDeletions + "-F" + ciagrFieldIndex + "-" + (numCiagrFields - 1);

  }



  public BedEntry (String bedEntryString) {
    
    try {
      
      StringTokenizer st = new StringTokenizer (bedEntryString, "\t");
      String token = "";
      int i = 0;
      while (st.hasMoreTokens()) {

	token = st.nextToken();
	if (i == 0) {
	  chromosome = token;
	  i++;
	} else if (i == 1) {
	  start = Integer.parseInt (token);
	  i++;
	} else if (i == 2) {
	  end = Integer.parseInt (token);
	  i++;
	} else if (i == 3) {
	  name = token;
	  i++;
	} else if (i == 4) {
	  score = token.equals("*")?-1:Integer.parseInt (token);
	  i++;
	} else if (i == 5) {
	  strand = token;
	  i++;
	}
      }

    } catch (Exception e) {
      System.out.println("ERROR: problem with bed entry: " + bedEntryString);
    }
  }

  public String getChromosome () {
    return (chromosome);
  }
  
  public int getStart () {
    return (start);
  }
  
  public int getEnd () {
    return (end);
  }
  
  public String getName () {
    return (name);
  }

  public void setName (String value) {
    name = value;
  }
  
  public int getScore () {
    return (score);
  }
  
  public String getStrand () {
    return (strand);
  }

  public boolean isStartEntry () {
    return ciagrFieldIndex == 0;
  }

  public boolean isEndEntry () {
    return ciagrFieldIndex == numCiagrFields;
  }

  public boolean sharesLeftExonBoundary () {
    return ciagrFieldIndex > 0;
  }

  public boolean sharesRightExonBoundary () {
    return ciagrFieldIndex < numCiagrFields;
  }

  
  /***********************************************************************************
   *
   *                              toString
   *


   ***********************************************************************************/
  
  public String toString (String alignmentId, String alignmentParam, String referenceVersionId) {

    if (debugLevel >= 2) {
      System.out.println("toString (String alignmentId, String alignmentParam, String referenceVersionId)");
      System.out.println(toString());
    }

    Pattern pattern = Pattern.compile("(/[SP][12])$");
    Matcher matcher = pattern.matcher(name);

    String fragmentId = matcher.replaceFirst(alignmentId + alignmentParam + "$1");

    return(referenceVersionId + "\t"  + start + "\t"  + end + "\t"  + fragmentId + "\t"  + score + "\t" + strand);
    
  }


  public String toString (String alignmentId, String referenceVersionId) {

    if (debugLevel >= 2) {
      System.out.println("toString (String alignmentId, String referenceVersionId)");
      System.out.println(toString());
    }
   
    Pattern pattern = Pattern.compile("(/[SP][12])$");
    Matcher matcher = pattern.matcher(name);

    String fragmentId = matcher.replaceAll(alignmentId + additionalParameter + "$1");

    return(referenceVersionId + "\t"  + start + "\t"  + end + "\t"  + fragmentId + "\t"  + score + "\t" + strand);
    
  }

  public String toStringSingleEnd (String alignmentId, String referenceVersionId) {

    if (debugLevel >= 2) {
      System.out.println("toString (String alignmentId, String referenceVersionId)");
      System.out.println(toString());
    }
   
    Pattern pattern = Pattern.compile("/[SP]([12])$");
    Matcher matcher = pattern.matcher(name);

    String fragmentId = matcher.replaceAll(alignmentId + additionalParameter + "/S$1");

    return(referenceVersionId + "\t"  + start + "\t"  + end + "\t"  + fragmentId + "\t"  + score + "\t" + strand);
    
  }
  

  public String toString () {

    if (debugLevel >= 2) {
      System.out.println("Simple toString");
    }

    return(chromosome + "\t"  + start + "\t"  + end + "\t"  + name + additionalParameter + "\t"  + score + "\t" + strand);
    
  }

  
}
