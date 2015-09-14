/**File: Exon.java 

Original Author: Sven Schuierer
Date: 20/12/2011

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
 *                              Class Exon
 *
 *   The Exon class is a mere wrapper class for string of the type
 *   <chromosome>/<start>/<end>/<strand> where only <chromosome>, <start>,
 *   and <end> determine the exon for equality, toString, and hashCode.
 *
 *   The Object methods just include get-methods, compareTo, equals, toString,
 *   and hashCode.
 *
 *   The relationship exon - count object is maintained via a HashSetTable.
 *
 ***********************************************************************************/

class Exon implements Comparable {

  private static final int debugLevel = 0;

  private static HashSetTable<Exon, String> countObjectTable = null;
  private static TreeMap<String, String>    countObjectIds = null;

  public static HashSetTable<Exon, String> getCountObjectTable () {
    return countObjectTable;
  }

  public static TreeMap<String, String> getCountObjectIds () {
    return countObjectIds;
  }

  private static int countObjectMapColumn = 2;
  public static int getCountObjectMapColumn () {
    return countObjectMapColumn;
  }
  
  private static Pattern versionPattern = Pattern.compile("[.][0-9]+$");
  

  /***********************************************************************************
   * 
   *                           loadCountObjectFile
   *
   ***********************************************************************************/

  public static void loadCountObjectFile (BufferedReader reader) throws IOException {

    if (reader == null) {
      countObjectTable = new HashSetTable<Exon, String> ();
      countObjectIds   = new TreeMap<String, String> ();
      return;
    }
    
    countObjectTable = new HashSetTable<Exon, String> (700 * 1000);
    countObjectIds   = new TreeMap<String, String> ();
    
    HashMap<String, String> exonIds = new HashMap<String, String> (500 * 1000);

    int countUnit = 500 * 1000;

    String line = reader.readLine();
    int lineNumber  = 1;
    while (line != null) {
      int fieldNumber = 0;
      StringTokenizer st = new StringTokenizer (line, "\t");

      String exonId = "";
      if (st.hasMoreTokens()) {
	exonId = st.nextToken ();	
      } else {
	throw new IOException ("No exon id found.");
      }
      fieldNumber++;

      /* Try to save space by reusing identical strings */
      if (exonIds.containsKey (exonId)) {
	exonId = exonIds.get (exonId);
      } else {
	exonIds.put (exonId, exonId);
      }
	
      if (! exonId.equals("Exon Id")) {
	String countObjectId = "";
	if (st.hasMoreTokens()) {
	  countObjectId = st.nextToken ();
	} else {
	  throw new IOException ("No count object for exon " + exonId + " found on loading.");
	}
	fieldNumber++;

	/* Check if there is a third column; if so, use the id in the third column as the countObjectId */
	if (st.hasMoreTokens()) {
	  countObjectId = st.nextToken ();
	  fieldNumber++;
	}

	if (fieldNumber < countObjectMapColumn) {
	  throw new IOException ("Not all lines have " + countObjectMapColumn + " columns.");
	} else if (fieldNumber > countObjectMapColumn) {
	  countObjectMapColumn = fieldNumber;
	}

	/* Try to save space by reusing identical strings */
	if (countObjectIds.containsKey (countObjectId)) {
	  countObjectId = countObjectIds.get (countObjectId);
	} else {
	  countObjectIds.put (countObjectId, countObjectId);
	}
	
	countObjectTable.add (new Exon (exonId, "/"), countObjectId);
      }

      if (lineNumber % countUnit == 0) {
	System.err.print (".");
      }
      
      line = reader.readLine();
      lineNumber++;
      
    }

    if (lineNumber > countUnit) {
      System.err.println (".");
    }

    reader.close ();

  }

  /***********************************************************************************
   * 
   *                           computeExonSetLength
   *
   ***********************************************************************************/

  public static int computeExonSetLength (HashSet<Exon> exonSet) {

    TreeSet<Exon> sortedExonSet = new TreeSet<Exon> (exonSet);
    int genomicLength = 0;

    Exon oldExon = null;
    int start = 0;
    int end = -1;
    for (Exon exon: sortedExonSet) {

      if (debugLevel >= 2) {
	System.out.println (exon);
      }
      
      if (end < exon.getChrStart () || (oldExon != null && ! exon.getChromosome ().equals(oldExon.getChromosome ()))) {
	
	if (oldExon != null && ! UtilLib.isTranscriptExon(oldExon.getChromosome ())) {

	  if (debugLevel >= 2) {
	    System.out.println ("Adding " + (end - start + 1));
	  }
	  
	  genomicLength = genomicLength + end - start + 1;
	}
	start = exon.getChrStart ();

	if (oldExon != null && ! exon.getChromosome ().equals(oldExon.getChromosome ())) {
	  end = -1;
	}
      }

      end = Math.max (end, exon.getChrEnd ());
      oldExon = exon;
      
    }
    
    if (oldExon != null && ! UtilLib.isTranscriptExon(oldExon.getChromosome ())) {
      if (debugLevel >= 2) {
	System.out.println ("Adding " + (end - start + 1) + " - done");
      }
      genomicLength = genomicLength + end - start + 1;
    }
      
    return genomicLength;
  }

  
  /***********************************************************************************
   * 
   *                         Object variables and methods
   *
   ***********************************************************************************/

  String geneId          = "";
  String chromosome      = "";
  int    chromosomeStart = 0;
  int    chromosomeEnd   = 0;
  String strand          = "";

  /* Since the strand may changed in BedRecord we keep the original strand separate which is
     important since we also use <chr>/<start>/<end>/<strand> as an identifier. */
  String originalStrand  = "";

  boolean isTranscriptExon = false;
  boolean isGenomicExon    = false;

  public Exon () {
  }


  public Exon (String exonString, String separator) {
    
    try {
      
      StringTokenizer st = new StringTokenizer (exonString, separator);
      String token = "";
      String chromosomeStartToken = "";
      String chromosomeEndToken   = "";
      
      int i = 0;
      while (st.hasMoreTokens()) {

	token = st.nextToken();
	if (i == 0) {
	  chromosome = token;
	  i++;
	} else if (i == 1) {
	  chromosomeStartToken = token;
	  i++;
	} else if (i == 2) {
	  chromosomeEndToken = token;
	  i++;
	} else if (i == 3) {
	  strand = token;
	  i++;
	} else if (i == 4) {
	  geneId               = chromosome;
	  chromosome           = chromosomeStartToken;
	  chromosomeStartToken = chromosomeEndToken;
	  chromosomeEndToken   = strand;
	  strand               = token;
	  i++;
	} else {
	  throw new IOException ("Too many fields: " + exonString);
	}
      }

      if (i < 3) {
	throw new IOException ("Less than three fields are separated by " + separator + " in exon: " + exonString);
      }

      if (i == 3) {
	strand = "+";
      }
      
      chromosomeStart = Integer.parseInt(chromosomeStartToken);
      if (chromosomeStart == 0) {
	chromosomeStart = 1;
      }
      
      chromosomeEnd   = Integer.parseInt(chromosomeEndToken);

      isTranscriptExon = UtilLib.isTranscriptExon(chromosome);
      /* Replace versions for "chromosomes" that are transcript identifiers - only occurs for tr-exons. */
      if (isTranscriptExon) {
	Matcher matcher = versionPattern.matcher(chromosome);
	chromosome = matcher.replaceAll("");
      }
      originalStrand = strand;

    } catch (Exception e) {
      System.err.println("ERROR: problem with exon: " + exonString + " - " + e);
    }
  }


  public Exon (String exonString, String separator, String referenceSequenceId) {
    
    this (exonString, separator);
    isGenomicExon = chromosome.equals (referenceSequenceId);

  }



  /***********************************************************************************
   * 
   *                          Basic object methods
   *
   ***********************************************************************************/

  public String getGeneId () {
    return geneId;
  }

  public String getChromosome () {
    return chromosome;
  }


  public boolean isTranscriptExon () {
    return isTranscriptExon;
  }


  public boolean isGenomicExon () {
    return isGenomicExon;
  }

  public int getChrStart () {
    return chromosomeStart;
  }

  public int getChrEnd () {
    return chromosomeEnd;
  }

  public String getStrand () {
    return strand;
  }

  public String getOriginalStrand () {
    return originalStrand;
  }

  public void setStrand (String value) {
    strand = value;
  }

  public void invertStrand () {
    if (strand.equals("+")) {
      strand = "-";
    } else {
      strand = "+";      
    }
  }

  public int getLength () {
    return chromosomeEnd - chromosomeStart + 1;
  }

  public boolean genomicCoordinatesEquals(Exon e) {
    return chromosome.equals(e.getChromosome()) && chromosomeStart == e.getChrStart () && chromosomeEnd == e.getChrEnd ();
  }

  
  public boolean equals (Object o) {

    if (o == null) {
      return false;
    }

    Exon e = (Exon) o;

    return chromosome.equals(e.getChromosome()) && chromosomeStart == e.getChrStart () && chromosomeEnd == e.getChrEnd () && originalStrand.equals(e.getOriginalStrand());
 
  }


  /***********************************************************************************
   * 
   *                          compareTo
   *
   ***********************************************************************************/

   public int compareTo (Object o) {

    Exon e = (Exon) o;

    if (! getChromosome ().equals(e.getChromosome ())) {
      return getChromosome ().compareTo(e.getChromosome ());
    }


    if (getChrStart () != e.getChrStart ()) {
      if (getChrStart () > e.getChrStart ()) {
	return 1;
      } else {
	return -1;
      }
    }


    if (getChrEnd () != e.getChrEnd ()) {
      if (getChrEnd () > e.getChrEnd ()) {
	return 1;
      } else {
	return -1;
      }
    }

    return getStrand ().compareTo(e.getStrand ());

  }


  public String toBedString (String sep) {
    
    return chromosome + sep + (chromosomeStart - 1) + sep + chromosomeEnd;
    
  }



  public String toString (String sep) {
    
    return chromosome + sep + chromosomeStart + sep + chromosomeEnd + sep + originalStrand;
    
  }


  public String toString () {
    
    return toString("/");
    
  }


  public int hashCode () {

    return toString().hashCode ();
    
  }
  
}
