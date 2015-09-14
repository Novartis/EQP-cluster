/**File: SamReader.java 

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
import java.util.regex.Matcher;
import java.util.regex.Pattern;


/***********************************************************************************
 *
 *                              Class SamReader
 *
 *   Reads a sam file line by line which is given by BufferedReader object that is
 *   passed to SamReader on construction. It constructs a SamRecord object for each
 *   line that does not start with a "@". SamRecords with the same fragment name are
 *   collected in the vector SamRecords. Once the fragment name changes, SamRecords 
 *   is sent to/processed by the processSamRecords method of the samProcessor object
 *   which is passed to the readSamFile method.
 *
 ***********************************************************************************/

public class SamReader {

  private Hashtable<String, Integer> processedReads = null;
  private BufferedReader reader = null;
  private boolean saveProcessedReads = false;
  private String filename = null;

  private boolean warningsOn    = false;
  private boolean printLines    = false;
  private boolean pairedEndOnly = false;
  private boolean useQueryName  = false;

  private int countUnit = 5 * 1000 * 1000;

  private int numReads       = 0;
  private int numMappedReads = 0;
    
  private int numFragments       = 0;
  private int numMappedFragments = 0;

  private Counter fragmentCounter = new Counter (9);

  SamReader (BufferedReader reader, boolean saveProcessedReads, String filename, boolean useQueryName) {
    
    this.reader             = reader;
    this.saveProcessedReads = saveProcessedReads;
    this.filename           = filename;
    this.useQueryName       = useQueryName;
      
  }


  SamReader (BufferedReader reader, boolean saveProcessedReads, String filename) {
    
    this.reader             = reader;
    this.saveProcessedReads = saveProcessedReads;
    this.filename           = filename;
      
  }


  SamReader (BufferedReader reader, String filename) {

    this (reader, false, filename);
        
  }
  

  SamReader (BufferedReader reader) {

    this (reader, false, null);
        
  }
  

  public void setWarningsOn (boolean value) {
    warningsOn = value;
  }

  public void setPrintLines (boolean value) {
    printLines = value;
  }

  public void setCountUnit (int value) {
    countUnit = value;
  }

  public void init () {
    if (saveProcessedReads) {
      processedReads = new Hashtable<String, Integer> (20 * 1000 * 1000);
    }
  }

  public int numProcessedReads () {
    return(processedReads.size());
  }

  public int getNumReads () {
    return numReads;
  }
  
  public int getNumMappedReads () {
    return numMappedReads;
  }

  public int getNumFragments () {
    return numFragments;
  }
  
  public int getNumMappedFragments () {
    return numMappedFragments;
  }
  

  /***********************************************************************************
   *
   *                             readSamFile
   *
   ***********************************************************************************/
  
  public void readSamFile (SamProcessor samProcessor) throws IOException, Exception {

    final int debugLevel = UtilLib.getDebugLevel ();

    String line = "";
    int lineNumber = 1;
    String oldQueryName = "";
    String oldMappedQueryName = "";
    String oldFragmentName = "";
    String oldMappedFragmentName = "";

    boolean mappedRecordFound = false;
    boolean counterChecked = false;

    Vector<SamRecord> samRecords = new Vector<SamRecord> (500);

    try {

      Pattern pattern = Pattern.compile("^F[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]$");

      if (debugLevel >= 1) {
	System.err.println ("Processing " + filename);
      }
    
      line = reader.readLine();
      lineNumber++;

      while (line != null) {

	if (debugLevel >= 3) {
	  System.err.println ("Reading line: " + line);
	  System.err.flush();
	}

	if (! line.startsWith("@") && ! line.equals("")) {
	
	  SamRecord samRecord = new SamRecord (line, lineNumber, filename);
	  if (! counterChecked) {
	    counterChecked = true;
	    Matcher match = pattern.matcher(samRecord.getFragmentName());
	    if (! match.find()) {
	      System.err.println("Changing fragment id format from original id format to: Fnnnnnnnnn.");
	      SamRecord.setCounter (fragmentCounter);
	    }
	  }

	  String queryName    = samRecord.getOriginalQueryName();
	  String fragmentName = samRecord.getOriginalFragmentName();

	  if (debugLevel >= 2) {
	    System.err.println ("fragmentName: " + fragmentName + ", oldFragmentName: " + oldFragmentName + " is mapped: " + samRecord.isMapped());
	  }

	  if (printLines) {
	    System.err.println (samRecord.toString());
	  }

	  /*  Consider only aligned reads */
	  if (samRecord.isMapped()) {

	    if (! mappedRecordFound) {
	      samProcessor.init (samRecord);
	      mappedRecordFound = true;
	    }
	    
	    if (debugLevel >= 2) {
	      System.err.println ("fragmentName: " + fragmentName + ", oldMappedFragmentName: " + oldMappedFragmentName);
	    }
	    
	    if ((useQueryName && ! queryName.equals(oldMappedQueryName) && oldMappedQueryName != "") ||
		(! fragmentName.equals(oldMappedFragmentName) && oldMappedFragmentName != "")) {
	      if (debugLevel >= 2) {
		System.err.println ("Processing SAM records: " + samRecords);
	      }
	      
	      samProcessor.processSamRecords (samRecords);
	      fragmentCounter.inc();
	      samProcessor.init (samRecord);
	      // samRecords.clear();
	      samRecords = new Vector<SamRecord> (500);
	    }

	    if (debugLevel >= 3) {
	      System.err.println ("Adding SAM record: " + samRecord);
	      System.err.flush();
	    }
	    samRecord.fixFragmentName ();
	    if (debugLevel >= 3) {
	      System.err.println ("Fragment name fixed");
	      System.err.flush();
	    }
	    
	    samRecords.add(samRecord);
	    
	    if (debugLevel >= 3) {
	      System.err.println ("SAM record added");
	      System.err.flush();
	    }
	    
	    if (! queryName.equals(oldMappedQueryName)) {
	      numMappedReads++;
	    }	      
	    if (! fragmentName.equals(oldMappedFragmentName)) {
	      numMappedFragments++;
	    }
	      
	    oldMappedQueryName    = queryName;
	    oldMappedFragmentName = fragmentName;
	      
	  } else if (! fragmentName.equals(oldFragmentName)) {
	    fragmentCounter.inc();
	  }

	  
	  if (! queryName.equals(oldQueryName)) {
	    numReads++;
	  }
	  if (! fragmentName.equals(oldFragmentName)) {
	    numFragments++;
	  }
	  
	  oldQueryName    = queryName;
	  oldFragmentName = fragmentName;
	  
	  if (lineNumber % countUnit == 0) {
	    System.err.print(".");
	  }
	  
	  lineNumber++;

	}
	
	line = reader.readLine();
	
      }

      /* process last set of sam records */
      if (samRecords.size () > 0) {
	if (debugLevel >= 1) {
	  System.err.println ("Processing last SAM Records set.");
	  System.err.flush();
	}
	samProcessor.processSamRecords (samRecords);
      }

      if (lineNumber >= countUnit) {
	System.err.println();
      }

      System.err.format (numMappedReads + " of " + numReads + " reads are mapped (%.2f%%).%n", numMappedReads * 100.0 / Math.max(1, numReads));
      System.err.format (numMappedFragments + " of " + numFragments + " fragments are mapped (%.2f%%).%n", numMappedFragments * 100.0 / Math.max(1, numFragments));
      System.err.flush();

      reader.close();

    }
    catch (IOException e) {
      throw new IOException ("IO ERROR: Problem in line: " + line + ", message: " + (e==null?"No error message":e.getMessage()) + "\n" + "SamRecords: " + samRecords);
    }
    catch (Exception e) {
      throw new Exception ("ERROR: Problem in line: " + line + ", message: " + (e==null?"No error message":e.getMessage()));
    }
    
  }

}
