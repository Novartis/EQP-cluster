/**File: CreateUniqueList.java 

Original Author: Sven Schuierer
Date: 24/11/2014

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
 *                              Class CreateUniqueList
 *
 *   Reads an input file line by line and output the unique lines
 *
 ***********************************************************************************/

public class CreateUniqueList {

  private static final int debugLevel = UtilLib.getDebugLevel ();
  private static int countUnit = 25 * 1000 * 1000;

  /***********************************************************************************/

   private static void printHelp () {
    System.out.print("CreateUniqueList\n" +                                              
    "USAGE: CreateUniqueList -i <input file> -o <unique string file>\n");
  }
                                     
  /***********************************************************************************/

  public static void main (String [] args) {

    String inputFilename = "-";
    String outputFilename = "-";
    
    Getopt g = new Getopt("CreateUniqueList.java", args, "o:i:h");
    
    int c;
    String arg = "";
    
    c = g.getopt();
    
    while (c  != -1) {
      switch(c) {
      case 'o':
	outputFilename = g.getOptarg();
	break;
      case 'i':
	inputFilename = g.getOptarg();
	break;
      case 'h':
	printHelp();
	System.exit(0);
	break;
      default:
	System.out.print("Error: getopt() returned unknown option: " + c + "\n");
      }
      c = g.getopt();
    }

    String line = "";
    int lineNumber = 0;
    try {
      
      /***********************************************************************************
       *
       *  Read input reference ids
       *
       ***********************************************************************************/

      System.err.println("Reading input file " + (inputFilename.equals("-")?"stdin":inputFilename));
      System.err.println("(. = " + countUnit + " entries.)");
      System.err.flush();
      BufferedReader inputReader = UtilLib.getBufferedReader (inputFilename);

      System.err.println("Writing reference names to " + (outputFilename.equals("-")?"stdout":("file:\n" + outputFilename)));
      PrintWriter outputWriter = UtilLib.getPrintWriter (outputFilename);
      
      line = inputReader.readLine();
      lineNumber = 0;
      HashSet<String> stringSet = new HashSet<String> (2 * 1000 * 1000);
      while (line != null) {

	if (! stringSet.contains (line)) {
	  stringSet.add (line);
	  outputWriter.println (line);
	}
	
	lineNumber++;
	if (lineNumber % countUnit == 0) {
	  System.err.print(".");
	}
	  
	line = inputReader.readLine();
	
      }

      inputReader.close();
      outputWriter.close ();

      if (lineNumber > countUnit) {
	System.err.println();
      }
      System.err.flush();
     
    }
    catch (IOException e) {
      System.err.println ("IO ERROR: Problem with line: " + line + ", message: " + (e==null?"No error message":e.getMessage()));
    }
    catch (Exception e) {
      System.err.println ("ERROR: Problem with line: " + line + ", message: " + (e==null?"No error message":e.getMessage()));
    }
  }
}
