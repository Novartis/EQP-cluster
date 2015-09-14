/**File: ShortenFastqEntries.java

Original Author: Sven Schuierer
Date: 13/01/2012

Classes :   ShortenFastqEntries

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
 *  Class ShortenFastqEntries
 *
 **********************************************************************************/

public class ShortenFastqEntries {
 
  /***********************************************************************************/
  
  private static void printHelp () {
    System.out.println("ShortenFastqEntries\n" +                                                 
      "   -- a program which outputs fastq entries with new identifiers of the form Fxxxxx/[12].\n" +               
      "\n" +
      "\n" +                                                              
      "usage: command line invocation: je ShortenFastqEntries <options>\n" +   
      "\n" +
      "Options:\n" + 
      "-f <file name>   -  name of the fastq input file with the reads (- for STDIN) [-]\n" + 
      "-F <file name>   -  name of the 2. fastq input file with the reads (- for STDIN) [-]\n" + 
      "-o <file name>   -  name of the output file (- for STDOUT) [].\n" +
      "-O <file name>   -  name of the 2. output file (- for STDOUT) [].\n" +
      "-l <read length> -  length of the reads in the output file.\n" +
      "-u               -  uncompress the fastq files\n" +
      "-h               -  help: display this information.");
  }
                                     
    

  /***********************************************************************************/

 
  public static void main (String [] args) {
    

    try {

      String inputFilename1   = "-";
      String outputFilename1  = "-";
      String inputFilename2   = "";
      String outputFilename2  = "";
 
      int readLength = 1;
      boolean compress = true;
      
      Getopt g = new Getopt("ShortenFastqEntries", args, "f:F:o:O:l:h");

      int c;
      String arg = "";


      c = g.getopt();
      
      while (c  != -1) {
	switch(c) {
	case 'f':
	  arg = g.getOptarg();
	  inputFilename1 = arg;
	  break;
	case 'F':
	  arg = g.getOptarg();
	  inputFilename2 = arg;
	  break;
	case 'o':
	  arg = g.getOptarg();
	  outputFilename1 = arg;
	  break;
	case 'O':
	  arg = g.getOptarg();
	  outputFilename2 = arg;
	  break;
	case 'l':
	  readLength = Integer.parseInt(g.getOptarg());
	  break;
	case 'u':
	  compress = false;
	  break;
	case 'h':
	  printHelp();
	  System.exit(0);
	  break;
	default:
	  System.err.print("Error: getopt() returned " + c + "\n");
	}
	c = g.getopt();
      }


      if (inputFilename1.equals("-") && inputFilename2.equals("-")) {
	throw new IOException ("ERROR: Only one input file can be set to std in.");
      }

      if (outputFilename1.equals("-") && outputFilename2.equals("-")) {
	throw new IOException ("ERROR: Only one output file can be set to std in.");
      }

      if (inputFilename1 != "" && outputFilename1 == "") {
	outputFilename1 = inputFilename1;
      }

      if (inputFilename2 != "" && outputFilename2 == "") {
	outputFilename2 = inputFilename2;
      }
      
      if (UtilLib.getFastqFilenameBase (inputFilename1) == null || UtilLib.getFastqFilenameBase (outputFilename1) == null) {
	throw new Exception ("Please specify the input and output filenames of the fastq files with the correct extension (.fq or .fastq).");
      }

      if (inputFilename2 != "" && UtilLib.getFastqFilenameBase (inputFilename2) == null || UtilLib.getFastqFilenameBase (outputFilename1) == null) {
	throw new Exception ("Please specify the input and output filenames of the 2. fastq files with the correct extension (.fq or .fastq).");
      }

      
      BufferedReader fastqReader1 = null;
      BufferedReader fastqReader2 = null;
      PrintWriter    fastqWriter1 = null;
      PrintWriter    fastqWriter2 = null;

      fastqReader1 = UtilLib.getBufferedReader (inputFilename1);
      if (inputFilename2 != "") {
	fastqReader2 = UtilLib.getBufferedReader (inputFilename2);
      }
      
      int countUnit = 1 * 1000 * 1000;

      String outputFilename1Base = UtilLib.getFastqFilenameBase (outputFilename1);
      String outputFilename2Base = UtilLib.getFastqFilenameBase (outputFilename2);
      
      String outputFilename1BaseSuffix = "";
      if (outputFilename1Base.endsWith("_1") || outputFilename1Base.endsWith("-1")) {
	outputFilename1BaseSuffix = "_1";
	outputFilename1Base = outputFilename1Base.substring(0, outputFilename1Base.length() - 2);
      } else if (outputFilename1Base.endsWith("_2") || outputFilename1Base.endsWith("-2")) {
	outputFilename1BaseSuffix = "_2";
	outputFilename1Base = outputFilename1Base.substring(0, outputFilename1Base.length() - 2);
      }

      String outputFilename2BaseSuffix = "";
      if (outputFilename2Base.endsWith("_1") || outputFilename2Base.endsWith("-1")) {
	outputFilename2BaseSuffix = "_1";
	outputFilename2Base = outputFilename2Base.substring(0, outputFilename2Base.length() - 2);
      } else if (outputFilename2Base.endsWith("_2") || outputFilename2Base.endsWith("-2")) {
	outputFilename2BaseSuffix = "_2";
	outputFilename2Base = outputFilename2Base.substring(0, outputFilename2Base.length() - 2);
      }

      if (outputFilename1BaseSuffix.equals(outputFilename2BaseSuffix)) {
	throw new IOException ("ERROR: the output file basenames have the same suffix: " + outputFilename1BaseSuffix);
      }
      
      String fastqSuffix = "fq.gz";
      if (! compress) {
	fastqSuffix = "fq";
      }

      FastqEntry fastqEntry1 = null;
      FastqEntry fastqEntry2 = null;

      int i = 0;
      
      System.err.println ("Shortening the reads of input files:\n" +
			  (inputFilename1=="-"?"stdin":inputFilename1) + " and \n" +
			  (inputFilename2=="-"?"stdin":inputFilename2) + "\n" +
			  "to length " + readLength + ".");
      try {

	String fastqOutputFilename1 = outputFilename1Base + "-L" + readLength + outputFilename1BaseSuffix + "." + fastqSuffix;
	System.err.println ("Writing to file " + fastqOutputFilename1 + " (line: " + i + ")");
	fastqWriter1 = UtilLib.getPrintWriter (fastqOutputFilename1);
		
	if (inputFilename2 != "") {
	  String fastqOutputFilename2 = outputFilename2Base + "-L" + readLength + outputFilename2BaseSuffix + "." + fastqSuffix;
	  System.err.println ("Writing to file " + fastqOutputFilename2 + " (line: " + i + ")");
	  fastqWriter2 = UtilLib.getPrintWriter (fastqOutputFilename2);
	}

	while (true) {

	  fastqEntry1 = new FastqEntry (fastqReader1);
	  fastqWriter1.println (fastqEntry1.toString (readLength));
	  if (inputFilename2 != "") {
	    fastqEntry2 = new FastqEntry (fastqReader2);
	    fastqWriter2.println (fastqEntry2.toString (readLength));
	  }
	      
	  i++;
	  if (i % countUnit == 0) {
	    System.err.print (".");
	  }

	}

	
      } catch (IOException e) {
	
	if (i >= countUnit) {
	  System.err.println ();
	}
	System.err.println (e==null?"Null error message":e.getMessage());
	
      }

      fastqWriter1.close ();
      fastqWriter2.close ();

    }
    catch (Exception e) {
      System.err.println(e==null?"Null error message":e.getMessage());
    }
  }
  
}
