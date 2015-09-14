/**File: FilterUnalignedFastqEntries

Original Author: Sven Schuierer
Date: 13/01/2012

Classes :   FilterUnalignedFastqEntries

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
 *  Class FilterUnalignedFastqEntries
 *
 **********************************************************************************/

public class FilterUnalignedFastqEntries {
  
 
  /***********************************************************************************/
  
  private static void printHelp () {
    System.out.println("FilterUnalignedFastqEntries\n" +                                                 
      "   -- a program which outputs only those fastq entries that occur in both input files.\n" +               
      "\n" +
      "\n" +                                                              
      "usage: command line invocation: je FilterUnalignedFastqEntries <options>\n" +   
      "\n" +
      "Options:\n" + 
      /* "-c <config file name> -  name of a config file (default: \"stringIndexSearch.cfg\")\n" + */
      "-C <common prefix> - Common prefix of all read ids which is removed.\n" +
      "   [<empty>]\n" +
      "-i <increment> -  increment of start position for reads to be ouput [5]\n" + 
      "-l <length> -     length of the reads to be output; output complete sequence\n" +
      "                  if negative [-1]\n" + 
      "-f <file name> -  name of the fastq input file with the reads (- for STDIN) [-]\n" + 
      "-w <file name> -  name of the read weight file with the number of alignments\n" +
      "                  for the reads\n" + 
      "-o <file name> -  name of the output file (- for STDOUT) [-].\n" +
      "-h             -  help: display this information.");
  }
                                     
    

  /***********************************************************************************/

 
  public static void main (String [] args) {
    

    try {

      String configFilename = "filterUnalignedFastqEntries.cfg";
      String inputFilename      = "-";
      String readWeightFilename = "";
      String outputFilename     = "-";
      String commonPrefix       = "";

      int inc    = 5;
      int length = -1;
      
      Getopt g = new Getopt("FilterUnalignedFastqEntries", args, "c:C:f:i:l:o:w:h");

      int c;
      String arg = "";


      c = g.getopt();
      
      while (c  != -1) {
	switch(c) {
	case 'c':
	  arg = g.getOptarg();
	  configFilename = arg;
	  break;
	case 'C':
	  commonPrefix = g.getOptarg();
	  break;
	case 'f':
	  arg = g.getOptarg();
	  inputFilename = arg;
	  break;
	case 'i':
	  inc = Integer.parseInt(g.getOptarg());
	  break;
	case 'l':
	  length = Integer.parseInt(g.getOptarg());
	  break;
	case 'o':
	  arg = g.getOptarg();
	  outputFilename = arg;
	  break;
	case 'w':
	  arg = g.getOptarg();
	  readWeightFilename = arg;
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

      /*
	try {
	loadConfigFile(configFilename);
	}
	catch (IOException e) {
	System.out.println ((e==null?"ERROR: could not load " + configFilename +
	" - null pointer exception":e.getMessage()));
	}
      */
      
      if (commonPrefix != "") {
	UtilLib.setCommonPrefix (commonPrefix);
      }

      if (! UtilLib.checkFastqFilename (inputFilename)  || ! UtilLib.checkFastqFilename (outputFilename)) {
	throw new Exception ("Please specify the input and output filenames of the Fastq files.");
      }

      BufferedReader fastqReader  = UtilLib.getBufferedReader (inputFilename);
      BufferedReader weightReader = UtilLib.getBufferedReader (readWeightFilename);

      PrintWriter writer = UtilLib.getPrintWriter (outputFilename);

      Hashtable<String, FragmentEntry> fragmentEntryTable = new Hashtable<String, FragmentEntry> ();
      
      int countUnit = 1000 * 1000;

      try {

	FastqEntry fastqEntry = new FastqEntry (fastqReader);
	
	String        line               = weightReader.readLine ();
	FragmentEntry fragmentEntry      = new FragmentEntry (line);
	String        modifiedFragmentId = UtilLib.modifyFragmentId (fragmentEntry.getFragmentName ());

	int i = 0;;
	while (line != null) {

	  if (modifiedFragmentId.equals(UtilLib.modifyFragmentId (fastqEntry.getFragmentId()))) {
	    line = weightReader.readLine ();
	    if (line != null) {
	      fragmentEntry      = new FragmentEntry (line);
	      modifiedFragmentId = UtilLib.modifyFragmentId (fragmentEntry.getFragmentName ());
	    }
	  } else {
	    if (length < 0) {
	      writer.println (fastqEntry.toString ());
	    } else {
	      writer.println (fastqEntry.toString(length, inc));
	    }
	  }

	  fastqEntry = new FastqEntry (fastqReader);

	  i++;
	  if (i % countUnit == 0) {
	    System.err.print (".");
	    // System.err.println ("FastqEntryTable size: " + fastqEntryTable.size ());
	  }
	}

	while (true) {
	  fastqEntry = new FastqEntry (fastqReader);

	  if (length < 0) {
	    writer.println (fastqEntry.toString ());
	  } else {
	    writer.println (fastqEntry.toString(length, inc));
	  }
	  
	}
	
      } catch (IOException e) {
	
	System.err.println ();
	System.err.println (e==null?"Null error message":e.getMessage());
	
      }

      writer.close ();

    }
    catch (Exception e) {
      System.err.println(e==null?"Null error message":e.getMessage());
    }
  }
  
}
