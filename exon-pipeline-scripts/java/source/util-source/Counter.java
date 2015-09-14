/**File: ConvertSamBed.java 

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
 *                              Class Counter
 *
 ***********************************************************************************/

public class Counter {

  private int counter;
  private int countDiv;
  private int curDigits;
  private int digits;
  private String prefix;

  public Counter (int d) {
    
    digits    = d;

    counter   = 1;
    countDiv  = 10;
    curDigits = 1;
    
    prefix = computePrefix ();
    
  }

  public String computePrefix (int digits) {
    String prefix = "";
    for (int i = 0; i < digits - curDigits; i++) {
      prefix = prefix + "0";
    }
    return (prefix);
  }
  
  public String computePrefix () {
    return (computePrefix(digits));
  }


  public void setDigits (int d) {
    digits = d;
  }

  public void inc () throws IOException {
    
    counter++;
    if (counter >= countDiv) {

      if (prefix.length() == 0) {
	throw new IOException ("Counter is larger than maximal value: " + (counter - 1));
      }
      
      countDiv = countDiv * 10;
      curDigits++;
      prefix = prefix.substring(0, prefix.length()-1);     
    }

  }


  public void dec () throws IOException {
    counter--;

    if (counter < 0) {
      System.err.println("Counter is smaller than 0.");
      throw new IOException ("Counter is less than 0.");
    }
    
    if (counter < countDiv / 10) {
      countDiv = countDiv / 10;
      curDigits--;
      prefix = prefix + "0";     
    }
  }

  public int getCount() {
    return (counter);
  }

  public String getZeroFilledCount () {
    return prefix + counter;
  }


  public String toString () {
    return getZeroFilledCount ();
  }



  /***********************************************************************************/

  public static void main (String [] args) {

    try {
      Counter counter = new Counter (6);
      
      for (int i = 0; i < 100; i++) {
	
	counter.inc();
	
	if (i % 4 == 0) {
	  System.out.println (counter);
	}
	
      }
      
      for (int i = 0; i < 98; i++) {
	
	counter.dec();
	
	if (i % 4 == 0) {
	  System.out.println (counter);
	}
	
      }
      
    }
    catch (IOException e) {
      System.out.println ((e==null?"Null message in Counter":e.getMessage ()));
    }

  }

}

