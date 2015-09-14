/**File:Interval.java 

Original Author: Sven Schuierer
Date: 30/09/2014

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
 *                              Class Interval
 *
 ***********************************************************************************/

public class Interval implements Comparable<Interval> {

  private static final int ODD_PRIME_NUMBER = 37;

  private static int firstTerm(int seed){
    return ODD_PRIME_NUMBER * seed;
  }

  public static int hash(int seed , int value) {
    return firstTerm(seed) + value;
  }

  private static boolean quantify = false;
  public static void setQuantify () {
    quantify = true;
  }
  
  /***********************************************************************************
   *
   *                         Object variables
   *
   ***********************************************************************************/

  private String chromosome = "";
  private int [] intervalArray = {0, 0};
  private HashSet<String> fragmentNameSet = null;

  
  /***********************************************************************************
   *
   *                           Constructors
   *
   ***********************************************************************************/

  Interval (String chromosome, int left, int right) {
    this.chromosome = chromosome;
    this.intervalArray[0] = left;
    this.intervalArray[1] = right;
  }


  Interval (String chromosome, int [] intervalArray) {
    this.chromosome = chromosome;
    this.intervalArray = intervalArray;
  }

  Interval (String chromosome, int [] intervalArray, String fragmentName) {
    this(chromosome, intervalArray);
    if (fragmentNameSet == null) {
      fragmentNameSet = new HashSet<String> (10);
    }
    fragmentNameSet.add(fragmentName);
  }

  
  /***********************************************************************************
   *
   *                           get functions
   *
   ***********************************************************************************/
  
  public String getChromosome () {
    return chromosome;
  }
  
  public int [] getIntervalArray () {
    return intervalArray;
  }

  public int getStart () {
    return intervalArray[0];
  }

  public int getEnd () {
    return intervalArray[1];
  }

  public HashSet<String> getFragmentNameSet () {
    return fragmentNameSet;
  }

  public void addFragmentNameSet (HashSet<String> fragmentNameSet) {
    if (this.fragmentNameSet != fragmentNameSet) {
      this.fragmentNameSet.addAll (fragmentNameSet);
    }
  }

  public void addFragmentNameSet (Interval interval) {
    addFragmentNameSet(interval.getFragmentNameSet ());
  }

  public void setFragmentNameSet (HashSet<String> fragmentNameSet) {
    this.fragmentNameSet = fragmentNameSet;
  }

  public void clearFragmentNameSet () {
    this.fragmentNameSet.clear();
  }

  
  /***********************************************************************************
   *
   *                           compareTo
   *
   ***********************************************************************************/
  
  public int compareTo (Interval interval) {

    String intervalChromosome = interval.getChromosome ();
    int chromosomeComp = chromosome.compareTo(intervalChromosome);
    if (chromosomeComp != 0) {
      return chromosomeComp;
    }

    int [] intervalArray = interval.getIntervalArray ();
    if (this.intervalArray[0] != intervalArray[0]) {
      return UtilLib.compareInt (this.intervalArray[0], intervalArray[0]);
    }

    int returnValue = UtilLib.compareInt (this.intervalArray[1], intervalArray[1]);
    if (returnValue == 0) {
      addFragmentNameSet (interval);
      interval.setFragmentNameSet (getFragmentNameSet());
    }
      
    return returnValue;
    
  }


  /***********************************************************************************
   *
   *                  hashCode, equals, and toString
   *
   ***********************************************************************************/
  
  public int hashCode () {
    return hash(Arrays.hashCode (intervalArray), chromosome.hashCode());
  }

  
  public boolean equals (Object o) {

    Interval interval = (Interval) o;
    int [] compIntArray = interval.getIntervalArray ();
    String compChrom = interval.getChromosome ();
    
    boolean isEqual = chromosome.equals(compChrom) && intervalArray[0] == compIntArray[0] && intervalArray[1] == compIntArray[1];

    if (isEqual && quantify) {
      addFragmentNameSet (interval.getFragmentNameSet());
    }

    return isEqual;
  }

  
  public String toString () {
    
    return chromosome + "\t" + intervalArray[0] + "\t" + intervalArray[1];
    
  }
}
