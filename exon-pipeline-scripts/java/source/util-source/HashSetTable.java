/**File: HashSetTable.java

Original Author: Sven Schuierer
Date: 26/01/2007

Filename    : $RCSfile: HashSetTable.java,v $

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



import java.util.*;


/***********************************************************************************
 *
 *
 *
 *
 ***********************************************************************************/

public class HashSetTable<E, F> extends Hashtable<E, HashSet<F>> {


  HashSetTable () {
    super ();
  }


  HashSetTable (int capacity) {
    super (capacity);
  }

  

  /**********************************************************************************/
  
  public HashSet<F> add (E key, F value) {

    HashSet<F> values = super.get(key);
    if (values == null) {
      values = new HashSet<F> ();
      super.put(key, values);
    }

    values.add(value);
    return values;
    
  }

  
  /**********************************************************************************/
  
  public void putValue (E key, F value) {

    HashSet<F> values = super.get(key);
    if (values == null) {
      values = new HashSet<F> ();
      super.put(key, values);
    }

    values.add(value);
    
  }

  /**********************************************************************************/
  
  public void putSet (E key, HashSet<F> valueSet) {

    HashSet<F> values = super.get(key);
    if (values == null) {
      values = new HashSet<F> ();
      super.put(key, values);
    }

    values.addAll(valueSet);
    
  }


  /**********************************************************************************/
  
  public void remove (E key, F value) {

    HashSet<F> values = super.get(key);
    if (values == null) {
      return;
    }

    values.remove(value);

    if (values.size() == 0) {
      super.remove(key);
    }
    
  }


  /**********************************************************************************/
  
  public void removeNotEmpty (E key, F value) {

    HashSet<F> values = super.get(key);
    
    if (values == null) {
      return;
    }

    if (! values.contains(value)) {
      return;
    }

    if (values.size () > 1) {
      values.remove(value);
    }

  }

  
  /**********************************************************************************/
  
  public HashSet<F> getSet (E key) {

    return super.get(key);
    
  }
  

  /**********************************************************************************/
  
  public Iterator getIterator (E key) {

    HashSet<F> values = super.get(key);
    
    if (values != null) {
      return values.iterator ();
    }
    
    return null;
    
  }


  /**********************************************************************************/
  
  public HashSet<E> getKeys () {

    HashSet<E> result = new HashSet<E>();

    Enumeration<E> e = super.keys();

    while (e.hasMoreElements()) {

      result.add(e.nextElement());
      
    }

    return result;
    
  }


  /**********************************************************************************/
  
  public HashSet<E> getKeys (F value) {

    HashSet<E> result = new HashSet<E>();

    Enumeration<E> e = super.keys();

    while (e.hasMoreElements()) {

      E key = e.nextElement();

      HashSet<F> values = super.get(key);
      if (values != null && values.contains(value)) {
	result.add(key);
      }
      
    }

    return result;
    
  }

  /**********************************************************************************/
  
  public HashSet<E> getKeysExcluding (F value) {

    HashSet<E> result = new HashSet<E>();

    Enumeration<E> e = super.keys();

    while (e.hasMoreElements()) {

      E key = e.nextElement();

      HashSet<F> values = super.get(key);
      if (values == null || ! values.contains(value)) {
	result.add(key);
      }
      
    }

    return result;
    
  }
  
  

  /*********************************************************************************/
  
  public String toString () {
    
    String returnString = "";

    Enumeration e = super.keys ();
      
    if (e != null) {
      while (e.hasMoreElements()) {
	Object key = e.nextElement();
	returnString = returnString + key.toString() + ": ";
	HashSet<F> values = super.get(key);
	if (values != null) {
	  Iterator it = values.iterator();
	  while (it.hasNext()) {
	    returnString = returnString + it.next().toString() + (it.hasNext()?", ":"");
	  }
	}
	returnString = returnString + (e.hasMoreElements()?"; ":"");
      }
    }
    
    
    return returnString;
    
  }

  
  
}
