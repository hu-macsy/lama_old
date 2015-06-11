/*
 * MyHashTable.java
 * 
 * Frame that displays the call graph.
 * 
 * Created: 2006-02-20 Thomas Brandes <thomas.brandes@scai.fraunhofer.de>
 * Changed:
 * 
 * $Id$
 * 
 * Copyright (C) 2006 Fraunhofer SCAI, Germany
 * 
 * All rights reserved
 *
 * http://www.scai.fhg.de/EP-CACHE/adaptor
 */
package adaptor.Calltree;

import org.apache.log4j.Logger;

/**
 * Own class to implement efficient hash table for objects
 * of type CTNode and CTEdge.
 *
 * @version $LastChangedRevision$
 * @author Dr. Thomas Brandes
 */

public class MyHashTable {

    /**
     * Logger for this class.
     */
    private static Logger logger = Logger.getLogger(MyHashTable.class);

    /**
     * Own class for entries of the hash table.
     *
     * @version $LastChangedRevision$
     * @author Thomas Brandes
     */
    private class HashTableEntry {

        /**
         * The name attribute for the entry.
         */
        private Object name;
        
        /**
         * The value attribute for the entry.
         */
        private int value;
        
        /**
         * Pointer to the next entry in the hash table.
         */
        private HashTableEntry next;

        /**
         * Constructor for a simple hash table entry.
         * 
         * @param theName is the name or key field
         * @param theValue is the value to be stored for name
         */
        HashTableEntry(Object theName, int theValue) {

            this.name = theName;
            this.value = theValue;
            this.next = null;

        }

        /**
         * Add an entry behind me.
         * 
         * @param entry will be appended
         */
        void appendEntry(HashTableEntry entry) {

            if (next == null) {
                
                next = entry;
                return;
            }

            next.appendEntry(entry);
        }

        /**
         * Searches a specific name in my chain.
         * 
         * @param searchName is the name to be search
         * @return the full entry with the searched name
         */
        HashTableEntry findEntry(Object searchName) {

            HashTableEntry entry = null;
            
            if (name.equals(searchName)) {

                entry = this;
                
            } else if (next != null) {
                
                entry = next.findEntry(searchName);
                
            }
            
            return entry;


        } // findEntry
        
        /**
         * This function returns for the collision list in the table
         * all nodes that equals the searched name.
         * 
         * @param searchName is the object we want to find
         * @return the number of entries found
         */
        int countEntries(Object searchName) {
            
            int n = 0;         
            
            HashTableEntry elem = this;
            
            while (elem != null) {
                
                if (elem.name.equals(searchName)) { 
                    
                    n++;
                }

                elem = elem.next;
                
            }

            return n;
        }

        /**
         * This routine sets for all entries of the collision list in the
         * array values the value of the objects where the name equals.
         * 
         * @param searchName is the name searched
         * @param values is an array whose size must be equal to countEntries(searchName)
         */
        
        public void setEntries(Object searchName, int[] values) {

            int n = 0;         
            
            HashTableEntry elem = this;
            
            while (elem != null) {
                
                if (elem.name.equals(searchName)) { 
                    
                    values[n++] = elem.value;
                }

                elem = elem.next;
                
            }
            
        }

    } // class HashTableEntry

    /**
     * This will be the siye of my hash table.
     */
    private int hashSize = 0;

    /**
     * Hash table that will be allocated for a certain size.
     * Conflict entries will be appended as a linked list.
     */
    private HashTableEntry[] table;

    /**
     * The constructor creates a new hash table of a certain size.
     * @param size specifies the size of the HashTable
     */
    public MyHashTable(int size) {

        hashSize = size;
        table = new HashTableEntry[hashSize];
        
        logger.info("Hash table of " + size + " entries generated");
    }

    /**
     * Add a new entry (string, int) in this hash table.
     * @param name is the string for which we store a key
     * @param val is the value associated with name
     */
    public void addEntry(Object name, int val) {

        int key = Math.abs(name.hashCode() % hashSize);
        
        HashTableEntry entry = new HashTableEntry(name, val);

        if (table[key] == null) {

            table[key] = entry;
            
        } else {
            
            table[key].appendEntry(entry);            
        }

    } // addEntry

    /**
     * Query the key for (Object name, int key) in the table.
     * 
     * @param name is the Object to be searched
     * @return key value for the string (-1 if not found)
     */
    int queryEntry(Object name) {

        final int notFoundVal = -1;
        
        int key = Math.abs(name.hashCode() % hashSize);
         
        int val = notFoundVal;

        HashTableEntry entry = table[key];

        if (entry != null) {
            
            entry = entry.findEntry(name);            
        }

        if (entry != null) {
             
            val = entry.value;
        }

        return val;

    } // queryEntry

    /**
     * This routine returns all for all entries in the Hash table
     * with a given name the corresponding key value.
     * 
     * @param name is the name we search for (name, value)
     * @return an array of values with (name, value) in the table
     */
    int[] queryEntries(Object name) {
        
        int[] values = null;
        
        int key = Math.abs(name.hashCode() % hashSize);
        
        HashTableEntry entry = table[key];
        
        if (entry != null) {
            
            int numberMatches = entry.countEntries(name);
            
            if (numberMatches > 0) {

                values = new int[numberMatches];
            
                entry.setEntries(name, values);
            }
            
        }
        
        return values;
    }
    
} // MyHashTable
