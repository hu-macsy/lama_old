/*
 * CounterTable.java
 *
 * Abstract table for the user selected counters.
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

package adaptor.EditPM;

import javax.swing.table.AbstractTableModel;

/**
 * The <code>CounterTable</code> contains all counters that will be enabled.
 *
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */

class CounterTable extends AbstractTableModel {
    
    /**
     * These are the headers of the columns in the counter table.
     */
    private static final String[] COLUMN_NAMES = { "Event", "Event Description", "Excl", "Scale", "Width", "Precision" };
    
    /**
     * This is the array of counters managed by this abstract table model.
     */
    private Counter[] myCounterList;
    
    /**
     * Number of entries in use of myCounterList.
     */
    private int numberCounters = 0;
    
    /**
     * Constructor of a counter table. The initial size
     * is zero counters.
     * 
     */
    public CounterTable() {
        
        final int initialSize = 256;
        
        myCounterList = new Counter[initialSize];
        numberCounters = 0;
        
    } 
    
   
    /**
     * Dynamic extension of the array of counters.
     * 
     * @param nr is the minimal size that myCounterList must have
     */
    private void extendCounterList(int nr) {
        
        int sizeList = myCounterList.length;
        
        while (nr >= sizeList) {
            
            // double number of Nodes
            
            Counter[] newList = new Counter[2*sizeList];
            
            for (int i = 0; i < sizeList; i++) {
                newList[i] = myCounterList[i];              
            }
           
            sizeList = 2*sizeList;
            myCounterList = newList;
            
        } // while not enough size
        
    } 
    
    /**
     * This routine adds a new counter to the counter table.
     * 
     * @param counter is the new counter
     */
    public void addCounter(Counter counter) {
        
        extendCounterList(numberCounters);
        
        myCounterList[numberCounters] = counter;
        numberCounters++;
    }
    
    /**
     * Adding a new counter at a certain position.
     * 
     * @param row is the position where to insert the counter
     * @param counter will be inserted.
     */
    public void addCounter(int row, Counter counter) {
        
        // move the counters above row one position up
        
        for (int i = numberCounters - 1; i >= row; i--) {
            
            myCounterList[i + 1] = myCounterList[i];
        }
        
        myCounterList[row] = counter;
        numberCounters++;
    }
     
    /**
     * Adding a new counter for a performance event.
     * 
     * @param event is the performance event.
     */
    public void addCounter(PerformanceEvent event) {
        
        addCounter(new Counter(event));
    }
    
    /**
     * {@inheritDoc}
     *
     * @see javax.swing.table.TableModel#getColumnCount()
     */
    public int getColumnCount() {
        return COLUMN_NAMES.length;
    }
    
    /**
     * {@inheritDoc}
     *
     * @see javax.swing.table.TableModel#getRowCount()
     */
    public int getRowCount() {
        
        return numberCounters;
    }
    
    /**
     * {@inheritDoc}
     *
     * @see javax.swing.table.TableModel#getColumnName(int)
     */
    public String getColumnName(int col) {
        
        return COLUMN_NAMES[col];
    }
    
    /**
     * {@inheritDoc}
     *
     * @see javax.swing.table.TableModel#getValueAt(int, int)
     */
    public Object getValueAt(int row, int col) {
        
        // row should be a legal index for CounterList
        
        return myCounterList[row].getCounterProperty(col);
        
    } // method getValueAt
    
    /**
     * {@inheritDoc}
     *
     * @see javax.swing.table.TableModel#getColumnClass(int)
     */
    public Class getColumnClass(int c) {
        
        return getValueAt(0, c).getClass();
    }
    
    /**
     * {@inheritDoc}
     *
     * @see javax.swing.table.TableModel#isCellEditable(int, int)
     */
    public boolean isCellEditable(int row, int col) {
        
        return true;
    }
    
    /**
     * {@inheritDoc}
     *
     * @see javax.swing.table.TableModel#setValueAt(java.lang.Object, int, int)
     */
    public void setValueAt(Object value, int row, int col) {
        
        fireTableCellUpdated(row, col);
        
        myCounterList[row].setCounterProperty(col, value);
        
    } // setValueAt
    
    
    /**
     * This routine is used to delete a user specified counter.
     * 
     * @param row specifies the row from which the counter will be removed
     * @return the deleted counter
     */
    public Counter deleteCounter(int row) {
        
        // illegal index with throw IndexOutOfBoundException
        
        // at first we save the deleted counter to return nit
        
        Counter deletedCounter = myCounterList[row];
        
        for (int i = row; i < numberCounters - 1; i++) {
              
            myCounterList[i] = myCounterList[i + 1];            
        }
        
        numberCounters--;
        
        return deletedCounter;
    }

    /**
     * Clears the counter table.
     * 
     */
    public void clear() {
        
        numberCounters = 0;
        
    }

    /**
     * Getter routine for a counter at a certain position.
     * 
     * @param row is the position for the counter
     * @return the wanted counter
     */
    public Counter getCounter(int row) {

        return myCounterList[row];
    }
    
} // class CounterTable
