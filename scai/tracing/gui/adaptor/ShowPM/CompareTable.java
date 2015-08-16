/*
 * CompareTable.java
 * 
 * CompareTabel is an abstract table to compare performance data from
 * different runs.
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

package adaptor.ShowPM;

/**
 * CompareTable is an abstract table to compare performance data from
 * different runs. Processor, thread and counter are fixed.
 *
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */
class CompareTable extends PMTableModel {

    /**
     * This is a constant for an unselected region or user counter.
     */
    private static final int NO_SELECTION = -1;
    
    /**
     * For the comparison table we do not support a popup menu at this time.
     */
    private static final String[] MENU_ITEMS = {};

    /**
     * These are all the data sets used for comparison.
     */
    private PMData[] myPMDataSets = null;
    
    /**
     * The processor is fixed for the compare table.
     */
    private int selectedProc;
    
    /**
     * The thread is fixed for the compare table.
     */
    private int selectedThread; // selected processor counter
    
    /**
     * The region counter is fixed or not selected at all.
     */
    private int selRegionCounter = NO_SELECTION;
    
    /**
     * The user counter is fixed or not selected at all. 
     */
    private int selUserCounter = NO_SELECTION;

    /**
     * This is the constructor for a new compare table. The
     * array may be larger than the number of data sets really
     * used.
     * 
     * @param dataSets are the PM data sets
     * @param numberSets is the number of used PM data sets
     */
    public CompareTable(PMData[] dataSets, int numberSets) {
    
        myPMDataSets = new PMData[numberSets];
    
        // we copy the the (pointer to the) data sets
    
        for (int i = 0; i < numberSets; i++) {
    
            myPMDataSets[i] = dataSets[i];
        }
    
        selectedProc = 0;
        selectedThread = 0;
        selRegionCounter = 0;
    
    }

    /**
     * {@inheritDoc}
     *
     * @see javax.swing.table.TableModel#getColumnCount()
     */
    public int getColumnCount() {

        // The first column is the for the region names
        
        int numberCols = 1;
        
        if (myPMDataSets != null) {
            
            // for each data set we have an additional column
            
            numberCols += myPMDataSets.length;
        }

        return numberCols;
    }

    /**
     * {@inheritDoc}
     *
     * @see javax.swing.table.TableModel#getRowCount()
     */
    public int getRowCount() {

        // number of rows is given by the number of regions

        return myPMDataSets[0].getNumberRegions();
    }

    /**
     * {@inheritDoc}
     *
     * @see javax.swing.table.TableModel#getColumnName(int)
     */
    public String getColumnName(int col) {

        String name = "Region";
        
        if (col > 0) {
            
            name = myPMDataSets[col - 1].getName();
        }

        return name; 

    }

    /**
     * {@inheritDoc}
     *
     * @see javax.swing.table.TableModel#getValueAt(int, int)
     */
    public Object getValueAt(int row, int col) {

        if (col == 0) {
            
            return myPMDataSets[0].getRegion(row).getName();            
        }

        long[] values = myPMDataSets[col - 1].getRegionCounterVals(selectedProc, selectedThread, row);
        
        if (selRegionCounter >= 0) {
            
            return new Long(values[selRegionCounter]);
            
        } else {
            // important: take the user counter from the first data set
            // so we have not to make sure that all sets have the same
            // user counters available
            
            UserCounter counter = myPMDataSets[0].getUserCounter(selUserCounter);
            
            return new Double(counter.getScaledValue(values));
        }
    }

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

        return false;

    }

    /**
     * Set a new selection for user counter.
     * 
     * @param counter is the index of the selected user counter
     * @param proc is the selected processor
     * @param thread is the selected thread
     */
    void setSelectedUC(int counter, int proc, int thread) {

        selectedThread = thread;
        selectedProc = proc;
        selUserCounter = counter;
        selRegionCounter = NO_SELECTION;
    }

    /**
     * Set a new selection for region counter.
     * 
     * @param counter is the index of the selected region counter
     * @param proc is the selected processor
     * @param thread is the selected thread
     */
    void setSelectedRC(int counter, int proc, int thread) {

        selectedThread = thread;
        selectedProc = proc;
        selUserCounter = NO_SELECTION;
        selRegionCounter = counter;

    }

    /**
     * This routine returns the counter description of the currently
     * selected user or region counter for this table.
     * 
     * @return description of the selected counter for the table
     */
    String getCounterDescription() {

        String description = "<no counter>";

        if (selUserCounter != NO_SELECTION) {
            
            description = myPMDataSets[0].getCounterHeader(true, selUserCounter);
            
        } else if (selRegionCounter != NO_SELECTION) {
            
            description = myPMDataSets[0].getCounterHeader(false, selRegionCounter);            
        }

        return description;

    } // getCounterDescription
    
    /**
     * {@inheritDoc}
     *
     * @see adaptor.ShowPM.PMTableModel#getPopupItems()
     */
    public String[] getPopupItems() {
        
        return MENU_ITEMS;
    }
    
    /**
     * {@inheritDoc}
     *
     * @see adaptor.ShowPM.PMTableModel#actionPopupItem(int, int, int)
     */
    public void actionPopupItem(int kind, int col, int row) {

         // never called as we have no MENU_ITEMS here
    }
    
    /**
     * {@inheritDoc}
     *
     * @see adaptor.ShowPM.PMTableModel#mousePressed(int, int)
     */
    public void mousePressed(int col, int row) {
        
    }

} // PMCompareTable
