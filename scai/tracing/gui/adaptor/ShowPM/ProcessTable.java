/*
 * TableDisplay.java
 * 
 * Utitilies realizes some help functions for strings.
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

import org.apache.log4j.Logger;

/**
 * The ProcessTable is a table that
 * shows the values of one counter for all regions and processors.
 *
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */
class ProcessTable extends PMTableModel {

    /**
     * These are the menu items for the popup menu of this class.
     */
    private static final String[] MENU_ITEMS = { "Show region", "Info region"};
    

    /**
     * Logger for this class.
     */
    private static Logger logger = Logger.getLogger(ProcessTable.class);

    /**
     * The performance data on which this process table works.
     */
    private PMData myPMData;
    
    /**
     * The interface to get and make selections.
     */
    private SelectedInterface myShowPM;

    /**
     * Creates a process table on which the routines will work.
     * 
     * @param inPM is the input PM data
     * @param main is the implementation of a selection interface
     */
    ProcessTable(PMData inPM, SelectedInterface main) {

        myPMData = inPM;
        myShowPM = main;
        
        logger.info("new ProcessTable has been created"); 
    }

    /**
     * {@inheritDoc}
     *
     * @see javax.swing.table.TableModel#getColumnCount()
     */
    public int getColumnCount() {

        return (myPMData.getNumberThreads() * myPMData.getNumberProcesses()) + 1;
    }

    /**
     * {@inheritDoc}
     *
     * @see javax.swing.table.TableModel#getRowCount()
     */
    public int getRowCount() {

        // number of rows is given by the number of regions

        return myPMData.getNumberRegions();
    }

    /**
     * {@inheritDoc}
     *
     * @see javax.swing.table.TableModel#getColumnName(int)
     */
    public String getColumnName(int col) {

        if (col == 0) {
            
            return "Region";
            
        } else {
            
            // col - 1 stands for the processor
            
            return myPMData.getProcessorString(col - 1);
        }


    }

    /**
     * {@inheritDoc}
     *
     * @see javax.swing.table.TableModel#getValueAt(int, int)
     */
    public Object getValueAt(int row, int col) {

        if (col == 0) {
            
            return myPMData.getRegion(row).getName();         
        }

        // col - 1 is index processor, row is the region
        
        long[] counterVals = myPMData.getRegionCounterVals(col - 1, row);
        
        int indexCounter = myShowPM.getSelectedRegionCounter();
        
        if (indexCounter >= 0) {
            
            // so the counter is really a region counter
            
            return new Long(counterVals[indexCounter]);
            
        } else {

            // so the counter is a user metric counter
            
            indexCounter = myShowPM.getSelectedUserCounter();
            
            UserCounter counter = myPMData.getUserCounter(indexCounter);
            
            return new Double(counter.getScaledValue(counterVals));
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
     
        if (kind == 0) {
            
            // item "Show region"
            
            myShowPM.showFile(myPMData, row);
            
        } else if (kind == 1) {
            
            myShowPM.infoFile(myPMData, row);
        }
    }
    
    /**
     * {@inheritDoc}
     *
     * @see adaptor.ShowPM.PMTableModel#mousePressed(int, int)
     */
    public void mousePressed(int col, int row) {

        logger.info("ProcessTable: row = " + row + ", col = " 
                    + col + " has been selected by mouse"); 

        if (col == 0) {
            
            actionPopupItem(0, col, row);

        }
        
    }

}
