/*
 * RegionTable.java
 *
 * Abstract Table for all Regions.
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

package adaptor.EditRIF;

import javax.swing.table.AbstractTableModel;

import adaptor.General.FileDescriptor;
import adaptor.General.RegionDescriptor;

/**
 * The class RegionTable defines an abstract table that is used to make
 * a display out of it.
 *
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */

class RegionTable extends AbstractTableModel {
    
    // data of the table, allocated later
    
    /**
     * Pointer to the editor of RIF data to call interface routines.
     */
    private EditRIFInterface myEditRIF;
    
    /**
     * This routine creates a region table. 
     * 
     * @param editor is the interface to access region descriptors
     */
    public RegionTable(EditRIFInterface editor) {
        
        myEditRIF = editor;
        
    } // SetRegionTable
    
    /**
     * {@inheritDoc}
     *
     * @see javax.swing.table.TableModel#getColumnCount()
     */
    public int getColumnCount() {
        
        return RegionDescriptor.REGION_ENTRIES.length;
    }
    
    /**
     * {@inheritDoc}
     *
     * @see javax.swing.table.TableModel#getRowCount()
     */
    public int getRowCount() {
                
        return myEditRIF.noRegions();
    
    }
    
    /**
     * {@inheritDoc}
     *
     * @see javax.swing.table.TableModel#getColumnName(int)
     */
    public String getColumnName(int col) {
        
        return RegionDescriptor.REGION_ENTRIES[col];
    }
    
    /**
     * {@inheritDoc}
     *
     * @see javax.swing.table.TableModel#getValueAt(int, int)
     */
    public Object getValueAt(int row, int col) {
        
        RegionDescriptor theRegion = myEditRIF.getRegion(row);
        
        return theRegion.getRegionObject(col);
    }
    
    /**
     * {@inheritDoc}
     *
     * JTable uses this method to determine the default renderer/
     * editor for each cell.  If we didn't implement this method,
     * then the last column would contain text ("true"/"false"),
     * rather than a check box.
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
        
        RegionDescriptor theRegion = myEditRIF.getRegion(row);
        
        return theRegion.isEditable(col);
        
    }
    
    /**
     * {@inheritDoc}
     *
     * @see javax.swing.table.TableModel#setValueAt(java.lang.Object, int, int)
     */
    public void setValueAt(Object value, int row, int col) {
        
        RegionDescriptor theRegion = myEditRIF.getRegion(row);
        
        theRegion.setRegionObject(value, col);
        fireTableCellUpdated(row, col);
        
        if (col == 6) {
            
            fireTableCellUpdated(row, 7);
        }
        
        // now we assume that the RIF file has been changed
        
        myEditRIF.setModified();
    }
    
    // get longest value of a column 
    
    /**
     * This routine returns for a column the maximal length
     * of characters that is used for the representation.
     * 
     * @param column for which we want to find the maximal length
     * @return the maximal number of chars used for the column
     */
    public Object longValue(int column) {
        
        int n = myEditRIF.noRegions(); // number of Regions
        
        int maxLength = 0;
        
        Object resobj = null;
        
        for (int i = 0; i < n; i++) {
            
            RegionDescriptor theRegion = myEditRIF.getRegion(i);
            
            Object obj = theRegion.getRegionObject(column);
            
            int len = obj.toString().length();
            
            if (len > maxLength) {
                resobj = obj;
                maxLength = len;
            }
        }
        
        return resobj;
        
    } // longValue
    
    /**
     * Return the file descriptor for the region in row 'row'.
     * 
     * @param row is the number of row
     * @return the file descriptor for the region in row
     */
    public FileDescriptor getFileDescriptor(int row) {
        
        RegionDescriptor theRegion = myEditRIF.getRegion(row);
        
        return theRegion.getFile();
    }
    
    /**
     * This routine returns the start line of the region in row.
     * 
     * @param regionId is the region id
     * @return the first line in the source code file
     */
    public int getLineStart(int regionId) {

        RegionDescriptor theRegion = myEditRIF.getRegion(regionId);
        
        int line = theRegion.getFirstLine();
        
        return line;
    }

    /**
     * This routine returns the last line of the region.
     * 
     * @param row is the region id row
     * @return the line number of last source code line
     */
    public int getLineStop(int row) {
        
        RegionDescriptor theRegion = myEditRIF.getRegion(row);
                
        int line = theRegion.getLastLine();
        
        return line;
    }
    
} // class RegionTable 
