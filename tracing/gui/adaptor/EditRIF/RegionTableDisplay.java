/*
 * RegionTableDisplay.java
 *
 * Panel to dipslay the region table
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

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.event.MouseEvent;

import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSplitPane;
import javax.swing.JTable;
import javax.swing.ListSelectionModel;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;

import org.apache.log4j.Logger;

import adaptor.Calltree.TableSorter;
import adaptor.General.FileDescriptor;
import adaptor.General.ShowFilePanel;

/**
 * This class is the panel for the display of the region table.
 *
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */
public class RegionTableDisplay extends JPanel {

    /**
     * Default dimension for the region tabel display.
     */
    private static final Dimension MY_DIMS = new Dimension(800, 300);
    
    /**
     * Logger for this class.
     */
    private static Logger logger = Logger.getLogger(RegionTableDisplay.class);
    
    /**
     * Anchor to the editor to which this panel belongs.
     */
    private EditRIFInterface myEditRIF;
    
    /**
     * The region table provides the abstract table routines needed for the JTable.
     */
    private RegionTable myRegionTable; // is the abstract table containing the data

    /**
     * This is the JTable for the regions.
     */
    private JTable myRegionJTable;
    
    /**
     * The region table can be sorted via this class.
     */
    private TableSorter mySortedRegionTable;
    
    /**
     * Panel to show files that contain a selected region.
     */
    private ShowFilePanel myShowFilePanel;

    /**
     * The region table display is constructed by the editor.
     * 
     * @param editor has information about all needed regions.
     */
    public RegionTableDisplay(EditRIFInterface editor) {

        super(new BorderLayout());

        myEditRIF = editor;
        
        myRegionTable = new RegionTable(myEditRIF);        

        // Tabelle mit Daten erzeugen

        mySortedRegionTable = new TableSorter(myRegionTable);

        myRegionJTable = new JTable(mySortedRegionTable) {

            // Implement table cell tool tips.

            public String getToolTipText(MouseEvent e) {

                String tip = "no tool tip available";
                
                java.awt.Point p = e.getPoint();

                int rowIndex = rowAtPoint(p);

                // problem : columns might have been moved

                int colIndex = columnAtPoint(p);
                
                int colIndex1 = convertColumnIndexToModel(colIndex);

                // colIndex == 3 would result in problems here

                if (colIndex1 == 3) {
                    
                    /* fileId */
                    
                    Object fileObj = getValueAt(rowIndex, colIndex);

                    int fileId = ((Integer) fileObj).intValue();

                    tip =  "File: " + myEditRIF.getFileDescriptor(fileId).getShortFileName();
                    
                } else if (colIndex1 == 0) {
                    
                    tip = "unique identification of region";
                    
                } else if (colIndex1 == 1) {
                    
                    /* Region Id */

                    tip = "name of the region";
                    
                } else if (colIndex1 == 2) {
                    
                    /* ClassName  */
                    
                    tip = "name of the group/class";
                    
                } else if (colIndex1 == 4) {
                    
                    /* Lines  */

                    tip = "<line_start>:<line_stop>";
                    
                } else if (colIndex1 == 5) {
                    
                    /* Kind  */

                    tip = "kind of the region";
                    
                } else if (colIndex1 == 6) {
                    
                    /* Counting */
                    
                    Object counting = getValueAt(rowIndex, colIndex);
                    
                    boolean val = ((Boolean) counting).booleanValue();
                    
                    if (val) {
                        tip = "click box to switch off event counting";
                    } else {
                        tip = "click box to switch on event counting";
                    }
                    
                } else if (colIndex1 == 8) {
                    
                    /* Data */
                    
                    Object counting = myRegionTable.getValueAt(rowIndex, 6);
                    
                    // can be wrong:  Object Counting = getValueAt (rowIndex, 6);
                    
                    boolean val = ((Boolean) counting).booleanValue();
                    
                    if (!val) {
                        
                        tip = "Data sampling only possible if counting is enabled";
                        
                    } else {
                        
                        Object data = getValueAt(rowIndex, colIndex);
                        val = ((Boolean) data).booleanValue();
                        
                        if (val) {
                            tip = "click box to switch off data sampling";
                        } else {
                            tip = "click box to switch on data sampling";
                        }

                    }

                }

                return tip;

            } // table cell tool tips

        };

        mySortedRegionTable.setTableHeader(myRegionJTable.getTableHeader()); // needed

        // make good column sizes

        // initColumnSizes (RT_Table);

        myRegionJTable.setPreferredScrollableViewportSize(MY_DIMS);

        //Create the scroll pane and add the RT_Table to it.

        JScrollPane scrollPane = new JScrollPane(myRegionJTable);

        // Set the selection model: only one single row can be selected

        myRegionJTable.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);

        ListSelectionModel rowSM = myRegionJTable.getSelectionModel();

        rowSM.addListSelectionListener(new ListSelectionListener() {

            public void valueChanged(ListSelectionEvent e) {

                if (e.getValueIsAdjusting()) {
                    
                    return;                    
                }

                ListSelectionModel lsm = (ListSelectionModel) e.getSource();

                if (lsm.isSelectionEmpty()) {

                    logger.info("No rows are selected.");

                } else {

                    int selectedRow = lsm.getMinSelectionIndex();

                    // selected row in RT_Table is not necessary the row in RT

                    selectedRow = mySortedRegionTable.modelIndex(selectedRow);

                    FileDescriptor fileDSP = myRegionTable.getFileDescriptor(selectedRow);
                    
                    int startLine = myRegionTable.getLineStart(selectedRow);
                    
                    int stopLine = myRegionTable.getLineStop(selectedRow);

                    if (stopLine == 0) {
                        
                        stopLine = startLine;                        
                    }

                    myShowFilePanel.setFileDescriptor(fileDSP, startLine, stopLine);

                }
            }
        }); // ListSelectionListener 

        //  add the Region Table Menu at the top of the frame

        add("North", new RIFMenuBar(myEditRIF));

        myShowFilePanel = new ShowFilePanel();

        // make a split pane of Table and View Pane 

        JSplitPane centerPane = new JSplitPane(JSplitPane.VERTICAL_SPLIT, scrollPane, myShowFilePanel);

        add("Center", centerPane);

    }
}

