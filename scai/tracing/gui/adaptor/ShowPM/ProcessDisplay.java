/*
 * ProcessDisplay.java
 *
 * Display for the performance data of a process table (counter is fixed).
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

import java.awt.BorderLayout;

import javax.swing.JFrame;
import javax.swing.JMenuBar;
import javax.swing.JPanel;

/**
 * The class ProcessDisplay provides a fram for the data of a PM processor table.
 * In this table data is fixed for a region or user counter, but values for all threads
 * and/or processors are available.
 *
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */
class ProcessDisplay extends JFrame
{

    /**
     * This variable contains the performance data that will be displayed in
     * the process table.
     */
    private PMData myPMData;

    /**
     * This display will provide the panel for the process table.
     */
    private TableDisplay myTableDisplay;

    /**
     * Constructor for setting up a new process display.
     *
     * @param inPM is the performance data
     * @param main is the main GUI, needed to set up process table
     */
    ProcessDisplay( PMData inPM, SelectedInterface main )
    {

        super( "PM Processor Frame: " + inPM.getName() );

        myPMData = inPM;

        JPanel viewPane = new JPanel();

        viewPane.setLayout( new BorderLayout() );

        ProcessTable myProcessTable = new ProcessTable( inPM, main );

        myTableDisplay = new TableDisplay( myProcessTable, "<counter value>" );

        JMenuBar menuBar = new JMenuBar();

        menuBar.add( new ChartMenu( "Charts", myTableDisplay ) );

        viewPane.add( "North", menuBar );
        viewPane.add( "Center", myTableDisplay );

        setContentPane( viewPane );

        pack();

    } // constructor PMProcessDisplay

    /**
     * Selection of a new counter for the process table. This
     * routine only changes the header of the table, selection
     * itself is handled over the PM data.
     *
     * @param isUser flag is true if it is a user counter
     * @param index is the index of the counter array
     */
    public void setCounter( boolean isUser, int index )
    {

        String description = myPMData.getCounterHeader( isUser, index );

        myTableDisplay.setTitle( description );

    } // setCounter


} // class ProcessDisplay
