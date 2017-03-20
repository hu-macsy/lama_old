/*
 * ShowPM.java
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

import java.awt.BorderLayout;

import javax.swing.JFrame;
import javax.swing.JMenuBar;
import javax.swing.JPanel;

/**
 * This class is a JFrame comparing performance values of different input files,
 * e.g. from different runs.
 *
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */
class CompareDisplay extends JFrame
{

    /**
     * This is the table used for the comparison of performance data.
     */
    private CompareTable myCompareTable;

    /**
     * This is the display used for the comparison table.
     */
    private TableDisplay myTableDisplay;

    /**
     * Constructor to create a display for comparison of performance
     * data sets.
     *
     * @param dataSets is an array of dataSets
     * @param number is the number of used items in dataSets
     * @param main is pointer back to the Calltree
     */
    CompareDisplay( PMData[] dataSets, int number, ShowPM main )
    {

        super( "PM Comparison Frame" );

        JPanel viewPane = new JPanel();

        viewPane.setLayout( new BorderLayout() );

        myCompareTable = new CompareTable( dataSets, number );

        myTableDisplay = new TableDisplay( myCompareTable, "<counter>" );

        JMenuBar menuBar = new JMenuBar();
        menuBar.add( new ChartMenu( "Charts", myTableDisplay ) );

        viewPane.add( "North", menuBar );
        viewPane.add( "Center", myTableDisplay );

        setContentPane( viewPane );

        pack();

    } // constructor PMCompareDisplay

    /**
     * This routine is used to change the counter, processor or thread for
     * the comparison table.
     *
     * @param isUser is true if indexCounter stands for a user counter.
     * @param indexCounter is the index of region or user counter.
     * @param indexProc is the selected processor.
     * @param indexThread is the selected thread.
     */
    void setCounter( boolean isUser, int indexCounter, int indexProc, int indexThread )
    {

        String title;

        if ( isUser )
        {

            myCompareTable.setSelectedUC( indexCounter, indexProc, indexThread );

        }
        else
        {

            myCompareTable.setSelectedRC( indexCounter, indexProc, indexThread );
        }

        title = myCompareTable.getCounterDescription();
        title = title + " (IP = " + indexProc + ", IT = " + indexThread + ")";

        myTableDisplay.setTitle( title );

    } // setCounter

} // class PMCompareDisplay
