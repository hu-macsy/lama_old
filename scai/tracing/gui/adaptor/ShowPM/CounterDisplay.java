/*
 * CounterDisplay.java
 *
 * This class displays the two panels for the region and the user
 * counters.
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

import javax.swing.JPanel;
import javax.swing.JSplitPane;

/**
 * This class realizess a new panel that combines the two panels for the region and the user counters.
 *
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */
public class CounterDisplay extends JPanel
{

    /**
     * Title for the user table.
     */
    private static final String USER_TABLE_TITLE = "User Counters";

    /**
     * Title for the region table.
     */
    private static final String REGION_TABLE_TITLE = "Region Counters";

    /**
     * Panel for the region table.
     */
    private TableDisplay myRegionTableDisplay;

    /**
     * Panel for the user table.
     */
    private TableDisplay myUserTableDisplay;

    /**
     * Pane that splits user and counter panel.
     */
    private JSplitPane centerPane;

    /**
     * Constructor to build panel with counter and region table.
     *
     * @param theData is the performance data
     * @param main is feedback to the main GUI for selection routines
     */
    public CounterDisplay( PMData theData, ShowPM main )
    {

        super( new BorderLayout() );

        myRegionTableDisplay = new TableDisplay( new RegionTable( theData, main ), REGION_TABLE_TITLE );
        myUserTableDisplay = new TableDisplay( new UserTable( theData, main ), USER_TABLE_TITLE );

        centerPane = new JSplitPane( JSplitPane.VERTICAL_SPLIT, myRegionTableDisplay, myUserTableDisplay );

        // CenterPane.setDividerLocation (0.8);

        // CenterPane.setOneTouchExpandable (true);

        add( "North", new PMMenuBar( myRegionTableDisplay, myUserTableDisplay, main ) );
        add( "Center", centerPane );
    }

    /**
     * This routine sets an extension to the title of the user and region table.
     *
     * @param extension is appended to the normal title
     */
    void setTitleExtension( String extension )
    {

        if ( extension.length() == 0 )
        {

            myRegionTableDisplay.setTitle( REGION_TABLE_TITLE );
            myUserTableDisplay.setTitle( USER_TABLE_TITLE );

        }
        else
        {

            myRegionTableDisplay.setTitle( REGION_TABLE_TITLE + " (" + extension + ")" );
            myUserTableDisplay.setTitle( USER_TABLE_TITLE + " (" + extension + ")" );
        }

    } // setTitleExtension

    /**
     * This routine allows to focus on one of the region or user table. A
     * small value focuses the user table, a larger value makes the region
     * table more visible.
     *
     * @param divider must be a double value greater tan  0.0 and less than 1.0
     */
    public void setView( double divider )
    {

        centerPane.setDividerLocation( divider );
        repaint();
    }


} // class CounterDisplay
