/*
 * PMMenuBar.java
 *
 * Menu bar for the ShowPM main frame window.
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

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;

import adaptor.General.AboutAdaptorFrame;

/**
 * This class realizes the menu bar for ShowPM. Routines from
 * the SelectedInterface can be called from the different menu
 * items.
 *
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */
class PMMenuBar extends JMenuBar implements ActionListener
{

    /**
     * Divider factor to show user and region display.
     */
    private static final double VIEW_BOTH = 0.5;

    /**
     * Divider factor to show only user display.
     */
    private static final double VIEW_USER = 0.001;

    /**
     * Divider factor to show only region display.
     */
    private static final double VIEW_REGION = 0.999;

    /**
     * Menu item to quit the application.
     */
    private JMenuItem quitItem = new JMenuItem( "Quit" );

    /**
     * Menu item to export an RIF file.
     */
    private JMenuItem rifItem = new JMenuItem( "Export RIF" );

    /**
     * Menu item to remove the last data set of comparison.
     */
    private JMenuItem rmCmpItem = new JMenuItem( "Remove comparison data" );

    /**
     * Menu item to read another data set for comparison.
     */
    private JMenuItem readCmpItem = new JMenuItem( "Read comparison data" );

    /**
     * Menu item to select the next process whose counter values will be
     * displayed.
     */
    private JMenuItem nextProcItem = new JMenuItem( "Next process" );

    /**
     * Menu item to select the next thread whose counter values will be
     * displayed.
     */
    private JMenuItem nextThreadItem = new JMenuItem( "Next thread" );

    /**
     * Focus the view on the region table.
     */
    private JMenuItem viewRegionItem = new JMenuItem( "Region" );

    /**
     * Menu item to focus the view on the table with user counter values.
     */
    private JMenuItem viewUserItem = new JMenuItem( "User" );

    /**
     * Menu item to show both tables in equal parts.
     */
    private JMenuItem viewBothItem = new JMenuItem( "Both" );

    /**
     * Menu item to display about info for this tool.
     */
    private JMenuItem aboutItem = new JMenuItem( "About" );

    /**
     * This frame will be generated/shown when "About" is selected.
     */
    private AboutAdaptorFrame myAbout = null;

    /**
     * This is the pointer to the main GUI for which we will call
     * routines when menu items have been selected.
     */
    private ShowPM myShowPM = null;

    /**
     * Constructor for a menu bar by region and user table.
     *
     * @param regionTable is the display for the region table.
     * @param userTable is the display for the user table.
     * @param editor is the callback to the ShowPM interface
     */
    PMMenuBar( TableDisplay regionTable, TableDisplay userTable, ShowPM editor )
    {

        this.myShowPM = editor;

        // define the items of the FileMenu

        quitItem.setToolTipText( "quit the PM Result Editor" );
        quitItem.addActionListener( this );

        // define the items of the FileMneu

        rifItem.setToolTipText( "export a region information file" );
        rifItem.addActionListener( this );

        // FileMenu

        JMenu fileMenu = new JMenu( "File" );

        fileMenu.add( quitItem );
        fileMenu.add( rifItem );

        add( fileMenu );

        // ViewMenu

        viewRegionItem.addActionListener( this );
        viewUserItem.addActionListener( this );
        viewBothItem.addActionListener( this );

        JMenu viewMenu = new JMenu( "View" );

        viewMenu.add( viewRegionItem );
        viewMenu.add( viewUserItem );
        viewMenu.add( viewBothItem );

        add( viewMenu );

        // ComparisonMenu

        JMenu comparisonMenu = new JMenu( "Comparison" );

        rmCmpItem.addActionListener( this );
        readCmpItem.addActionListener( this );

        comparisonMenu.add( readCmpItem );
        comparisonMenu.add( rmCmpItem );

        add( comparisonMenu );

        // ChartMenu for Regions

        add( new ChartMenu( "Region Charts", regionTable ) );

        // ChartMenu for User Counters

        add( new ChartMenu( "User Charts", userTable ) );

        // SelectMenu

        JMenu selectMenu = new JMenu( "Select" );

        nextProcItem.addActionListener( this );
        nextThreadItem.addActionListener( this );

        selectMenu.add( nextProcItem );
        selectMenu.add( nextThreadItem );

        add( selectMenu );

        JMenu helpMenu = new JMenu( "Help" );

        aboutItem.addActionListener( this );

        helpMenu.add( aboutItem );

        add( helpMenu );

    } // PMMenuBar

    /**
     * {@inheritDoc}
     *
     * @see java.awt.event.ActionListener#actionPerformed(java.awt.event.ActionEvent)
     */
    public void actionPerformed( ActionEvent e )
    {

        Object actionSource = ( Object ) e.getSource();

        if ( !( actionSource instanceof JMenuItem ) )
        {

            return;
        }

        JMenuItem actionItem = ( JMenuItem ) actionSource;

        String cmd = actionItem.getText();

        if ( cmd.equals( "Quit" ) )
        {

            System.exit( 0 );

        }
        else if ( actionItem.equals( rifItem ) )
        {

            myShowPM.exportRIF();

        }
        else if ( actionItem.equals( nextProcItem ) )
        {

            myShowPM.nextProcess();

        }
        else if ( actionItem.equals( nextThreadItem ) )
        {

            myShowPM.nextThread();

        }
        else if ( actionItem.equals( viewRegionItem ) )
        {

            myShowPM.setView( VIEW_REGION );

        }
        else if ( actionItem.equals( viewUserItem ) )
        {

            myShowPM.setView( VIEW_USER );

        }
        else if ( actionItem.equals( viewBothItem ) )
        {

            myShowPM.setView( VIEW_BOTH );

        }
        else if ( actionItem.equals( readCmpItem ) )
        {

            myShowPM.readCompareData();

        }
        else if ( actionItem.equals( rmCmpItem ) )
        {

            myShowPM.removeCompareData();

        }
        else if ( actionItem.equals( aboutItem ) )
        {

            if ( myAbout == null )
            {

                String info = "Show Results of Performance Monitoring";

                myAbout = new AboutAdaptorFrame( "ShowPM", info );
            }

            myAbout.setVisible( true );

        }
    } // actionPerformed
}
