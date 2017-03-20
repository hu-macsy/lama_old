/*
 * RIFMenuBar.java
 *
 * Menu bar class for the RIF editor.
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

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;

import org.apache.log4j.Logger;

import adaptor.General.AboutAdaptorFrame;
import adaptor.General.RegionDescriptor;

/**
 * The class RIFMenuBar provides the menues for the RIF Editor.
 *
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */
class RIFMenuBar extends JMenuBar implements ActionListener
{

    /**
     * Logger for class RIFMenuBar.
     */
    private static Logger logger = Logger.getLogger( RIFMenuBar.class );

    /**
     * Pointer to the RIF editor to call its routines.
     */
    private EditRIFInterface myEditRIF;

    // these are the menu items for which actions are specified

    /**
     * Menu item to quit the application.
     */
    private JMenuItem quitItem = new JMenuItem( "Quit" );

    /**
     * Menu item to save the current RIF file.
     */
    private JMenuItem saveItem = new JMenuItem( "Save" );

    /**
     * Menu item to save the current RIF data to a new file.
     */
    private JMenuItem saveAsItem = new JMenuItem( "Save As ..." );

    /**
     * Menu item to get more information about this application.
     */
    private JMenuItem aboutItem = new JMenuItem( "About" );

    /**
     * The AboutFrame to be shown when selected.
     */
    private AboutAdaptorFrame myAbout = null;

    /***********************************************************************************************
     * * constructor for RIFMenuBar * *
     **********************************************************************************************/

    /**
     * Constructor for the RIF menu bar.
     *
     * @param editor is the editor to which the menu will be anchored
     */
    public RIFMenuBar( EditRIFInterface editor )
    {

        // define the items of the FileMneu

        myEditRIF = editor;

        quitItem.setToolTipText( "quit the RIF Editor (without save)" );
        quitItem.addActionListener( this );

        saveItem.setToolTipText( "save the RIF file with all changes" );
        saveItem.addActionListener( this );

        saveAsItem.setToolTipText( "save the RIF file" );
        saveAsItem.addActionListener( this );

        JMenu fileMenu = new JMenu( "File" );

        fileMenu.add( saveItem );
        fileMenu.add( saveAsItem );
        fileMenu.add( quitItem );

        // Help Menu

        aboutItem.addActionListener( this );

        JMenu helpMenu = new JMenu( "Help" );
        helpMenu.add( aboutItem );

        JMenu enableMenu = new JMenu( "Enable" );
        enableMenu.setToolTipText( "enable event counting for certain regions" );

        JMenu disableMenu = new JMenu( "Disable" );
        disableMenu.setToolTipText( "disable event counting for certain regions" );

        JMenuItem allOnItem = new JMenuItem( "Enable all" );
        JMenuItem allOffItem = new JMenuItem( "Disable all" );

        enableMenu.add( allOnItem );
        disableMenu.add( allOffItem );
        helpMenu.add( aboutItem );

        // for each region kind we add a menu item to enable or disable this kind of region

        for ( int k = 0; k < RegionDescriptor.RegionKind.values().length; k++ )
        {

            String kindString = RegionDescriptor.RegionKind.values()[k].toString();

            JMenuItem kindItem = new JMenuItem( "Enable " + kindString );
            kindItem.addActionListener( this );
            enableMenu.add( kindItem );

            kindItem = new JMenuItem( "Disable " + kindString );
            kindItem.addActionListener( this );
            disableMenu.add( kindItem );

        }

        allOnItem.addActionListener( this );
        allOffItem.addActionListener( this );

        // make some mnemonics

        fileMenu.setMnemonic( 'F' );
        enableMenu.setMnemonic( 'E' );
        disableMenu.setMnemonic( 'D' );
        helpMenu.setMnemonic( 'H' );

        // now add all menus to the menu bar

        add( fileMenu );
        add( enableMenu );
        add( disableMenu );
        add( helpMenu );

    } // RIFMenuBar

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

        String cmd = ( ( JMenuItem ) actionSource ).getText();

        logger.info( "perform action for menu item " + cmd );

        if ( actionSource == quitItem )
        {
            myEditRIF.quitRIF();
        }

        if ( actionSource == saveItem )
        {
            myEditRIF.save();
        }

        if ( actionSource == saveAsItem )
        {
            myEditRIF.saveAs();
        }

        if ( actionSource == aboutItem )
        {

            if ( myAbout == null )
            {

                String info = "Editor for Region Information File";

                myAbout = new AboutAdaptorFrame( "EditRIF", info );
            }

            myAbout.setVisible( true );
        }

        Boolean val = null;

        if ( cmd.startsWith( "Enable " ) )
        {
            val = Boolean.TRUE;
        }
        else if ( cmd.startsWith( "Disable " ) )
        {
            val = Boolean.FALSE;
        }

        if ( val != null )
        {

            // found Enable or Disable

            String subCommand = cmd.substring( cmd.indexOf( " " ) + 1 );

            // logger.info(cmd + " sub=" + sub_cmd + " val=" + val);

            int noRegions = myEditRIF.noRegions();

            for ( int i = 0; i < noRegions; i++ )
            {

                RegionDescriptor theRegion = myEditRIF.getRegion( i );

                if ( subCommand.equals( "all" ) || subCommand.equals( theRegion.getRegionKind() ) )
                {

                    theRegion.setRegionObject( val, 6 );
                }

            }

            myEditRIF.setModified();

        }

    } // actionPerformed
}
