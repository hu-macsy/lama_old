/*
 * EditRIF.java
 *
 * Main class for EditRIF: GUI to edit RIF files.
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


import java.awt.event.WindowEvent;
import java.awt.event.WindowListener;

import java.io.File;

import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JOptionPane;

import org.apache.log4j.Logger;

import adaptor.General.FileDescriptor;
import adaptor.General.OptionsAdaptor;
import adaptor.General.RIFData;
import adaptor.General.RegionDescriptor;

/**
 * The class EditRIF implements an editor for RIF files used by the ADAPTOR tool environment.
 *
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */

public class EditRIF implements EditRIFInterface, WindowListener
{

    /**
     * The logger variable for this class.
     */
    private static Logger logger = Logger.getLogger( EditRIF.class );

    /**
     * The data of the RIF file for which this editor stands.
     */
    private RIFData theRIFData = null;

    /**
     * indicates whether RIF data has been modified.
     */
    private boolean rifChanged = false;

    /**
     * The main frame of the RIF Editor.
     */
    private JFrame myEditorFrame;

    /**
     * The display of the region table.
     */
    private RegionTableDisplay myRegionDisplay; // accessed by the RIF Menu

    /**
     * Constructor for the RIF Editor.
     *
     * @param filename is the name of the RIF file
     */
    public EditRIF( String filename )
    {

        theRIFData = new RIFData( filename );

        // Schedule a job for the event-dispatching thread:
        // creating and showing this application's GUI.

        JFrame.setDefaultLookAndFeelDecorated( false );

        // create and set up the window

        myEditorFrame = new JFrame( "Region Information File: " + filename );

        // might be better: call also QuitRIF in this case

        myEditorFrame.addWindowListener( this );

        // frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

        myRegionDisplay = new RegionTableDisplay( this );

        logger.info( "created region table display." );

        myRegionDisplay.setOpaque( true );
        myEditorFrame.setContentPane( myRegionDisplay );
        myEditorFrame.pack();
        myEditorFrame.setVisible( true );
    }

    /**
     * {@inheritDoc}
     *
     * @see adaptor.EditRIF.EditRIFInterface#quitRIF()
     */
    public void quitRIF()
    {

        // make sure that a modified RIF file is saved

        if ( rifChanged )
        {

            String title = "Save Region Information File";

            String msg = "Do you want to save file " + theRIFData.getFilename() + "?";

            int response = JOptionPane.showConfirmDialog( null, msg, title, JOptionPane.YES_NO_OPTION, JOptionPane.ERROR_MESSAGE );

            if ( response == 0 )
            {

                theRIFData.writeRIF();
            }

        }

        System.exit( 0 );
    }

    /**
     * {@inheritDoc}
     *
     * @see adaptor.EditRIF.EditRIFInterface#getFileDescriptor(int)
     */
    public FileDescriptor getFileDescriptor( int fileId )
    {

        return theRIFData.getFile( fileId );

    }

    /**
     * {@inheritDoc}
     *
     * @see adaptor.EditRIF.EditRIFInterface#getRegions()
     */
    public RegionDescriptor getRegion( int index )
    {

        return theRIFData.getRegion( index );

    }

    /**
     * {@inheritDoc}
     *
     * @see adaptor.EditRIF.EditRIFInterface#noRegions()
     */
    public int noRegions()
    {

        return theRIFData.noRegions();
    }

    /**
     * Static main routine that can be called.
     *
     * @param arguments are the command line arguments
     */
    public static void main( String[] arguments )
    {

        final int errorExitCode = -1;

        // describe this program for parsing arguments

        OptionsAdaptor theOptions = new OptionsAdaptor( "EditRIF", "Editor for Region Information File" );

        // parse the arguments

        String[] inputFiles = null;

        try
        {

            inputFiles = theOptions.parseArguments( arguments );

        }
        catch ( IllegalArgumentException e )
        {

            logger.error( e.getMessage() );
            theOptions.printUsage();
            System.exit( errorExitCode );

        }

        logger.info( "main program started (command line arguments have been evaluated)" );

        // we create an instance of the editor

        if ( inputFiles.length == 0 )
        {

            // take the default RIF file

            new EditRIF( "RIF" );

        }
        else
        {

            for ( int i = 0; i < inputFiles.length; i++ )
            {

                new EditRIF( inputFiles[i] );

            }
        }

        logger.info( "main terminates, GUI thread for EditRIF remains" );

    } // main

    /**
     * {@inheritDoc}
     *
     * @see adaptor.EditRIF.EditRIFInterface#save()
     */
    public void save()
    {

        theRIFData.writeRIF();

    }

    /**
     * {@inheritDoc}
     *
     * @see adaptor.EditRIF.EditRIFInterface#saveAs()
     */
    public void saveAs()
    {

        String currentDir = System.getProperty( "user.dir" );

        JFileChooser chooser = new JFileChooser( currentDir );

        chooser.setDialogTitle( "Select file to save:" );

        int result = chooser.showSaveDialog( myEditorFrame );

        if ( result == JFileChooser.CANCEL_OPTION )
        {

            return;
        }

        File selectedFile = chooser.getSelectedFile();

        // make sure that that really a file has been selected

        if ( selectedFile == null )
        {

            return;
        }

        theRIFData.writeRIF( selectedFile.getPath() );

    }

    /**
     * {@inheritDoc}
     *
     * @see adaptor.EditRIF.EditRIFInterface#setModified()
     */
    public void setModified()
    {

        rifChanged = true;

        // repaint also the RegionTableDisplay

        myRegionDisplay.repaint();

    }

    /**
     * {@inheritDoc}
     *
     * @see java.awt.event.WindowListener#windowOpened(java.awt.event.WindowEvent)
     */
    public void windowOpened( WindowEvent arg0 )
    {

        // logger.info(arg0);

    }

    /**
     * {@inheritDoc}
     *
     * @see java.awt.event.WindowListener#windowClosing(java.awt.event.WindowEvent)
     */
    public void windowClosing( WindowEvent arg0 )
    {

        // logger.info(arg0);

        quitRIF();

    }

    /**
     * {@inheritDoc}
     *
     * @see java.awt.event.WindowListener#windowClosed(java.awt.event.WindowEvent)
     */
    public void windowClosed( WindowEvent arg0 )
    {

        // logger.info(arg0);

    }

    /**
     * {@inheritDoc}
     *
     * @see java.awt.event.WindowListener#windowIconified(java.awt.event.WindowEvent)
     */
    public void windowIconified( WindowEvent arg0 )
    {

        // logger.info(arg0);

    }

    /**
     * {@inheritDoc}
     *
     * @see java.awt.event.WindowListener#windowDeiconified(java.awt.event.WindowEvent)
     */
    public void windowDeiconified( WindowEvent arg0 )
    {

        // logger.info(arg0);

    }

    /**
     * {@inheritDoc}
     *
     * @see java.awt.event.WindowListener#windowActivated(java.awt.event.WindowEvent)
     */
    public void windowActivated( WindowEvent arg0 )
    {

        // logger.info(arg0);

    }

    /**
     * {@inheritDoc}
     *
     * @see java.awt.event.WindowListener#windowDeactivated(java.awt.event.WindowEvent)
     */
    public void windowDeactivated( WindowEvent arg0 )
    {

        // logger.info(arg0);

    }

} // class EditRIF

