/*
 * EditPM.java
 *
 * Main class for EditPM: GUI to configure performance measurement.
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

package adaptor.EditPM;

import java.io.BufferedReader;
import java.io.EOFException;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;

import java.awt.event.WindowEvent;
import java.awt.event.WindowListener;

import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JOptionPane;

import org.apache.log4j.Logger;

import adaptor.General.OptionsAdaptor;

/**
 * The class EditPM implements an editor for configuration files to specify
 * performance events.
 *
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */

/**
 * TODO: brandes: Enter comment!
 *
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */
public class EditPM implements WindowListener
{

    /**
     * This string is used for save dialogs.
     */
    private static final String SAVE_STRING = "save";

    /**
     * This is the logger variable for this call to use log4j.
     */
    private static Logger logger = Logger.getLogger( EditPM.class );

    /**
     * This is the main frame for this GUI and gets the WindowListener.
     * It will be allocated with this class constructor.
     */
    private JFrame myMainFrame = null;

    /**
     * This variable contains all the counter that will be defined
     * for this PM configuration. A counter counts basic, composed, or
     * derived events and has some other properties for presentations.
     */
    private CounterTable      myCounterTable = new CounterTable();

    /**
     * This variable contains the frame for the counter table. It allows
     * to edit the different counters, e.g. delete and copy of counters,
     * change the properties of the counters.
     */
    private CounterTablePanel myCTFrame;

    /**
     * This frame is used to define new composed performance events.
     */
    private ComposedFrame myComposedFrame = null;

    /**
     * This frame is used to define new derived performance events.
     */
    private DerivedFrame myDerivedFrame = null;

    /**
     * The PM configuration contains all information about supported
     * events, but not the counters itself.
     */
    private PMConfiguration myPMConfiguration = new PMConfiguration();

    /**
     * Constructor for a Performance Monitoring Editor that reads
     * supported events from the executable, reads configuration files
     * and writes an output configuration file.
     *
     * @param configFiles are the name of the configuration files
     * @param executable is the name of the executable (can be null)
     * @param outfile is the name of the output configuration file
     */
    public EditPM( String [] configFiles, String executable, String outfile )
    {

        if ( executable != null )
        {

            // run the executable with help flag to find out supported events

            parseExecutable( executable );

            // set the name of the executable in the configuration

            myPMConfiguration.setExecutable( executable );

        }

        // now we read all configuration files

        for ( int i = 0; i < configFiles.length; i++ )
        {

            myPMConfiguration.readConfigFile( configFiles[i], myCounterTable );

        }

        if ( outfile == null )
        {

            // no output file, so we make a new GUI

            makeGUI();

            logger.info( "EditPM terminates, GUI thread remains running" );

        }
        else
        {

            // write configuration directly in the new output file

            myPMConfiguration.writeConfigFile( outfile, myCounterTable );

            logger.info( "output file " + outfile + " written, main terminates" );

        }

    }


    /**
     * This routine adds a new counter for a given performance event.
     * The counter will be added to the counter table, the frame for
     * the table will be updated. It is called as an action for the
     * button belonging to the performance event.
     *
     * @param newEvent is the event for which the counter will be defined.
     */
    private void addCounter( PerformanceEvent newEvent )
    {

        logger.info( "added counter for " + newEvent.getName() );

        myCounterTable.addCounter( newEvent );

        // actualize the frame for the counter table

        myCTFrame.update();

    } // addCounter


    /**
     * This routine is called when an event button has been pressed. The selected
     * event will be set in the composed or derived frame if open. Otherwise a new
     * counter for this event will be added to the counter table.
     *
     * @param pressedPerformanceEvent is the event for which the event button has been pressed
     */
    public void actionEvent( PerformanceEvent pressedPerformanceEvent )
    {

        // when the ComposedEventFrame is visible select there the event

        if ( myComposedFrame.isShowing() )
        {

            if ( pressedPerformanceEvent instanceof BasicEvent )
            {

                myComposedFrame.setEvent( ( BasicEvent ) pressedPerformanceEvent );

            }

            if ( pressedPerformanceEvent instanceof ComposedEvent )
            {
                myComposedFrame.editEvent( ( ComposedEvent ) pressedPerformanceEvent );
            }

        }
        else if ( myDerivedFrame.isShowing() )
        {

            // when the DerivedEventFrame is visible select there the event

            if ( pressedPerformanceEvent instanceof BasicEvent )
            {
                myDerivedFrame.setEvent( ( BasicEvent ) pressedPerformanceEvent );
            }

            if ( pressedPerformanceEvent instanceof DerivedEvent )
            {
                myDerivedFrame.editEvent( ( DerivedEvent ) pressedPerformanceEvent );
            }
        }
        else
        {

            // none of the frames for composed / derived counters is open

            addCounter( pressedPerformanceEvent );

        }

    }

    /**
     * Show the frame with buttons for all supported performance events.
     *
     */
    public void showEventFrame()
    {

        // the event frame is available via the PM configuration
        // but the frame needs my GUI for button actions on the frame

        myPMConfiguration.showEventFrame( this );
    }


    /**
     * Show the frame to define a composed performance event (metric).
     *
     */
    public void showComposedFrame()
    {

        myDerivedFrame.setVisible( false );
        myComposedFrame.setVisible( true );

    }

    /**
     * Show the frame to define a derived performance event (metric).
     *
     */
    public void showDerivedFrame()
    {

        myComposedFrame.setVisible( false );
        myDerivedFrame.setVisible( true );

    }

    /**
     * This routine pos up an open dialog t select a file for reading
     * or writing.
     *
     * @param windowTitle is the title of the dialog wind
     * @param openFlag is true for open a file, otherwise it will be written
     * @return the selected file name (null for cancel)
     */
    private String getFileName( String windowTitle, boolean openFlag )
    {

        String fileName = null;

        String thePWD = System.getProperty( "user.dir" );
        JFileChooser chooser = new JFileChooser( thePWD );
        int result;

        chooser.setDialogTitle( "Select config file to " + windowTitle + ":" );

        if ( openFlag )
        {
            result = chooser.showOpenDialog( myMainFrame );
        }
        else
        {
            result = chooser.showSaveDialog( myMainFrame );
        }

        if ( result == JFileChooser.CANCEL_OPTION )
        {
            return fileName;
        }

        File selectedFile = chooser.getSelectedFile();

        if ( selectedFile == null )
        {

            return fileName;
        }

        fileName = selectedFile.getPath();

        // check whether PWD appears

        if ( fileName.indexOf( thePWD ) >= 0 )
        {

            fileName = fileName.substring( thePWD.length() + 1 );

        }

        return fileName;
    }

    /***********************************************************************************************
     * * clear: make original configuration * *
     **********************************************************************************************/

    /**
     * This routine is called to reset the counter table and to set default
     * values for all extra features.
     *
     */
    public void clear()
    {

        // remove all counters of Counter Table

        myCTFrame.clear();

        // reset all extra values

        myPMConfiguration.clear();

    }

    /**
     * This routine will be called to save the configuration in a file. If
     * we have a file name in the configuration, we take it.
     *
     */
    public void save()
    {

        String configFile = myPMConfiguration.getConfigurationFile();

        if ( configFile == null )
        {

            // no filename available so we call saveAS

            saveAs();

        }
        else
        {

            myPMConfiguration.writeConfigFile( configFile, myCounterTable );
        }

    }

    /**
     * This routine will be called to save the configuration in a
     * new file whose name has to be specified.
     *
     */
    public void saveAs()
    {

        String fileName = getFileName( SAVE_STRING, false );

        if ( fileName == null )
        {

            return;
        }

        myPMConfiguration.writeConfigFile( fileName, myCounterTable );

    }

    /**
     * This routine is called to read a new configuration file. The
     * name of the file will be given by a dialow.
     *
     */
    public void load()
    {

        String configFile = getFileName( SAVE_STRING, false );

        if ( configFile == null )
        {

            return;
        }

        clear();

        myPMConfiguration.readConfigFile( configFile, myCounterTable );

    }

    /**
     * This routine is called to read a configuration file and to
     * add the performance events and counters to the existing features.
     * Certain value might be overwritten.
     *
     */
    public void add()
    {

        String fileName = getFileName( SAVE_STRING, false );

        if ( fileName == null )
        {

            return;
        }

        // Note: we do not clear the counter table

        myPMConfiguration.readConfigFile( fileName, myCounterTable );

    }

    /**
     *
     * This routine will be called to quit the PM editor. A dialog might
     * appear to save the last configuration.
     *
     * TODO: keep track of changes to get a useful value for the hasChanged.
     *
     */
    public void quit()
    {

        boolean hasChanged = true;

        if ( hasChanged )
        {

            String configFile = myPMConfiguration.getConfigurationFile();

            String title = "Save PM Configuration";

            String msg = "Would you like to save the configuration";

            if ( configFile == null )
            {

                msg += "?";

            }
            else
            {

                msg += " " + configFile + "?";
            }

            int response = JOptionPane.showConfirmDialog( null, msg, title, JOptionPane.YES_NO_OPTION, JOptionPane.ERROR_MESSAGE );

            if ( response == 0 )
            {

                // so we have to da SAVE

                if ( configFile == null )
                {

                    saveAs();

                }
                else
                {

                    save();

                }

            }
        }

        System.exit( 0 );

    }

    /**
     * This routine will start the instrumented executable to find out
     * which performance counters are supported.
     *
     * @param executable is the name of the ADAPTOR instrumented executable
     */
    void parseExecutable( String executable )
    {

        // would be helpful to test whether it is really an executable

        String [] cmd = { executable, "-pm", "-help" };

        try
        {

            Runtime theRuntime = Runtime.getRuntime();
            Process theProcess = theRuntime.exec( cmd );

            logger.info( theProcess + " has started" );

            try
            {

                // myPMConfiguration.init();

                boolean eof = false;

                BufferedReader in = new BufferedReader( new InputStreamReader( theProcess.getInputStream() ) );

                while ( !eof )
                {

                    String inputLine = in.readLine();

                    if ( inputLine == null )
                    {

                        eof = true;

                    }
                    else
                    {

                        myPMConfiguration.parseExecutableLine( inputLine );

                    }

                } // while not eof

            }
            catch ( EOFException e )
            {
                logger.info( "EOF reached (piped input from " + executable + ")" );
            }

        }
        catch ( IOException e )
        {

            logger.error( "Error: " + e.getMessage(), e );
        }

        int noEvents = myPMConfiguration.getNumberEvents();

        if ( noEvents == 0 )
        {

            String command = "";

            for ( int i = 0; i < cmd.length; i++ )
            {

                command += " " + cmd[i];
            }

            logger.error( "No events found for this command: " + command );
            System.exit( 0 );

        }
        else
        {

            logger.info( "Number of known events: " + noEvents );
        }
    }

    /**
     * This routine builds up all elements of the graphical user interface.
     *
     */
    private void makeGUI()
    {

        //Schedule a job for the event-dispatching thread:
        //creating and showing this application's GUI.

        JFrame.setDefaultLookAndFeelDecorated( true );

        // create and set up the window

        myMainFrame = new JFrame( "Performance Monitoring - Configuration" );

        myMainFrame.addWindowListener( this );

        myCTFrame = new CounterTablePanel( myCounterTable );

        myCTFrame.add( "North", new PMMenuBar( this ) );

        myCTFrame.add( "South", myPMConfiguration.getPanel() );

        myCTFrame.setOpaque( true );

        // Frame for composed events (do not show now)

        myComposedFrame = new ComposedFrame( myPMConfiguration );

        // Frame for derived events (do not show now)

        myDerivedFrame = new DerivedFrame( myPMConfiguration );

        myMainFrame.setContentPane( myCTFrame );
        myMainFrame.pack();
        myMainFrame.setVisible( true );

    }

    /**
     * This class contains a main routine so that it can be called.
     *
     * @param args These are the command line arguments
     */

    public static void main( String[] args )
    {

        final int errorExitCode = -1;

        String helpText = "Editor for Performance Monitoring";

        // EditPM [executable [ configFile]]

        OptionsAdaptor theOptions = new OptionsAdaptor( "EditPM", helpText );

        theOptions.addValueOption( "-e", "executable (to ask for supported events" );
        theOptions.addValueOption( "-o", "outputfile (no GUI)" );

        String[] configFiles = null;

        try
        {

            configFiles = theOptions.parseArguments( args );

        }
        catch ( IllegalArgumentException e )
        {

            logger.error( e.getMessage() );
            theOptions.printUsage();
            System.exit( errorExitCode );

        }

        String executable    = theOptions.getOptionVal( "-e" );
        String outfile       = theOptions.getOptionVal( "-o" );

        logger.info( "main program started (command line arguments have been evaluated)" );

        // we create an instance of the editor

        new EditPM( configFiles, executable, outfile );

        logger.info( "main thread terminates" );

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
        quit();

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


    /**
     * TODO: brandes: Enter comment!
     *
     */
    public void deleteCounter()
    {

        myCTFrame.deleteCounter();

    }


    /**
     * TODO: brandes: Enter comment!
     *
     */
    public void copyCounter()
    {

        myCTFrame.copyCounter();

    }


    /**
     * TODO: brandes: Enter comment!
     *
     */
    public void pasteCounter()
    {

        myCTFrame.pasteCounter();

    }


    /**
     * This method deletes all counters from the counter table.
     *
     */
    public void clearCounters()
    {

        myCTFrame.clear();

    }
}
