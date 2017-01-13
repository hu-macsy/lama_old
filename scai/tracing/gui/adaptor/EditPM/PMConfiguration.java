/*
 * PMConfiguration.java
 *
 * Class contains all configuration data for PM measurement
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
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import javax.swing.JOptionPane;
import javax.swing.JPanel;

import org.apache.log4j.Logger;

import adaptor.General.ArraySelection;

/**
 * The class PMConfiguration contains all data needed to
 * describe the current configuration of a Performance
 * Monitoring (PM) experiment.
 *
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */

/**
 * TODO: adaptor: Enter comment!
 *
 * @version $LastChangedRevision$
 * @author adaptor
 */
public class PMConfiguration
{

    /**
     * This string identifies comments in the configuration file.
     */
    private static final String COMMENT_SIGN = "#";

    /**
     * This string is used at the begin of each comment line.
     */
    private static final String COMMENT_START = COMMENT_SIGN + " ";

    /**
     * Logger variable for this class.
     */
    private static Logger logger = Logger.getLogger( PMConfiguration.class );

    /**
     * This is the default kind ("CLASS") of how regions are grouped.
     */
    private static final int GROUP_ITEM_DEFAULT = 2;

    /**
     * These are all possibilities of how regions can grouped.
     * ALL groups all regions to one group. FILE groups regions
     * by the file in which they apper. CLASS groups regions by
     * identical class names. SINGLE make every region to a single
     * group.
     */
    private static final String[] GROUP_ITEMS = { "ALL", "FILE", "CLASS", "SINGLE" };

    /**
     * This is the number of supported performance events.
     */
    private int noPerformanceEvents = 0;


    /**
     * This array contains all supported performance events. It may be allocated larger than number of supported events.
     */
    private PerformanceEvent[] myPerformanceEvents;

    /**
     * Panel needed for all values that are not related to performance counter.
     */
    private ExtraPanel myExtraPanel = null;

    /**
     * Frame that contains buttons for each performance event.
     */
    private EventFrame myEventFrame = null;

    /**
     * This is the name of the executable from which the performance data
     * will be collected.
     */
    private String myExecutable = null;

    /**
     * This is the name of the configuration file.
     */
    private String myConfigFile = null;

    /**
     * This is the maximal profiling depth where depth < 0 stands
     * for unlimited profiling.
     */
    private int profilingDepth = -1;

    /**
     * Flag that indicates if data profiling is enabled.
     */
    private boolean dataProfiling = false;

    /**
     * Value for the sampling rate in case of data profiling.
     */
    private int rateSampling = 0;

    /**
     * Maximal number of samples in case of data profiling.
     */
    private int sizeSampling = 0;

    /**
     * Flag that indiciates if calltree data is collected at
     * runtime.
     */
    private boolean doCalltree = false;

    /**
     * Flag that indicates if calltree data for data profiling
     * is collected at runtime.
     */
    private boolean doCalldata = false;

    /**
     * Flag that indicates if tracing is enabled.
     */
    private boolean doTracing  = false;

    /**
     * The current selection of how regions are grouped.
     */
    private int groupIndex = GROUP_ITEM_DEFAULT;

    /**
     * Name of the RIF file that is used for the regions.
     */
    private String rifFile = null;

    /**
     * Name of the output file for performance data.
     */
    private String outFile = null;

    /**
     * Searching for a performance event with a given name.
     *
     * @param name is the name of an event to be searched
     * @return the PerformanceEvent with the given name
     */
    public PerformanceEvent findSupportedEvent( String name )
    {

        PerformanceEvent foundEvent = null;

        for ( int i = 0; i < noPerformanceEvents; i++ )
        {

            PerformanceEvent event = myPerformanceEvents[i];

            if ( event.hasName( name ) )
            {

                foundEvent = event;
                break; // stops further search
            }
        }

        // we will return an event != null if found

        return foundEvent;
    }

    /**
     * This routine returns the index for a given performance event.
     *
     * @param event is the performance event
     * @return index of event in my event array
     */
    public int getEventIndex( PerformanceEvent event )
    {


        for ( int i = 0; i < noPerformanceEvents; i++ )
        {


            if ( event == myPerformanceEvents[i] )
            {

                return i;
            }
        }

        return ArraySelection.NO_SELECTION;

    } // getEventIndex

    /**
     * This routine adds a new performance event to my known
     * events.
     *
     * @param event is the new event.
     */
    private void addEvent( PerformanceEvent event )
    {

        if ( myPerformanceEvents == null )
        {

            final int defaultSize = 64;

            myPerformanceEvents = new PerformanceEvent[defaultSize];
            noPerformanceEvents = 0;

        }

        int len = myPerformanceEvents.length;

        // make sure that array is not too small

        if ( noPerformanceEvents >= len )
        {

            final int addedSize = 32;

            PerformanceEvent[] helpArray;

            helpArray = new PerformanceEvent[noPerformanceEvents + addedSize];

            for ( int i = 0; i < myPerformanceEvents.length; i++ )
            {

                helpArray[i] = myPerformanceEvents[i];
            }

            for ( int i = 0; i < addedSize; i++ )
            {

                helpArray[noPerformanceEvents + i] = null;
            }

            myPerformanceEvents = helpArray;
        }

        myPerformanceEvents[noPerformanceEvents++] = event;

        // add the event to the event frame if it already exsits

    }

    /**
     * This routine adds a performance event to the array of events.
     *
     * @param event is the new event
     */
    public void addNewEvent( PerformanceEvent event )
    {

        String name = event.getName();

        name = name.toUpperCase();
        name = name.trim();

        if ( name.length() == 0 )
        {

            throw new RuntimeException( "no name specified for new event" );
        }

        if ( name.indexOf( " " ) >= 0 )
        {

            throw new RuntimeException( "new event name must not contain blanks" );
        }


        if ( findSupportedEvent( name ) != null )
        {

            throw new RuntimeException( name + " is already defined event" );
        }

        // add the new Event to the Event Table

        addEvent( event );

        // add the new Event to the Event Frame

        if ( myEventFrame != null )
        {

            myEventFrame.addEvent( event );
        }

    } // addNewEvent

    /**
     * Asking for the total number of known events in the configuration.
     *
     * @return number of known events
     */
    public int getNumberEvents()
    {

        return noPerformanceEvents;

    }

    /**
     * Show the frame with all events.
     *
     * @param editor is needed for feedback
     */
    public void showEventFrame( EditPM editor )
    {

        if ( myEventFrame == null )
        {

            // we have to generate a new frame for the events

            myEventFrame = new EventFrame( editor, myPerformanceEvents );

        }

        myEventFrame.pack();
        myEventFrame.setVisible( true );

    }

    /**
     * This routine looks for the basic event with the given name.
     *
     * @param name is the name of the event
     * @return the basic event
     */
    public PerformanceEvent getBasicEvent( String name )
    {

        PerformanceEvent event = findSupportedEvent( name );

        if ( event == null )
        {

            throw new RuntimeException( name + " is unknown event" );

        }
        else  if ( event instanceof BasicEvent )
        {

            return event;
        }

        throw new RuntimeException( name + " is not basic event" );

    } // GetBasicEvent

    /**
     * This routine returns an array with all strings of the basic events.
     *
     * @return array with name strings
     */
    public String[] getBasicEventStrings()
    {

        String[] choices = new String[noPerformanceEvents];

        for ( int i = 0; i < noPerformanceEvents; i++ )
        {

            choices[i] = myPerformanceEvents[i].getName();
        }

        return choices;
    }

    /**
     * This routine is used to write the counter table in a configuration file.
     *
     * @param fileName is the name of the output file
     * @param counterTable is the counter table to write
     */
    public void writeConfigFile( String fileName, CounterTable counterTable )
    {

        try
        {

            int i;
            int numberComposedEvents = 0;
            int numberDerivedEvents = 0;

            PerformanceEvent event;

            FileWriter file = new FileWriter( fileName );
            BufferedWriter buff = new BufferedWriter( file );

            // make a nice header line to describe this file

            buff.write( COMMENT_START + "performance monitoring configuration file" );
            buff.newLine();
            buff.write( COMMENT_START + "generated by EditPM" );
            buff.newLine();
            buff.newLine();
            buff.write( COMMENT_START + "these are all supported basic performance events" );
            buff.newLine();

            for ( i = 0; i < noPerformanceEvents; i++ )
            {

                event = myPerformanceEvents[i];

                if ( event instanceof ComposedEvent )
                {
                    numberComposedEvents++;
                }
                else if ( event instanceof DerivedEvent )
                {
                    numberDerivedEvents++;
                }
                else if ( event instanceof BasicEvent )
                {
                    buff.write( COMMENT_START + event.getName() + " : " + event.getDescription() );
                    buff.newLine();
                }

            }

            buff.newLine();

            if ( numberComposedEvents > 0 )
            {
                buff.write( COMMENT_START + numberComposedEvents + " composed events" );
                buff.newLine();

                for ( i = 0; i < noPerformanceEvents; i++ )
                {
                    event = myPerformanceEvents[i];

                    if ( event instanceof ComposedEvent )
                    {
                        buff.write( event.getName() + " = " + event.getDescription() );
                        buff.newLine();
                    } // ComposedEvent
                } // for

                buff.newLine();
            } // writing all composed events

            if ( numberDerivedEvents > 0 )
            {
                buff.write( COMMENT_START + numberDerivedEvents + " derived events" );
                buff.newLine();

                for ( i = 0; i < noPerformanceEvents; i++ )
                {
                    event = myPerformanceEvents[i];

                    if ( event instanceof DerivedEvent )
                    {
                        buff.write( event.getName() + " = " + event.getDescription() );
                        buff.newLine();
                    } // DerivedEvent
                } // for

                buff.newLine();
            } // writing all derived events

            // now we write the table of counters

            for ( i = 0; i < counterTable.getRowCount(); i++ )
            {

                Counter c = counterTable.getCounter( i );
                buff.write( COMMENT_START + c.getDescription() );
                buff.newLine();
                buff.write( c.makeLine() );
                buff.newLine();
                buff.newLine();

            } // for loop

            writeConfigExtras( buff );

            buff.close();

        }
        catch ( IOException e )
        {

            String dialogTitle = "File Output Error";
            String dialogMsg = e.getMessage();

            JOptionPane.showMessageDialog( null, dialogMsg, dialogTitle, JOptionPane.ERROR_MESSAGE );

        } // catch

        // we take the file name as the default name for the configuration file

        myConfigFile = fileName;

        if ( myExtraPanel != null )
        {

            myExtraPanel.showValues();

        }

    } // WritePMFile

    /**
     * Additional routine to write all configuration data except the
     * counters in a buffered output file.
     *
     * @param buff is the output file
     * @throws IOException in case of any IO problem
     */
    private void writeConfigExtras( BufferedWriter buff ) throws IOException
    {

        // now we write the extra items

        buff.write( COMMENT_START + "additional entries" );
        buff.newLine();

        if ( profilingDepth < 0 )
        {

            buff.write( COMMENT_START + "DEPTH not set for profiling" );

        }
        else
        {

            buff.write( "DEPTH " + profilingDepth );

        }

        buff.newLine();

        if ( !dataProfiling )
        {

            buff.write( COMMENT_START + "no data profiling specified" );

        }
        else
        {

            buff.write( "SAMPLING " + rateSampling + " " + sizeSampling );
        }

        buff.newLine();

        if ( !doCalltree )
        {

            buff.write( COMMENT_START + "calltree not enabled here" );

        }
        else
        {

            buff.write( "CALLTREE on" );
        }

        buff.newLine();

        if ( !doCalldata )
        {

            buff.write( COMMENT_START + "calltree not enabled here" );

        }
        else
        {

            buff.write( "CALLDATA on" );
        }

        buff.newLine();

        if ( !doTracing )
        {

            buff.write( COMMENT_START + "trace not enabled here" );

        }
        else
        {

            buff.write( "TRACE on" );
        }

        buff.newLine();

        if ( rifFile == null )
        {

            buff.write( COMMENT_START + "no RIF file specified (default)" );

        }
        else
        {

            buff.write( "RIF " + rifFile );
        }

        buff.newLine();

        if ( outFile == null )
        {

            buff.write( COMMENT_START + "no output file specified (default)" );

        }
        else
        {

            buff.write( "OUTFILE " + outFile );
        }

        buff.newLine();


        buff.write( "GROUP " + GROUP_ITEMS[groupIndex] );
        buff.newLine();

    }

    /**
     * This routine parses a line from the configuration file and defines
     * the corresponding performance event, performance counter or other
     * additional data.
     *
     * @param line is the input line read from the configuration file
     * @param counterTable is the table of counters to which we add the counter
     */
    public void parseConfigLine( String line, CounterTable counterTable )
    {

        // split the entry in separate items to make parsing easier

        String[] lineItems = line.split( " +" );

        String firstItem = "";  // stands for an empty line

        if ( lineItems.length > 0 )
        {

            firstItem = lineItems[0];
        }

        if ( firstItem.length() == 0 )
        {

            return;
        }

        if ( firstItem.startsWith( COMMENT_SIGN ) )
        {

            // without executable we take the description of performance
            // events from the comment lines
            // # <counter_name> <counter_mode> : <counter_description>

            boolean isDef = false;

            for ( int i = 1; i < lineItems.length; i++ )
            {

                if ( lineItems[i].equals( ":" ) )
                {

                    isDef = true;
                }
            }

            // now we act on this line if it contains definition of basic event

            if ( isDef )
            {

                String eventName = lineItems[1];
                String eventDescription = line.substring( line.indexOf( ":" ) + 2 );

                PerformanceEvent event = findSupportedEvent( eventName );

                if ( event == null )
                {

                    addEvent( new BasicEvent( false, eventName, eventDescription ) );

                    logger.info( "Event, name=" + eventName + ", descr=" + eventDescription );

                }
                else if ( !eventDescription.equals( event.getDescription() ) )
                {

                    // make sure that we have the same event

                    logger.error( "Event " + eventName + " mismatches other event" );
                    logger.error( "new description: " + eventDescription );
                    logger.error( "old description: " + event.getDescription() );

                }

            }

        }
        else if ( firstItem.equals( "SAMPLING" ) )
        {

            // SAMPLING rate size

            dataProfiling = true;

            rateSampling = Integer.parseInt( lineItems[1] );
            sizeSampling = Integer.parseInt( lineItems[2] );

        }
        else if ( firstItem.equals( "DEPTH" ) )
        {

            // DEPTH <profile_depth>

            profilingDepth = Integer.parseInt( lineItems[1] );

        }
        else if ( firstItem.equals( "TRACE" ) )
        {

            // TRACE

            doTracing = true;

        }
        else if ( firstItem.equals( "CALLTREE" ) )
        {

            // CALLTREE

            doCalltree = true;

        }
        else if ( firstItem.equals( "CALLDATA" ) )
        {

            // CALLTREE

            doCalldata = true;

        }
        else if ( firstItem.equals( "GROUP" ) )
        {

            // find the index position in GroupItems

            for ( int i = 0; i < GROUP_ITEMS.length; i++ )
            {

                if ( lineItems[1].equals( GROUP_ITEMS[i] ) )
                {

                    groupIndex = i;
                }
            }

        }
        else if ( firstItem.equals( "RIF" ) )
        {

            rifFile = lineItems[1];

        }
        else if ( firstItem.equals( "OUTFILE" ) )
        {

            outFile = lineItems[1];

        }
        else if ( ( lineItems.length > 1 ) && lineItems[1].equals( "=" ) )
        {

            // check for new event definition <event> = ....

            if ( lineItems.length < 4 )
            {

                throw new IllegalArgumentException( "illegal counter definition" );
            }

            PerformanceEvent event1;
            PerformanceEvent event2;

            String name = lineItems[0];
            String op = lineItems[3];

            event1 = getBasicEvent( lineItems[2] );

            if ( lineItems.length == 5 )
            {

                event2 = getBasicEvent( lineItems[4] );
                ComposedEvent compEvent = new ComposedEvent( name, event1, op, event2 );
                addNewEvent( compEvent );

            }
            else
            {

                DerivedEvent derEvent = new DerivedEvent( name, event1, op );
                addNewEvent( derEvent );
            }

        }
        else
        {

            // now we think it will be a counter line

            Counter newCounter = new Counter( lineItems, this );

            counterTable.addCounter( newCounter );

        }

    } // ParsePMLine

    /***********************************************************************************************
     * * ReadPMFile (file_name) : read a configuration file * *
     **********************************************************************************************/

    /**
     * This routine is used to read a configuration file.
     *
     * @param configFile is the name of the configuration file to read.
     * @param counterTable is the table to which the counters of the configuration file are added.
     */
    public void readConfigFile( String configFile, CounterTable counterTable )
    {

        if ( configFile == null )
        {

            logger.info( "no configuration file available for reading" );
            return;

        }

        try
        {

            int line = 0;
            boolean eof = false;

            FileReader file = new FileReader( configFile );
            BufferedReader buff = new BufferedReader( file );

            while ( !eof )
            {

                String inputLine = buff.readLine();
                line++;

                if ( inputLine == null )
                {

                    eof = true;

                }
                else
                {

                    // handle error in lines separately

                    try
                    {

                        parseConfigLine( inputLine, counterTable );

                    }
                    catch ( IllegalArgumentException e )
                    {

                        String title = "Input error at line " + line;
                        String msg = e.getMessage() + "\n" + "Want to continue?";
                        int response = JOptionPane
                                       .showConfirmDialog( null, msg, title, JOptionPane.YES_NO_OPTION, JOptionPane.ERROR_MESSAGE );

                        if ( response != 0 )
                        {

                            throw new IOException( "stop due to error" );
                        }

                    }
                }

            } // while

            buff.close();

        }
        catch ( IOException e )
        {

            logger.error( e.getMessage(), e );

            String title = "File Input Error";
            String msg = e.getMessage();

            JOptionPane.showMessageDialog( null, msg, title, JOptionPane.ERROR_MESSAGE );

        } // catch

        // successful read, so we set the configuration file if not done yet
        // otherwise we have to assume that another file will be added

        if ( myConfigFile == null )
        {

            myConfigFile = configFile;

        }

        if ( myExtraPanel != null )
        {

            myExtraPanel.showValues();

        }

    }

    /**
     * Getter routine for the name of the configuraton file.
     *
     * @return name of the current config file (can be null)
     */
    public String getConfigFileName()
    {

        return myConfigFile;
    }

    /**
     * Getter routine for the name of the executable.
     *
     * @return name of the current executable (can be null)
     */
    public String getExecutable()
    {

        return myExecutable;
    }

    /**
     * This routine parses a line that will be printed when the -pm help is printed
     * for an executable (e.g. a.out -pm -help).
     *
     * @param line is the line to be parsed
     * @throws IOException when the line cannot be identified
     */
    public void parseExecutableLine( String line ) throws IOException
    {

        // split the entry in separate itmes

        // logger.info("parsing help line items: " + line);

        String[] items = line.split( " +" );

        // for (i=0; i<Items.length; i++)
        // logger.info ("Item " + i + "=<" + Items[i] + ">");

        // identification of a counter line

        if ( items.length < 5 )
        {
            return;
        }

        String eventId = items[3];

        String kindId = items[2];

        boolean isEventLine = kindId.equals( "HPM" ) || kindId.equals( "ADP" );

        if ( !isEventLine )
        {

            return;
        }

        if ( items[4].equals( "=" ) )
        {

            // definiton of an event, e.g. TOT_CYC = PAPI_TOT_CYC

            if ( items.length == 6 )
            {

                // that is ADAPTOR_NAME = HPM_NAME

                PerformanceEvent event = findSupportedEvent( items[5] );

                // add the ADAPTOR name for this event

                event.adaptorName = eventId;
            }

            return;
        }

        // so it is a basic event

        String longName = line.substring( line.indexOf( "(" ) + 1, line.length() - 1 );

        boolean isADAPTOR = true;

        if ( kindId.equals( "HPM" ) )
        {
            isADAPTOR = false;
        }

        PerformanceEvent event = new BasicEvent( isADAPTOR, eventId, longName );

        addEvent( event );

    }

    /**
     * The the name of the executable for which PM data stands.
     *
     * @param name is the name of the executable
     */
    public void setExecutable( String name )
    {

        myExecutable = name;

    }

    /**
     * Getter routine for the configuration file.
     *
     * @return name of the configuration file
     */
    public String getConfigurationFile()
    {

        return myConfigFile;
    }

    /**
     * Getter routine for the panel for the extra configuration.
     *
     * @return the panel
     */
    public JPanel getPanel()
    {

        if ( myExtraPanel == null )
        {

            myExtraPanel = new ExtraPanel( this );

        }

        return myExtraPanel;

    }

    /**
     * Clear all values of the current PM configuration.
     *
     */
    public void clear()
    {

        dataProfiling = false;

        rateSampling = 0;
        sizeSampling = 0;

        doCalltree = false;
        doCalldata = false;
        doTracing  = false;

        groupIndex = GROUP_ITEM_DEFAULT;

        rifFile = null;

        outFile = null;

        // if the GUI is available we show the new values

        if ( myExtraPanel != null )
        {

            myExtraPanel.showValues();

        }

    }

    /**
     * Getter routine for the sampling rate.
     *
     * @return the current sampling rate
     */
    public int getSamplingRate()
    {

        return rateSampling;
    }

    /**
     * Getter routine for the sampling size.
     *
     * @return the current size of sampling
     */
    public int getSamplingSize()
    {

        return sizeSampling;
    }

    /**
     * Getter routine for data profiling flag.
     *
     * @return true if data profiling is enabled
     */
    public boolean doDataProfling()
    {

        return dataProfiling;
    }

    /**
     * Getter routine for the group index.
     *
     * @return the current group setting as an index
     */
    public int getGroupIndex()
    {

        return groupIndex;
    }

    /**
     * Setter routine for the RIF filename.
     *
     * @param name is the name of the RIF file
     */
    public void setRIF( String name )
    {

        rifFile = name;

    }

    /**
     * Setter routine for the tracing.
     *
     * @param enabled is the value to be set for tracing
     * (true for enabled, false for disabled).
     */
    public void setTracing( boolean enabled )
    {

        doTracing = enabled;

    }

    /**
     * Setter routine for the sampling size.
     *
     * @param size is the sampling size
     */
    public void setSamplingSize( int size )
    {

        sizeSampling = size;

    }

    /**
     * Setter routine for the sampling rate.
     *
     * @param rate is the sampling rate
     */
    public void setSamplingRate( int rate )
    {

        rateSampling = rate;

    }

    /**
     * Setter routine for the group index.
     *
     * @param index is the new group index.
     */
    public void setGroupIndex( int index )
    {

        groupIndex = index;

    }

    /**
     * Seter routine for the name of the output file.
     *
     * @param text is the name of the output file
     */
    public void setOutFile( String text )
    {

        outFile = text;

    }

    /**
     * Getter routine for the calltree flag.
     *
     * @return value of doCalltree
     */
    public boolean doCalltree()
    {

        return doCalltree;
    }

    /**
     * Getter routine for the calldata flag.
     *
     * @return value of doCalldata
     */
    public boolean doCalldata()
    {

        return doCalldata;
    }

    /**
     * Setter routine for calldata.
     *
     * @param enabled is the new flag for calldata
     */
    public void setCalldata( boolean enabled )
    {

        doCalldata = enabled;

    }

    /**
     * Setter routine for calltree.
     *
     * @param enabled is the new flag for calltree
     */
    public void setCalltree( boolean enabled )
    {

        doCalltree = enabled;

    }

    /**
     * Getter routine for outFile.
     *
     * @return outfile name
     */
    public String getOutFile()
    {

        return outFile;
    }

    /**
     * Getter routine for doTracing.
     *
     * @return true if tracing is enabled
     */
    public boolean doTracing()
    {

        return doTracing;
    }

    /**
     * Getter routine for the name of the RIF file. Null
     * indicates that no RIF file should be used.
     *
     * @return name of the RIF file
     */
    public String getRIF()
    {

        return rifFile;

    }

    /**
     * Enable or disable data profiling.
     *
     * @param enabled value is used to enable or disable data profiling
     */
    public void setDataProfling( boolean enabled )
    {

        dataProfiling = enabled;

    }

    public void setProfilingDepth( int depth )
    {

        profilingDepth = depth;

    }

    public int getProfilingDepth()
    {

        return profilingDepth;
    }
}
