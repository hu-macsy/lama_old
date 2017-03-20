/*
 * CTDataFile.java
 *
 * Data for a call tree that comes from the ADAPTOR performance runtime system.
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

package adaptor.Calltree;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import org.apache.log4j.Logger;

import adaptor.General.Addr2Line;
import adaptor.General.CallNode;
import adaptor.General.CountedEvent;
import adaptor.General.CounterMetric;
import adaptor.General.FileTable;
import adaptor.General.FileDescriptor;
import adaptor.General.RegionTable;
import adaptor.General.RegionDescriptor;

/**
 * CTData build from CT output file of ADAPTOR run.
 *
 * @version $LastChangedRevision$
 * @author Dr. Thomas Brandes
 */

public class CTDataFile extends CTData
{

    /**
     * Suffix to be used for calltree data files.
     */
    public static final String SUFFIX = "ct";

    /**
     * Suffix to be used for calltree data files.
     */
    public static final String PM_SUFFIX = "pms";

    /**
     * The following token is used to identify the command (name of the executable).
     */
    private static final String CMD_TOKEN = "cmd";

    /**
     * The following token is used to identify the process id.
     */
    private static final String PID_TOKEN = "pid";

    /**
     * Logger for this class.
     */
    private static Logger logger = Logger.getLogger( CTDataFile.class );

    /**
     * This is the process id for the calltree data (not used).
     */
    private int pid;

    /**
     * This is the name of the executable from which the performance data
     * has begin generated.
     */
    private String myExecutable;

    /**
     * Last selection for caller node.
     */
    private CTNode actualCaller = null;

    /**
     * Last selection for called node.
     */
    private CTNode actualCallee = null;

    /**
     * Object for the table with all file descriptors.
     */
    private FileTable myFileTable = null;

    /**
     * Object for the table with all region descriptors.
     */
    private RegionTable myRegionTable = null;


    /**
     * Construct CTData from an input file. The format of this file
     * must satisfy the CT output file format.
     *
     * @param filename is the input file
     */
    CTDataFile( String filename )
    {

        logger.info( "start constructor of CTData from file " + filename );

        myRegionTable = new RegionTable();
        myFileTable   = new FileTable();

        int lineNr = 0;

        try
        {

            String inputLine;

            boolean stop;

            int lastFileId = 0;

            int noEvents = 0; // counts 'event' entries in CT file

            // open the PM result file

            FileReader file = new FileReader( filename );
            BufferedReader buff = new BufferedReader( file );

            stop = false;

            while ( !stop )
            {

                inputLine = buff.readLine();
                lineNr++;

                // logger.debug("read line " + LineNr + ": " + Line);

                stop = ( inputLine == null );

                if ( !stop )
                {

                    String[] inputItems = inputLine.split( " +" );

                    String first = inputItems[0];

                    // only for debug mode
                    // logger.info ("line " + LineNr + ": " + Line + " : "
                    //                      + Items.length + " items");

                    boolean error = false;

                    if ( first.equals( PID_TOKEN ) )
                    {

                        setPid( inputItems[1] );

                    }
                    else if ( first.equals( CMD_TOKEN ) )
                    {

                        setExecutable( inputItems[1] );

                    }
                    else if ( first.equals( "part" ) )
                    {

                        setPart( inputItems[1] );

                    }
                    else if ( first.equals( "event" ) )
                    {

                        setEvent( noEvents++, inputItems );

                    }
                    else if ( first.equals( "define" ) )
                    {

                        setMetric( inputItems );

                    }
                    else if ( first.equals( "events" ) )
                    {

                        dummy( "events seen" );

                    }
                    else if ( first.equals( "totals" ) )
                    {

                        setTotals( inputItems );

                    }
                    else if ( first.equals( "#" ) )
                    {

                        dummy( "comment line" );

                    }
                    else if ( first.length() == 0 )
                    {

                        dummy( "empty line" );

                    }
                    else if ( Character.isDigit( first.charAt( 0 ) ) )
                    {

                        setCallCosts( inputItems );

                    }
                    else if ( first.startsWith( "-" ) )
                    {

                        setCallCosts( inputItems );

                    }
                    else if ( first.equals( "fl" ) || first.equals( "cfl" ) )
                    {

                        // [c]fl  <number>  [ <file_name> ]

                        lastFileId = Integer.parseInt( inputItems[1] );

                        FileDescriptor dspFile;

                        if ( inputItems.length > 2 )
                        {

                            dspFile = myFileTable.getFile( lastFileId, inputItems[2] );

                        }
                        else
                        {

                            dspFile = myFileTable.getFile( lastFileId );
                        }

                        if ( dspFile == null )
                        {

                            logger.error( "could not get file by id " + lastFileId + " (" + inputItems.length + " items)" );

                            error = true;
                        }

                    }
                    else if ( first.equals( "fn" ) )
                    {

                        setRegionNode( inputItems, myFileTable.getFile( lastFileId ), false ); // not called

                    }
                    else if ( first.equals( "cfn" ) )
                    {

                        setRegionNode( inputItems, myFileTable.getFile( lastFileId ), true ); // is callee

                    }
                    else if ( first.equals( "calls" ) )
                    {

                        setNoRunCalls( inputItems );

                    }
                    else
                    {

                        logger.error( filename + ", line = " + lineNr + " : line not identified" );

                    }

                    if ( error )
                    {

                        logger.error( "serious error appeared in line " + lineNr );
                        logger.error( inputLine );
                        System.exit( 0 );

                    }

                }

            } // while

            buff.close();

        }
        catch ( FileNotFoundException e )
        {

            logger.error( e.getMessage() );
            System.exit( 0 );

        }
        catch ( IOException e )
        {

            logger.error( "Error(IO) in " + filename + ", Line " + lineNr + ": " + e.toString() );
            System.exit( 0 );

        }
        catch ( RuntimeException e )
        {

            logger.error( "Error(Runtime) in " + filename + ", Line " + lineNr + ": " + e.toString() );
            System.exit( 0 );

        }


        // new calculation of inclusive costs

        // calcInclCosts ();

    } // constructor CTDataFile

    /**
     * Dummy routine needed to avoid empty statement.
     *
     * @param msg is dummy argument
     */
    private void dummy( String msg )
    {

    }

    /**
     * This routine sets the process id.
     *
     * @param v is a String containing the integer for pid
     */
    private void setPid( String v )
    {

        pid = Integer.parseInt( v );

    }

    /**
     * Setter routine to set the name of the executable.
     *
     * @param cmd is the name of the executable
     */
    private void setExecutable( String cmd )
    {

        myExecutable = cmd;
    }

    /**
     * Set name of the part.
     *
     * @param v is name of the part
     */
    private void setPart( String v )
    {

        // currently ignored

    }

    /**
     * This routine defines a new event.
     *
     * @param eventIndex is the index for the event
     * @param items contains name of the event and description
     */
    private void setEvent( int eventIndex, String[] items )
    {

        String eventName;
        String eventDescription;

        int noItems = items.length;

        eventName = items[1];

        eventDescription = "";

        for ( int i = 2; i < noItems; i++ )
        {

            eventDescription = eventDescription + items[i];

            if ( i + 1 < noItems )
            {

                eventDescription = eventDescription + " ";
            }

        } // for all other items

        // logger.info("Event :" + Event + ": = <" + Event_Description + ">");

        CounterData.defineEvent( eventIndex, eventName, eventDescription );

    } // setEvent

    /**
     * Define a new metric specified by the input items.
     * Here are some examples:
     *
     * <ul>
     * <li>
     * <code>define WALL_TIME 3.563528e-10 WALL_TICKS</code>
     * <li>
     * <code>define BUS_RATIO / 1.0 BUS 1.0 WALL_TIME</code>
     * <li>
     * <code>define MFLOPS / 1.0e-6 FP_INS_C 1.0 WALL_TIME</code>
     * </ul>
     *
     * @param inputItems are the blank separated items of the input line
     */
    private void setMetric( String[] inputItems )
    {

        String metricName = inputItems[1];

        if ( inputItems.length == 4 )
        {

            double rate = Double.parseDouble( inputItems[2] );

            CountedEvent event = CounterData.getCountedEvent( inputItems[3] );

            if ( event == null )
            {

                logger.info( "set this simple metric : event " + inputItems[3] );

                return;
            }

            CounterMetric newMetric = new CounterMetric( metricName, rate, event );

            CounterData.defineMetric( newMetric );

        }
        else
        {

            char mode = inputItems[2].charAt( 0 );

            double rate1 = Double.parseDouble( inputItems[3] );
            double rate2 = Double.parseDouble( inputItems[5] );

            CountedEvent event1 = CounterData.getCountedEvent( inputItems[4] );

            if ( event1 == null )
            {

                logger.info( "set_derived: event " + inputItems[4] + " not counted" );
            }

            CountedEvent event2 = CounterData.getCountedEvent( inputItems[6] );

            if ( event2 == null )
            {

                logger.info( "set_derived: event " + inputItems[6] + " not counted" );
            }

            if ( ( event1 == null ) || ( event2 == null ) )
            {

                return;
            }

            CounterMetric newMetric = new CounterMetric( metricName, mode, rate1, event1, rate2, event2 );
            CounterData.defineMetric( newMetric );

        }

        logger.info( "new derived event " + metricName + " defined" );

    } // setMetric

    /**
     * This routines make a region to the current node.
     *
     * @param inputItems are the items parsed
     * @param file is the file descriptor belonging to the region
     * @return the region descriptor for the defined region
     */
    private RegionDescriptor defineRegion( String[] inputItems, FileDescriptor file )
    {

        // fn  <number> <region_name> <region_kind> <line_start> <line_stop> <class_name>

        int nr = Integer.parseInt( inputItems[1] );

        RegionDescriptor region;

        if ( inputItems.length > 2 )
        {

            // the first two items "c[fn] <number>" are irrelevant

            int regionId = Integer.parseInt( inputItems[2] );

            String regionName = inputItems[3];

            int regionKind = Integer.parseInt( inputItems[4] );

            int startLine = 0;

            if ( inputItems.length > 5 )
            {

                startLine = Integer.parseInt( inputItems[5] );
            }

            int stopLine = startLine;

            if ( inputItems.length > 6 )
            {

                stopLine = Integer.parseInt( inputItems[6] );
            }

            // avoid empty class names as they might cause problems

            String className = "__";

            if ( inputItems.length > 7 )
            {

                className = inputItems[7];
            }

            region = new RegionDescriptor( regionId, regionName, className, regionKind, file, startLine, stopLine );

            // if there is an existing descriptor we will take this one

            region = myRegionTable.defineRegion( nr, region );

        }
        else
        {

            // the region has already been defined in the CT file
            // Attention getRegion (nr) can be another region

            region = getNode( nr ).getRegion();

        }

        return region;

    }

    /**
     * This routines make a region to the current node.
     *
     * @param inputItems are the items parsed
     * @param file is the file descriptor belonging to the region
     * @param calleeFlag if true sets actualCallee otherwise actualCaller
     */
    private void setRegionNode( String[] inputItems, FileDescriptor file, boolean calleeFlag )
    {

        // fn  <number> <region_name> <region_kind> <line_start> <line_stop> <class_name>

        int nr = Integer.parseInt( inputItems[1] );

        RegionDescriptor region = defineRegion( inputItems, file );

        logger.debug( "setRegionNode, nr = " + nr + " is " + region.getName() );

        CTNode node = getNode( nr );

        if ( node == null )
        {

            node = new CTNode( new CallNode( region ) );

            addCTNode( nr, node );

        }

        if ( calleeFlag )
        {

            actualCallee = node;

        }
        else
        {

            actualCaller = node;
        }

    } // setRegionNode

    /**
     * This routine sets the number of runtime calls
     *
     * @param inputItems is { "calls", "15" }
     */

    private void setNoRunCalls( String[] inputItems )
    {

        String val = inputItems[1]; // calls 135

        int calls = Integer.parseInt( val );

        logger.debug( "add calls for " + actualCaller.getName() + " -> "
                      + actualCallee.getName() + " = " + val );

        // get edge with actual caller and actual callee

        CTEdge edge = getEdge( actualCaller, actualCallee );

        edge.addNoRunCalls( calls );

        logger.debug( "add calls for " + edge.getSource().getName() + " -> "
                      + edge.getTarget().getName() + " = " + val );

    }

    /**
     * Set call costs for the current selected node pair.
     *
     * @param costItems is line costs1 costs2 ... costsn
     */
    private void setCallCosts( String[] costItems )
    {

        // important: first value stands for the line

        long[] costs = new long[costItems.length - 1];

        for ( int i = 1; i < costItems.length; i++ )
        {

            costs[i - 1] = Long.parseLong( costItems[i] );

        }

        // make sure that we really have correct number of values

        CounterData.checkCounters( costItems.length - 1, "set values" );

        if ( actualCaller == null )
        {

            logger.fatal( "no actual caller set for set values" );
            System.exit( 0 );

        }

        if ( actualCallee == null )
        {

            actualCaller.addExclusiveCosts( costs );
            actualCaller.addInclusiveCosts( costs );

        }
        else
        {

            /* debug code:

             logger.info("add call costs " + actual_caller.getName () +
             " -> " + actual_callee.getName ());

             for (int j=0; j<Values.length; j++)
             logger.info ("  call costs counter " + j + " = " + Values[j]);

             */

            CTEdge edge = getEdge( actualCaller, actualCallee );

            edge.addCallCosts( costs );

        }

        // now we have to reset the current entries

        actualCaller = null;
        actualCallee = null;

    } // setCallCosts

    /**
     * Set values for the toal costs of the event counters, e.g.
     * <code>totals: 1610613052 4688257947 11636557472</code>.
     *
     * @param totalItems are the item strings of the line
     */

    private void setTotals( String[] totalItems )
    {

        // important: first Item is "totals:"

        // make sure that Items.length - 1 == NoCounters

        long[] totals = new long[totalItems.length - 1];

        for ( int i = 1; i < totalItems.length; i++ )
        {

            totals[i - 1] = Long.parseLong( totalItems[i] );

        }

        CounterData.setTotals( totals );

    } // setTotals

    private static void outNode( BufferedWriter buff, int index, CTNode node, boolean full, boolean callee ) throws IOException
    {

        RegionDescriptor region = node.getRegion();
        FileDescriptor fdsp = region.getFile();

        if ( callee )
        {
            buff.write( "c" );
        }

        buff.write( "fl " + fdsp.getFileId() );

        if ( full )
        {

            buff.write( " " + fdsp.getLongFileName() );
        }

        buff.newLine();

        if ( callee )
        {
            buff.write( "c" );
        }

        buff.write( "fn " + index );

        if ( full )
        {

            buff.write( " " + region.getRegionId() );
            buff.write( " " + region.getName() + " " + region.getKind() );
            buff.write( " " + region.getFirstLine() + " " + region.getLastLine() );
            buff.write( " " + region.getClassName() );
        }

        buff.newLine();

    }

    /**
     * Output of the calltree data into a file. The kind of output
     * is determined by the suffix of the filename.
     *
     * @param outFile is the name of the output file.
     */
    void output( String outFile )
    {

        String fileName = outFile;

        int pointIndex = outFile.lastIndexOf( "." );

        logger.info( "output of " + outFile + ", pointIndex = " + pointIndex );

        String suffix;

        if ( pointIndex < 0 )
        {

            suffix = "ct";
            fileName += "." + suffix;

        }
        else
        {

            suffix  = outFile.substring( pointIndex + 1 );
        }

        logger.info( "suffix is <" + suffix + ">" );

        if ( suffix.equals( "rif" ) )
        {

            exportRIF( fileName );

        }
        else if ( suffix.equals( PM_SUFFIX ) )
        {

            exportPMS( fileName );

        }
        else
        {

            exportCT( fileName );
        }

        logger.info( "Output is available in file " + fileName );
    }

    /**
     * Translate region binary addresses to region names + file ids.
     * Be careful: nodes are not always regions
     */
    void addr2line()
    {

        Addr2Line.addr2line( myRegionTable, myFileTable, myExecutable );
    }
    /**
     * This routine exports the visible PM data to an output
     * file that can be read by the visualization tool ShowPM.
     *
     * @param filename is the name for the output file of the PMS data
     */
    void exportPMS( String filename )
    {

        CTNode node;
        int i;

        try
        {

            FileWriter file = new FileWriter( filename );
            BufferedWriter buff = new BufferedWriter( file );

            buff.write( "COMMAND=" + myExecutable );
            buff.newLine();

            buff.write( "NP=1" );
            buff.newLine();

            buff.write( "NT=1" );
            buff.newLine();

            int counterRegions = countNodeRegions();

            buff.write( "NR=" + counterRegions );
            buff.newLine();

            for ( i = 0; i < getNumberNodes(); i++ )
            {

                node = getNode( i );

                if ( node != null )
                {

                    RegionDescriptor region = node.getRegion();

                    if ( region != null )
                    {

                        // PMS has own representation for regions

                        buff.write( region.getAttributeString() );
                        buff.newLine();

                    }
                }

            } // for all Regions

            CounterData.writeCounters( buff );

            for ( i = 0; i < getNumberNodes(); i++ )
            {

                node = getNode( i );

                if ( node != null )
                {

                    if ( node.getRegion() != null )
                    {

                        node.writeCounterVals( buff );
                    }
                }

            } // for all Regions

            buff.close();

        }
        catch ( IOException e )
        {

            logger.error( "IO Exception  -- " + e.toString() );
        }

        logger.info( ">>>>> PMS output done" );

    }

    /**
     * Output of the calltree data into a File.
     *
     * @param outFile is the name of the output file for the CT data.
     */
    void exportCT( String outFile )
    {

        try
        {

            FileWriter file = new FileWriter( outFile );
            BufferedWriter buff = new BufferedWriter( file );

            buff.write( "# compressed PM output data" );
            buff.newLine();
            buff.write( PID_TOKEN + " " + pid );
            buff.newLine();
            buff.write( CMD_TOKEN + " " + myExecutable );
            buff.newLine();
            buff.write( "part 1" );
            buff.newLine();

            // Ouptut of all Events + Metrics

            buff.newLine();
            buff.write( "# output of counters and metrics" );
            buff.newLine();

            CounterData.output( buff );

            // Output of all nodes / regions

            buff.newLine();
            buff.write( "# output of all regions and their costs" );
            buff.newLine();

            for ( int i = 0; i < getNumberNodes(); i++ )
            {

                CTNode node = getNode( i );

                if ( node != null )
                {

                    // we should be sure that it is a region

                    outNode( buff, i, node, true, false );

                    long[] counterVals = node.getExclusiveCosts();

                    buff.write( Long.toString( node.getRuntimeCalls() ) );

                    for ( int k = 0; k < counterVals.length; k++ )
                    {

                        buff.write( " " + counterVals[k] );
                    }

                    buff.newLine();

                }

            } // for all Regions

            buff.newLine();
            buff.write( "# output of all calls and their costs" );
            buff.newLine();

            for ( int i = 0; i < getNumberEdges(); i++ )
            {

                CTEdge edge = getEdge( i );

                CTNode source = edge.getSource();
                CTNode target = edge.getTarget();

                // out of caller and callee node

                outNode( buff, getNodeIndex( source ), source, false, false );
                outNode( buff, getNodeIndex( target ), target, false, true );

                buff.write( "calls " + edge.getRuntimeCalls() + " 0" );
                buff.newLine();

                long[] counterVals = edge.getCallCosts();

                // the first entry is the line id of the call (but we make no difference)

                buff.write( "1" );

                for ( int k = 0; k < counterVals.length; k++ )

                {
                    buff.write( " " + counterVals[k] );
                }

                buff.newLine();

            }

            // Output of the total counter values for verification

            CounterData.outTotals( buff );

            buff.close();

        }
        catch ( IOException e )
        {

            logger.error( "IO Exception  -- " + e.toString() );
        }
    }

} // class CTDataFile
