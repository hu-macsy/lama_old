/*
 * PMData.java
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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

import org.apache.log4j.Logger;

import adaptor.General.Addr2Line;
import adaptor.General.ArraySelection;
import adaptor.General.FileDescriptor;
import adaptor.General.FileTable;
import adaptor.General.RIFData;
import adaptor.General.RegionDescriptor;
import adaptor.General.RegionTable;
import adaptor.General.Utilities;

/**
 * The class PMData contains all data of an PM output file. There might be many instances of this
 * class containing performance data from different runs.
 *
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */

public class PMData
{

    /**
     * This is the split pattern used for splitting entries of an input line.
     * Entries are splitted due to a non-empty set of blanks.
     */
    private static final String SPLIT_STRING = " +";

    /**
     * Logger variable for this class.
     */
    private static Logger logger = Logger.getLogger( PMData.class );

    /**
     * Array with all region counters.
     */
    private RegionCounter[] myRegionCounters;

    /**
     * Array with all user counters.
     */
    private UserCounter[] myUserCounters;

    /**
     * Object with all regions. Region descriptors can be directly
     * accessed by the index.
     */
    private RegionTable myRegionTable;

    /**
     * Array with file descriptor for each region.
     */
    private FileTable myFileTable;

    /**
     * This is the name of the file where the PM data comes from.
     * The name is used also for identification in certail tables.
     */
    private String myFileName; // data is identified by the name

    /**
     * This is the name of the executable for which performance
     * data has been collected.
     */
    private String myExecutable = "";

    /**
     * myPMCounterValues [numberProcesses][numberThreads] [numberRegions][numberRegionCounters]
     * contains the values of the counted events for all processors, threads, regions and counters.
     */

    private long[][][][] myPMCounterValues;

    /**
     * This is the number of processors for which performance data is available.
     * For single processes we have the value 1.
     */
    private int numberProcesses;

    /**
     * This is the number of threads for which performance data is available.
     */
    private int numberThreads;

    /**
     * This is the number of regions.
     */
    private int numberRegions; // number of regions

    /***********************************************************************************************
     * constructor PMData (String PM_filename) * - reads all the performance results from a file *
     **********************************************************************************************/

    /**
     * Constructor for performance data by the name of the file with PM data. File will be
     * read and the PM data constructed.
     *
     * @param inputFileName is the name of the input file read for performance data.
     */
    PMData( String inputFileName )
    {

        int lineNumber = 0;

        logger.info( "read PMData from file " + inputFileName );

        try
        {

            String inputLine;

            int indexProc;
            int indexThread;
            int indexRegion;
            int indexCounter;

            // open the PM result file

            FileReader file = new FileReader( inputFileName );
            BufferedReader buff = new BufferedReader( file );

            /***************************************************************************************
             * COMMAND=<executable>                                                                *
             **************************************************************************************/

            inputLine = buff.readLine();
            lineNumber++;
            myExecutable = Utilities.readString( inputLine, "COMMAND" );
            logger.info( "COMMAND = " + myExecutable );

            /***************************************************************************************
             * NP=<number_of_processes> * NT=<number_of_threads> *
             **************************************************************************************/

            // read from line: NP=<number_of_processes>
            inputLine = buff.readLine();
            lineNumber++;
            numberProcesses = Utilities.readInt( inputLine, "NP" );
            logger.info( "NP = " + numberProcesses );

            // read from line NT=<number_of_threds>

            inputLine = buff.readLine();
            lineNumber++;
            numberThreads = Utilities.readInt( inputLine, "NT" );
            logger.info( "NT = " + numberThreads );

            /***************************************************************************************
             * NR=<number_of_regions> *
             **************************************************************************************/

            // read from line NR=<number_of_regions>
            inputLine = buff.readLine();
            lineNumber++;
            numberRegions = Utilities.readInt( inputLine, "NR" );
            logger.info( "NR = " + numberRegions );

            myRegionTable = new RegionTable( numberRegions );

            // Attention: number of (different) files are not known

            myFileTable   = new FileTable( numberRegions );

            for ( indexRegion = 0; indexRegion < numberRegions; indexRegion++ )
            {

                inputLine = buff.readLine();

                lineNumber++;

                if ( inputLine == null )
                {

                    throw new IOException( onlyErrorMessage( "regions", indexRegion, numberRegions ) );

                }

                lineNumber++;

                /***********************************************************************************
                 * region_id region_name class kind file line_start line_stop                      *
                 **********************************************************************************/

                String[] items = inputLine.split( SPLIT_STRING );

                if ( items.length < 7 )
                {

                    String msg = "no region line, expected: id name class kind file start stop";

                    logger.error( msg );

                    throw new IOException( msg );
                }

                // get a new or an old file descriptor for the filename = items[4]

                FileDescriptor fileDSP = myFileTable.getFile( items[4] );

                // regionId = indexRegion would be another possibility

                int    regionId   = Integer.parseInt( items[0] );
                String regionName = items[1];                   //  can be "?"
                String className  = items[2];                   //  can be "?"
                int    regionKind = Integer.parseInt( items[3] );
                int    lineStart  = Integer.parseInt( items[5] );
                int    lineStop   = Integer.parseInt( items[6] );

                RegionDescriptor rdsp = new RegionDescriptor( regionId, regionName, className, regionKind, fileDSP, lineStart, lineStop );

                myRegionTable.defineRegion( indexRegion, rdsp );
            }

            // read now line:  NC=<number_of_region_counters>

            inputLine = buff.readLine();

            lineNumber++;

            int numberRegionCounters = Utilities.readInt( inputLine, "NC" );

            logger.info( "NC = " + numberRegionCounters );

            myRegionCounters = new RegionCounter[numberRegionCounters];

            for ( indexCounter = 0; indexCounter < numberRegionCounters; indexCounter++ )
            {

                inputLine = buff.readLine();

                lineNumber++;

                if ( inputLine == null )
                {

                    throw new IOException( onlyErrorMessage( "region counters", indexCounter, numberRegionCounters ) );

                }

                // region counter is: <counter_name> <counter_mode>

                String[] items = inputLine.split( SPLIT_STRING );

                myRegionCounters[indexCounter] = new RegionCounter( items[0], indexCounter, items[1] );

            }

            /***************************************************************************************
             * NU=<number_of_user_counters> *
             **************************************************************************************/

            inputLine = buff.readLine();
            lineNumber++;

            int numberUserCounters = Utilities.readInt( inputLine, "NU" );

            logger.info( "NU = " + numberUserCounters );

            myUserCounters = new UserCounter[numberUserCounters];

            for ( indexCounter = 0; indexCounter < numberUserCounters; indexCounter++ )
            {

                inputLine = buff.readLine();

                if ( inputLine == null )
                {

                    throw new IOException( onlyErrorMessage( "user counters", indexCounter, numberUserCounters ) );

                }

                lineNumber++;

                String[] items = inputLine.split( SPLIT_STRING );

                if ( items.length == 7 )
                {

                    String name  = items[0];
                    String mode  = items[1];
                    String scale = items[2];

                    // Items [3] == precision is ignored
                    // Items [4] == width is ignored

                    double rate = Double.parseDouble( items[5] );
                    int index = Integer.parseInt( items[6] );


                    myUserCounters[indexCounter] = new UserCounter( name, mode, scale, rate, myRegionCounters[index] );

                }
                else if ( items.length == 10 )
                {


                    String name  = items[0];
                    String mode  = items[1];
                    String scale = items[2];

                    // Items [3] == precision is ignored
                    // Items [4] == width is ignored

                    String op = items[5];

                    double rate1 = Double.parseDouble( items[6] );
                    int index1   = Integer.parseInt( items[7] );

                    double rate2 = Double.parseDouble( items[8] );
                    int index2   = Integer.parseInt( items[9] );

                    myUserCounters[indexCounter] = new UserCounter( name, mode, scale, op,
                            rate1, myRegionCounters[index1],
                            rate2, myRegionCounters[index2] );
                }
                else
                {

                    String msg = "illegal user counter line (7 or 10 items expected)";

                    logger.error( msg );

                    throw new IOException( msg );

                }

            }

            // now we read the info about the user counters

            myPMCounterValues = new long[numberProcesses][numberThreads][numberRegions][numberRegionCounters];

            for ( indexProc = 0; indexProc < numberProcesses; indexProc++ )
            {

                for ( indexThread = 0; indexThread < numberThreads; indexThread++ )
                {

                    for ( indexRegion = 0; indexRegion < numberRegions; indexRegion++ )
                    {

                        inputLine = buff.readLine();
                        lineNumber++;

                        String[] items = inputLine.split( SPLIT_STRING );

                        if ( items.length == numberRegionCounters )
                        {

                            // we have really one value for each region counter

                            for ( indexCounter = 0; indexCounter < numberRegionCounters; indexCounter++ )
                            {

                                myPMCounterValues[indexProc][indexThread][indexRegion][indexCounter] = Long.parseLong( items[indexCounter] );

                            }

                        }
                        else
                        {

                            throw new IOException( onlyErrorMessage( "region counter values", items.length, numberRegionCounters ) );

                        }
                    }

                }

            }

            buff.close();

        }
        catch ( FileNotFoundException e )
        {

            logger.error( e.getMessage() );
            System.exit( 0 );

        }
        catch ( IOException e )
        {

            logger.error( "PM data error in line " + lineNumber + " of file " + inputFileName );
            logger.fatal( e.toString() );
            System.exit( 0 );

        }

        myFileName = inputFileName;

        addr2line();

        logger.info( "PMData has been construced from file " + inputFileName );

    } // constructor

    /**
     * Error message if only a certain number of values has been read.
     *
     * @param kind is a string for the kind of item
     * @param realValue is the value of read items
     * @param expectedValue is the value of expected items
     * @return the error message
     */
    private String onlyErrorMessage( String kind, int realValue, int expectedValue )
    {

        String msg = "read only " + realValue + " of " + expectedValue + " expected " + kind;

        logger.error( msg );

        return msg;

    }

    /**
     * This routine returns a name for this performance data. It is the file name
     * but where the path will not appear any more.
     *
     * @return String for the name of the performance data
     */
    String getName()
    {

        File pmFile = new File( myFileName );

        return pmFile.getName();
    }

    /**
     * Get the number of processes for which performance data is available.
     *
     * @return integer value for number of processes
     */
    public int getNumberProcesses()
    {

        return numberProcesses;
    }

    /**
     * Get the number of threads for which performance data is available.
     *
     * @return integer value for number of threads
     */
    public int getNumberThreads()
    {

        return numberThreads;
    }

    /**
     * Get the number of user counters that have been defined.
     *
     * @return integer value for number of user counters
     */

    public int getNumberUserCounters()
    {

        return myUserCounters.length;
    }

    /**
     * Get the number of region counters (for each region we have this
     * number of counted values).
     *
     * @return integer value for number of processes
     */

    public int getNumberRegionCounters()
    {

        return myRegionCounters.length;

    }

    /**
     * Get the number of regions.
     *
     * @return integer value for number of regions
     */

    public int getNumberRegions()
    {

        return numberRegions;

    }

    /**
     * Getter routine for the name of the executable.
     *
     * @return the name of the executable
     */
    public String getExecutable()
    {

        return myExecutable;
    }

    /**
     * Get the descriptor of a region.
     *
     * @param index is the index of the region
     * @return the descriptor of the region
     */
    public RegionDescriptor getRegion( int index )
    {

        return myRegionTable.getRegion( index );
    }

    /**
     * Get the file descriptor for a region.
     *
     * @param index is the file identification as index of file table
     * @return the file descriptor for this region
     */
    public FileDescriptor getFileDescriptor( int index )
    {

        return myFileTable.getFile( index );

    }

    /**
     * This routine translates region addresses to region names.
     * We know that we have a region address if the filename is the executable.
     */

    public void addr2line()
    {

        logger.info( "addr2line: find addresses of executable " + myExecutable );

        Addr2Line.addr2line( myRegionTable, myFileTable, myExecutable );
    }

    /***********************************************************************************************
     * compressPM - remove entries with all counter values equal 0 *
     **********************************************************************************************/

    /**
     * This routine removes the performance data for all regions that have only zero counter values.
     * This is likely to happen when using RIF files with regions that are never called at all.
     */
    void compressPM()
    {

        logger.info( "compressPM for " + numberRegions + " regions" );

        int numberUsedRegions; // number of regions really used

        numberUsedRegions = 0;

        for ( int region = 0; region < numberRegions; region++ )
        {

            boolean used = !isZeroRegion( region );

            if ( used )
            {

                // compression works fine as NR_used <= IR

                inheritRegionVals( numberUsedRegions, myPMCounterValues, region );

                // ATTENTION: defineRegion is not the right routine here as it does
                //            not allow to override regions in the region table

                myRegionTable.setRegion( numberUsedRegions, myRegionTable.getRegion( region ) );
                myFileTable.defineFile( numberUsedRegions, myFileTable.getFile( region ) );

                numberUsedRegions++;

            } // used region

        } // for all regions

        if ( numberUsedRegions < numberRegions )
        {

            logger.info( ( numberRegions - numberUsedRegions ) + " regions ignored (zero values)" );

        }

        logger.info( "compressPM for " + numberRegions + " regions, used " +
                     numberUsedRegions + " regions" );

        numberRegions = numberUsedRegions;

    } // compressPM

    /**
     * This routine returns the name of a counter and will be
     * used as header in the tables.
     *
     * @param isUser is true if user counter
     * @param selCounter is the index of the counter
     * @return string for the header
     */
    public String getCounterHeader( boolean isUser, int selCounter )
    {

        String header;

        if ( isUser )
        {

            header = myUserCounters[selCounter].getHeader();

        }
        else
        {

            header = myRegionCounters[selCounter].getHeader();
        }

        return header;

    } // getCounterDescription

    /**
     * This routine gets the value of a user counter for given
     * processor and region and counter number.
     *
     * @param process is the index of the process
     * @param thread is the index of the thread
     * @param region is the index of the region
     * @param counter is the index of the counter
     * @return a double value with the value of the user counter
     */
    double getUserVal( int process, int thread, int region, int counter )
    {

        UserCounter uC = myUserCounters[counter];

        long[] counterVals = myPMCounterValues[process][thread][region];

        return uC.getScaledValue( counterVals );

    }

    /**
     * Direct access to a user counter by its index position.
     *
     * @param counter is the index of the user counter
     * @return a user counter at the given position.
     */
    UserCounter getUserCounter( int counter )
    {

        return myUserCounters[counter];
    }

    /**
     * Match of regions of this performance data set with another set of performance data.
     *
     * @param newPM is the other performance data set against which this set is matched
     * @return "" in case of success, otherwise match problems
     */
    String matchRegions( PMData newPM )
    {

        // seems to be okay now

        if ( numberRegions != newPM.numberRegions )
        {

            return "PMData: new set has " + newPM.numberRegions + " regions instead of " + numberRegions;

        }

        for ( int region = 0; region < numberRegions; region++ )
        {

            // make sure that we have the same regions

            String regionName1 = myRegionTable.getRegion( region ).getName();
            String regionName2 = newPM.myRegionTable.getRegion( region ).getName();

            if ( !regionName1.equals( regionName2 ) )
            {

                return "PMData: region " + region + " " + regionName1 + " not " + regionName2;

            }

        } // for all regions

        return "";

    } // matchRegions

    /**
     * This routine sets all counter values of a region to zero.
     *
     * @param region is the index of the region
     */
    private void setRegionZero( int region )
    {

        int noCounters = myRegionCounters.length;

        for ( int indexProc = 0; indexProc < numberProcesses; indexProc++ )
        {

            for ( int indexThread = 0; indexThread < numberThreads; indexThread++ )
            {

                for ( int indexCounter = 0; indexCounter < noCounters; indexCounter++ )
                {

                    myPMCounterValues[indexProc][indexThread][region][indexCounter] = 0;
                }
            }
        }
    }

    /**
     * This routine inherits for a region the values of another counter value array.
     *
     * @param region is the region for which counter values to set
     * @param vals is the array with all counter values
     * @param oldIndex is the region from which we inherit the values
     */
    private void inheritRegionVals( int region, long [][][][] vals, int oldIndex )
    {

        int noCounters = myRegionCounters.length;

        for ( int indexProc = 0; indexProc < numberProcesses; indexProc++ )
        {

            for ( int indexThread = 0; indexThread < numberThreads; indexThread++ )
            {

                for ( int indexCounter = 0; indexCounter < noCounters; indexCounter++ )
                {

                    long val = vals[indexProc][indexThread][oldIndex][indexCounter];

                    myPMCounterValues[indexProc][indexThread][region][indexCounter] = val;


                }
            }
        }
    }

    /**
     * Query method to see whether all values of a region are zero.
     *
     * @param region is the index of the region we want to check
     * @return true if all counter values of region are zero
     */
    private boolean isZeroRegion( int region )
    {

        int noCounters = myRegionCounters.length;

        for ( int indexProc = 0; indexProc < numberProcesses; indexProc++ )
        {

            for ( int indexThread = 0; indexThread < numberThreads; indexThread++ )
            {

                for ( int indexCounter = 0; indexCounter < noCounters; indexCounter++ )
                {

                    if ( myPMCounterValues[indexProc][indexThread][region][indexCounter] != 0 )
                    {

                        return false;
                    }

                }
            }
        }

        return true;
    }

    /**
     * This routine aligns the regions of this performance data with the regions of another
     * performance data set. Afterwards both sets have the same table of regions.
     *
     * @param newPM is the other performance data set.
     */
    void alignRegions( PMData newPM )
    {

        // set already the new number of regions and save old one

        int oldNumberRegions = numberRegions;
        int newNumberRegions = 0; // count new filled regions

        numberRegions = newPM.numberRegions;

        // save the old performance data

        long[][][][] oldPMVals = myPMCounterValues;

        int noCounters = myRegionCounters.length;

        myPMCounterValues = new long[numberProcesses][numberThreads][numberRegions][noCounters];

        // make new RegionTable and FileTable

        RegionTable oldRegionTable = myRegionTable;

        myRegionTable = new RegionTable( numberRegions );
        myFileTable   = new FileTable( numberRegions );

        // now create new performance data

        for ( int indexRegion = 0; indexRegion < numberRegions; indexRegion++ )
        {

            // try to find region in

            RegionDescriptor regionDSP = newPM.myRegionTable.getRegion( indexRegion );

            myRegionTable.defineRegion( indexRegion, regionDSP );

            myFileTable.defineFile( indexRegion, newPM.myFileTable.getFile( indexRegion ) );

            int newRegion = ArraySelection.NO_SELECTION;

            for ( int region1 = 0; region1 < oldNumberRegions; region1++ )
            {

                RegionDescriptor regionDSP1 = oldRegionTable.getRegion( region1 );

                if ( regionDSP1.getName().equals( regionDSP.getName() ) )
                {

                    newRegion = region1;
                }

            }

            if ( newRegion == ArraySelection.NO_SELECTION )
            {

                // this data set has not a region so we set default counter values with 0

                newNumberRegions++;

                setRegionZero( indexRegion );

            }
            else
            {

                // take over the old values

                inheritRegionVals( indexRegion, oldPMVals, newRegion );

            }

        } // for all new regions

        logger.info( "alignment of regions: " + newNumberRegions
                     + " regions inserted, " + ( oldNumberRegions - ( numberRegions - newNumberRegions ) )
                     + " regions removed" );

    } // alignRegions

    /**
     * This routine gets all the values of the counted events for a given
     * region and a given process-thread pair.
     *
     * @param selectedProc is the index of the process
     * @param selectedThread is the index of the thread
     * @param indexRegion is the index of the region
     * @return the values of the counted events
     */
    public long[] getRegionCounterVals( int selectedProc, int selectedThread, int indexRegion )
    {

        return myPMCounterValues[selectedProc][selectedThread][indexRegion];

    }

    /**
     * This routine gets all the values of the counted events for a given
     * region and a given processor.
     *
     * @param indexProcessor is the index of the processor
     * @param indexRegion is the index of the region
     * @return an array with the values of the counted events
     */
    public long[] getRegionCounterVals( int indexProcessor, int indexRegion )
    {

        int iP = indexProcessor / numberThreads;

        int iT = indexProcessor - iP * numberThreads;

        return myPMCounterValues[iP][iT][indexRegion];

    }

    /**
     * This routine returns a string specifying the process
     * and thread of a given processor.
     *
     * @param indexProcessor is the index of the processor (starting with 0)
     * @return a string like "P=1,T=2"
     */
    public String getProcessorString( int indexProcessor )
    {

        String procString = "";

        // thread and process is coded in one value

        int indexProcess = indexProcessor / numberThreads;

        if ( numberProcesses > 1 )
        {

            procString += "P=" + indexProcess;

        }

        int indexThread = indexProcessor - indexProcess * numberThreads;

        if ( numberThreads > 1 )
        {

            // make a separator if necessary

            if ( procString.length() != 0 )
            {

                procString += ", ";
            }

            procString += "T=" + indexThread;

        }

        if ( procString.length() == 0 )
        {

            procString = getName();
        }

        return procString;
    }


    /**
     * Write an RIF file from the current PM data.
     *
     * @param filename is the name of the output file
     *
     */
    public void writeRIF( String filename )
    {

        // Step 1: generate RIF data from the PM data

        RIFData rif = new RIFData( myFileTable, myRegionTable );

        // Step 2: write the RIF data

        rif.writeRIF( filename );

    }

} // class PMData

