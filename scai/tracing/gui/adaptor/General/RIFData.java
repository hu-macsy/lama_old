/*
 * RIFData.java
 *
 * Class for all data of an RIF file.
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

package adaptor.General;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import org.apache.log4j.Logger;

/**
 * class RIFData contains all data of an RIF file. It uses
 * the classes FileDescriptor, RegionDescriptor, and CallDescriptor.
 *
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */
public class RIFData
{

    /**
     * Suffix that should be used for region information files.
     */
    public static final String SUFFIX = "rif";

    /**
     * Logger for the class RIFData.
     */
    private static Logger logger = Logger.getLogger( RIFData.class );

    /**
     * This object contains all file descriptors of the RIF file.
     */
    private FileTable  myFileData = null;

    /**
     * All region descriptors of the RIF file will be stored in this object.
     * But careful: myRegions[i] contains the region with region id i+1
     * as we do not have any region 0 at this time.
     */
    private RegionTable myRegionData = null;

    /**
     * All call descriptors of the RIF file.
     */
    private CallDescriptor[] myCalls = null;

    /**
     * This is the name of the RIF file.
     */
    private String filenameRIF = "";

    /**TODO: adaptor: enter comment!
     * constructor for RIF data read from a file.
     *
     * @param filename is the name of the RIF file
     */
    public RIFData( String filename )
    {

        this.filenameRIF = filename;

        int lineNr = 0;

        try
        {

            String currentLine;
            int noFiles;
            int noRegions;
            int noCalls;
            int i;

            // open the RIF file

            FileReader file = new FileReader( filename );
            BufferedReader buff = new BufferedReader( file );

            // first line:  Files=<NoFiles>

            currentLine = buff.readLine();
            lineNr++;

            noFiles = Utilities.readInt( currentLine, "Files" );

            logger.info( "Files = " + noFiles );

            myFileData = new FileTable( noFiles );

            for ( i = 0; i < noFiles; i++ )
            {

                currentLine = buff.readLine();
                lineNr++;

                try
                {

                    myFileData.defineFile( i, new FileDescriptor( currentLine ) );

                }
                catch ( IOException e )
                {

                    String msg = "Line " + lineNr
                                 + " (File " + ( i + 1 )
                                 + " of " + noFiles + ") : "
                                 + e.getMessage();

                    logger.error( msg );
                }

            }

            // for regions:  Regions=<NoRegions>

            currentLine = buff.readLine();
            lineNr++;
            noRegions = Utilities.readInt( currentLine, "Regions" );

            logger.info( "Regions = " + noRegions );

            myRegionData = new RegionTable( noRegions );

            for ( i = 0; i < noRegions; i++ )
            {

                currentLine = buff.readLine();
                lineNr++;

                try
                {

                    RegionDescriptor region = makeRegionDescriptor( currentLine );

                    myRegionData.defineRegion( i, region );

                }
                catch ( RuntimeException e )
                {

                    String msg = "Line " + lineNr
                                 + " (Region " + ( i + 1 )
                                 + " of " + noRegions + ") : "
                                 + e.getMessage();

                    logger.error( msg );
                }
            }

            // for calls:  Calls=<NoCalls>

            currentLine = buff.readLine();
            lineNr++;

            try
            {

                noCalls = Utilities.readInt( currentLine, "Calls" );

            }
            catch ( RuntimeException e )
            {

                logger.error( e.getMessage() );
                noCalls = 0;
            }

            logger.info( "Calls = " + noCalls );

            myCalls = new CallDescriptor[noCalls];

            for ( i = 0; i < noCalls; i++ )
            {

                currentLine = buff.readLine();
                lineNr++;

                try
                {

                    myCalls[i] = new CallDescriptor( currentLine, myRegionData );

                }
                catch ( RuntimeException e )
                {

                    String msg = "Error in line " + lineNr + ": reading call " + ( i + 1 );

                    logger.error( msg, e );
                }

            }

            buff.close();

        }
        catch ( FileNotFoundException e )
        {

            // we make a fatal error here

            logger.fatal( e.getMessage(), e );

            System.exit( 0 );

        }
        catch ( IOException e )
        {

            // add the filename to the error message

            String msg = "RIFData of file " + filename
                         + " : " + e.getMessage();
            // we make a fatal error for any IOException

            logger.fatal( msg, e );

            System.exit( 0 );
        }

        logger.info( "RIFData from file " + filename + " read." );

    } // RIFData


    /**
     * RIF data might also be generated from file and region descriptors.
     *
     * @param fileTable is an array with all file descriptors
     * @param regionTable is an array with all region descriptors
     */
    public RIFData( FileDescriptor[] fileTable, RegionDescriptor[] regionTable )
    {

        myFileData   = new FileTable( fileTable );
        myRegionData = new RegionTable( regionTable );
        myCalls      = new CallDescriptor[0];

    }

    /**
     * RIF data might also be generated from file and region descriptors.
     *
     * @param fileTable is an array with all file descriptors
     * @param regionTable is an array with all region descriptors
     */
    public RIFData( FileTable fileTable, RegionTable regionTable )
    {

        myFileData   = fileTable;
        myRegionData = regionTable;
        myCalls      = new CallDescriptor[0];

    }

    /**
     * Help routine to make a region descriptor from an input line of the RIF file.
     *
     * @param inputLine is an input line from the RIF file
     * @throws IOException in case of IO problems.
     * @return a new region descriptor for data of input line.
     */

    private RegionDescriptor makeRegionDescriptor( String inputLine ) throws IOException
    {

        String inputVal;

        int pos;

        // may be possible that we have not read a line

        if ( inputLine == null )
        {

            logger.error( "expected line for region" );
            throw new IOException( "not enough lines for regions" );
        }

        // region=15 file=3 lines=211:255 class=CLASS name=REGION kind=0
        // enable=1 calldepth=-1 dataprof=0

        int regionId = Utilities.readInt( inputLine, "region" );

        int fileId = Utilities.readInt( inputLine, "file" );

        FileDescriptor file = myFileData.getFile( fileId );

        inputVal = Utilities.readString( inputLine, "lines" );

        pos = inputVal.indexOf( RegionDescriptor.LINE_SEPARATOR );
        int lineStart = Integer.parseInt( inputVal.substring( 0, pos ) );
        int lineStop = Integer.parseInt( inputVal.substring( pos + 1 ) );

        String className = Utilities.readString( inputLine, "class" );
        String regionName = Utilities.readString( inputLine, "name" );

        int regionKind = Utilities.readInt( inputLine, "kind" );

        RegionDescriptor region = new RegionDescriptor( regionId, regionName, className, regionKind, file, lineStart, lineStop );

        int enable = Utilities.readInt( inputLine, "enable", 1 );
        int depth  = Utilities.readInt( inputLine, "depth", -1 );

        region.setProfEnabled( enable != 0 );
        region.setNestProfiling( depth );
        region.setDataEnabled( enable > 1 );

        return region;

    } // RegionDescriptor

    /**
     * Get the name of the file from which the RIF data comes.
     *
     * @return array of all file descriptors
     */
    public String getFilename()
    {

        return filenameRIF;
    }

    /**
     * Get the file descriptor from the table at a certain position.
     *
     * @param index is index position
     * @return the file descriptor at the given index position.
     */
    public FileDescriptor getFile( int index )
    {

        return myFileData.getFile( index );
    }

    /**
     * This routine returns the maximal index of a defined file descriptor.
     *
     * @return highest index of a file descriptor.
     */
    public int noFiles()
    {

        return myFileData.noFiles();
    }
    /**
     * Getter routine for a region at a given position in RIF data.
     *
     * @param index is the index position
     * @return array of all region descriptors
     */
    public RegionDescriptor getRegion( int index )
    {

        return myRegionData.getRegion( index );
    }

    /**
     * This routine returns the maximal index of a defined region descriptor.
     *
     * @return highest index of a region descriptor.
     */
    public int noRegions()
    {

        return myRegionData.noRegions();
    }

    /**
     * get all calls of RIF data.
     *
     * @return array of all call descriptors
     */
    public CallDescriptor[] getCalls()
    {

        return myCalls;
    }

    /**
     * write the RIF data back to the file.
     */
    public void writeRIF()
    {

        logger.info( "write RIFData back to file " + filenameRIF );

        try
        {

            String currentLine;
            int i;

            int noFiles   = myFileData.noFiles();
            int noRegions = myRegionData.noRegions();
            int noCalls   = myCalls.length;

            FileWriter file = new FileWriter( filenameRIF );
            BufferedWriter buff = new BufferedWriter( file );

            currentLine = "Files=" + noFiles;
            buff.write( currentLine );
            buff.newLine();

            for ( i = 0; i < noFiles; i++ )
            {
                FileDescriptor f = myFileData.getFile( i );

                if ( f == null )
                {
                    currentLine = FileDescriptor.makeDummyString( i );
                }
                else
                {
                    currentLine = f.makeFileString();
                }

                buff.write( currentLine );
                buff.newLine();
            }

            currentLine = "Regions=" + noRegions;
            buff.write( currentLine );
            buff.newLine();

            for ( i = 0; i < noRegions; i++ )
            {
                currentLine = myRegionData.getRegion( i ).makeRegionString();
                buff.write( currentLine );
                buff.newLine();
            }

            currentLine = "Calls=" + noCalls;
            buff.write( currentLine );
            buff.newLine();

            for ( i = 0; i < noCalls; i++ )
            {
                currentLine = myCalls[i].makeCallString();
                buff.write( currentLine );
                buff.newLine();
            }

            buff.close();

        }
        catch ( IOException e )
        {

            logger.fatal( e.getMessage(), e );
        }

        logger.info( "RIFData back written back to file " + filenameRIF );

    } // writeRIF

    /**
     * write the RIF data back to the file.
     * @param filename is the name of the file that will be written.
     */
    public void writeRIF( String filename )
    {

        filenameRIF = filename;
        writeRIF();

    }

    /**
     * Getter routine for the RegionData of the RIF data.
     *
     * @return the RegionData of the RIF data
     */
    public RegionTable getRegionData()
    {

        return myRegionData;
    }

    /**
     * Getter routine for the FileData of the RIF data.
     *
     * @return the FileData of the RIF data
     */
    public FileTable getFileData()
    {

        return myFileData;
    }

} // class RIFData

