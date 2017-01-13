/*
 * Addr2Line.java
 *
 * Utility class that is used to translate binary addresses in the executable
 * to routine names with information about file location.
 *
 * Created: 2007-10-20 Thomas Brandes <thomas.brandes@scai.fraunhofer.de>
 * Changed:
 *
 * $Id$
 *
 * Copyright (C) 2007 Fraunhofer SCAI, Germany
 *
 * All rights reserved
 *
 * http://www.scai.fhg.de/EP-CACHE/adaptor
 */

package adaptor.General;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import org.apache.log4j.Logger;

/***
 * Utility class that is used to translate binary addresses in the executable
 * to routine names with information about file location.
 *
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */
public final class Addr2Line
{

    /**
     * Use this name as input file for "addrline".
     */
    private static final String INFILE_NAME = "addr2line.in";

    /**
     * Use this name as output file for "addrline".
     */
    private static final String OUTFILE_NAME = "addr2line.out";

    /**
     * Logger variable for this class.
     */
    private static Logger logger = Logger.getLogger( Addr2Line.class );

    /**
     * Variable will contain the name of the program executable.
     */
    private static String theExecutable = "";

    /**
     * overwrite default constructor.
     */
    private Addr2Line()
    {
    }

    /**
     * This routine is used to set the name of the executable. The routine
     * might search it in the path and append a suffix.
     *
     * @param name is the (short) name of the executable.
     * @return true if the executable has been found.
     */
    private static boolean setExecutable( String name )
    {

        theExecutable = name;

        File f = new File( theExecutable );

        // ready if we can access the executable

        if ( f.exists() )
        {

            return true;
        }

        // on Windows we give it a try with the suffix exe

        String osName = System.getProperty( "os.name" );

        logger.info( "os.name = " + osName );

        // It might be useful to add the suffix .exe for Windows

        if ( System.getProperty( "os.name" ).startsWith( "Win" ) )
        {

            logger.info( "Executable " + name + " not found, try " + name + ".exe" );

            theExecutable = name + ".exe";

            f = new File( theExecutable );

            if ( f.exists() )
            {

                return true;
            }
        }

        return false;
    }


    /**

     */

    /**
     * This routine translates region addresses to region names.
     *
     * @param addresses is an array of addresses
     * @param filenames will contain the name of the source file for location addresses[i]
     * @param lines [i] will contain the line number for location addresses[i]
     */
    private static void translate( String[] addresses, String[] filenames, int[] lines )
    {

        // ToDo: find the executable if it is not in the current dir

        try
        {

            BufferedWriter out = new BufferedWriter( new FileWriter( "addr2line.in" ) );

            for ( String addr : addresses )
            {

                out.write( addr );
                out.newLine();

            }

            out.close();

        }
        catch ( IOException e )
        {

            logger.error( e.getMessage() );

        }

        logger.info( "have written " + addresses.length + " records in file " + INFILE_NAME );

        // run the executable addr2line

        try
        {

            // option -e <executable>
            // option -f Show function names
            // -C        Demangle function names

            String cmd = "addr2line -e " + theExecutable + " -f -C";

            Runtime myRuntime = Runtime.getRuntime();

            logger.info( "exec: " + cmd );

            Process layoutProcess = myRuntime.exec( cmd );

            int rc = ProcessUtility.handleProcess( layoutProcess, INFILE_NAME, OUTFILE_NAME );

            logger.info( "addr2line ready, result = " + rc );

        }
        catch ( IOException e )
        {

            logger.error( "exec has exception: " + e.getMessage() );

        }
        catch ( RuntimeException e )
        {

            logger.error( "addr2line terminated with runtime exception: " + e.getMessage() );
        }

        // read the output file with names and files

        try
        {

            FileReader file = new FileReader( OUTFILE_NAME );

            BufferedReader buff = new BufferedReader( file );

            for ( int i = 0; i < addresses.length; i++ )
            {

                String line = buff.readLine();

                // line is the name of the region

                addresses[i] = line;

                line = buff.readLine();

                // <filename>:<line_nr>

                int pos = line.lastIndexOf( ":" );

                filenames[i] = line.substring( 0, pos );

                lines[i] = Integer.parseInt( line.substring( pos + 1 ) );

                logger.debug( "line2: File = " + filenames[i] + " line = " + lines[i] );
            }

            buff.close();

        }
        catch ( IOException e )
        {

            logger.error( e.getMessage() );

        }

    }

    /**
     * This routine translates region addresses to region names.
     * We know that we have a region address if the filename is the executable.
     *
     * @param myRegionTable is the table with all regions.
     * @param myFileTable is the table of all files that might be updated with the new source files.
     * @param executable is the name of the executable where we will find the binary addresses.
     */
    public static void addr2line( RegionTable myRegionTable, FileTable myFileTable, String executable )
    {

        int i;
        int k;

        String executableName = new File( executable ).getName();

        int noRegions = 0; // counter for address translations needed

        // Step 1: count the number of translations that are needed

        for ( i = 0; i < myRegionTable.noRegions(); i++ )
        {

            RegionDescriptor region = myRegionTable.getRegion( i );

            if ( region == null )
            {

                continue;
            }

            if ( !region.getFile().getShortFileName().equals( executableName ) )
            {

                continue;
            }

            if ( region.getRegionKind().equals( RegionDescriptor.RegionKind.CYGFUNC ) )
            {

                noRegions += 1;

            }
            else if ( region.getRegionKind().equals( RegionDescriptor.RegionKind.CYGUSER ) )
            {

                noRegions += 2;
            }
        }

        logger.info( "addr2line will get " + noRegions + " entries to translate" );

        // we are ready if there are no __cyg profiled functions

        if ( noRegions == 0 )
        {

            return;
        }

        String[] regionNames = new String[noRegions];

        k = 0;

        final String hexPrefix = "0x";

        for ( i = 0; i < myRegionTable.noRegions(); i++ )
        {

            RegionDescriptor region = myRegionTable.getRegion( i );

            if ( region == null )
            {

                continue;
            }

            if ( !region.getFile().getShortFileName().equals( executableName ) )
            {

                continue;
            }

            if ( region.getRegionKind().equals( RegionDescriptor.RegionKind.CYGFUNC ) )
            {

                regionNames[k++] = hexPrefix + Integer.toHexString( region.getFirstLine() );

            }
            else if ( region.getRegionKind().equals( RegionDescriptor.RegionKind.CYGUSER ) )
            {

                regionNames[k++] = hexPrefix + Integer.toHexString( region.getFirstLine() );
                regionNames[k++] = hexPrefix + Integer.toHexString( region.getLastLine() );
            }

        }

        // target arrays for address translation get same size as input array region_names

        String[] fileNames = new String[noRegions];

        int[] lines = new int[noRegions];

        // generate a file with all region addresses

        boolean done = setExecutable( executable );

        if ( !done )
        {

            logger.error( "could not find the executable " + executable );
            return;
        }

        try
        {

            translate( regionNames, fileNames, lines );

        }
        catch ( RuntimeException e )
        {

            logger.error( "Addr2Line fails: " + e.getMessage() );
            return;
        }

        k = 0;

        for ( i = 0; i < myRegionTable.noRegions(); i++ )
        {

            RegionDescriptor region = myRegionTable.getRegion( i );

            if ( region == null )
            {

                continue;
            }

            if ( !region.getFile().getShortFileName().equals( executableName ) )
            {

                continue;
            }

            if ( region.getRegionKind().equals( RegionDescriptor.RegionKind.CYGFUNC ) )
            {

                region.setName( regionNames[k] );

                FileDescriptor file = myFileTable.getFile( fileNames[k] );

                region.setFileInfo( file, lines[k], lines[k] );

                logger.info( "cyg func region " + k + " = " + regionNames[k] + " in " +
                             fileNames[k] + "(fileId = " + file.getFileId() + "):" + lines[k] );

                k++;

            }
            else if ( region.getRegionKind().equals( RegionDescriptor.RegionKind.CYGUSER ) )
            {

                // Attention: we do not update the name of the user region

                // FileDescriptor file = new FileDescriptor(i, fileNames[i]);

                // myFileTable.defineFile(i, file);

                FileDescriptor file = myFileTable.getFile( fileNames[k] );

                region.setFileInfo( file, lines[k], lines[k + 1] );

                logger.info( "cyg user region " + region.getName() + ": " + fileNames[k] +
                             "(fileId = " + file.getFileId() + "):" + lines[k] + "-" + lines[k + 1] );

                k++;
                k++;

            }
        }

    }

} // class Addr2Line


