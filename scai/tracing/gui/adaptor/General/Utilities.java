/*
 * Utilities.java
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

package adaptor.General;

import java.io.IOException;

import org.apache.log4j.Logger;

/***
 * Utilities is a class of methods used in other routines.
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */
public final class Utilities
{

    /**
     * Logger variable for this utility class.
     */
    private static Logger logger = Logger.getLogger( Utilities.class );

    /**
     * overwrite default constructor.
     */
    private Utilities()
    {
    }

    /**
     * checks for "name=value" in inputLine.
     *
     * @param inputLine will be searched for entry 'name=...'
     * @param name is the reserved word
     * @return value specified for name
     * @throws IOException if entry "name=..." is not in the inputLine
     */
    public static String readString( String inputLine, String name ) throws IOException
    {

        if ( inputLine == null )
        {

            String msg = "cannot find entry <" + name + "> on null input line";

            logger.debug( msg );

            throw new IOException( msg );

        }

        int pos1 = inputLine.indexOf( name + "=" );

        if ( pos1 < 0 )
        {

            String msg = "entry <" + name + "> not available";

            // this is not very serious, some entries might not be available

            // logger.debug(msg);

            throw new IOException( msg );

        }

        pos1 += name.length() + 1;

        int pos2 = inputLine.indexOf( " ", pos1 );

        if ( pos2 < 0 )
        {

            pos2 = inputLine.length();

        }

        return inputLine.substring( pos1, pos2 );

    } // ReadString

    /**
     * checks for "name=val" in Line where val is an int.
     *
     * @param inputLine xxx
     * @param name xxx
     * @return int val (may be 0 if line is not of this kind)
     * @throws IOException if val is not an integer
     */
    public static int readInt( String inputLine, String name ) throws IOException
    {

        // note that readString might throw already IOException

        String value = readString( inputLine, name );

        int val = 0;

        try
        {

            val = Integer.parseInt( value );

        }
        catch ( NumberFormatException e )
        {

            // this is a serious error, report it here

            String msg = "int expected here: " + name + "=" + value;
            msg += "\n" + "inputLine: " + inputLine;
            logger.error( msg );

            throw new IOException( "illegal input line" );

        }

        return val;

    }

    /**
     * checks for "name=val" in Line where val is an int.
     *
     * @param inputLine is line
     * @param name is searched valued
     * @param defaultValue to be returned by exception
     * @return val may be defaultValue if line is not of this kind
     */
    public static int readInt( String inputLine, String name, int defaultValue )
    {

        int outValue;

        try
        {

            outValue = readInt( inputLine, name );

        }
        catch ( IOException e )
        {

            outValue = defaultValue;

        }

        return outValue;

    }

    /**
     * Help routine to convert boolean value to integer value.
     *
     * @param boolVal is the input value
     * @return 1 for true, 0 for false
     */
    public static int bool2int( boolean boolVal )
    {

        int val = 0;

        if ( boolVal )
        {

            val = 1;
        }

        return val;
    }

    /**
     * Help routine to convert integer to boolean (0 for false, otherwise
     * it is true).
     *
     * @param intVal is an integer value
     * @return false if intVal equal zero, true otherwise
     */
    public static boolean int2bool( int intVal )
    {

        boolean val = ( intVal != 0 );

        return val;
    }


} // class Utilities

