/*
 * ProcessUtility.java
 *
 * ProcessUtility realizes a help routine for the execution of processes.
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
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.OutputStream;

import org.apache.log4j.Logger;

/***
 * Utilities is a class of methods used in other routines.
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */
public final class ProcessUtility
{

    /**
     * Logger variable for this class.
     */
    private static Logger logger = Logger.getLogger( ProcessUtility.class );


    /**
     * Constant for eofChar.
     */
    private static final int EOF_CHAR = -1;

    /**
     * This help class defines a thread that provides an input file
     * as an output stream. It becomes necessary for processes that work
     * line to line on the input file.
     *
     * @version $LastChangedRevision$
     * @author ThomasBrandes
     */
    private static class InputProvider extends Thread
    {

        /**
         * myFileName is the name of the input file and will be set by
         * the constructor.
         */

        private String myFileName = "";

        /**
         * myOutputStream will be the output stream for the process and
         * is set by the constructor.
         */
        private OutputStream myOutputStream = null;

        /**
         * Constructor for an instance that provides an input file.
         *
         * @param inFile is the name of the input file.
         * @param output is the stream that will get the content of the input file.
         */
        InputProvider( String inFile, OutputStream output )
        {

            myFileName = inFile;
            myOutputStream = output;

        }

        /**
         * {@inheritDoc}
         *
         * @see java.lang.Thread#run()
         */
        public void run()
        {

            try
            {

                boolean eof = false;

                BufferedReader inbuff = new BufferedReader( new FileReader( myFileName ) );

                DataOutputStream out = new DataOutputStream( myOutputStream );

                logger.info( "InputProvider starts piping" );

                while ( !eof )
                {

                    int nextChar = inbuff.read();

                    eof = ( nextChar == EOF_CHAR );

                    if ( !eof )
                    {

                        out.write( nextChar );

                    }

                }

                out.close();

                logger.info( "InputProvider finishes" );

            }
            catch ( IOException e )
            {

                logger.error( "IOException for input thread", e );

            }

        }
    }

    /**
     * overwrite default constructor.
     */
    private ProcessUtility()
    {
    }

    /**
     * This routine provides the layout process with input data and collects the output data.
     *
     * @param layoutProcess
     *            is the created process
     * @param input
     *            is the name of the inputfile
     * @param output
     *            is the name of the outputfile
     * @return 0 for successful completion
     */
    static int handleProcess( Process layoutProcess, String input, String output )
    {

        int rc = 0; // returnCode (0 for success)

        // pipe "infile" into the output stream of process P
        // pipe input stream of process P into "outfile"

        try
        {

            InputProvider inputThread = new InputProvider( input, layoutProcess.getOutputStream() );

            inputThread.start();

            boolean eof = false;

            // set for P : > output

            FileOutputStream outFile = new FileOutputStream( output );

            DataInputStream in = new DataInputStream( layoutProcess.getInputStream() );

            eof = false;

            while ( !eof )
            {

                int nextChar = in.read();

                eof = ( nextChar == EOF_CHAR );

                if ( !eof )
                {

                    outFile.write( nextChar );
                }

            } // while not eof

            outFile.close();

            logger.info( "Output of process P is closed" );

            if ( inputThread.isAlive() )
            {

                logger.info( "Input thread for P alive, wait for it" );

                inputThread.wait();
            }

            logger.info( "Input thread for P finished" );

            rc = layoutProcess.waitFor();

            logger.info( "P finished with code " + rc );

        }
        catch ( IOException e )
        {

            logger.error( "Error during I/O handling dot process: " + e.getMessage() );

            rc = 1;

        }
        catch ( InterruptedException e )
        {

            logger.error( "process interrupted: " + e.getMessage() );

            rc = 1;

        }

        return rc;

    } // handleProcess

}