/*
 * LayoutThread.java
 * 
 * Thread to make a graph layout.
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

import java.io.InputStream;
import java.io.OutputStream;
import java.io.FileOutputStream;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;

import java.io.FileWriter;

import java.io.BufferedWriter;
import java.io.BufferedReader;

import org.apache.log4j.Logger;

/**
 * The LayoutThread is a thread that runs the graph
 * layout program dot. As this might take a long while
 * it is important to cancel such a thread.
 *
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */
class LayoutThread extends Thread {

    /**
     * This is the name of the program that makes the layout
     * of graphs.
     */
    private static final String LAYOUT_CMD = "dot";
    
    /**
     * Logger for this class.
     */
    private static Logger logger = Logger.getLogger(LayoutThread.class);

    /**
     * This is the name of the input file.
     */
    private String inputName;
    
    /**
     * This is the name of the output file.
     */
    private String outputName;
    
    /**
     * This is the suffix for the layout output file. The
     * program dot can generate different formats.
     *  
     */
    private String suffix;
    
    /**
     * Pointer back to the Calltree editor.
     */
    private CTInterface myCT;
    
    /**
     * If this variable is true the output of the layout
     * program will be displayed after successful termination
     * of the layout process.
     */
    private boolean display;
    
    /**
     * This is the process that runs the 'dot' program. 
     * It will be destroyed if this thread is canceled.
     */
    private Process myLayoutProcess;

    /**
     * This constructor sets up an usual layout thread
     * to generate a dot file.
     * 
     * @param inputFile is the file containing the input data
     * @param main Pointer back to the editor to display the graph
     */
    LayoutThread(String inputFile, CTInterface main) {

        this.inputName = inputFile;
        this.outputName = "graph_out.dot";
        this.suffix = "";
        this.display = true;
        this.myCT = main;

    } // LayoutThread for show Graph with GRAPPA

    /**
    /**
     * This constructor sets up a layout thread
     * to generate a file in any format.
     * 
     * @param input is the file containing the input data
     * @param output is the name of the output file
     * @param suf is the suffix of output file used for conversion
     */
    LayoutThread(String input, String output, String suf) {

        this.inputName = input;
        this.outputName = output;
        this.suffix = suf;
        this.display = false;
 
    } // LayoutThread for export of the Graph

    /***********************************************************************************************
     * * void setupProcess (String suffix) * *
     **********************************************************************************************/

    /**
     * This routine sets up the layout process. Input and output
     * will be piped later.
     * 
     * @param suffix is an optional argument for specific conversions
     * @return pointer to the created process
     */
    private static Process setupProcess(String suffix) {

        String cmd = LAYOUT_CMD;

        Process layoutProcess = null;

        if (suffix.length() > 0) {
            
            // if there is a suffix we have to add a flag
            
            cmd = cmd + " -T" + suffix;
          
        }

        logger.info("setup dot process: " + cmd);

        try {

            Runtime myRuntime = Runtime.getRuntime();

            layoutProcess = myRuntime.exec(cmd);

        } catch (IOException e) {

            logger.error("exec has exception: " + e.getMessage());

        }

        logger.info("setup dot process done: P = " + layoutProcess);

        return layoutProcess;

    } // setupProcess

    /**
     * This routine provides the layout process with input data and collects
     * the output data.
     * 
     * @param layoutProcess is the created process
     * @param input is the name of the inputfile
     * @param output is the name of the outputfile
     * @return 0 for successful completion
     */
    static int handleProcess(Process layoutProcess, String input, String output) {

        final int eofChar = -1;
        
        int rc = 0; // returnCode (0 for success)

        // pipe "infile" into the output stream of process P
        // pipe input stream of process P into "outfile"

        try {

            // set for P : > output

            FileOutputStream outFile = new FileOutputStream(output);
            FileInputStream  infile  = new FileInputStream(input);

            OutputStream out = layoutProcess.getOutputStream();
            InputStream  in  = layoutProcess.getInputStream();

            final int BufferSize = 2048;

            byte[] buffer = new byte[BufferSize];
            int noOfBytes = 0;

            while ((noOfBytes = infile.read(buffer)) != -1) {
                out.write(buffer, 0, noOfBytes);
            }

            out.close();

            while ((noOfBytes = in.read(buffer)) != -1) {
                outFile.write(buffer, 0, noOfBytes);
            }

            outFile.close();

            rc = layoutProcess.waitFor();

            logger.info("P finished with code " + rc);

            // set P back to null so it will no more killed

            layoutProcess = null;

        } catch (IOException e) {

            logger.error("Error during I/O handling dot process: " + e.getMessage());

            rc = 1;

        } catch (InterruptedException e) {
            
            logger.error("process interrupted: " + e.getMessage());

            rc = 1;
            
        }

        return rc;

    } // handleProcess

    protected void readReplace(String filename, 
            String searchPattern1, String replacementPattern1,
            String searchPattern2, String replacementPattern2) throws IOException {
        String line;
        StringBuffer buffer = new StringBuffer();
        FileInputStream fileInputStream = new FileInputStream(filename);
        BufferedReader reader = new BufferedReader(
                new InputStreamReader(fileInputStream));
        while((line = reader.readLine()) != null) {
            String line1 = line.replaceAll(searchPattern1, replacementPattern1);
            String line2 = line1.replaceAll(searchPattern2, replacementPattern2);
            buffer.append(line2 + "\n");
        }
        reader.close();
        BufferedWriter out = new BufferedWriter(new FileWriter(filename));
        out.write(buffer.toString());
        out.close();
    }

    /**
     * {@inheritDoc}
     *
     * @see java.lang.Runnable#run()
     */
    public void run() {

        // generates the graph with the current selected properties

        myLayoutProcess = setupProcess(suffix);

        if (myLayoutProcess != null) {
            
            int rc = handleProcess(myLayoutProcess, inputName, outputName);

            // new version of dot requires update of output
            // widht=3.14 -> width="3.14", height=1.56 -> height="1.56"

            // sed -e "s/width=\([0-9]*\.[0-9]*\)/width=\"\1\"/g"      \
            //     -e "s/height=\([0-9]*\.[0-9]*\)/height=\"\1\"/g"      \
            //     -i $1

            String searchPattern1 = "width=(\\d*\\.\\d*)";
            String replacePattern1 = "width=\"$1\""; 
            String searchPattern2 = "height=([0-9]*\\.[0-9]*)";
            String replacePattern2 = "height=\"$1\""; 
 
            if ((rc == 0) && display) {
                
                // replace only for layout files, not for exported files like png, jpg, ...

                try {
                    readReplace( outputName, searchPattern1, replacePattern1,
                                             searchPattern2, replacePattern2 );
                } catch (IOException e) {
                    e.printStackTrace();
                }

               // in case of error return

                myCT.showGraphFile(outputName);
            }
        }
        
    } // run

    /**
     * This routine destroys the layout process and interrupts
     * the layout thread.
     * 
     */
    public void cancel() {

        logger.info("cancel layout thread");

        if (myLayoutProcess != null) {

            myLayoutProcess.destroy();
            logger.info("have killed layout process");

        }

        // now we interrupt the thread execution 

        interrupt();

    }

}
