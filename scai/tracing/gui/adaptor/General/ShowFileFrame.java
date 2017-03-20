/*
 * ShowFileFrame.java
 *
 * Frame to display the content of file
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

import javax.swing.JFrame;

import org.apache.log4j.Logger;


/**
 * Frame to display a file with line numbers.
 *
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */
public class ShowFileFrame extends JFrame
{

    /**
     * Logger for this class.
     */
    private static Logger logger = Logger.getLogger( ShowFileFrame.class );

    /**
     * Variable for the panel to show the file.
     */
    private ShowFilePanel myFilePanel;

    /**
     * This is the constructor for this class and has no special arguments.
     *
     */
    public ShowFileFrame()
    {

        super( "File Display" );

        myFilePanel = new ShowFilePanel();
        setContentPane( myFilePanel );
        pack();

    } // constructor ShowFileFrame

    /**
     * This routine shows a file in this frame.
     *
     * @param fileDSP is the descriptor of the file to display.
     * @param start is the first line to highlight
     * @param stop is the last line to highlight
     */
    public void showFile( FileDescriptor fileDSP, int start, int stop )
    {

        myFilePanel.setFileDescriptor( fileDSP, start, stop );

        logger.info( "show File " + fileDSP.getShortFileName() + ", lines = " + start + ":" + stop );

        setVisible( true );
    }

} // class ShowFileFrame

