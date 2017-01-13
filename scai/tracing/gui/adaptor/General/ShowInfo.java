/*
 * ShowInfo.java
 *
 * Utility class to display text information in a frame.
 *
 * Created: 2007-01-02 Thomas Brandes <thomas.brandes@scai.fraunhofer.de>
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

import java.awt.Font;

import javax.swing.JFrame;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;

/***
 * Utility class to display text information in a frame.
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */
public final class ShowInfo
{

    /**
     * Constant value for the initial width of text window.
     */
    private static final int TEXT_WIDTH = 50;

    /**
     * Constant value for the initial hight of text window.
     */
    private static final int TEXT_HEIGHT = 20;

    /**
     * Info frame to display infos about nodes in the CallTree.
     */
    private static JFrame myInfoFrame = null;

    /**
     * Text area used for the InfoFrame.
     */
    private static JTextArea textArea = null;

    /**
     * overwrite default constructor.
     */
    private ShowInfo()
    {
    }

    /**
     * This routine displays an information text in the associated information frame.
     *
     * @param infoText
     *            is the text to be displayed
     */
    public static void showInfoText( String infoText )
    {

        final int fontSize = 12;

        // we initialize InfoFrame and textArea the first time and reuse it later.

        if ( myInfoFrame == null )
        {

            myInfoFrame = new JFrame( "INFO" );

            textArea = new JTextArea( TEXT_HEIGHT, TEXT_WIDTH );
            textArea.setFont( new Font( "Courier", Font.PLAIN, fontSize ) );
            textArea.setText( "Info" );

            myInfoFrame.getContentPane().add( new JScrollPane( textArea ) );
            myInfoFrame.pack();
        }

        textArea.setText( infoText );
        myInfoFrame.setVisible( true );
    }

} // class ShowInfo


