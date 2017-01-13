/*
 * ShowFilePanel.java
 *
 * ShowFilePanel is a panel to display the content of a file
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

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Font;

import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;

import org.apache.log4j.Logger;

/**
 * Panel to show files and highlight a range of lines in it.
 *
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */
public class ShowFilePanel extends JPanel
{

    /**
     * Predefined value for number of text columns.
     */
    private static final int TEXT_COLUMNS = 50;

    /**
     * Predefined value for number of text rows.
     */
    private static final int TEXT_ROWS = 20;

    /**
     * Predefined value for the text font size.
     */
    private static final int TEXT_FONT_SIZE = 12;

    /**
     * Logger for this class.
     */
    private static Logger logger = Logger.getLogger( ShowFilePanel.class );

    /**
     * The fileLabel is a label for the selected file.
     */
    private JLabel fileLabel;    // Label for the selected file

    /**
     * The textArea contains the text of the selected region.
     */
    private JTextArea textArea;

    /**
     *
     * This constructor creates a panel to show files with default value.
     *
     */
    public ShowFilePanel()
    {

        textArea = new JTextArea( TEXT_ROWS, TEXT_COLUMNS );
        textArea.setFont( new Font( "Courier", Font.PLAIN, TEXT_FONT_SIZE ) );
        textArea.setText( "display of the selected region" );

        fileLabel = new JLabel( "<no file>" );

        setLayout( new BorderLayout() );

        add( "North", fileLabel );
        add( "Center", new JScrollPane( textArea ) );

    } // ShowFilePanel

    /**
     *
     * This routines sets the file to show and the first and last line of the selection.
     *
     * @param file is the descriptor of the file to show
     * @param startLine is the first line of the selection
     * @param stopLine ist the last line of the selection
     */
    public void setFileDescriptor( FileDescriptor file, int startLine, int stopLine )
    {

        try
        {

            textArea.setText( file.getContent() );

            fileLabel.setText( file.getShortFileName() );

            int startPos = file.getLineIndex( startLine - 1 );
            int stopPos = file.getLineIndex( stopLine );

            // textArea.select (start_pos, stop_pos);

            textArea.setSelectedTextColor( Color.red );
            textArea.setSelectionStart( startPos );
            textArea.setSelectionEnd( stopPos );

            logger.info( "setFileDescriptor, selected " + startLine + ":" + stopLine + " is pos " + startPos + ":" + stopPos );

        }
        catch ( IOException exception )
        {

            textArea.setText( "" );
            fileLabel.setText( exception.getMessage() );
        }

    } // setFileDescriptor

} // class ShowFilePanel
