/*
 * CallEnvironmentFrame.java
 *
 * Frame to define the call environment for the call tree
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

import java.awt.Color;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JPanel;
import javax.swing.border.Border;

import org.apache.log4j.Logger;

/**
 * Panel to define the call environment for the call tree.
 *
 * @version $LastChangedRevision$
 * @author Dr. Thomas Brandes
 */
public class CallEnvironmentPanel extends JPanel implements ActionListener
{

    /**
     * The generated serial version ID.
     */
    private static final long serialVersionUID = 8453734186917099850L;

    /**
     * Default calling depth.
     */
    private static final int DEFAULT_CALLING_DEPTH = 2;

    /**
     * Default called depth.
     */
    private static final int DEFAULT_CALLED_DEPTH = 1;

    /**
     * Number of pixel for the gaps in the layout.
     */
    private static final int GAP_SIZE = 10;

    /**
     * Logger for this class.
     */
    private static Logger logger = Logger.getLogger( CallEnvironmentPanel.class );

    /**
     * Array of entries to appear in menu for calling depth.
     */
    private static int[] callingDepthItems = { 0, 1, 2, 3, 5, 10 };

    /**
     * Array of entries to appear in menu for called depth.
     */
    private static int[] calledDepthItems = { 0, 1, 2, 3, 5, 10 };

    /**
     * Current value used to mark called routines to a certain depth.
     */
    private int callingDepth = DEFAULT_CALLING_DEPTH;

    /**
     * Current value used to mark routines calling me to a certain depth.
     */
    private int calledDepth = DEFAULT_CALLED_DEPTH; // mark routines calling me

    /**
     * Combo box for selection of calling depth.
     */
    private JComboBox myCallingBox;

    /**
     * combo box for selection of called depth.
     */
    private JComboBox myCalledByBox;

    /**
     * Button to set default values for CallEnvironment.
     */
    private JButton defaultButton;

    /**
     * Constructor to make a panel for selection of call environment.
     *
     */
    public CallEnvironmentPanel()
    {

        super();

        Border myBorder = BorderFactory.createLineBorder( Color.BLACK );

        myBorder = BorderFactory.createTitledBorder( myBorder, "Call Environment" );

        setBorder( myBorder );

        setLayout( new GridLayout( 0, 1, GAP_SIZE, GAP_SIZE ) );

        // Row 1: Calling Depth

        myCallingBox = new JComboBox();

        myCallingBox.addActionListener( this );

        for ( int i = 0; i < callingDepthItems.length; i++ )
        {

            myCallingBox.addItem( "Calling depth : " + callingDepthItems[i] );
        }

        add( myCallingBox );

        // Row 2: Called Depth

        myCalledByBox = new JComboBox();

        myCalledByBox.addActionListener( this );

        for ( int i = 0; i < calledDepthItems.length; i++ )
        {

            myCalledByBox.addItem( "Called by depth : " + calledDepthItems[i] );

        }

        add( myCalledByBox );

        // take default values to set index of the combo boxes
        // be careful: call(ing)Depth might be overwritten due to actions

        defaultEnvSelection();

        // Row 3: default Button

        defaultButton = new JButton( "Default" );

        defaultButton.addActionListener( this );

        add( defaultButton );

        logger.info( "CallEnvironmentPanel created, calling = " + callingDepth + ", called = " + calledDepth );

    }

    /**
     * search value in an array.
     * @param val is the value to be search
     * @param array is the array where we look for val
     * @return the index in the array, -1 if not found.
     */
    static int getIndex( int val, int[] array )
    {

        final int errorIndex = -1;

        int len = array.length;

        for ( int i = 0; i < len; i++ )
        {

            if ( array[i] == val )
            {
                return i;
            }
        }

        return errorIndex;

    } // getIndex

    /**
     * Gets the value for calling depth.
     *
     * @return current setting of calling depth
     */
    public int getCallingDepth()
    {

        return callingDepth;
    }

    /**
     * Gets the value for called depth.
     *
     * @return current setting of called depth
     */
    public int getCalledDepth()
    {

        return calledDepth;
    }

    /**
     * set current selection to actual values.
     */
    private void defaultEnvSelection()
    {

        callingDepth = DEFAULT_CALLING_DEPTH;
        calledDepth  = DEFAULT_CALLED_DEPTH;

        setComboBoxes();

    }

    /**
     * Sets the actual values of callingDepth and calledDepth
     * as the current selection in the combo boxes.
     */
    private void setComboBoxes()
    {

        myCallingBox.setSelectedIndex( getIndex( callingDepth, callingDepthItems ) );
        myCalledByBox.setSelectedIndex( getIndex( calledDepth, calledDepthItems ) );

    }

    /**
     * {@inheritDoc}
     *
     * @see java.awt.event.ActionListener#actionPerformed(java.awt.event.ActionEvent)
     */
    public void actionPerformed( ActionEvent e )
    {

        Object actionSource = e.getSource();

        logger.info( "action " + e );

        if ( actionSource instanceof JButton )
        {

            JButton actionButton = ( JButton ) actionSource;

            if ( actionButton.equals( defaultButton ) )
            {

                defaultEnvSelection();

            }

        } // Source instanceof Button

        if ( actionSource instanceof JComboBox )
        {

            JComboBox actionComboBox = ( JComboBox ) actionSource;

            int selectedIndex = actionComboBox.getSelectedIndex();

            if ( actionComboBox == myCallingBox )
            {

                callingDepth = callingDepthItems[selectedIndex];

                logger.info( "set callingDepth = " + callingDepth );
            }

            if ( actionComboBox == myCalledByBox )
            {

                calledDepth = calledDepthItems[selectedIndex];

                logger.info( "set calledDepth = " + calledDepth );
            }

        }

    } // actionPerformed

} // CallEnvironmentPanel