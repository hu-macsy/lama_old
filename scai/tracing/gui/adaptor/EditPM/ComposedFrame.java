/*
 * EditRIFInterface.java
 *
 * Interface that has to be implemented by the RIF Editor.
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

package adaptor.EditPM;

import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.KeyStroke;
import javax.swing.SwingConstants;

import org.apache.log4j.Logger;

/**
 * The class ComposedFrame provides a frame to define composed events.
 *
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */
public class ComposedFrame extends JFrame implements ActionListener
{

    /**
     * Size of the text field for the name of the composed event.
     */
    private static final int TEXT_FIELD_SIZE = 20;

    /**
     * Logger for the class ComposedFrame.
     */
    private static Logger logger = Logger.getLogger( ComposedFrame.class );

    /**
     * This text field holds the name of the composed event.
     */
    private JTextField name;

    /**
     * This combo box is used to select the first basic event.
     */
    private JComboBox event1;

    /**
     * This combo box is used to select the second basic event.
     */
    private JComboBox event2;

    /**
     * This menu bar is used to select the binary operation.
     */
    private JMenuBar opMenuBar;

    /**
     * This label contains the currently selected string for the
     * binary operation.
     */
    private JLabel opLabel;

    /**
     * The switch button switches the two basic events.
     */
    private JButton mySwitchButton;

    /**
     * The define button defines a new composed event with the
     * current selections.
     */
    private JButton myDefineButton;

    /**
     * Pointer back to the configuration to get info about supported events.
     */
    private PMConfiguration myPMConfiguration;

    /**
     * The variable switchPos is 0 or 1 and indicates whether event1 or event2 can be set next.
     */
    private int switchPos = 0;

    /**
     * Constructor sets up the frame to define composed events.
     *
     * @param config is pointer to the PM configuration
     */
    ComposedFrame( PMConfiguration config )
    {

        super( "Composed Event" );

        myPMConfiguration = config;

        JPanel pane = new JPanel();

        pane.setLayout( new GridLayout( 5, 2, 10, 10 ) );

        // choose color for the pane that fits color of ComposedEvents

        pane.setBackground( ComposedEvent.getEventColor() );

        JLabel label;

        name = new JTextField( TEXT_FIELD_SIZE );

        mySwitchButton = new JButton( "Switch" );
        mySwitchButton.setToolTipText( "switches Event 1 and Event 2" );
        myDefineButton = new JButton( "Define" );

        String[] events = myPMConfiguration.getBasicEventStrings();

        event1 = new JComboBox();
        event1.setBackground( BasicEvent.getEventColor() );
        event2 = new JComboBox();
        event2.setBackground( BasicEvent.getEventColor() );

        myDefineButton.addActionListener( this );
        mySwitchButton.addActionListener( this );

        // we make a menu for the possible operations to set ToolTips

        opMenuBar = new JMenuBar();
        opMenuBar.setBorderPainted( true );

        JMenu opMenu = new JMenu( "Op" );

        String[] binaryOperations = { "+", "-", "*", "/", "#", "~" };

        for ( int i = 0; i < binaryOperations.length; i++ )
        {

            String binop = binaryOperations[i];

            JMenuItem opMenuItem = new JMenuItem( binop );

            opMenuItem.setToolTipText( "Event1 " + binop + " Event2" );
            opMenuItem.addActionListener( this );
            opMenuItem.setAccelerator( KeyStroke.getKeyStroke( binop ) );
            opMenu.add( opMenuItem );
        }

        opMenuBar.add( opMenu );

        for ( int i = 0; i < events.length; i++ )
        {

            event1.addItem( events[i] );
            event2.addItem( events[i] );

        }

        label = new JLabel( "New Event", SwingConstants.RIGHT );
        pane.add( label );
        pane.add( name );

        label = new JLabel( "Event 1", SwingConstants.RIGHT );
        pane.add( label );
        pane.add( event1 );

        opLabel = new JLabel( "+", SwingConstants.LEFT );
        pane.add( opMenuBar );
        pane.add( opLabel );

        label = new JLabel( "Event 2", SwingConstants.RIGHT );
        pane.add( label );
        pane.add( event2 );

        pane.add( myDefineButton );
        pane.add( mySwitchButton );

        setContentPane( pane );
        pack();

        switchPos = 0;

        logger.info( "ComposedFrame constructted" );

    } // ComposedFrame

    /**
     * {@inheritDoc}
     *
     * @see java.awt.event.ActionListener#actionPerformed(java.awt.event.ActionEvent)
     */
    public void actionPerformed( ActionEvent e )
    {

        Object actionSource = e.getSource();

        if ( actionSource instanceof JMenuItem )
        {

            JMenuItem item = ( JMenuItem ) actionSource;

            logger.debug( "action on menu item " + item.getText() );

            opLabel.setText( item.getText() );
            opLabel.setToolTipText( item.getToolTipText() );

        }

        if ( actionSource instanceof JButton )
        {

            JButton actionButton;

            actionButton = ( JButton ) actionSource;

            logger.debug( "action on button " + actionButton.getText() );

            if ( actionButton == mySwitchButton )
            {

                int index1 = event1.getSelectedIndex();
                event1.setSelectedIndex( event2.getSelectedIndex() );
                event2.setSelectedIndex( index1 );
            }

            if ( actionButton == myDefineButton )
            {

                PerformanceEvent p1;
                PerformanceEvent p2;

                String newName;

                boolean error = false;

                p1 = myPMConfiguration.findSupportedEvent( ( String ) event1.getSelectedItem() );
                p2 = myPMConfiguration.findSupportedEvent( ( String ) event2.getSelectedItem() );

                newName = name.getText();

                // make a preliminary definition

                ComposedEvent event = new ComposedEvent( newName, p1, opLabel.getText(), p2 );

                // try to add it to the table of events

                try
                {

                    myPMConfiguration.addNewEvent( event );

                }
                catch ( RuntimeException ex )
                {

                    // show error message

                    JOptionPane.showMessageDialog( null, ex.getMessage(), "Error for composed event", JOptionPane.ERROR_MESSAGE );

                    error = true; // in case of error we do not close

                } // catch

                if ( !error )
                {

                    // hide this frame if composed event has been defined successfully

                    setVisible( false );

                }
            }
        }

    } // actionPerformed

    // method is called when Event Button has been selected while
    // this ComposedFrame is open

    /**
     * This routine sets a basic event in the frame for composed events.
     * It becomes one of the events for the binary operation.
     *
     * @param event is the basic event taken for definition
     */
    public void setEvent( BasicEvent event )
    {

        int index = myPMConfiguration.getEventIndex( event );

        if ( index >= 0 )
        {

            // so we have a valid index

            if ( switchPos == 0 )
            {

                event1.setSelectedIndex( index );

            }
            else
            {

                event2.setSelectedIndex( index );
            }

            switchPos = 1 - switchPos;
        }
    }

    /**
     * This routine takes an already defined composed event to
     * fill the components of this editor.
     *
     * @param event is an already defined composed event.
     */
    public void editEvent( ComposedEvent event )
    {

        event1.setSelectedIndex( myPMConfiguration.getEventIndex( event.getP1() ) );
        event2.setSelectedIndex( myPMConfiguration.getEventIndex( event.getP2() ) );
        opLabel.setText( event.getOp() );
        name.setText( event.getName() );

    }
}

