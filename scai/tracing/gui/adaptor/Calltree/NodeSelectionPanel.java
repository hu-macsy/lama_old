/*
 * NodeSelectionFrame.java
 *
 * Frame that is used to select nodes for the calltree.
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
import java.awt.TextField;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.border.Border;

import org.apache.log4j.Logger;

/**
 * Frame that is used to select nodes for the calltree.
 *
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */
public class NodeSelectionPanel extends JPanel implements ActionListener
{

    /**
     * This array contains short strings describing the different kind of groups that
     * are supported for the node selection.
     */
    private static final String[] GROUP_STRINGS = { "none", "class", "file", "dir-1", "dir-2" };

    /**
     * Logger for this class.
     */
    private static Logger logger = Logger.getLogger( NodeSelectionPanel.class );

    /**
     * The group for the nodes can be selected by a Combo Box.
     */
    private JComboBox groupBox = null;

    /**
     * The only root flag can be set by a CheckBox.
     */
    private JCheckBox rootBox = null;

    /**
     * The leaf root flag can be set by a CheckBox. Value of this box will be
     * taken when OK or Apply is selected.
     */
    private JCheckBox leafBox = null;

    /**
     * Label field to show the input search pattern.
     */
    private JLabel myInputLabel;

    /**
     * TextField that contains the pattern for which nodes are searched.
     */
    private TextField mySearchPattern;

    /**
     * CheckBox to select if environment of selected nodes will be enabled.
     */
    private JCheckBox myEnvironmentBox;

    /**
     * CheckBox to enable a search pattern.
     */
    private JCheckBox myPatternBox;

    /**
     * This is the reference to the Calltree functionality.
     */
    private CTInterface myCalltree = null;

    /**
     * Constructor for setting up the node selection frame.
     *
     * @param calltree is the reference to the Calltree
     * @param group is the kind of group that should be used
     * @param rootFlag is the flag to select only root nodes
     * @param leafFlag is the flag to select only leaf nodes
     */
    public NodeSelectionPanel( CTInterface calltree, int group, boolean rootFlag, boolean leafFlag )
    {

        super();

        this.myCalltree = calltree;

        setLayout( new GridLayout( 0, 1, 5, 5 ) );

        // SELECTION 1: the Group Combo Box

        groupBox = new JComboBox();
        add( groupBox );

        for ( int i = 0; i < GROUP_STRINGS.length; i++ )
        {
            groupBox.addItem( "Group = " + GROUP_STRINGS[i] );
        }

        groupBox.setSelectedIndex( group );


        // SELECTION 2: the Root/Leaf Only Flags

        rootBox = new JCheckBox( "root only", null, rootFlag );
        rootBox.setSelected( rootFlag );
        add( rootBox );

        leafBox = new JCheckBox( "leaf only", null, leafFlag );
        leafBox.setSelected( leafFlag );
        add( leafBox );

        // SELECTION 4: the Call Environment Flag

        myEnvironmentBox = new JCheckBox( "+ call environemnt", null, false );
        add( myEnvironmentBox );

        // CheckBox: search pattern (enable/disable)

        myPatternBox = new JCheckBox( "search pattern", null, false );
        myPatternBox.addActionListener( this );

        add( myPatternBox );

        // Label: Input search pattern

        myInputLabel = new JLabel( "Input search pattern:" );
        myInputLabel.setVisible( false );
        add( myInputLabel );

        // Text Field for Input Pattern

        final int fieldSize = 20;

        mySearchPattern = new TextField( fieldSize );
        mySearchPattern.setText( "*" );
        mySearchPattern.setEditable( true );
        mySearchPattern.addActionListener( this );
        mySearchPattern.setVisible( false );

        add( mySearchPattern );

        Border myBorder = BorderFactory.createLineBorder( Color.BLACK );
        myBorder = BorderFactory.createTitledBorder( myBorder, "Node Selection" );
        setBorder( myBorder );
    }

    /**
     * This routine return the current search pattern. The UNIX
     * style is translated to the JAVA style for pattern.
     *
     * @return pattern string as used in Java
     */
    public String getPattern()
    {

        String outputPattern = null;

        if ( myPatternBox.isSelected() )
        {

            String inputPattern = mySearchPattern.getText();
            outputPattern = "";

            for ( int i = 0; i < inputPattern.length(); i++ )
            {

                char c = inputPattern.charAt( i );

                if ( c == '*' )
                {

                    outputPattern += "\\p{Graph}*";

                }
                else if ( c == '.' )
                {

                    outputPattern += "\\.";

                }
                else
                {

                    outputPattern += c;
                }
            }

        }

        return outputPattern;

    } // getPattern


    /**
     * {@inheritDoc}
     *
     * @see java.awt.event.ActionListener#actionPerformed(java.awt.event.ActionEvent)
     */
    public void actionPerformed( ActionEvent e )
    {

        Object actionSource = e.getSource();

        if ( actionSource instanceof JButton )
        {

            // OK, Cancel, Apply

            JButton actionButton = ( JButton ) actionSource;

            logger.info( "actionButton pressed " + actionButton.getText() );

        }
        else if ( actionSource instanceof JCheckBox )
        {

            JCheckBox actionBox = ( JCheckBox ) actionSource;

            if ( actionBox.equals( myPatternBox ) )
            {

                boolean isSelected = myPatternBox.isSelected();

                mySearchPattern.setVisible( isSelected );
                myInputLabel.setVisible( isSelected );

            }

        }
        else if ( actionSource.equals( mySearchPattern ) )
        {

            // <CR> key in TextField

            logger.info( "Search pattern is : " + mySearchPattern.getText() );

            myCalltree.recalc();
        }

    }

    /**
     * This routine returns the array with shorts strings for the different
     * kinds of how to build groups.
     *
     * @return array with short strings for the group kinds
     */
    public static String[] getGroupItems()
    {

        return GROUP_STRINGS;

    } // getGroupItems

    /**
     * Get the actual setting of the root flag.
     *
     * @return boolean value of the root flag.
     */
    public boolean getRootFlag()
    {

        return rootBox.isSelected();
    }

    /**
     * Get the actual setting of the leaf flag.
     *
     * @return boolean value of the leaf flag.
     */
    public boolean getLeafFlag()
    {

        return leafBox.isSelected();
    }

    /**
     * Get the current selection of the group.
     *
     * @return the group index as an int value
     */
    public int getSelectedGroup()
    {

        return groupBox.getSelectedIndex();
    }

    /**
     * Get the current flag whether neighbors are enabled.
     *
     * @return current setting of the flag for call environment
     */
    public boolean getNeighborsEnabled()
    {

        return myEnvironmentBox.isSelected();
    }

} // NodeSelectionFrame
