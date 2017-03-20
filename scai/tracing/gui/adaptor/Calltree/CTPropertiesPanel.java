/*
 * CTPropertiesPanel.java
 *
 * Panel to select properties for the Visualization of Calltree Data.
 *
 * Created: 2008-01-02 Thomas Brandes <thomas.brandes@scai.fraunhofer.de>
 * Changed:
 *
 * $Id$
 *
 * Copyright (C) 2008 Fraunhofer SCAI, Germany
 *
 * All rights reserved
 *
 * http://www.scai.fhg.de/EP-CACHE/adaptor
 */

package adaptor.Calltree;

import java.util.Hashtable;
import java.util.Map;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JPanel;
import javax.swing.border.Border;

import org.apache.log4j.Logger;

import adaptor.General.PropertyMenu;


/**
 * Panel to select properties for the Visualization of Calltree Data.
 *
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */
public class CTPropertiesPanel extends JPanel implements ActionListener
{

    /**
     * We give the CheckBoxes a name to identify header and title. We
     * take this string to separate header:title.
     */
    private static final String BOX_NAME_SEPARATOR = ":";

    /**
     * Each boxItem gets a label title=value.
     */
    private static final String ITEM_ASSIGN = "=";

    /**
     * Tooltip for the show menu.
     */
    private static final String SHOW_TOOL_TIP = "Select the info to be displayed in the calltree.";

    /**
     * Logger for this class.
     */
    private static Logger logger = Logger.getLogger( CTPropertiesPanel.class );

    /**
     * We need a callback object to redraw the graph after changing properties.
     */
    private CTInterface myCallback;

    /**
     * Global variable that will be set to false if changes of the combo box
     * should not make any action.
     */
    private boolean comboBoxChangedAction = true;

    /**
     * We need a mapping from a property menu to the combo box that
     * belongs to the property menu. It allows the update of all
     * properties on the panel.
     */
    private Map<String, JComboBox> allBoxes;

    /**
     * Button that will redraw the Calltree.
     */
    private JButton applyButton = new JButton( "Apply" );

    /**
     * Help routine to create a new menu with different items and tool tip
     * and adds this class as ActionListener for all entries.
     *
     * @param groupName is the group to which the box belongs.
     * @param title is the title of the new combo box
     * @param tip is the tool tip for this combo box
     * @param items are the items of this combo box
     * @param selection is the current selection for this combo box
     * @return a new menu
     */
    private JComboBox newSelectionBox( String groupName, String title, String tip, String[] items, String selection )
    {

        // CounterMenu

        String boxTitle = groupName + BOX_NAME_SEPARATOR + title;

        JComboBox box = new JComboBox();

        for ( String item : items )
        {

            box.addItem( title + ITEM_ASSIGN + item );
        }

        box.setToolTipText( tip );

        box.setName( boxTitle );

        box.setSelectedItem( title + ITEM_ASSIGN + selection );

        box.addActionListener( this );

        allBoxes.put( boxTitle, box );

        return box;

    }

    /**
     * Help routine to make a JComboBox from a PropertyMenu.
     *
     * @param groupName is the group name to which the propMenu belongs.
     * @param propMenu is the PropertyMenu
     * @return a JMenu with items of the PropertyMenu.
     */
    private JComboBox newMenu( String groupName, PropertyMenu propMenu )
    {

        String[] menuItems = propMenu.getItems();
        String   title     = propMenu.getTitle();
        String   help      = propMenu.getHelp();
        String   selection = propMenu.getSelection();

        return newSelectionBox( groupName, title, help, menuItems, selection );

    }

    /**
     * Help routine to make a panel with buttons for all Show entries.
     *
     * @return a panel with a button for each show entry.
     */
    private JPanel showButtonPanel()
    {

        final String title = "Show";

        String[] menuItems = CTProperties.getShowItems();

        JPanel buttonPanel = new JPanel();

        buttonPanel.setLayout( new GridLayout( 1, 0, 3, 3 ) );

        for ( int i = 0; i < menuItems.length; i++ )
        {

            JButton button = new JButton( menuItems[i] );
            button.setToolTipText( SHOW_TOOL_TIP );
            button.addActionListener( this );

            buttonPanel.add( button );
        }

        Border myBorder = BorderFactory.createLineBorder( Color.BLACK );
        buttonPanel.setBorder( BorderFactory.createTitledBorder( myBorder, title ) );

        return buttonPanel;
    }

    /**
     * Help routine to make panels.
     *
     * @param header is the global header name for which property menus
     * are created.
     * @return a panel with the different entries.
     */
    private JPanel newHeaderPanel( String header )
    {

        String[] titles = PropertyMenu.getAllTitles( header );

        JPanel mainPanel = new JPanel();

        mainPanel.setLayout( new GridLayout( titles.length, 1, 2, 2 ) );

        Border myBorder = BorderFactory.createLineBorder( Color.BLACK );
        mainPanel.setBorder( BorderFactory.createTitledBorder( myBorder, header ) );

        for ( int i = 0; i < titles.length; i++ )
        {

            PropertyMenu menu = PropertyMenu.findMenu( header, titles[i] );

            mainPanel.add( newMenu( header, menu ) );

        }

        return mainPanel;
    }

    /**
     * This help routine allows to set all property values in the
     * corresponding combo boxes.
     *
     */
    private void setCurrentSelections()
    {

        // disable actions on the combo boxes when system sets selections

        comboBoxChangedAction = false;

        String[] menuHeaders = PropertyMenu.getAllHeaders();

        for ( int i = 0; i < menuHeaders.length; i++ )
        {

            String header = menuHeaders[i];

            String[] menuTitles = PropertyMenu.getAllTitles( header );

            for ( int j = 0; j < menuTitles.length; j++ )
            {

                String title = menuTitles[j];

                PropertyMenu menu = PropertyMenu.findMenu( header, title );

                String boxTitle = header + BOX_NAME_SEPARATOR + title;

                JComboBox box = allBoxes.get( boxTitle );

                if ( box == null )
                {

                    logger.error( "could not find JComboBox for " + boxTitle );
                }

                String selectedItem = title + ITEM_ASSIGN + menu.getSelection();

                box.setSelectedItem( selectedItem );

                logger.debug( "JComboBox " + boxTitle + " -> " + selectedItem );

            }

        }

        // now enable the combo box action again so that user changes apply

        comboBoxChangedAction = true;

    }

    /**
     * Constructor for the panel to manage all Calltree properties.
     *
     * @param main is used for callbacks to routines on the Calltree
     * @param selectionPanels is an array of panels for the node selection
     */
    public CTPropertiesPanel( CTInterface main, JPanel[] selectionPanels )
    {

        super();

        allBoxes = new Hashtable<String, JComboBox>();

        myCallback = main;

        setLayout( new BorderLayout() );

        // North : Panel of Show Buttons

        add( "North", showButtonPanel() );

        JPanel mainPanels = new JPanel();

        mainPanels.setLayout( new GridBagLayout() );

        GridBagConstraints constraints = new GridBagConstraints();

        // general setting for the constraints

        // fill = HORIZONTAL guarantees that all panels in one column
        // (e.g. for Node + Edge) have the same size
        // default for (width, height) = (1, 1)

        constraints.fill  = GridBagConstraints.HORIZONTAL;
        constraints.gridheight = 1;
        constraints.gridwidth  = 1;
        constraints.anchor = GridBagConstraints.NORTHWEST;

        int maxY = 0;

        constraints.gridx = 0;

        for ( int i = 0; i < selectionPanels.length; i++ )
        {

            JPanel panel = selectionPanels[i];

            if ( panel != null )
            {

                constraints.gridy = maxY++;
                mainPanels.add( panel, constraints );
            }
        }

        JPanel headerPanel;

        headerPanel = newHeaderPanel( CTProperties.COSTS_MENU );
        constraints.gridx = 1;
        constraints.gridy = 0;
        mainPanels.add( headerPanel, constraints );

        headerPanel = newHeaderPanel( CTProperties.NODE_MENU );
        constraints.gridx = 2;
        constraints.gridy = 0;
        mainPanels.add( headerPanel, constraints );

        headerPanel = newHeaderPanel( CTProperties.EDGE_MENU );
        constraints.gridx = 2;
        constraints.gridy = 1;
        mainPanels.add( headerPanel, constraints );

        headerPanel = newHeaderPanel( CTProperties.SHAPE_MENU );
        constraints.gridx = 3;
        constraints.gridy = 0;
        constraints.gridheight = maxY;
        mainPanels.add( headerPanel, constraints );

        add( "Center", mainPanels );

        applyButton = new JButton( "Apply" );

        applyButton.addActionListener( this );

        add( "South", applyButton );

        setCurrentSelections();

    }

    /**
     * {@inheritDoc}
     *
     * @see java.awt.event.ActionListener#actionPerformed(java.awt.event.ActionEvent)
     */
    public void actionPerformed( ActionEvent e )
    {

        Object actionSource = e.getSource();

        String command = e.getActionCommand();

        if ( actionSource instanceof JComboBox )
        {

            if ( !comboBoxChangedAction )
            {

                return;
            }

            JComboBox box = ( JComboBox ) actionSource;

            String selection = ( String ) box.getSelectedItem();

            String[] sel = selection.split( ITEM_ASSIGN );

            String key = sel[0];
            String val = sel[1];

            logger.info( command + " action JComboBox " + box.getName() + " selection, key = " + key + ", val = " + val );

            String name = box.getName();

            sel = name.split( BOX_NAME_SEPARATOR );

            String header = sel[0];
            String title  = sel[1];

            logger.info( "header = " + header + ", titel = " + title );

            // last argument null will not redraw the graph now

            CTProperties.setProperty( val, title, header );

            // Attention: changing one property can imply changes of other properties

            setCurrentSelections();


        }
        else if ( actionSource instanceof JButton )
        {

            JButton button = ( JButton ) actionSource;

            if ( button == applyButton )
            {

                myCallback.recalc();

            }
            else
            {

                // it must be one of the Show buttons

                String buttonText = button.getText();

                logger.info( buttonText + " has been pressed" );

                CTProperties.setShowDefaults( buttonText );

                setCurrentSelections();
            }

        }

    } // Source was JCheckBox

} // CTPropertiesPanel
