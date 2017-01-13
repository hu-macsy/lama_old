/*
 * ChartMenu.java
 *
 * Menu for TableDisplay to generate charts.
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

package adaptor.ShowPM;

import java.awt.Component;
import java.awt.Container;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JMenu;
import javax.swing.JMenuItem;
import javax.swing.JPopupMenu;

import org.jfree.util.TableOrder;

import adaptor.General.PropertyMenu;

/**
 * The class ChartMenu provides entries for making charts.
 * It can be allocated for a TableDisplay as methods of this
 * class will be called.
 *
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */
class ChartMenu extends JMenu implements ActionListener
{

    /**
     * The chart menu is anchored to a table display from
     * which we will create the relevant values for the chart.
     */
    private TableDisplay myTableDisplay;

    /**
     * Charts need a category name (not used here).
     */
    private final String myCategoryName = null;

    /**
     * Charts need a value name (not used here).
     */
    private final String myValueName = null;

    /**
     * Menu to select the orientation of the Charts. Items for
     * this menu are HORIZONTAL and VERTICAL. The orientation will
     * be set globally for the Charts properties.
     */
    private JMenu menuOrientation = new JMenu( "Set Orientation" );

    /**
     * Menu to select labels for Pie Charts.
     */
    private JMenu menuLabel = new JMenu( "Set Labels (Pie)" );


    /**
     * Menu item to generate a Bar Chart.
     */
    private JMenuItem menuBarItem = new JMenuItem( "Bar Chart" );

    /**
     * Menu item to generate a Stacked Bar Chart.
     */
    private JMenuItem menuStBarItem = new JMenuItem( "Stacked Bar Chart" );

    /**
     * Menu item to generate a Statistical Bar Chart.
     */
    private JMenuItem menuStatBarItem = new JMenuItem( "Statistical Bar Chart" );

    /**
     * Menu item to generate a Multiple Pie Chart (for each performance event).
     */
    private JMenuItem menuPieItem = new JMenuItem( "Multiple Pie Charts" );

    /**
     * Menu item to generate a Multiple Pie Chart (for each region).
     */
    private JMenuItem menuPieItem1 = new JMenuItem( "Multiple Pie Charts (Regions)" );

    /**
     * Constructor for a new ChartMenu.
     *
     * @param title becomes the title of this menu
     * @param newTD ist the TableDisplay to which the menu will be anchored
     */
    ChartMenu( String title, TableDisplay newTD )
    {

        super( title );

        myTableDisplay = newTD;

        int i;

        JMenuItem numberItem;

        DatasetProperties.makePropertyMenus();

        // Property Menus

        String[] headers = PropertyMenu.getAllHeaders();

        for ( i = 0; i < headers.length; i++ )
        {

            String [] titles = PropertyMenu.getAllTitles( headers[i] );

            for ( int j = 0; j < titles.length; j++ )
            {

                PropertyMenu propMenu = PropertyMenu.findMenu( headers[i], titles[j] );

                JMenu menu = new JMenu( titles[j] );

                String[] items = propMenu.getItems();

                for ( int k = 0; k < items.length; k++ )
                {

                    JMenuItem menuItem = new JMenuItem( items[k] );
                    menuItem.addActionListener( this );
                    menu.add( menuItem );
                }

                add( menu );
            }
        }

        // Orientation Menu

        String[] itemsOrientation = ChartProperties.ITEMS_ORIENTATION;

        for ( i = 0; i < itemsOrientation.length; i++ )
        {

            numberItem = new JMenuItem( itemsOrientation[i] );
            numberItem.addActionListener( this );
            menuOrientation.add( numberItem );
        }

        // Labels Menu

        String[] itemsLabel = ChartProperties.ITEMS_LABEL;

        for ( i = 0; i < itemsLabel.length; i++ )
        {

            numberItem = new JMenuItem( itemsLabel[i] );
            numberItem.addActionListener( this );
            menuLabel.add( numberItem );
        }

        menuBarItem.addActionListener( this );
        menuStBarItem.addActionListener( this );
        menuStatBarItem.addActionListener( this );
        menuPieItem.addActionListener( this );
        menuPieItem1.addActionListener( this );

        add( menuOrientation );
        add( menuLabel );
        add( menuBarItem );
        add( menuStBarItem );
        add( menuStatBarItem );
        add( menuPieItem );
        add( menuPieItem1 );

    } // end construcutor ChartMenu

    /**
     * {@inheritDoc}
     *
     * @see java.awt.event.ActionListener#actionPerformed(java.awt.event.ActionEvent)
     */
    public void actionPerformed( ActionEvent e )
    {

        Object actionSource = ( Object ) e.getSource();

        if ( !( actionSource instanceof JMenuItem ) )
        {

            return;
        }

        String cmd = ( ( JMenuItem ) actionSource ).getText();

        String title = cmd;

        Container parent = ( ( JMenuItem ) actionSource ).getParent();

        if ( parent instanceof JPopupMenu )
        {

            Component invoker = ( ( JPopupMenu ) parent ).getInvoker();

            if ( invoker instanceof JMenu )
            {

                title = ( ( JMenu ) invoker ).getText();
            }
        }

        // some commands belong to the DatasetProperties

        DatasetProperties.setProperty( title, cmd );

        // some commands belong to the ChartProperties

        ChartProperties.setProperty( cmd );

        title = myTableDisplay.getTitle();

        if ( cmd.equals( menuBarItem.getText() ) )
        {

            ChartProperties.makeBarChart( title, myCategoryName, myValueName, myTableDisplay.createDataset() );

        } // action for Horizontal Bar Chart

        if ( cmd.equals( menuStBarItem.getText() ) )
        {

            ChartProperties.makeStackedBarChart( title, myCategoryName, myValueName, myTableDisplay.createDataset() );

        } // action for Stacked Bar Chart

        if ( cmd.equals( menuStatBarItem.getText() ) )
        {

            ChartProperties.makeStatBarChart( title, myCategoryName, myValueName, myTableDisplay.createStatDataset() );

        } // action for Horizontal Bar Chart

        if ( cmd.equals( menuPieItem.getText() ) )
        {

            ChartProperties.makePieChart( title, myCategoryName, myValueName, myTableDisplay.createDataset(), TableOrder.BY_ROW );

        } //

        if ( cmd.equals( menuPieItem1.getText() ) )
        {

            ChartProperties.makePieChart( title, myCategoryName, myValueName, myTableDisplay.createDataset(), TableOrder.BY_COLUMN );

        } // action for Horizontal Bar Chart

    } // actionPerformed

} // class ChartMenu
