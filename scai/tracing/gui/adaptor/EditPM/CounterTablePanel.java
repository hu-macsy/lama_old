/*
 * CounterTablePanel.java
 *
 * Panel for the the table that contains the user selected counters.
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

import java.awt.BorderLayout;
import java.awt.Dimension;

import javax.swing.DefaultCellEditor;
import javax.swing.JComboBox;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.table.DefaultTableCellRenderer;
import javax.swing.table.TableColumn;

import org.apache.log4j.Logger;

/**
 * The class CounterTablePanel defines a panel containing the table with the defined counters.
 *
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */
public class CounterTablePanel extends JPanel
{

    /**
     * Logger for this class.
     */
    private static Logger logger = Logger.getLogger( CounterTablePanel.class );

    /**
     * This is the preferred size for the counter table.
     */
    private static final Dimension COUNTER_TABLE_DIMENSION = new Dimension( 800, 300 );

    /**
     * Display for the counter table.
     */
    private JTable myCounterJTable;

    /**
     * This counter is used for delete and paste.
     */
    private Counter mySavedCounter = null;

    /**
     * This is the underlying table model for the specified counters.
     */
    private CounterTable myCounterTable;

    /**
     * This constructor makes a panel for the counter table.
     *
     * @param table is the counter table
     */
    public CounterTablePanel( CounterTable table )
    {

        super( new BorderLayout() );

        myCounterTable = table;

        myCounterJTable = new JTable( myCounterTable );

        myCounterJTable.setPreferredScrollableViewportSize( COUNTER_TABLE_DIMENSION );

        // Create the scroll pane and add the CTable to it.

        JScrollPane scrollPane = new JScrollPane( myCounterJTable );

        setUpScaleBox( myCounterJTable, myCounterJTable.getColumnModel().getColumn( 3 ) );

        // JSplitPane hpane = new JSplitPane (JSplitPane.HORIZONTAL_SPLIT,
        // scrollEventMenu, scrollPane);

        add( "Center", scrollPane );

        // Attention: north and south can be filled up by caller

    } // constructor CounterTableDisplay

    /**
     * This routine makes a scale box for a certain column in the counter table.
     *
     * @param table is the JTable for which a box is generated
     * @param scaleColumn is the colum for which the scale box will be enabled
     */
    public void setUpScaleBox( JTable table, TableColumn scaleColumn )
    {

        // Set up the editor for the scale cells.

        JComboBox comboBox = new JComboBox();

        comboBox.addItem( "m" );
        comboBox.addItem( "%" );
        comboBox.addItem( "" );
        comboBox.addItem( "k" );
        comboBox.addItem( "M" );
        comboBox.addItem( "G" );

        scaleColumn.setCellEditor( new DefaultCellEditor( comboBox ) );

        // Set up tool tips for the scale cells.

        DefaultTableCellRenderer renderer = new DefaultTableCellRenderer();
        renderer.setToolTipText( "Click for combo box" );
        scaleColumn.setCellRenderer( renderer );

    } // method setUpScaleBox

    /**
     * This methode deletes the counter in the counter table
     * at the currently selected row.
     *
     */
    public void deleteCounter()
    {

        int row = myCounterJTable.getSelectedRow();

        logger.info( "remove selected row (is " + row + ")" );

        if ( ( row >= 0 ) && ( row < myCounterTable.getRowCount() ) )
        {

            // so this is a valid row that can be deleted

            mySavedCounter = myCounterTable.deleteCounter( row );

            myCounterJTable.revalidate();
        }

    } // deleteCounter

    /**
     * This routine copies the counter of the currently selected
     * row into the clipboard.
     *
     */
    void copyCounter()
    {

        int row = myCounterJTable.getSelectedRow();

        if ( ( row >= 0 ) && ( row < myCounterTable.getRowCount() ) )
        {

            mySavedCounter = myCounterTable.getCounter( row );

        }

    } // copyCounter

    /**
     * This routine pastes the counter from the Clipboard
     * to the currently selected row.
     *
     */
    void pasteCounter()
    {

        int row = myCounterJTable.getSelectedRow();

        row = Math.max( 0, row );
        row = Math.min( row, myCounterTable.getRowCount() );

        myCounterTable.addCounter( row, mySavedCounter.copyCounter() );

        myCounterJTable.revalidate();

    } // pasteCounter

    /**
     * This method deletes all counters from the counter table.
     *
     */
    public void clear()
    {

        myCounterTable.clear();

        // size of table has changed, so do a repaint

        myCounterJTable.revalidate();

    } // remove all

    /**
     * This method makes sure that the counter table will be displayed
     * correctly after any updates.
     *
     */
    public void update()
    {

        // we need to revalidate the table

        myCounterJTable.revalidate();

    }

} // class CounterTableDisplay
