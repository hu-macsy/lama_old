/*
 * TableDisplay.java
 *
 * Table display of an abstract table model that allow selection and hiding of columns.
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

import java.awt.BorderLayout;
import java.awt.Point;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.InputEvent;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;

import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JPopupMenu;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.event.TableModelEvent;
import javax.swing.event.TableModelListener;
import javax.swing.table.TableColumn;
import javax.swing.table.TableColumnModel;

import org.apache.log4j.Logger;

import org.jfree.data.category.CategoryDataset;
import org.jfree.data.category.DefaultCategoryDataset;
import org.jfree.data.statistics.DefaultStatisticalCategoryDataset;

import adaptor.Calltree.TableSorter;

/**
 * TableDisplay is a panel that displays abstract tables used for ShowPM.
 * Generally, it supports the abstract table model used for ShowPM, in
 * particular these are the region table, user table, process table and
 * comparison table. The table displays itself provides a menu that allows
 * to hide selected columns and to show them again.
 *
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */
/**
 * TODO: brandes: Enter comment!
 *
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */
class TableDisplay extends JPanel implements MouseListener, ActionListener, TableModelListener
{

    /**
     * This is the factor used to scale to percent.
     */
    private static final double PERCENTAGE_FACTOR = 0.01;

    /**
     * Logger variable for this class.
     */
    private static Logger logger = Logger.getLogger( TableDisplay.class );

    /**
     * This is the instance of abstract table that the display shows currently.
     */
    private PMTableModel myTableModel = null;

    /**
     * This is the sorted table view of the abstract table.
     */
    private TableSorter myTableSorter = null;

    /**
     * This is JTable used for the table display and becomes the central
     * part of the table display.
     */
    private JTable myJTable = null;

    /**
     * This is the title of the display. Getter and setter are provided.
     */
    private JLabel titleLabel = null;

    /**
     * This is the menu entry that is used to show hided columns again.
     */
    private JMenuItem noneItem = new JMenuItem( "NONE" );

    /**
     * This is the menu entry that is used to hide the selected columns.
     */
    private JMenuItem selItem = new JMenuItem( "selected col" );

    /**
     * This is the menu entry for the popup Menu to hide the selected column.
     */
    private JMenuItem hideItem = new JMenuItem( "Hide column" );

    /**
     * This variable counts the removed columns.
     */
    private int noHiddenColumns = 0;

    /**
     * This variable keeps track of all hidden columns.
     */
    private TableColumn[] hiddenColumns = null;

    /**
     * This variable contains the popup menu. It will only be built once
     * and then reused.
     */
    private JPopupMenu popup = null;

    /**
     * Constructor for a table display by the abstract table model
     * and the title.
     *
     * @param model is the abstract table model to display
     * @param title is the title for the display
     */
    TableDisplay( PMTableModel model, String title )
    {

        super( new BorderLayout() );

        myTableModel = model;

        model.addTableModelListener( this );

        myTableSorter = new TableSorter( model );
        myJTable = new JTable( myTableSorter );

        // selection of single columns and rows

        myJTable.setRowSelectionAllowed( false );

        // add the mouse Listener to the JTable

        myJTable.addMouseListener( this );

        myTableSorter.setTableHeader( myJTable.getTableHeader() );

        titleLabel = new JLabel( title );

        TableColumnModel colModel = myJTable.getColumnModel();
        colModel.setColumnSelectionAllowed( true );

        // logger.info("Table has " + M.getColumnCount() + " columns");

        JMenuBar menuBar = new JMenuBar();
        JMenu theHideMenu = new JMenu( "Hide" );

        noneItem.addActionListener( this );
        selItem.addActionListener( this );

        theHideMenu.add( noneItem );
        theHideMenu.add( selItem );

        menuBar.add( titleLabel );
        menuBar.add( theHideMenu );

        add( "North", menuBar );
        add( "Center", new JScrollPane( myJTable ) );

    }

    /**
     * This routine build the popup menu used for the table cells on right
     * mouse button click.
     *
     */
    private void makePopup()
    {

        popup = new JPopupMenu();

        hideItem.addActionListener( this );
        popup.add( hideItem );

        String [] myItems = myTableModel.getPopupItems();

        for ( int i = 0; i < myItems.length; i++ )
        {

            JMenuItem newItem = new JMenuItem( myItems[i] );
            newItem.addActionListener( this );
            popup.add( newItem );

        }

    } // makePopup

    /**
     * Help routine to make a string "row = r, col = c".
     *
     * @param row is the row
     * @param col is the column
     * @return a string containing both values
     */
    private String makePosString( int row, int col )
    {

        return "row = " + row + ", col = " + col;

    }
    /**
     * This routine will show the popup menu if there is a
     * selection available.
     *
     * @param e is used to get the position of the menu
     */
    private void showPopup( MouseEvent e )
    {

        int row = myJTable.getSelectedRow();
        int col = myJTable.getSelectedColumn();

        if ( ( row >= 0 ) && ( col >= 0 ) )
        {

            // so we are sure that there is really a selection

            if ( popup == null )
            {

                makePopup();
            }

            Point mousePoint = e.getPoint();

            popup.show( this, ( int ) mousePoint.getX(), ( int ) mousePoint.getY() );

        }
        else
        {

            logger.info( "no selection, no popup" );
        }
    }

    /**
     * This routine is called when the mouse has been pressed
     * on a selected cell.
     *
     */
    private void doMouseSelection()
    {

        int row = myJTable.getSelectedRow();
        int col = myJTable.getSelectedColumn();

        logger.info( "mouse pressed on selected " + makePosString( row, col ) );

        if ( ( row >= 0 ) && ( col >= 0 ) )
        {

            // there was really a selection before

            row = myTableSorter.modelIndex( row );
            col = myJTable.convertColumnIndexToModel( col );

            logger.info( "mouse pressed on model " + makePosString( row, col ) );

            // we pass the event to the different table realizations.

            myTableModel.mousePressed( col, row );
        }

    }

    /**
     * This routine is called when an item of the popup menu
     * has been selected.
     *
     * @param kind is the index number of menu entry in the popu menu
     */
    private void doActionSelection( int kind )
    {

        int row = myJTable.getSelectedRow();
        int col = myJTable.getSelectedColumn();

        logger.info( "action " + kind + " on selected " + makePosString( row, col ) );

        if ( ( row >= 0 ) && ( col >= 0 ) )
        {

            row = myTableSorter.modelIndex( row );
            col = myJTable.convertColumnIndexToModel( col );

            logger.info( "action " + kind + " on model " + makePosString( row, col ) );

            // we pass the event to the different table realizations.

            myTableModel.actionPopupItem( kind, col, row );

        }

    }

    /**
     * Right button of mouse will pop up the popup menu.
     *
     * {@inheritDoc}
     *
     * @see java.awt.event.MouseListener#mousePressed(java.awt.event.MouseEvent)
     */
    public void mousePressed( MouseEvent e )
    {

        if ( e.getButton() != MouseEvent.BUTTON3 )
        {

            return;
        }

        int modifiers = e.getModifiers();

        int flag = modifiers & ( InputEvent.BUTTON2_MASK | InputEvent.BUTTON3_MASK );

        if ( flag != 0 && flag == modifiers )
        {

            showPopup( e );

        } // mouse selection

    }

    /**
     * {@inheritDoc}
     *
     * @see java.awt.event.MouseListener#mouseReleased(java.awt.event.MouseEvent)
     */
    public void mouseReleased( MouseEvent e )
    {

    }

    /**
     * {@inheritDoc}
     *
     * @see java.awt.event.MouseListener#mouseClicked(java.awt.event.MouseEvent)
     */
    public void mouseClicked( MouseEvent e )
    {

        // logger.info(e);

        int modifiers = e.getModifiers();
        int clicks    = e.getClickCount();

        if ( ( modifiers == InputEvent.BUTTON1_MASK ) && ( clicks == 2 ) )
        {

            doMouseSelection();

        }
    }

    /**
     * {@inheritDoc}
     *
     * @see java.awt.event.MouseListener#mouseEntered(java.awt.event.MouseEvent)
     */
    public void mouseEntered( MouseEvent e )
    {
    }

    /**
     * {@inheritDoc}
     *
     * @see java.awt.event.MouseListener#mouseExited(java.awt.event.MouseEvent)
     */
    public void mouseExited( MouseEvent e )
    {

    }

    /**
     * This routine hides a column. It will be removed from the model and saved
     * so that it can be diplayed again later.
     *
     * @param col is the column to hide
     */
    private void hideColumn( int col )
    {

        TableColumnModel model = myJTable.getColumnModel();

        TableColumn column = model.getColumn( col );

        logger.info( "hide current col = " + col + ", model index = " + column.getModelIndex() );

        model.removeColumn( column );

        // stack the hidden column to make it reappear later

        hiddenColumns[noHiddenColumns++] = column;

    }

    /**
     * This routine adds back the hidden columns.
     *
     */
    private void showHiddenColumns()
    {

        int i;

        TableColumnModel model = myJTable.getColumnModel();

        for ( i = 0; i < noHiddenColumns; i++ )
        {

            model.addColumn( hiddenColumns[i] );
        }


        noHiddenColumns = 0;
    }

    /**
     * Help routine needed to rehide columns again if the
     * whole table has changed.
     *
     */
    private void rehideColumns()
    {

        TableColumnModel model = myJTable.getColumnModel();

        // it is very likely that all columns have been rebuilt

        for ( int i = 0; i < noHiddenColumns; i++ )
        {

            TableColumn oldColumn = hiddenColumns[i];

            int col = oldColumn.getModelIndex();

            for ( int j = 0;  j < model.getColumnCount(); j++ )
            {

                TableColumn newColumn = model.getColumn( j );

                if ( col == newColumn.getModelIndex() )
                {

                    // this is the colum we have to remove again

                    logger.info( "rehide model column " + col + ", is at current col " + col );

                    model.removeColumn( newColumn );

                    hiddenColumns[i] = newColumn;

                }
            }

        }
    }

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


        JMenuItem itemSource = ( JMenuItem ) actionSource;

        if ( itemSource.equals( noneItem ) )
        {

            logger.info( "get back " + noHiddenColumns + " hidden columns" );

            showHiddenColumns();

        } // action for NONE

        if ( itemSource.equals( selItem ) || itemSource.equals( hideItem ) )
        {

            int[] selectedColumns = myJTable.getSelectedColumns();

            if ( selectedColumns.length == 0 )
            {

                return;
            }


            // make sure that we do not hide the columns with region names

            for ( int col = 0; col < selectedColumns.length; col++ )
            {

                int mcol = selectedColumns[col];

                mcol = myJTable.convertColumnIndexToModel( mcol );

                if ( mcol == 0 )
                {

                    JOptionPane.showMessageDialog( null, "Region column cannot be hidden", "ShowPM Hide ERROR", JOptionPane.ERROR_MESSAGE );

                    return;
                }
            }

            TableColumnModel model = myJTable.getColumnModel();

            if ( hiddenColumns == null )
            {

                hiddenColumns = new TableColumn[model.getColumnCount()];

            }

            // now we hide all selected columns

            for ( int col = selectedColumns.length - 1; col >= 0; col-- )
            {

                hideColumn( selectedColumns[col] );

            } // for loop

        } // action for hiding columns

        String [] myItems = myTableModel.getPopupItems();

        for ( int i = 0; i < myItems.length; i++ )
        {

            if ( myItems[i].equals( itemSource.getText() ) )
            {

                doActionSelection( i );

            }

        }


    } // actionPerformed

    /**
     * Setter routine for the title of the Table Display.
     *
     * @param newTitle is the string to be set for title
     */
    public void setTitle( String newTitle )
    {

        titleLabel.setText( newTitle );

    }

    /**
     * Getter routine for the title of the Table Display.
     *
     * @return the title of the display as a string
     */
    public String getTitle()
    {

        return titleLabel.getText();

    }

    /**
     * This routine returns the current selected row.
     *
     * @return the index of the selected row
     */
    public int getSelectedRow()
    {

        int row = myJTable.getSelectedRow();

        // be careful: modelIndex (-1) throws exception

        if ( row >= 0 )
        {

            row = myTableSorter.modelIndex( row );

        }

        return row;

    }

    /**
     * This routine returns the selected column with its model index.
     *
     * @return integer value for the model index of the column
     */
    public int getSelectedColumn()
    {

        int col = myJTable.getSelectedColumn();
        col = myJTable.convertColumnIndexToModel( col );
        return col;

    }

    /**
     * This routine returns all selected columns by their model indexes
     * (and not by their positions).
     *
     * @return an array with all model indexes of the selected columns.
     */
    public int[] getSelectedColumns()
    {

        int[] cols = myJTable.getSelectedColumns();

        int i;

        String val = "table columns = ";

        for ( i = 0; i < cols.length; i++ )
        {

            val = val + cols[i] + " ";
        }

        // logger.info (val);

        for ( i = 0; i < cols.length; i++ )
        {

            cols[i] = myJTable.convertColumnIndexToModel( cols[i] );
        }

        val = "model columns ";

        for ( i = 0; i < cols.length; i++ )
        {

            val = val + cols[i] + " ";
        }

        // logger.info(val);

        return cols;

    }

    /**
     * This method returns the selected columns
     * (note: might also return column with regions).
     *
     * @return integer array with indexes of all selected columns
     */
    public int[] getColumnSelection()
    {

        int[] columns;
        int numberColumns; // number of selected columns
        int col; // loop variable

        columns = myJTable.getSelectedColumns();

        numberColumns = 0;

        if ( columns != null )
        {

            numberColumns = columns.length;
        }


        // if we have no selected columns we take all columns
        // that are not hidden

        if ( numberColumns == 0 )
        {

            numberColumns = myJTable.getColumnCount();

            columns = new int[numberColumns];

            for ( col = 0; col < numberColumns; col++ )
            {

                columns[col] = col;
            }


        }

        // now we have to map back to the columns of the model

        for ( col = 0; col < numberColumns; col++ )
        {

            int col1 = myJTable.convertColumnIndexToModel( columns[col] );

            columns[col] = col1;

            /* was useful for debugging:

             logger.info("column " + col + " " +
             + " is model column " + col1 + " " +
             Sorted_Table.getColumnName (col1));

             */

        } // for loop

        return columns;

    }

    // the following method gets the data of selected columns

    /**
     * This methode creates a two-dimensional array with double values
     * of the abstract table. Only selected columns are taken for the values.
     *
     * @param rowCount is the number of row
     * @param columns contains the selected column indexes
     * @param colCount is the number of columns, can be less than Columns.length
     *
     * @return a two-dimensional array with table values converted to double
     */
    double[][] getTableValues( int rowCount, int[] columns, int colCount )
    {

        double[][] values = new double[rowCount][colCount];

        for ( int row = 0; row < rowCount; row++ )
        {

            for ( int col = 0; col < colCount; col++ )
            {

                Object value = myTableSorter.getValueAt( row, columns[col] );

                if ( value instanceof Double )
                {

                    values[row][col] = ( ( Double ) value ).doubleValue();

                }
                else if ( value instanceof Long )
                {

                    values[row][col] = ( ( Long ) value ).longValue();

                }
                else
                {

                    values[row][col] = 0.0;
                }

            }
        }

        return values;

    }

    /**
     * This routine creates a category data set from the table.
     *
     * @return data set that is used to make a chart of it
     */
    public CategoryDataset createDataset()
    {

        double depth = DatasetProperties.getSelectedDepth();

        int[] columns;
        int numberColumns; // number of selected coluns
        int row;
        int col; // loop variables

        // we create the dataset only for the selected columns

        columns = getColumnSelection();

        // we do not consider the column with the regions

        numberColumns = 0;

        for ( col = 0; col < columns.length; col++ )
        {

            if ( columns[col] > 0 )
            {

                columns[numberColumns++] = columns[col];
            }

        }

        if ( numberColumns == 0 )
        {

            return null;
        }

        int numberRows = myTableSorter.getRowCount();

        int entries = DatasetProperties.getSelectedEntries();

        // entries > 0 : restrict the number of regions

        if ( entries > 0 )
        {

            // so the number of regions should be restriccted

            numberRows = Math.min( numberRows, entries );

        }

        final DefaultCategoryDataset dataset = new DefaultCategoryDataset();

        double[][] values = getTableValues( numberRows, columns, numberColumns );

        // now we build the sum of values to limit entries

        double[] columnSums = new double[numberColumns];

        for ( col = 0; col < numberColumns; col++ )
        {

            columnSums[col] = 0.0;

            for ( row = 0; row < numberRows; row++ )
            {

                columnSums[col] += values[row][col];
            }

        }

        // now we create the data set

        for ( row = 0; row < numberRows; row++ )
        {

            boolean isRowVisible = false;

            // step 1 : check if we take the row (at least one value big enough)

            for ( col = 0; col < numberColumns; col++ )
            {

                double dvalue = values[row][col];

                double threshold = columnSums[col] * ( depth * PERCENTAGE_FACTOR );

                if ( dvalue >= threshold )
                {

                    isRowVisible = true;

                    break;  // no more checks necessary
                }

            }

            if ( isRowVisible )
            {

                for ( col = 0; col < numberColumns; col++ )
                {

                    String series;
                    String category;
                    double dvalue;

                    // each region becomes its own category

                    category = ( String ) myTableSorter.getValueAt( row, 0 );
                    series = myTableSorter.getColumnName( columns[col] );

                    dvalue = values[row][col];

                    dataset.addValue( dvalue, series, category );
                }
            }


        } // for all rows

        return dataset;

    }

    /**
     * This method generates a dataset by extracting the data of the displayed table.
     *
     * @return a CatagoryDataset that can be used for making charts
     */
    public CategoryDataset createStatDataset()
    {

        int entries = DatasetProperties.getSelectedEntries();
        double depth = DatasetProperties.getSelectedDepth();

        int[] columns;
        int numberColumns; // number of selected coluns
        int row;
        int col; // loop variables

        // we create the dataset only for the selected columns

        columns = getColumnSelection();

        // we do not consider the column with the regions

        numberColumns = 0;

        for ( col = 0; col < columns.length; col++ )
        {

            if ( columns[col] > 0 )
            {

                columns[numberColumns++] = columns[col];
            }
        }

        if ( numberColumns == 0 )
        {

            return null;
        }


        int numberRows = myTableSorter.getRowCount();

        // entries > 0 : restrict the number of regions

        if ( entries > 0 )
        {

            numberRows = Math.min( numberRows, entries );
        }

        final DefaultStatisticalCategoryDataset dataset = new DefaultStatisticalCategoryDataset();

        double[][] values = getTableValues( numberRows, columns, numberColumns );

        // now we build the average of the column values for each region

        double[] valueAverages = new double[numberRows];

        double valueAverageAll = 0.0;

        for ( row = 0; row < numberRows; row++ )
        {

            double avg = 0.0;

            for ( col = 0; col < numberColumns; col++ )
            {

                avg += values[row][col];

            }

            valueAverages[row] = avg / numberColumns;

            valueAverageAll += avg;
        }

        valueAverageAll = valueAverageAll / numberColumns;

        // now we build the standard deviation

        double[] valueStdDeviation = new double[numberRows];

        for ( row = 0; row < numberRows; row++ )
        {

            double avg = valueAverages[row];

            double stddev = 0.0;

            for ( col = 0; col < numberColumns; col++ )
            {

                double v = values[row][col];

                stddev += ( avg - v ) * ( avg - v );

            }

            stddev = Math.sqrt( stddev / numberColumns );

            valueStdDeviation[row] = stddev;

        }

        // now we create the data set

        for ( row = 0; row < numberRows; row++ )
        {

            String category;

            double dvalue = valueAverages[row];

            // each region becomes its own category

            category = ( String ) myTableSorter.getValueAt( row, 0 );

            /* might be useful for debugging

             logger.info("value of region " + category + " = " +
             dvalue + " check against " +
             value_avg_all + " (depth " +
             value_avg_all * (depth / 100));
             */

            if ( dvalue >= valueAverageAll * ( depth * PERCENTAGE_FACTOR ) )
            {

                dataset.add( dvalue, valueStdDeviation[row], "", category );
            }
        }

        return dataset;

    }

    /**
     * {@inheritDoc}
     *
     * @see javax.swing.event.TableModelListener#tableChanged(javax.swing.event.TableModelEvent)
     */
    public void tableChanged( TableModelEvent event )
    {

        // This routine is called if the table has been changed

        logger.info( "call of tableChanged" );
        logger.info( event );

        int firstRow = event.getFirstRow();

        logger.info( "first row = " + event.getFirstRow() + ", last "
                     + makePosString( event.getLastRow(), event.getColumn() ) );

        if ( firstRow < 0 )
        {

            rehideColumns();

        }

    }

} // class TableDisplay
