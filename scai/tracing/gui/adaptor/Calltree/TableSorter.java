/*
 * TableSorter.java
 *
 * Sorting of abstract table model.
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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import java.awt.Color;
import java.awt.Component;
import java.awt.Graphics;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;

import javax.swing.Icon;
import javax.swing.JLabel;
import javax.swing.JTable;
import javax.swing.event.TableModelEvent;
import javax.swing.event.TableModelListener;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.JTableHeader;
import javax.swing.table.TableCellRenderer;
import javax.swing.table.TableColumnModel;
import javax.swing.table.TableModel;

/**
 * TableSorter is a decorator for TableModels; adding sorting
 * functionality to a supplied TableModel. TableSorter does
 * not store or copy the data in its TableModel; instead it maintains
 * a map from the row indexes of the view to the row indexes of the
 * model. As requests are made of the sorter (like getValueAt(row, col))
 * they are passed to the underlying model after the row numbers
 * have been translated via the internal mapping array. This way,
 * the TableSorter appears to hold another copy of the table
 * with the rows in a different order.
 * <p/>
 * TableSorter registers itself as a listener to the underlying model,
 * just as the JTable itself would. Events recieved from the model
 * are examined, sometimes manipulated (typically widened), and then
 * passed on to the TableSorter's listeners (typically the JTable).
 * If a change to the model has invalidated the order of TableSorter's
 * rows, a note of this is made and the sorter will resort the
 * rows the next time a value is requested.
 * <p/>
 * When the tableHeader property is set, either by using the
 * setTableHeader() method or the two argument constructor, the
 * table header may be used as a complete UI for TableSorter.
 * The default renderer of the tableHeader is decorated with a renderer
 * that indicates the sorting status of each column. In addition,
 * a mouse listener is installed with the following behavior:
 * <ul>
 * <li>
 * Mouse-click: Clears the sorting status of all other columns
 * and advances the sorting status of that column through three
 * values: {NOT_SORTED, ASCENDING, DESCENDING} (then back to
 * NOT_SORTED again).
 * <li>
 * SHIFT-mouse-click: Clears the sorting status of all other columns
 * and cycles the sorting status of the column through the same
 * three values, in the opposite order: {NOT_SORTED, DESCENDING, ASCENDING}.
 * <li>
 * CONTROL-mouse-click and CONTROL-SHIFT-mouse-click: as above except
 * that the changes to the column do not cancel the statuses of columns
 * that are already sorting - giving a way to initiate a compound
 * sort.
 * </ul>
 * <p/>
 * This is a long overdue rewrite of a class of the same name that
 * first appeared in the swing table demos in 1997.
 *
 * @author Philip Milne
 * @author Brendon McLean
 * @author Dan van Enckevort
 * @author Parwinder Sekhon
 * @version 2.0 02/27/04
 */


public class TableSorter extends AbstractTableModel
{


    /**
     * Constant value for desceding.
     */
    public static final int DESCENDING = -1;

    /**
     * Constand value for not sorted.
     */
    public static final int NOT_SORTED = 0;

    /**
     * Constand value for ascending order.
     */
    public static final int ASCENDING = 1;

    private static final Directive EMPTY_DIRECTIVE = new Directive( -1, NOT_SORTED );

    public static final Comparator COMPARABLE_COMAPRATOR = new Comparator()
    {

        public int compare( Object o1, Object o2 )
        {
            return ( ( Comparable ) o1 ).compareTo( o2 );
        }
    };

    public static final Comparator LEXICAL_COMPARATOR = new Comparator()
    {

        public int compare( Object o1, Object o2 )
        {
            return o1.toString().compareTo( o2.toString() );
        }
    };

    /**
     * This is the underlying table model for the sorted table.
     */
    protected TableModel tableModel;

    private Row[] viewToModel;
    private int[] modelToView;

    private JTableHeader tableHeader;
    private MouseListener mouseListener;
    private TableModelListener tableModelListener;
    private Map columnComparators = new HashMap();
    private List sortingColumns = new ArrayList();

    /**
     * Constructor for this class.
     *
     */
    public TableSorter()
    {
        this.mouseListener = new MouseHandler();
        this.tableModelListener = new TableModelHandler();
    }

    /**
     * Constructor where table model is already available.
     *
     * @param tableModel is the abstract table model to be sorted
     */
    public TableSorter( TableModel tableModel )
    {
        this();
        setTableModel( tableModel );
    }

    /**
     * Constructor with table model and table header.
     *
     * @param model is the model
     * @param header is the header
     */
    public TableSorter( TableModel model, JTableHeader header )
    {
        this();
        setTableHeader( header );
        setTableModel( model );
    }

    /**
     * Routine to reset the sorting state.
     *
     */
    private void clearSortingState()
    {
        viewToModel = null;
        modelToView = null;
    }

    /**
     * Getter routine for the table model.
     *
     * @return the table model of this table.
     */
    public TableModel getTableModel()
    {
        return tableModel;
    }

    /**
     * Setter routine for the table model.
     *
     * @param newTableModel is the new table model that will be sorted
     */
    public void setTableModel( TableModel newTableModel )
    {

        if ( this.tableModel != null )
        {

            this.tableModel.removeTableModelListener( tableModelListener );
        }

        this.tableModel = newTableModel;

        if ( this.tableModel != null )
        {

            this.tableModel.addTableModelListener( tableModelListener );
        }

        clearSortingState();
        fireTableStructureChanged();
    }

    /**
     * Getter routine for the table header.
     *
     * @return the tableHeader
     */
    public JTableHeader getTableHeader()
    {

        return tableHeader;
    }

    /**
     * Setter routine for the table header.
     *
     * @param newTableHeader is the new table header
     */
    public void setTableHeader( JTableHeader newTableHeader )
    {

        if ( this.tableHeader != null )
        {
            this.tableHeader.removeMouseListener( mouseListener );
            TableCellRenderer defaultRenderer = this.tableHeader.getDefaultRenderer();

            if ( defaultRenderer instanceof SortableHeaderRenderer )
            {
                this.tableHeader.setDefaultRenderer( ( ( SortableHeaderRenderer ) defaultRenderer ).tableCellRenderer );
            }
        }

        this.tableHeader = newTableHeader;

        if ( this.tableHeader != null )
        {
            this.tableHeader.addMouseListener( mouseListener );
            this.tableHeader.setDefaultRenderer( new SortableHeaderRenderer( this.tableHeader.getDefaultRenderer() ) );
        }
    }

    /**
     * Ask wheter sorting is enabled at a given time.
     *
     * @return true if number of sorted columns is not zero
     */
    public boolean isSorting()
    {

        return sortingColumns.size() != 0;
    }

    /**
     * Asking for the directive of a column.
     *
     * @param column is the index of the column.
     * @return directive for this column
     */
    private Directive getDirective( int column )
    {
        for ( int i = 0; i < sortingColumns.size(); i++ )
        {
            Directive directive = ( Directive ) sortingColumns.get( i );

            if ( directive.column == column )
            {
                return directive;
            }
        }

        return EMPTY_DIRECTIVE;
    }

    public int getSortingStatus( int column )
    {
        return getDirective( column ).direction;
    }

    private void sortingStatusChanged()
    {
        clearSortingState();
        fireTableDataChanged();

        if ( tableHeader != null )
        {
            tableHeader.repaint();
        }
    }

    public void setSortingStatus( int column, int status )
    {
        Directive directive = getDirective( column );

        if ( directive != EMPTY_DIRECTIVE )
        {
            sortingColumns.remove( directive );
        }

        if ( status != NOT_SORTED )
        {
            sortingColumns.add( new Directive( column, status ) );
        }

        sortingStatusChanged();
    }

    protected Icon getHeaderRendererIcon( int column, int size )
    {
        Directive directive = getDirective( column );

        if ( directive == EMPTY_DIRECTIVE )
        {
            return null;
        }

        return new Arrow( directive.direction == DESCENDING, size, sortingColumns.indexOf( directive ) );
    }

    /**
     * This routine cancels the sorting.
     *
     */
    private void cancelSorting()
    {
        sortingColumns.clear();
        sortingStatusChanged();
    }

    public void setColumnComparator( Class type, Comparator comparator )
    {
        if ( comparator == null )
        {
            columnComparators.remove( type );
        }
        else
        {
            columnComparators.put( type, comparator );
        }
    }

    protected Comparator getComparator( int column )
    {
        Class columnType = tableModel.getColumnClass( column );
        Comparator comparator = ( Comparator ) columnComparators.get( columnType );

        if ( comparator != null )
        {
            return comparator;
        }

        if ( Comparable.class.isAssignableFrom( columnType ) )
        {
            return COMPARABLE_COMAPRATOR;
        }

        return LEXICAL_COMPARATOR;
    }

    private Row[] getViewToModel()
    {
        if ( viewToModel == null )
        {
            int tableModelRowCount = tableModel.getRowCount();
            viewToModel = new Row[tableModelRowCount];

            for ( int row = 0; row < tableModelRowCount; row++ )
            {
                viewToModel[row] = new Row( row );
            }

            if ( isSorting() )
            {
                Arrays.sort( viewToModel );
            }
        }

        return viewToModel;
    }

    public int modelIndex( int viewIndex )
    {
        return getViewToModel()[viewIndex].modelIndex;
    }

    private int[] getModelToView()
    {
        if ( modelToView == null )
        {
            int n = getViewToModel().length;
            modelToView = new int[n];

            for ( int i = 0; i < n; i++ )
            {
                modelToView[modelIndex( i )] = i;
            }
        }

        return modelToView;
    }

    // TableModel interface methods

    /**
     * {@inheritDoc}
     *
     * @see javax.swing.table.TableModel#getRowCount()
     */
    public int getRowCount()
    {

        int count = 0;

        if ( tableModel != null )
        {

            count = tableModel.getRowCount();
        }

        return count;
    }

    /**
     * {@inheritDoc}
     *
     * @see javax.swing.table.TableModel#getColumnCount()
     */
    public int getColumnCount()
    {

        int count = 0;

        if ( tableModel != null )
        {

            count = tableModel.getColumnCount();
        }

        return count;

    }

    public String getColumnName( int column )
    {

        return tableModel.getColumnName( column );
    }

    public Class getColumnClass( int column )
    {

        return tableModel.getColumnClass( column );
    }

    /**
     * {@inheritDoc}
     *
     * @see javax.swing.table.TableModel#isCellEditable(int, int)
     */
    public boolean isCellEditable( int row, int column )
    {
        return tableModel.isCellEditable( modelIndex( row ), column );
    }

    /**
     * {@inheritDoc}
     *
     * @see javax.swing.table.TableModel#getValueAt(int, int)
     */
    public Object getValueAt( int row, int column )
    {
        return tableModel.getValueAt( modelIndex( row ), column );
    }

    /**
     * {@inheritDoc}
     *
     * @see javax.swing.table.TableModel#setValueAt(java.lang.Object, int, int)
     */
    public void setValueAt( Object aValue, int row, int column )
    {
        tableModel.setValueAt( aValue, modelIndex( row ), column );
    }

    // Helper classes

    private class Row implements Comparable
    {

        private int modelIndex;

        public Row( int index )
        {
            this.modelIndex = index;
        }

        public int compareTo( Object o )
        {

            final int lessValue = -1;

            int row1 = modelIndex;
            int row2 = ( ( Row ) o ).modelIndex;

            for ( Iterator it = sortingColumns.iterator(); it.hasNext(); )
            {
                Directive directive = ( Directive ) it.next();
                int column = directive.column;
                Object o1 = tableModel.getValueAt( row1, column );
                Object o2 = tableModel.getValueAt( row2, column );

                int comparison = 0;

                // Define null less than everything, except null.
                if ( o1 == null && o2 == null )
                {
                    comparison = 0;
                }
                else if ( o1 == null )
                {
                    comparison = lessValue;
                }
                else if ( o2 == null )
                {
                    comparison = 1;
                }
                else
                {
                    comparison = getComparator( column ).compare( o1, o2 );
                }

                if ( comparison != 0 )
                {

                    if ( directive.direction == DESCENDING )
                    {

                        comparison = - comparison;

                    }

                    return comparison;
                }
            }

            return 0;
        }
    }

    private class TableModelHandler implements TableModelListener
    {

        public void tableChanged( TableModelEvent e )
        {
            // If we're not sorting by anything, just pass the event along.
            if ( !isSorting() )
            {
                clearSortingState();
                fireTableChanged( e );
                return;
            }

            // If the table structure has changed, cancel the sorting; the
            // sorting columns may have been either moved or deleted from
            // the model.
            if ( e.getFirstRow() == TableModelEvent.HEADER_ROW )
            {
                cancelSorting();
                fireTableChanged( e );
                return;
            }

            // We can map a cell event through to the view without widening
            // when the following conditions apply:
            //
            // a) all the changes are on one row (e.getFirstRow() == e.getLastRow()) and,
            // b) all the changes are in one column (column != TableModelEvent.ALL_COLUMNS) and,
            // c) we are not sorting on that column (getSortingStatus(column) == NOT_SORTED) and,
            // d) a reverse lookup will not trigger a sort (modelToView != null)
            //
            // Note: INSERT and DELETE events fail this test as they have column == ALL_COLUMNS.
            //
            // The last check, for (modelToView != null) is to see if modelToView
            // is already allocated. If we don't do this check; sorting can become
            // a performance bottleneck for applications where cells
            // change rapidly in different parts of the table. If cells
            // change alternately in the sorting column and then outside of
            // it this class can end up re-sorting on alternate cell updates -
            // which can be a performance problem for large tables. The last
            // clause avoids this problem.
            int column = e.getColumn();

            if ( e.getFirstRow() == e.getLastRow() && column != TableModelEvent.ALL_COLUMNS && getSortingStatus( column ) == NOT_SORTED
                    && modelToView != null )
            {
                int viewIndex = getModelToView()[e.getFirstRow()];
                fireTableChanged( new TableModelEvent( TableSorter.this, viewIndex, viewIndex, column, e.getType() ) );
                return;
            }

            // Something has happened to the data that may have invalidated the row order.
            clearSortingState();
            fireTableDataChanged();
            return;
        }
    }

    private class MouseHandler extends MouseAdapter
    {

        private static final int NO_SELECTION = -1;

        /**
         * {@inheritDoc}
         *
         * @see java.awt.event.MouseListener#mouseClicked(java.awt.event.MouseEvent)
         */
        public void mouseClicked( MouseEvent e )
        {

            JTableHeader h = ( JTableHeader ) e.getSource();
            TableColumnModel columnModel = h.getColumnModel();
            int viewColumn = columnModel.getColumnIndexAtX( e.getX() );
            int column = columnModel.getColumn( viewColumn ).getModelIndex();

            if ( column != NO_SELECTION )
            {
                int status = getSortingStatus( column );

                if ( !e.isControlDown() )
                {
                    cancelSorting();
                }

                // Cycle the sorting states through {NOT_SORTED, ASCENDING, DESCENDING} or
                // {NOT_SORTED, DESCENDING, ASCENDING} depending on whether shift is pressed.
                status = status + ( e.isShiftDown() ? NO_SELECTION : 1 );
                status = ( status + 4 ) % 3 - 1; // signed mod, returning {-1, 0, 1}
                setSortingStatus( column, status );
            }
        }
    }

    /**
     * The class Arrow implements Icon and is used to indicate sorting
     * up or down.
     *
     * @version $LastChangedRevision$
     * @author Thomas Brandes
     */
    private static class Arrow implements Icon
    {

        /**
         * Flag to indicate up or down.
         */
        private boolean descending;

        /**
         * Value for the size is used for the width and
         * height of this icon.
         */
        private int size;

        /**
         * Value for the priority.
         */
        private int priority;

        /**
         * Constructor for the arrow.
         *
         * @param valDescending is flag for descending
         * @param valSize becomes the size of he arrow
         * @param valPriority the value for priority
         */
        public Arrow( boolean valDescending, int valSize, int valPriority )
        {

            this.descending = valDescending;
            this.size = valSize;
            this.priority = valPriority;
        }

        /**
         * {@inheritDoc}
         *
         * @see javax.swing.Icon#paintIcon(java.awt.Component, java.awt.Graphics, int, int)
         */
        public void paintIcon( Component c, Graphics g, int x, int y )
        {

            Color color = Color.GRAY;

            if ( c != null )
            {

                color = c.getBackground();
            }

            // In a compound sort, make each succesive triangle 20%
            // smaller than the previous one.

            int dx = ( int ) ( size / 2 * Math.pow( 0.8, priority ) );

            int dy = descending ? dx : -dx;

            // Align icon (roughly) with font baseline.

            // y = y + 5 * size / 6 + (descending ? -dy : 0);

            y = y + 5 * size / 6;

            if ( descending )
            {

                y -= dy;
            }

            int shift = 1;

            if ( !descending )
            {

                shift = -shift;
            }

            // int shift = descending ? 1 : -1;

            g.translate( x, y );

            // Right diagonal.
            g.setColor( color.darker() );
            g.drawLine( dx / 2, dy, 0, 0 );
            g.drawLine( dx / 2, dy + shift, 0, shift );

            // Left diagonal.
            g.setColor( color.brighter() );
            g.drawLine( dx / 2, dy, dx, 0 );
            g.drawLine( dx / 2, dy + shift, dx, shift );

            // Horizontal line.
            if ( descending )
            {
                g.setColor( color.darker().darker() );
            }
            else
            {
                g.setColor( color.brighter().brighter() );
            }

            g.drawLine( dx, 0, 0, 0 );

            g.setColor( color );
            g.translate( -x, -y );
        }

        /**
         * {@inheritDoc}
         *
         * @see javax.swing.Icon#getIconWidth()
         */
        public int getIconWidth()
        {
            return size;
        }

        /**
         * {@inheritDoc}
         *
         * @see javax.swing.Icon#getIconHeight()
         */
        public int getIconHeight()
        {
            return size;
        }
    }

    /**
     * Help class for rendering.
     *
     * @version $LastChangedRevision$
     * @author Thomas Brandes
     */
    private class SortableHeaderRenderer implements TableCellRenderer
    {

        private TableCellRenderer tableCellRenderer;

        /**
         * TODO: brandes: enter comment!
         *
         * @param tableCellRenderer
         */
        public SortableHeaderRenderer( TableCellRenderer tableCellRenderer )
        {
            this.tableCellRenderer = tableCellRenderer;
        }

        /**
         * {@inheritDoc}
         *
         * @see javax.swing.table.TableCellRenderer#getTableCellRendererComponent(javax.swing.JTable, java.lang.Object, boolean, boolean, int, int)
         */
        public Component getTableCellRendererComponent( JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int column )
        {
            Component c = tableCellRenderer.getTableCellRendererComponent( table, value, isSelected, hasFocus, row, column );

            if ( c instanceof JLabel )
            {
                JLabel l = ( JLabel ) c;
                l.setHorizontalTextPosition( JLabel.LEFT );
                int modelColumn = table.convertColumnIndexToModel( column );
                l.setIcon( getHeaderRendererIcon( modelColumn, l.getFont().getSize() ) );
            }

            return c;
        }
    }

    /**
     * Private class for a directive that consists of integer values
     * for column and direction.
     *
     * @version $LastChangedRevision$
     * @author Thomas Brandes
     */
    private static class Directive
    {

        /**
         * Index of the column.
         */
        private int column;

        /**
         * Value for the direction.
         */
        private int direction;

        /**
         * Constructor for a directive.
         *
         * @param theColumn becomes local value for this instance.
         * @param theDirection becomes local value for this instance.
         */
        public Directive( int theColumn, int theDirection )
        {

            this.column = theColumn;
            this.direction = theDirection;
        }
    }
}
