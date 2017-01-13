/*
 * UserTable.java
 *
 * File contains class PMUserTable.
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

import org.apache.log4j.Logger;



/**
 * This class defines an abstact table for the performance events as
 * specified by the user. Each user counter is a rated basic, derived or compound
 * performance event. The values of the user table are fixed for the actually
 * selected process and thread.
 *
 * @version $LastChangedRevision$
 * @author Dr. Thomas Brandes
 */
class UserTable extends PMTableModel
{

    /**
     * Value for cell update with unknown row.
     */
    private static final int NO_ROW_VAL = -1;

    /**
     * These are the items for the popup menu of UserTable.
     */
    private static final String[] MENU_ITEMS = { "Show region", "Info region", "Select counter", "Info Counter", "Scale +", "Scale -" };

    /**
     * Logger for this class.
     */
    private static Logger logger = Logger.getLogger( UserTable.class );

    /**
     * This variable is the performance data from which the user counter
     * values will be displayed.
     */
    private PMData myPMData;

    /**
     * Pointer back to the implementation of ShowPM selection facilities.
     */
    private SelectedInterface myShowPM;

    /**
     * The user table will be created by providing performance data
     * and an implementation for the SelectedInterface.
     *
     * @param newPM is the performance data used for the user table
     * @param main is implementation of SelectedInterface
     */
    public UserTable( PMData newPM, SelectedInterface main )
    {

        myPMData = newPM;
        myShowPM = main;

    } // constructor PMUserTable

    /**
     * {@inheritDoc}
     *
     * @see javax.swing.table.TableModel#getColumnCount()
     */
    public int getColumnCount()
    {

        return myPMData.getNumberUserCounters() + 1;
    }

    /**
     * {@inheritDoc}
     *
     * @see javax.swing.table.TableModel#getRowCount()
     */
    public int getRowCount()
    {

        // number of rows is given by the number of regions

        return myPMData.getNumberRegions();
    }

    /**
     * {@inheritDoc}
     *
     * @see javax.swing.table.TableModel#getColumnName(int)
     */
    public String getColumnName( int col )
    {

        String name = null;

        if ( col == 0 )
        {

            name = "Region";

        }
        else
        {

            name = myPMData.getCounterHeader( true, col - 1 );
        }

        return name;
    }

    /**
     * {@inheritDoc}
     *
     * @see javax.swing.table.TableModel#getValueAt(int, int)
     */
    public Object getValueAt( int row, int col )
    {

        Object result;

        if ( col == 0 )
        {

            result = myPMData.getRegion( row ).getName();

        }
        else
        {

            int sP = myShowPM.getSelectedProcess();
            int sT = myShowPM.getSelectedThread();

            result = new Double( myPMData.getUserVal( sP, sT, row, col - 1 ) );
        }

        return result;
    }

    /**
     * {@inheritDoc}
     *
     * @see javax.swing.table.TableModel#getColumnClass(int)
     */
    public Class getColumnClass( int c )
    {

        return getValueAt( 0, c ).getClass();
    }

    /**
     * {@inheritDoc}
     *
     * @see javax.swing.table.TableModel#isCellEditable(int, int)
     */
    public boolean isCellEditable( int row, int col )
    {

        return false;

    }

    /**
     * {@inheritDoc}
     *
     * @see adaptor.ShowPM.PMTableModel#getPopupItems()
     */
    public String[] getPopupItems()
    {

        return MENU_ITEMS;
    }

    /**
     * {@inheritDoc}
     *
     * @see adaptor.ShowPM.PMTableModel#actionPopupItem(int, int, int)
     */
    public void actionPopupItem( int kind, int col, int row )
    {

        if ( kind == 0 )
        {

            // item "Show region"

            myShowPM.showFile( myPMData, row );

        }
        else if ( kind == 1 )
        {

            myShowPM.infoFile( myPMData, row );

        }
        else if ( col < 1 )
        {

            // all other actions are for counters, but the column with
            // the regions has been selected here

            String msg = MENU_ITEMS[kind];

            msg += " can only be selected for columns with counter";

            myShowPM.showErrorMessage( msg );

        }
        else if ( kind == 2 )
        {

            // item "Select counter"

            myShowPM.counterSelection( true, col - 1 );

        }
        else if ( kind == 3 )
        {

            // item "Info Counter"

            UserCounter c = myPMData.getUserCounter( col - 1 );

            String info = "counter info: " + c.getHeader() + " "
                          + c.getDescription() + " " + c.getFormula( 0 );

            // TODO: window that shows info about the user counter

            logger.info( info );

        }
        else if ( kind == 4 )
        {

            // item "Scale +"

            UserCounter c = myPMData.getUserCounter( col - 1 );

            c.scaleUp();

            this.fireTableDataChanged();

            this.fireTableCellUpdated( NO_ROW_VAL, col );

            // myShowPM.refresh();

        }
        else if ( kind == 5 )
        {

            // item "Scale -"

            UserCounter c = myPMData.getUserCounter( col - 1 );

            c.scaleDown();

            this.fireTableDataChanged();

            this.fireTableCellUpdated( NO_ROW_VAL, col );

            // myShowPM.refresh();

        }


    }

    /**
     * {@inheritDoc}
     *
     * @see adaptor.ShowPM.PMTableModel#mousePressed(int, int)
     */
    public void mousePressed( int  col, int row )
    {

        if ( col == 0 )
        {

            actionPopupItem( 0, col, row );

        }
        else
        {

            actionPopupItem( 1, col, row );
        }

    }

} // PMUserTable
