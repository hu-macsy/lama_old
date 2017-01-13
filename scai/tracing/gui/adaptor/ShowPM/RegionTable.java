/*
 * RegionTable.java
 *
 * Utitilies realizes some help functions for strings.
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


/**
 * PMRegionTable realizes a PMTableModel and is used for the presentation of
 * the counted events for the regions. The values are fixed for the actual
 * selected process and thread.
 *
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 *
 */
class RegionTable extends PMTableModel
{

    /**
     * These are the items for the popup menu of UserTable.
     */
    private static final String[] MENU_ITEMS = { "Show region", "Info region", "Select counter" };

    /**
     * This variable is the performance data from which the user counter
     * values will be displayed.
     */

    private PMData myPMData;

    /**
     * Pointer back to the implementation of ShowPM selection facilities.
     */
    private SelectedInterface mySelection;

    /**
     * The region table will be created by providing performance data
     * and an implementation for the SelectedInterface.
     *
     * @param newPM is the performance data used for the user table
     * @param main is implementation of SelectedInterface
     */
    RegionTable( PMData newPM, SelectedInterface main )
    {

        myPMData = newPM;
        mySelection = main;

    } // constructor PMRegionTable

    /**
     * {@inheritDoc}
     *
     * @see javax.swing.table.TableModel#getColumnCount()
     */
    public int getColumnCount()
    {

        return myPMData.getNumberRegionCounters() + 1;
    }

    /**
     * {@inheritDoc}
     *
     * @see javax.swing.table.TableModel#getRowCount()
     */
    public int getRowCount()
    {

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

            name = myPMData.getCounterHeader( false, col - 1 );
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

        int selThread = mySelection.getSelectedThread();
        int selProcess = mySelection.getSelectedProcess();

        long[] threadPM = myPMData.getRegionCounterVals( selProcess, selThread, row );

        if ( col == 0 )
        {

            return myPMData.getRegion( row ).getName();

        }
        else
        {

            return new Long( threadPM[col - 1] );
        }
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

            mySelection.showFile( myPMData, row );

        }
        else if ( kind == 1 )
        {

            mySelection.infoFile( myPMData, row );

        }
        else if ( col < 1 )
        {

            // all other actions are for counters, but the column with
            // the regions has been selected here

            String msg = MENU_ITEMS[kind];

            msg += " can only be selected for columns with counter";

            mySelection.showErrorMessage( msg );

        }
        else if ( kind == 1 )
        {

            // item "Select counter"

            mySelection.counterSelection( false, col - 1 );

        }

    }

    /**
     * {@inheritDoc}
     *
     * @see adaptor.ShowPM.PMTableModel#mousePressed(int, int)
     */
    public void mousePressed( int col, int row )
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

} // PMRegionTable
