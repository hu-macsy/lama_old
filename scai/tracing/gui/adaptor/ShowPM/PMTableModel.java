/*
 * PMTableModel.java
 *
 * The PMTableModel extends the abstract table model of Java to
 * support additional features needed for ShowPM.
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

import javax.swing.table.AbstractTableModel;

/**
 * PMTableModel extends the abstract table model of Java to
 * support additional features needed for ShowPM. This includes
 * especially the popup menu and related call backs.
 *
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */
abstract class PMTableModel extends AbstractTableModel
{

    /**
     * This routine will return all supported items of the popup
     * menu used for the cells of the table.
     *
     * @return all names of items for the popup menu
     */
    public abstract String[] getPopupItems();

    /**
     * New routine that will be called if the mouse has been used
     * for selection of a value in this table.
     *
     * @param col is the selected column
     * @param row is the selected row
     */
    public abstract void mousePressed( int col, int row );
    /**
     * New routine that will be called if a popup item has been
     * selected for a cell of the table.
     *
     * @param kind is the item of the popup menu
     * @param col is the selected column
     * @param row is the selected row
     */
    public abstract void actionPopupItem( int kind, int col, int row );

}