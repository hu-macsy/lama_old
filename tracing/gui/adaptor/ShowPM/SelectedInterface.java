/*
 * SelectedInterface.java
 * 
 * This is the inferface for the selection of PM data to show.
 * 
 * Created: 05.04.2006 brandes <email>
 * Changed:
 * 
 * $Id$
 *
 * Copyright (C) 2006 Fraunhofer SCAI, Germany
 * 
 * All rights reserved
 * 
 * http://www.rcenviroment.de/
 */

package adaptor.ShowPM;

/**
 * This is the inferface for the selection of PM data to show.
 *
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */
public interface SelectedInterface {

    /**
     * Gets the currently selected process.
     * 
     * @return the selected process id
     */
    int getSelectedProcess();

    /**
     * Gets the currently selected thread.
     * 
     * @return index of the selected thread
     */
    int getSelectedThread();

    /**
     * Gets the currently selected user counter.
     * 
     * @return user counter, -1 if not selected
     */
    int getSelectedUserCounter();

    /**
     * Get the currently selected region counter.
     * 
     * @return region counter, -1 if not selected
     */
    int getSelectedRegionCounter();
    
    /**
     * Shows the region for the id where region data is
     * available in PMData.
     * 
     * @param data contains information about PM data
     * @param region is the region that should be shown.
     */
    void showFile(PMData data, int region);
    
    /**
     * Info about the region for the id where region data is
     * available in PMData.
     * 
     * @param data contains information about PM data
     * @param regionIndex is the index of the region that should be shown.
     */
    void infoFile(PMData data, int regionIndex);
    
    /**
     * Select a new user or region counter.
     * 
     * @param isUser true if it is a user counter, false for region counter
     * @param counter is the counter index
     */
    void counterSelection(boolean isUser, int counter);
    
    /**
     * Shows in an appropriate way an error message.
     * 
     * @param msg is a string containing the error message.
     */
    void showErrorMessage(String msg);

}