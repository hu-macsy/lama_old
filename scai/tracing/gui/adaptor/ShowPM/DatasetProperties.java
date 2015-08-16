/*
 * DatasetProperties.java
 * 
 * Class containing some global properties how to make data sets.
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

import adaptor.General.PropertyMenu;

/**
 * Class containing some global properties how to make data sets.
 *
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */
public final class DatasetProperties {

    /**
     * This is the default value for the number of maximal entries taken in the dataset.
     */
    private static final int ENTRY_ITEMS_DEFAULT = 12;

    /**
     * This is the title for the depth menu (value in percentage that
     * is the threshold whether data comes into the data set or not).
     */
    private static final String DEPTH_TITLE = "Set Depth";
    
    /**
     * This is the title of the menu to select the maximal number of
     * entries.
     */
    private static final String ENTRY_TITLE = "Set Max Regions";
    
    /**
     * Logger for this class.
     */
    private static Logger logger = Logger.getLogger(DatasetProperties.class);

    /**
     * The array depthItems is a list of possible values for Chart depth.
     */
    private static final double[] DEPTH_ITEMS = { 10.0, 5.0, 2.0, 1.0, 0.5, 0.2, 0.1, 0.0 };

    /**
     * Property menu for the depth.
     */
    private static PropertyMenu propDepth = new PropertyMenu(DEPTH_TITLE, DEPTH_TITLE, "", DEPTH_ITEMS, 1.0);
    
    /**
     * This is the list of possible values for maximal entries.
     */
    private static final int[] ENTRY_ITEMS = { 0, 6, 8, 10, 12, 16, 20, 24, 28, 32 };

    /**
     * Property menu for the maximal number of entries.
     */
    private static PropertyMenu propEntry = new PropertyMenu(ENTRY_TITLE, ENTRY_TITLE, "", ENTRY_ITEMS, ENTRY_ITEMS_DEFAULT);

    /**
     * Overrides the default constructor. 
     */
    private DatasetProperties() {
        
    }
    
    /**
     * This routine sets up all the property menus needed. It should
     * be called in any case to make sure that the whole class will be
     * initialized.
     * 
     */
    static void makePropertyMenus() {
        
        logger.info("initialization routine");
    }
    
    /**
     * This routine is called to change the item of a 
     * property menu belonging to this class.
     * 
     * @param title is the title of the menu
     * @param cmd is the selected item of the menu
     */
    static void setProperty(String title, String cmd) {
        
        logger.info("setProperty " + cmd + " for " + title);
        
        if (title.equals(propDepth.getTitle())) {
            
            propDepth.setItem(cmd);
            propEntry.setNone();
            
            
        } else if (title.equals(propEntry.getTitle())) {
            
            propEntry.setItem(cmd);
            propDepth.setNone();
        }
        
    }
        
    /**
     * This routine returns the maximal number of entries as it
     * is currently selected.
     * 
     * @return integer value for the number of selected entries.
     */
    static int getSelectedEntries() {
      
        int entries = 0;
        
        Object obj = propEntry.getSelectedObject();
        
        if (obj instanceof Integer) {
            
            entries = ((Integer) obj).intValue();
            
        }
        
        return entries;
    }
    
    /**
     * This routine returns the value of depth as it is
     * currently selected.
     * 
     * @return double value in percentage for the threshold.
     */
    static double getSelectedDepth() {
  
        double depth = 0.0;
        
        Object obj = propDepth.getSelectedObject();
        
        if (obj instanceof Double) {
            
            depth = ((Double) obj).doubleValue();
            
        }
        
        return depth;
    }

}