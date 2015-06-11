/*
 * PropertyMenu.java
 * 
 * This class used for property menus containing a title, a number of
 * items and a current selection.
 * 
 * Created: 04.05.2006 Thomas Brandes <Thomas.Brandes@scai.fraunhofer.de>
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
 
package adaptor.General;

import java.util.List;
import java.util.Vector;

import org.apache.log4j.Logger;


/**
 * This class used for property menus containing a title, a help info,
 * a number of items and a current selection.
 *
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */
public class PropertyMenu {
    
    /**
     * This is constant value that indicates that no element
     * is selected at all.
     */
    public static final int NO_SELECTION = -1;

    /**
     * Logger for this class. Here it is very helpful to give 
     * error messages about wrong menus in the application.
     */
    private static Logger logger = Logger.getLogger(PropertyMenu.class);
    
    /**
     * We collect all defined property menus in a list to have access
     * to all defined menus.
     */
    private static List<PropertyMenu> allMenus = new Vector<PropertyMenu>();
    
    /**
     * This is the header of this property menu. Menus with the same
     * header can be grouped.
     */
    private String myHeader;
    
    /**
     * This is the title of the property menu.
     */
    private String myTitle;
    
    /**
     * This string explains the content of the property menu.
     */
    private String myHelp;
    
    /**
     * This array contains all possible entries for this property menu.
     */
    private Object[] myItems;
    
    /**
     * These are the string represensations of the items. 
     */
    private String[] myItemStrings;
    
    /**
     * This is the current selection, i.e. myItems[myCurrentSelection]
     * is the selected value.
     */
    private int myCurrentSelection = NO_SELECTION;
    
    /**
     * Constructor for a property menu by header, title, possible items and 
     * the default selection.
     * 
     * @param header is the group to which the new property menu belongs
     * @param title is the title of the property menu
     * @param help is a tip what the property menu is about.
     * @param theItems These are the objects of all possible selections
     *        for the property.
     * @param defaultItem This is the default object for the property.
     */
    public PropertyMenu(String header, String title, String help, Object[] theItems, Object defaultItem) {
        
        myHeader = header;
        
        myTitle = title;
        
        myItems = theItems;
        
        myHelp = help;
        
        myItemStrings = new String [myItems.length];
        
        for (int i = 0; i < myItems.length; i++) {
            
            myItemStrings[i] = myItems[i].toString();
        }
        
        if (defaultItem != null) {
            
            myCurrentSelection = ArraySelection.findIndex(myItems, defaultItem);
            
            if (myCurrentSelection == ArraySelection.NO_SELECTION) {
                
                logger.error("constructor: default val " + defaultItem + " not in items");
                
                // So we take the first array value as default
                
                myCurrentSelection = NO_SELECTION;
            }           
        }
        
        allMenus.add(this);
        
    }
    
    /**
     * Constructor for a property menu without the default item.
     * 
     * @param header is the group to which the new property menu belongs
     * @param title is the title of the property menu
     * @param help is a tip what the property menu is about.
     * @param theItems These are the string representations of all possible selections
     *        for the property.
     */

    public PropertyMenu(String header, String title, String help, Object[] theItems) {
    
        this(header, title, help, theItems, null);
        
    }
 
    /**
     * Constructor for a property menu that contains double values.
     * 
     * @param header is the group to which the new property menu belongs
     * @param title is the title of the property menu
     * @param help is a tip what the property menu is about.
     * @param items These are all possible values that can be selected
     *        for the property.
     * @param defaultItem This is the default value for the property.
     */

    public PropertyMenu(String header, String title, String help, double[] items, double defaultItem) {
    
        this(header, title, help, ArraySelection.makeObjects(items), new Double(defaultItem));
                
    }
 
    /**
     * Constructor for a property menu that contains integer values.
     * 
     * @param header is the group to which the new property menu belongs
     * @param title is the title of the property menu
     * @param help is a tip what the property menu is about.
     * @param items These are all possible values that can be selected
     *        for the property.
     * @param defaultItem This is the default value for the property.
     */
    
    public PropertyMenu(String header, String title, String help, int[] items, int defaultItem) {

        this(header, title, help, ArraySelection.makeObjects(items), new Integer(defaultItem));
      
    }

    /**
     * This routine searches for the position of a given item in the array of 
     * possible items.
     * 
     * @param item is the searched item
     * @param msg is the name of the routine used for error messages
     * @return the position of the item.
     */
    private int getItemPos(Object item, String msg) {
        
        int pos;
        
        if (item == null) {
            
            // this is legal to make a no selection
            
            return NO_SELECTION;
            
        } else if (item instanceof String) {
            
            pos = ArraySelection.findIndex(myItemStrings, item);
            
        } else {
            
            pos = ArraySelection.findIndex(myItems, item);            
        }

        if (pos < 0) {
            
            logger.error(msg + ": val " + item + " not in items of " + myTitle + "(" + myHeader + ")");
            
            pos = myCurrentSelection;
        }

        return pos;
    }
   
    /**
     * This routine sets an item as the current selection of the property menu.
     * 
     * @param item is the object that should be selected now.
     */
    public void setItem(Object item) {
        
        myCurrentSelection = getItemPos(item, "setItem");
    }
    
    /**
     * This routine sets the current selection of the property menu to the last
     * item in the list.
     * 
     */
    public void setLastItem() {
        
        myCurrentSelection = myItems.length - 1;
        
    }
    
    /**
     * This routine deletes the current selection and there will be no selection
     * for the property.
     * 
     */
    public void setNone() {
        
        myCurrentSelection = NO_SELECTION;
    }
   
    /**
     * This routine sets an item as the current selection of the property menu.
     * 
     * @param item is the item to be set as the current selection.
     * @return true if the selection has changed.
     */
    public boolean changeItem(Object item) {
        
        int oldPos = myCurrentSelection;
        
        myCurrentSelection = getItemPos(item, "changeItem");
        
        return (oldPos != myCurrentSelection);                
    }
    
    /**
     * Asking whether in this property menu a given item is the current selection.
     * 
     * @param item is the item we ask for
     * @return true if item is the current selection.
     */
    public boolean isSelected(Object item) {
        
        int pos = getItemPos(item, "isSelected");

        return (pos == myCurrentSelection);
        
    }

    /**
     * Getter routine for the current selection.
     * 
     * @return the string representation of the current selection.
     */
    public String getSelection() {
        
        if (myCurrentSelection == NO_SELECTION) {
            
            return "";
            
        } else {

            return myItemStrings[myCurrentSelection];

        }             
    }

    /**
     * Getter routine for the current selection.
     * 
     * @return the object of the current selection.
     */

    public Object getSelectedObject() {
        
        if (myCurrentSelection == NO_SELECTION) {
            
            return null;
            
        } else {

            return myItems[myCurrentSelection];

        }        
    }
    
    /**
     * Getter routine for the info about the property menu.
     * 
     * @return the info as a String for this property menu.
     */
    public String getHelp() {
        
        if (myHelp.equals("")) {
            
            return "No help available for " + myTitle;
            
        }
        
        return myHelp;
        
    }
    
    /**
     * Getter routine for the title of this menu.
     * 
     * @return the title of this menu.
     */
    public String getTitle() {
        
        return myTitle;
        
    }

    /**
     * Getter routine for the items of the property menu.
     * 
     * @return an array of String for the items.
     */
    public String[] getItems() {
      
        return myItemStrings;
    }
    
    /**
     * This routine returns the headers of all property menus that
     * have been defined. Double values are eliminated.
     * 
     * @return an array with String representation of all headers.
     */
    public static String[] getAllHeaders() {
        
        String [] headers = new String [allMenus.size()];
        
        int n = 0;
        
        for (int i = 0; i < headers.length; i++) {
            
            // make sure that the same header does not exist
            
            boolean found = false;
            
            String header = ((PropertyMenu) allMenus.get(i)).myHeader;
            
            for (int j = 0; j < n; j++) {
                
                if (headers[j].equals(header)) {
                    
                    found = true;
                }
            }
            
            if (!found) {
                
                headers[n++] = header;
            }
        }
        
        String [] headers1 = new String [n];
        
        for (int i = 0; i < n; i++) {
            
            headers1[i] = headers[i];
        }
        
        return headers1;
    }
    
    /**
     * This routine returns all titles of the property menus for a given header.
     * 
     * @param header is the header for which title menus are searched
     * @return an array of string representation with all titles
     */
    public static String[] getAllTitles(String header) {
        
        String [] titles = new String [allMenus.size()];
        
        int n = 0;
        
        for (int i = 0; i < titles.length; i++) {
            
            // make sure that the same header does not exist
            
            PropertyMenu prop = (PropertyMenu) allMenus.get(i);
            
            if (prop.myHeader.equals(header)) {
                
                titles[n++] = prop.myTitle;
            }
        }
        
        String [] titles1 = new String [n];
        
        for (int i = 0; i < n; i++) {
            
            titles1[i] = titles[i];
        }
        
        return titles1;

    }

    /**
     * This routine searches a property menu with given header name and
     * given title.
     * 
     * @param header is the name of the header
     * @param title is the name of the title
     * @return a property menu with the given header and the given title
     */
    public static PropertyMenu findMenu(String header, String title) {
        
        int n = allMenus.size();
        
        for (int i = 0; i < n; i++) {
            
            // make sure that the same header does not exist
            
            PropertyMenu prop = (PropertyMenu) allMenus.get(i);
            
            if (prop.myHeader.equals(header)) {
                
                if (prop.myTitle.equals(title)) {
                    
                    return prop;
                }
            }
        }
        
        logger.error("findMenu: could not find menu with header " + header
                     + ", title = " + title);
        
        return null;

    }
    
}
