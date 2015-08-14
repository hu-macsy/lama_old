/*
 * ArraySelection.java
 * 
 * TODO: brandes Enter comment!
 * 
 * Created: 21.04.2006 brandes <email>
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


/**
 * This final class contains routines for selection in array.
 *
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */
public final class ArraySelection {

    /**
     * This is constant value that indicates that no array element
     * is selected at all.
     */
    public static final int NO_SELECTION = -1;
    
    /**
     * Private constructor to hide the default constructor.
     * 
     */
    private ArraySelection() {        
    }
    
    /**
     * This rouine is used to find an array position for scaling item.
     * 
     * @param items is the array of strings
     * @param item is the searched item
     * @return the index position of the item (-1 if not found)
     */
    public static int findIndex(Object [] items, Object item) {
        
        int pos = NO_SELECTION;
        
        for (int i = 0; i < items.length; i++) {
            
            if (item.equals(items[i])) {
                
                pos = i;
                break;
            }
        }
        
        return pos;
        
    }

    /**
     * The routine getIndex tries to find the index of a String in an 
     * array of strings. It is used for find positions of selected items of
     * menus.
     * 
     * @param item is the name that is search
     * @param items is the array of strings
     * @return the index of item in the array items
     * @throws NoSuchFieldException if item has not been found
     */

    public static int getIndex(Object[] items, Object item) throws NoSuchFieldException {
                
        for (int i = 0; i < items.length; i++) {
            
            if (item.equals(items[i])) {
                
                return i;
            }
        }
        
        throw new NoSuchFieldException("getIndex failed for " + item);
        
    }
    
    /**
     * This routine make an array of Objects from an array of integer
     * values.
     * 
     * @param val is an integer array.
     * @return is an array with the same size of Integer objects.
     */
    public static Object[] makeObjects(int[] val) {
        
        Object[] values = new Object[val.length];
        
        for (int i = 0; i < val.length; i++) {
            
            values[i] = new Integer(val[i]);
        }
        
        return values;
    }
 
    /**
     * This routine make an array of Objects from an array of double
     * values.
     * 
     * @param val is a double array.
     * @return an array with the same size of Double objects.
     */
    public static Object[] makeObjects(double[] val) {
        
        Object[] values = new Object[val.length];
        
        for (int i = 0; i < val.length; i++) {
            
            values[i] = new Double(val[i]);
        }
        
        return values;
    }


}
