/*
 * Counter.java
 *
 * Interface that has to be implemented by the RIF Editor.
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

/**
 * The <code>Counter</code> specifies the event to be counted and printed.
 * 
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */

class Counter {
    
    /**
     * This is the performance event that is counted.
     */
    private PerformanceEvent event;
    
    /**
     * Exclusive or inclusive counting.
     */
    private boolean exclusive = false;
    
    /**
     * User specified value for scaling. 
     */
    private String scale = "";
    
    /**
     * The width is the number of relevant digits.
     */
    private int width = 10;
    
    /**
     * The precision is the number of digits.
     */
    private int precision = 0;
    
    // construct a Counter by the PerformanceEvent to be counted
    
    /**
     * Make a counter for a performance event with default values.
     * 
     * @param countedEvent is the event that will be counted
     */
    public Counter(PerformanceEvent countedEvent) {
    
        event = countedEvent;
    }
    
    /**
     * This routine make a counter definition from a certain set of
     * input items. These items are usually from a line in the configuration
     * file.
     * 
     * @param counterItems are name mode width precision
     * @param config is pointer to supported events
     */
    public Counter(String[] counterItems, PMConfiguration config) {
        
        String eventName = counterItems[0];
        
        event = config.findSupportedEvent(eventName);
        
        // check for <rate><name>, e.g. mWALL_TIME
        
        if (event == null) {
            
            scale = eventName.substring(0, 1);
            
            eventName = eventName.substring(1);
            
            event = config.findSupportedEvent(eventName);
        }
        
        if (event == null) {
            
            // oops, may be we could not find any basic performance event
            
            throw new IllegalArgumentException(eventName + " unknown event");
        }
        
        if (counterItems.length > 1) {
          
            // mode: E for EXCLUSIVE or N for NETTO

            String mode = counterItems[1];
            
            exclusive = mode.startsWith("E") || mode.startsWith("N");
            
        }
        
        try {
            
            if (counterItems.length > 2) {
                
                width = Integer.parseInt(counterItems[2]);                
            }

        } catch (NumberFormatException e) {
            
            throw new IllegalArgumentException(counterItems[2] + ": width must be integer");
        }
        
        try {
            
            if (counterItems.length > 3) {

                precision = Integer.parseInt(counterItems[3]);
            }

        } catch (NumberFormatException e) {
            
            throw new IllegalArgumentException(counterItems[3] + ": precision must be integer");
        }
        
    } // constructor Counter with strings
    
    /**
     * This routine returns a duplcate of this counter. Only the mode
     * flag is switched.
     * 
     * @return a copy of this counter
     */
    public Counter copyCounter() {
        
        Counter theCopy = new Counter(event);
        theCopy.precision = precision;
        theCopy.width = width;
        theCopy.scale = scale;
        theCopy.exclusive = !exclusive;
        
        return theCopy;
        
    } // copyCounter
    
    /**
     * This routine returns the property of a counter.
     * 
     * @param pos is a value between 0 and 5 (six properties)
     * @return property of the counter as an object
     */
    public Object getCounterProperty(int pos) {

        // Counter has six properties

        Object property = null;

        switch (pos) {

        case 0:
            property = event.getName();
            break;

        case 1:
            property = event.getDescription();
            break;

        case 2:
            property = Boolean.FALSE;
            if (exclusive) {
                property = Boolean.TRUE;
            }
            break;

        case 3:
            property = scale;
            break;
            
        case 4:
            property = new Integer(width);
            break;
            
        case 5:
            property = new Integer(precision);
            break;

        default:
            property = null;
        }

        return property;

    } // method getCounterProperty
    
    /**
     * This routine can be used to set a new value for the counter.
     * 
     * @param pos is the position of the property
     * @param value is the new property
     */
    public void setCounterProperty(int pos, Object value) {
        
        if (pos == 2) {
            
            exclusive = ((Boolean) value).booleanValue();            
        }

        if (pos == 3) {

            scale = (String) value;
        }

        if (pos == 4) {
            
            width = ((Integer) value).intValue();            
        }

        if (pos == 5) {

            precision = ((Integer) value).intValue();
        }

    }
    
    /**
     * This routine generates for this counter a line that is
     * printed in the PM configuration file.
     * 
     * @return the line to be printed for this counter;
     */
    public String makeLine() {
        
        String exclString;
        
        if (exclusive) {
            
            exclString = "EXCL";            
            
        } else {

            exclString = "INCL";
        }
      
        return scale + event.getName() + " " + exclString + " " + width + " " + precision;
        
    } // makeLine
    
    /**
     * This routine delivers a string that represents a description of the counter
     * used in comment lines for the PM configuration file.
     * 
     * @return String with the description
     */
    public String getDescription() {
        
        String mode = null;
        
        if (exclusive) {
            
            mode = "(exclusive)";    
            
        } else {
            
            mode = "(inclusive)";            
        }
        
        return event.getName() +  " " + mode + " : " + event.getDescription();
        
    }
    
} // class Counter

