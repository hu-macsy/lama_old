/*
 * UserCounter.java
 * 
 * Class UserCounter stands for a user defined metric on counted events.
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

import adaptor.General.CounterMetric;


/**
 * UserCounter is a metric that combines one or two counted events and rates the values.
 *
 * @version $LastChangedRevision$
 * @author Dr. Thomas Brandes
 */
public class UserCounter extends CounterMetric {

    /**
     * A user counter has the additional attribute whether it is inclusvie
     * or exclusive .
     */
    private boolean isExclusive = false; 

    /**
     * Constructor for a user counter.
     * 
     * @param uname the name of the user counter/metric
     * @param umode mode of the counter
     * @param scaling used for scaling
     * @param rate rating value of the counted events
     * @param counter is the RegionCounter that is used for the metric
     */
    UserCounter(String uname, String umode, String scaling, double rate, RegionCounter counter) {

        super(uname, rate, counter);

        setMode(umode);
        setScaling(scaling);
        
    }

    // user counter for two composed region counters

    /**
     * Constructor for a binary user counter.
     * 
     * @param uname the name of the user counter/metric
     * @param umode mode of the counter
     * @param scaling used for scaling
     * @param op binary operator
     * @param rate1 rating for first counter
     * @param counter1 first region counter
     * @param rate2 rating for second counter
     * @param counter2 rating for second counter
     */
    UserCounter(String uname, String umode, String scaling, String op,
                double rate1, RegionCounter counter1, 
                double rate2, RegionCounter counter2) {
        
        super(uname, op.charAt(0), rate1, counter1, rate2, counter2);
        
        setMode(umode);
        setScaling(scaling);
 
    }
    
    /**
     * This routine returns a description of the user counter needed
     * for the tables.
     * 
     * @return a string with name and mode of the user counter
     */
    String getHeader() {

        String header;

        char mode;
        
        if (isExclusive) {
            
            mode = 'O'; 
        
        } else {
            
            mode = '+';
        }

        header = getName();

        if (!header.equals("CALLS")) {
            
            header += " " + mode;          
        }

        String scaling = getScaling();
        
        if (!scaling.equals("-")) {
            
            header += " (" + scaling + ")";
        }
        
        return header;
    }

    /**
     * This routine converts the string for the counter mode in a boolean flag.
     * 
     * @param umode is a String describing the counter mode
     */
    private void setMode(String umode) {
        
        char mode = ' ';
        
        if (umode != null) {
            
            if (umode.length() > 0) {
                
                mode = umode.charAt(0);
                
            }
           
        }
        
        isExclusive = (mode == 'E') || (mode == 'e') || (mode == 'O');
    
    }

} // class UserCounter
