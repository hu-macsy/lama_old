/*
 * RegionCounter.java
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

import adaptor.General.CountedEvent;


/**
 * A region counter is a counted event; it has as an additional attribute the
 * mode for inclusive or exclusive counting.
 *
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */
class RegionCounter extends CountedEvent {

    /**
     * The mode is either 'I' for inclusive or 'E' for exclusive.
     */
    private char mode; 

    /**
     * Constructor for a new region counter.
     * 
     * @param uname is the name of the counter.
     * @param index specifies the index in the counted values.
     * @param umode is either "I" or "E"
     */
    public RegionCounter(String uname, int index, String umode) {

        super(uname, "", index);
        
        mode = umode.charAt(0);
        
    }

    /**
     * This routine returns a header string for this region counter.
     * 
     * @return string with information to be used as header.
     */
    public String getHeader() {

        String description;

        if (mode == 'I') {
            
            mode = '+';
            
        } else if (mode == 'E') {
            
            mode = 'O';
        }

        description = getName();

        if (!description.equals("CALLS")) {
            
            description += " " + mode;
        }
        
        return description;
    }

    /**
     * This comparison function returns true if two region counters
     * count the same thing (not same values).
     * 
     * @param counter is the other region counter
     * @return 
     */
    public boolean sameRegionCounter(RegionCounter counter) {
       
        boolean same = (mode == counter.mode);
        
        if (same) {
            
            same = this.getName().equals(counter.getName());
            
        }

        return same;
        
    }

} // class RegionCounter
