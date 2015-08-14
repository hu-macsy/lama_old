/*
 * BasicEvent.java
 *
 * Basic performance event.
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

import java.awt.Color;

/**
 * The class BasicEvent is a basic performance event. 
 *
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */
public class BasicEvent extends PerformanceEvent {

    /**
     * Every basic event gets yellow color.
     */
    private static final Color BASIC_EVENT_COLOR = Color.YELLOW;

    /**
     * The long name describes the performance event in more detail.
     */
    private String longName = null;

    /**
     * Constructor for a new basic event.
     * 
     * @param isAdaptor indicates that this event is supported by ADAPTOR
     * @param sname is the short name of the basic event
     * @param lname is the description of the event
     */
    public BasicEvent(boolean isAdaptor, String sname, String lname) {
    
        super(isAdaptor, sname);
    
        longName = lname;
    
    } // constructor BasicEvent

    /**
     * This routine gets a color for this kind of event.
     * 
     * @return the color for this event
     */
    public static Color getEventColor() {

        return BASIC_EVENT_COLOR;

    }

    /**
     * {@inheritDoc}
     *
     * @see adaptor.EditPM.PerformanceEvent#getColor()
     */
    public Color getColor() {

        Color c = BASIC_EVENT_COLOR;

        if (isADAPTORCounter()) {

            c = c.darker();
        }

        return c;
    }

    /**
     * {@inheritDoc}
     *
     * @see adaptor.EditPM.PerformanceEvent#getDescription()
     */
    public String getDescription() {
        
        return longName;
    }
}
