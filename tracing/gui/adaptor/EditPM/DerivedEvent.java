/*
 * EditRIFInterface.java
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

import java.awt.Color;

/**
 * The class DerivedEvent derives a basic event for the walltime or CPU time.
 *
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */
public class DerivedEvent extends PerformanceEvent {

    /**
     * myEventColor is a light blue.
     */
    private static final Color MY_EVENT_COLOR = Color.CYAN;
    
    /**
     * Kind of derived event, can be ' (for walltime) or " (for CPU time).
     */
    private String op;

    /**
     * This basic performance event will be derived.
     */
    private BasicEvent myBasicEvent;
    
    /**
     * Constructor for a new derived performance event.
     * 
     * @param name is the name of the new event
     * @param e1 is the performance event that will be derived
     * @param cop is the kind of derivation (wall or cpu time)
     */
    public DerivedEvent(String name, PerformanceEvent e1, String cop) {
                
        super(true, name);
        
        myBasicEvent = (BasicEvent) e1;
        op = cop;
    
    } // DerivedEvent

    /**
     * Getter for the color of a derived event.
     * 
     * @return the color unique for the derived event.
     */
    static Color getEventColor() {
        
        return MY_EVENT_COLOR;
    }

    /**
     * {@inheritDoc}
     *
     * @see adaptor.EditPM.PerformanceEvent#getColor()
     */
    public Color getColor() {
        
        return MY_EVENT_COLOR;
    }
    
    /**
     * {@inheritDoc}
     *
     * @see adaptor.EditPM.PerformanceEvent#getDescription()
     */
    public String getDescription() {
        
        return myBasicEvent.getName() + " " + op;
        
    }

    /**
     * Getter routine for the basic event that is derived.
     * 
     * @return the event that is derived by thise event.
     */
    public PerformanceEvent getEvent() {
        
        return myBasicEvent;
    }
    
    /**
     * Getter routine for the kind of derivation.
     * 
     * @return either ' or "
     */
    public String getKind() {
        
        return op;
    }
}
