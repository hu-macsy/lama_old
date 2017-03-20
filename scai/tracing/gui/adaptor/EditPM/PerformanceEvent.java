/*
 * PerformanceEvent.java
 *
 * PerformanceEvent is an abstract class used for BasicEvent, Composed and Derived Event.
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
 * The class PerformanceEvent.
 *
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */

public abstract class PerformanceEvent
{

    /**
     * adaptorName is an additional name for an event that
     * is not system-dependent.
     */
    protected String adaptorName = null; // ADAPTOR name of the event

    /**
     * The systemName is a unique name for an event that might
     * be machine-specific.
     */
    private String systemName = null; // HPM name of the event

    /**
     * Constructor for a performance event.
     *
     * @param isADAPTOR is true if event is not system but ADAPTOR specific
     * @param sname is the name of the event
     */
    public PerformanceEvent( boolean isADAPTOR, String sname )
    {

        if ( isADAPTOR )
        {
            adaptorName = sname;
        }
        else
        {
            systemName = sname;
        }

    } // constructor PerformanceEvent

    /**
     * Every performance event has an color to be used in
     * the frames for events.
     *
     * @return color for the performance event.
     */
    abstract Color getColor();

    /**
     * getDescription returns a string that describes the performance event.
     * It is not implemented here as the kind of how to build it depends
     * on the kind of performance event.
     *
     * @return string with long description of the performance event.
     */
    abstract String getDescription();

    /**
     * This routine returns the name of the pereformance event. If
     * a more general name is available (adaptor name) this one is
     * taken.
     *
     * @return name of the performance event
     */
    public String getName()
    {

        String name = systemName;

        // if adaptorName is available we take this as the short name

        if ( adaptorName != null )
        {

            name = adaptorName;
        }

        return name;

    }

    /**
     * This routine returns true if it is a general ADAPTOR event.
     *
     * @return true if this counter is ADAPTOR specific and not machine specific.
     */
    public boolean isADAPTORCounter()
    {

        return ( adaptorName != null );
    }

    /**
     * This routine checks whether this event has a certain name.
     *
     * @param name that should be my name
     * @return true if this event has the given name
     */
    boolean hasName( String name )
    {

        boolean found = false;

        if ( name.equals( adaptorName ) )
        {
            found = true;
        }
        else if ( name.equals( systemName ) )
        {
            found = true;
        }

        return found;
    }

} // class

