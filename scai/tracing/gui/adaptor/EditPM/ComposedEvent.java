/*
 * ComposedEvent.java
 *
 * File contains class for a composed performance event.
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
 * The class ComposedEvent combines two basic events.
 *
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */
public class ComposedEvent extends PerformanceEvent
{

    /**
     * Color for composed events to distinguish them from others.
     */
    private static final Color COMPOSED_EVENT_COLOR = Color.pink;

    /**
     * This is the pointer to the first basic event.
     */
    private BasicEvent p1 = null;

    /**
     * This is the pointer to the second basic event.
     */
    private BasicEvent p2 = null;

    /**
     * The op string stands for the binary operation between
     * the two basic events.
     */
    private String op = "";

    /**
     * A composed event has a name and combines two basic events
     * to some kind of metric.
     *
     * @param name is the name of the composed event
     * @param e1 is the first basic event
     * @param cop is the binary operator
     * @param e2 is the second basic event
     */
    public ComposedEvent( String name, PerformanceEvent e1, String cop, PerformanceEvent e2 )
    {

        // name is always taken as ADAPTOR specific

        super( true, name );

        adaptorName = name;
        p1 = ( BasicEvent ) e1;
        p2 = ( BasicEvent ) e2;
        op = cop;

    } // ComposedEvent

    /**
     * This routine gives the color for composed events.
     *
     * @return the color for a composed event
     */
    public static Color getEventColor()
    {

        return COMPOSED_EVENT_COLOR;
    }

    /**
     * {@inheritDoc}
     *
     * @see adaptor.EditPM.PerformanceEvent#getColor()
     */
    public Color getColor()
    {

        return COMPOSED_EVENT_COLOR;
    }

    /**
     * {@inheritDoc}
     *
     * @see adaptor.EditPM.PerformanceEvent#getDescription()
     */
    public String getDescription()
    {

        return p1.getName() + " " + op + " " + p2.getName();
    }

    /**
     * This routine gets the binary operator for the composed event.
     *
     * @return Returns the op.
     */
    public String getOp()
    {
        return op;
    }


    /**
     * The first basic event.
     * @return Returns the p1.
     */
    public BasicEvent getP1()
    {
        return p1;
    }


    /**
     * The second basic event.
     * @return Returns the p2.
     */
    public BasicEvent getP2()
    {
        return p2;
    }

}
