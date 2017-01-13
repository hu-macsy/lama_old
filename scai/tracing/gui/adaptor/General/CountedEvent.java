/*
 * CountedEvent.java
 *
 * Class is used for a counted event during the performance analysis.
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

package adaptor.General;

import org.apache.log4j.Logger;

/**
 * CountedEvent is a class that stands for a counted event with its total value
 * (summarized over the whole data file). The value of a counted event is always
 * of type long.
 *
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */
public class CountedEvent
{

    /**
     * Logger for this class.
     */
    private static Logger logger = Logger.getLogger( CountedEvent.class );

    /**
     * Name of the counted event.
     */
    private String myName = null;

    /**
     * Detailled description of the event.
     */
    private String myDescription = null;

    /**
     * This is the total value of the event (especially used for
     * relative values).
     */
    private long myTotalValue = ( long ) 0;

    /**
     * We set the correct scaling before it will be used the first time.
     */
    private String myScaling = null;

    /**
     * Each counted event gets an index so that we know which
     * position we have to take in the array of counter values.
     * CounterData.Counters[Index] == this. The default value -1
     * is used to see that index has not been set yet.
     */
    private int myIndex = ArraySelection.NO_SELECTION; /*  */

    /**
     * Constructor for a new counted event.
     *
     * @param theName of the event
     * @param theDescription is a detailed description
     * @param theIndex is index, 0 <= index < # events
     */
    public CountedEvent( String theName, String theDescription, int theIndex )
    {

        myName        = theName;
        myDescription = theDescription;
        myIndex       = theIndex;
    }

    /**
     * getter for the total value.
     *
     * @return the total value
     */
    public long getTotalValue()
    {

        return myTotalValue;
    }

    /**
     * Setter routine for the total value.
     *
     * @param value
     *            is the value to set
     */
    public void setTotalValue( long value )
    {

        myTotalValue = value;

    }

    /**
     * This routine determines for this event the scaling that is used for auto scaling.
     */
    private void setScaling()
    {

        myScaling = Scaling.getAutoScaleOrder( ( double ) myTotalValue, false );

        logger.info( "Counter " + myName + " has total " + myTotalValue + ", autoscale = " + myScaling );

    }

    /**
     * Adder routine for the total value.
     *
     * @param value is the value to add
     */
    public void addTotalValue( long value )
    {

        myTotalValue += value;
    }

    /**
     * This routine gets the counter value for the event. It uses
     * the index of the counted event to pick up the right value.
     *
     * @param allCounterVals is the array with counter values for all counted events
     * @return the corresponding value for this event.
     */
    public long getCounterValue( long[] allCounterVals )
    {

        return allCounterVals[myIndex];
    }

    /**
     * From the array of counter values we determine for this
     * counter the relative costs.
     *
     * @param allCounterVals are the counter values
     * @return the relative costs (between 0.0 and 1.0)
     */
    public double getRelativeValue( long[] allCounterVals )
    {

        long value = allCounterVals[myIndex];

        double val = 0.0;  // used if total is zero

        if ( myTotalValue > 0 )
        {

            val = ( double ) value / ( double ) myTotalValue;

        }

        return val;

    }
    /**
     * Get the index position of this event.
     *
     * @return index value, 0 <= index < # events
     */
    public int getIndex()
    {

        return myIndex;
    }

    /**
     * Getter routine for the name.
     *
     * @return the name of this event
     */
    public String getName()
    {

        return myName;
    }

    /**
     * Getter routine for the description.
     *
     * @return the detailed description of this event
     */
    public String getDescription()
    {

        return myDescription;
    }

    /**
     * {@inheritDoc}
     *
     * @see java.lang.Object#toString()
     */
    public String toString()
    {

        return myName;
    }

    /**
     * This function returns the value for this event as
     * a string where the value is scaled and rounded. The
     * auto-scale or global-scale factor is appended.
     *
     * @param counterVals are the values of all counted events.
     * @return the value for this event as a scaled string.
     */
    public String getAbsValueString( long[] counterVals )
    {

        double val = getCounterValue( counterVals );

        if ( myScaling == null )
        {

            setScaling();

        }

        return Scaling.getValueString( val, myScaling );

    }

    /**
     * This function returns the relative value for this event as
     * a string where the value is scaled and rounded. Scaling
     * is done in percent.
     *
     * @param counterVals are the values of all counted events.
     * @return the value for this event as a percent string.
     */
    public String getRelValueString( long[] counterVals )
    {

        return Scaling.getPercentString( getRelativeValue( counterVals ) );
    }
}

