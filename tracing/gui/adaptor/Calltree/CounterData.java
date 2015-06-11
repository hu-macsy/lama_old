/*
 * CounterData.java
 * 
 * This class contains all data related to the counted events and the metrics.
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

package adaptor.Calltree;

import java.io.BufferedWriter;
import java.io.IOException;

import org.apache.log4j.Logger;

import adaptor.General.CountedEvent;
import adaptor.General.CounterMetric;

/**
 * This class contains all data related to the counted events and the metrics.
 * It is completely static so it can be used globally.
 *
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */
public final class CounterData {

    /**
     * Logger for this class.
     */
    private static Logger logger = Logger.getLogger(CalltreePM.class);

    /**
     * Initial size for number of counters.
     */
    private static final int MAX_COUNTERS = 20;
    
    /**
     * Initial size for number of metrics.
     */
    private static final int MAX_METRICS = 20;

    /**
     * Number of events that are currently defined.
     */
    private static int noEvents = 0;
    
    /**
     * Number of metrics that are currently defined.
     */
    private static int noMetrics = 0;

    /**
     * This is the array of all events that have been counted by the PM.
     */
    private static CountedEvent[] theCountedEvents = new CountedEvent[MAX_COUNTERS];
    
    /**
     * This is an array of all metrics that can be used for the counted events.
     */
    private static CounterMetric[] theCounterMetrics = new CounterMetric[MAX_METRICS];

    /**
     * Private constructor that is never used.
     * 
     */
    private CounterData() {
        
    }
    
    /**
     * This routine defines a new event with the corresponding index. The calls
     * must be in ascending order starting with 0.
     * 
     * @param index is the index for the new event
     * @param eventName is the name of the event
     * @param eventDescription is the detailled description of the event
     * @return true if call of the routine was without error
     */
    public static boolean defineEvent(int index, String eventName, String eventDescription) {
        
        boolean okay = false;
        
        if (index < noEvents) {
            
            // okay if the name fits with name already defined before
            
            okay = eventName.equals(theCountedEvents[index].getName());
            
            if (!okay) {
                
                logger.error("Event conflict for event " + index + ", defined as " + eventName + ", "
                             + "but was " + theCountedEvents[index].getName() + " before");
            }
            
        } else if (index > noEvents) {
            
            // not okay 
            
            logger.error("CounterData: define Event " + index + ", but must be " + noEvents);
            
        } else {
            
            okay = true;
            
            // make sure that array for Counters is big enough
            
            if (noEvents + 1 >= theCountedEvents.length) {
                
                int size = theCountedEvents.length;
                
                CountedEvent[] newCounters = new CountedEvent[2*size];
                
                for (int i = 0; i < size; i++) {
                    
                    newCounters[i] = theCountedEvents[i];
                    newCounters[i + size] = null;
                }
                
                logger.info("doubled size of array for counters, new size = " + 2*size);
                
                theCountedEvents = newCounters;
            }
            
            
            theCountedEvents[index] = new CountedEvent(eventName, eventDescription, index);
            
            noEvents++;
        }
        
        return okay;
        
    } // defineEvent

    /**
     * This routine tries to find a metric with the given name.
     * 
     * @param searchedMetricName is the metric
     * @return the metric if found otherwise null
     */
    private static CounterMetric findMetric(String searchedMetricName) {

        CounterMetric resultMetric = null;
        
        for (int i = 0; i < noMetrics; i++) {

            CounterMetric metric = theCounterMetrics[i];

            if (searchedMetricName.equals(metric.getName())) {
                
                resultMetric = metric;
                break;
            }
        }

        return resultMetric;

    } // getMetric

    /**
     * This routine adds a new CounterMetric to the counter data.
     * 
     * @param newMetric is the new CounterMetric to set
     */
    public static void defineMetric(CounterMetric newMetric) {

        // check whether the metric has already been defined

        CounterMetric oldMetric = findMetric(newMetric.getName());

        if (oldMetric != null) {

            // we should make sure that we have the same metrics

            String newFormula = newMetric.getFormula(0);
            String oldFormula = oldMetric.getFormula(0);

            if (newFormula.equals(oldFormula)) {
              
                return;
            }

            logger.error("defineMetric mismatch for metric " + oldMetric.getName());
            logger.error("   was: " + newFormula);
            logger.error("   now: " + oldFormula);
            System.exit(0);

        }

        if (noMetrics + 1 >= MAX_METRICS) {

            logger.error("maximal number of metrics (derived events) reached");
            System.exit(0);

        }

        theCounterMetrics[noMetrics++] = newMetric;

    } // defineMetric

    /**
     * This routine is used to verify that the number of counted events
     * has the right number.
     * 
     * @param n is the assumed number of events
     * @param msg is a message to print in case of error
     */
    public static void checkCounters(int n, String msg) {

        if (n != noEvents) {


            logger.error("Counters mismatch (have " + noEvents + ", but only " + n + " values): " + msg);

            System.exit(0);

        }

    } // CheckCounters

    /**
     * This routine sets the total values for all counted events.
     * 
     * @param totals are the total values
     */
    public static void setTotals(long[] totals) {

        checkCounters(totals.length, "set totals");

        for (int i = 0; i < totals.length; i++) {

            theCountedEvents[i].addTotalValue(totals[i]);
        }


    } // setTotals

    /**
     * This routine searches for a counted event with a given name.
     * 
     * @param name is the name of the searched event
     * @return pointer to the counted event (can be null)
     */
    static CountedEvent getCountedEvent(String name) {
    
        for (int i = 0; i < noEvents; i++) {
    
            CountedEvent event = theCountedEvents[i];
    
            if (name.equals(event.getName())) {
    
                return event;        
            }
        }
    
        return null; // not found 
    
    } // getCountedEvent

    /**
     * Getter routine for the number of known events.
     * 
     * @return the number of known events
     */
    public static int numberEvents() {
        
        return noEvents;
        
    }

    /**
     * This routine returns an array of strings with all the names of the 
     * counted events.
     * 
     * @return array with names of counted events
     */
    public static String[] getCounterItems() {

        String[] theCounterStrings = new String[noEvents];

        for (int i = 0; i < noEvents; i++) {

            theCounterStrings[i] = theCountedEvents[i].getName();

        }

        return theCounterStrings;

    } // getCounterItems

    /**
     * Get the name of a counted event.
     * 
     * @param counterIndex is the index
     * @return the name as a String
     */
    public static String getEventName(int counterIndex) {
    
        return theCountedEvents[counterIndex].getName();
    
    } // getEventName

    /**
     * Getter routine for the number of known metrics.
     * 
     * @return the number of known metrics
     */
    public static int numberMetrics() {
        
        return noMetrics;
        
    }

    /**
     * This routine returns an array of strings with all the names of the counter
     * metrics that have been defined.
     * 
     * @return array with the names of counter metrics
     */
    public static String[] getMetricItems() {

        String[] theMetricStrings = new String[noMetrics];

        int i;

        for (i = 0; i < noMetrics; i++) {

            theMetricStrings[i] = theCounterMetrics[i].getName();

        }

        return theMetricStrings;

    } // getMetricItems

    /**
     * Get the name of a counter metric.
     * 
     * @param metricIndex is the index
     * @return the name as a String
     */
    public static String getMetricName(int metricIndex) {

        return theCounterMetrics[metricIndex].getName();

    } // getMetricName

    /**
     * Get the description of an event.
     * 
     * @param counterIndex is the index of the counted event
     * @return String containing a description.
     */
    public static String getEventDescription(int counterIndex) {

        return theCountedEvents[counterIndex].getDescription();

    } // getEventDescription

    /**
     * Get the description of a metric.
     * 
     * @param metricIndex is the index of the metric
     * @return String containing a description.
     */
    public static String getMetricDescription(int metricIndex) {

        return theCounterMetrics[metricIndex].getDescription();

    } // getMetricDescription

    /**
     * This routine returns the formula used for a metric.
     * 
     * @param metricIndex is the metric
     * @param offset is an offset in the counter table 
     * @return string containing the formula for the metric
     */
    public static String getMetricFormula(int metricIndex, int offset) {

        return theCounterMetrics[metricIndex].getFormula(offset);

    } // getMetricFormula 

    /**
     * This routine returns the relative value for a counter value compared
     * against the total value of the counter.
     * 
     * @param counterIndex is the index of the counter
     * @param value is the value
     * @return the relative value (between 0.0 and 1.0)
     */
    static double getRelativeValue(int counterIndex, long value) {

        long total = theCountedEvents[counterIndex].getTotalValue();

        double val = 0.0;  // used if total is zero
        
        if (total > 0) {
            
            val = (double) value / (double) total;
            
        }

        return val;

    }

    /**
     * This routine returns the value  of the counter metric as a string.
     * 
     * @param metricIndex is the index of the metric
     * @param values are the values of the counted events
     * @return the metric value as a string
     */
    static String getMetricValueString(int metricIndex, long[] values) {

        return theCounterMetrics[metricIndex].getValueString(values);
    }

    /**
     * This routine returns the value  of the counter metric as a double value.
     * 
     * @param metricIndex is the index of the metric
     * @param values are the values of the counted events
     * @return the metric value
     */    
    static double getMetricValue(int metricIndex, long[] values) {

        return theCounterMetrics[metricIndex].getValue(values);
    }
    
    /**
     * This routine returns the relative value of the counter metric as a string.
     * 
     * @param metricIndex is the index of the metric
     * @param values are the values of the counted events
     * @return the relative value
     */
    static double getRelativeMetric(int metricIndex, long[] values) {

        return theCounterMetrics[metricIndex].getRelativeValue(values);
    }

    /**
     * This routine writes the description for all counted events and for
     * the metrics in a buffered output file. The format is used for .pms files.
     * 
     * @param buff is the buffer where to write
     * @throws IOException in case of error
     */
    static void writeCounters(BufferedWriter buff) throws IOException {

        buff.write("NC=" + 2 * noEvents);
        buff.newLine();

        // all INCL counters at first 

        for (int i = 0; i < noEvents; i++) {

            buff.write(getEventName(i));
            buff.write(" INCL");
            buff.newLine();
        }

        // now all EXCL counters 

        for (int i = 0; i < noEvents; i++) {

            buff.write(getEventName(i));
            buff.write(" EXCL");
            buff.newLine();
        }

        buff.write("NU=" + 2 * noMetrics);
        buff.newLine();

        for (int i = 0; i < noMetrics; i++) {
            
            buff.write(getMetricName(i));
            buff.write(" INCL");
            buff.write(" - 10 3 ");
            buff.write(getMetricFormula(i, 0));
            buff.newLine();
            buff.write(getMetricName(i));
            buff.write(" EXCL");
            buff.write(" - 10 3 ");
            buff.write(getMetricFormula(i, noEvents));
            buff.newLine();
        }

    } // writeCounters

    /**
     * This routine returns all events that are counted.
     * 
     * @return array with all counted events
     */
    public static CountedEvent[] getCountedEvents() {
        
        CountedEvent[] events = new CountedEvent[noEvents];
        
        for (int i = 0; i < noEvents; i++) {
            
            events[i] = theCountedEvents[i];
        }
        
        return events;
    }
    
    /**
     * This routine returns all defined metrics.
     * 
     * @return an array with all defined metrics
     */
    public static CounterMetric[] getCounterMetrics() {
        
        CounterMetric[] metrics = new CounterMetric[noMetrics];
        
        for (int i = 0; i < noMetrics; i++) {
            
            metrics[i] = theCounterMetrics[i];
        }
        
        return metrics;
    }
    
    /**
     * 
     * This routine writes the infos about events and metrics in a file.
     * 
     * @param buff is the buffered output file
     * 
     * @throws IOException can be thrown by writing into the buffer
     */
    public static void output(BufferedWriter buff) throws IOException {

        buff.write("events");

        for (int i = 0; i < noEvents; i++) {

            buff.write(" " + theCountedEvents[i].getName());

        }

        buff.newLine();

        for (int i = 0; i < noEvents; i++) {

            CountedEvent event = theCountedEvents[i];

            buff.write("event " + event.getName() + " " + event.getDescription());
            buff.newLine();

        }

        buff.write("# metrics");

        for (int i = 0; i < noMetrics; i++) {

            buff.write(" " + theCounterMetrics[i].getName());

        }

        buff.newLine();

        for (int i = 0; i < noMetrics; i++) {

            CounterMetric metric = theCounterMetrics[i];
            buff.write("define " + metric.getName() + " " + metric.getFullDescription());
            buff.newLine();

        }
    }

    /**
     * output of total values in calltree PM file.
     * 
     * @param buff is the buffered write file
     * 
     * @throws IOException in case of IO problems
     */
    public static void outTotals(BufferedWriter buff) throws IOException {

        buff.write("totals");

        for (int i = 0; i < noEvents; i++) {

            buff.write(" " + theCountedEvents[i].getTotalValue());

        }

        buff.newLine();

    }
    
} // class CounterData
