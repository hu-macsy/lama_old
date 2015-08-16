/*
 * CounterMetric.java
 * 
 * Class defines a metric for performance event counters.
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
 * The class CounterMetric defines a performance metric for one or two
 * performance events.
 *
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */
public class CounterMetric {

    /**
     * The following string identifies undefined counters/metrics.
     */
    private static final String UNDEFINED = "<undefined>";

    /**
     * Logger for this class.
     */
    private static Logger logger = Logger.getLogger(CounterMetric.class);

    /**
     * name of the Metric.
     */
    private String myName = "noName";

    /**
     * mode specifies the binary operation for the two performance events.
     */
    private char myMode = '+';

    /**
     * This is first counted performance event.
     */
    private CountedEvent event1 = null;
    
    /**
     * The second event is not mandatory (can be null).
     */
    private CountedEvent event2 = null;

    /**
     * Rating for the first performance event.
     */
    private double rate1;
    
    /**
     * Rating for the seconnd performance event.
     */
    private double rate2;

    /**
     * Scaling is an additional factor for the rating of the value.
     * It can be "M", "K", "m" or "%".
     */
    private String scaling = null;
    
    /**
     * 
     */
    private double total; // reference value 

    /**
     * Constructor of a metric for one performance event.
     * 
     * @param theName becomes the name of the metric
     * @param rate is the rating of the first event
     * @param event is the first performance event
     */ 
    public CounterMetric(String theName, double rate, CountedEvent event) {
    
        myName = theName;
    
        myMode = ' ';
        event1 = event;
        event2 = null;
    
        rate1 = rate;
        rate2 = 0.0;
    
    }

    /**
     * Constructor of a metric for two performance events.
     * 
     * @param name theName becomes the name of the metric
     * @param mode specifies the binary operatior
     * @param r1 rating of first event
     * @param cev1 first performance event
     * @param r2 rating of second event
     * @param cev2 second performance event
     */
    public CounterMetric(String name, char mode, double r1, CountedEvent cev1, double r2, CountedEvent cev2) {
    
        myName = name;
    
        this.myMode = mode;
    
        event1 = cev1;
        event2 = cev2;
    
        rate1 = r1;
        rate2 = r2;
    
    } // constructor CounterMetric

    /**
     * compute a metric value by counter values.
     * 
     * @param longVal1 is the long value of the first performance event.
     * @param longVal2 is the long value of the second performance event.
     * @return the rated value
     */
    private double computeValue(long longVal1, long longVal2) {

        double val = 0.0;
        
        if (event1 != null) {
            
            val = (double) longVal1;
            val *= rate1;
            
            if (event2 != null) {
                
                double val2 = (double) longVal2;
                val2 *= rate2;
                
                // logger.debug("computeValue, mode = " + mode + ", " + val + " " + val2);
                
                switch (myMode) {
                
                case '/' : 
                    val = val / val2;
                    break;
                
                case '*' : 
                    val = val * val2;
                    break;
                
                case '+' : 
                    val = val + val2;
                    break;
                
                case '-' : 
                    val = val - val2;
                    break;
                    
                case '#' :
                    val = val / (val + val2);
                    break;
                
                case '~' :
                    val = (val - val2) / val;
                    break;
                
                default:
                    
                    logger.error("illegal op (" + myMode + ")");
                }
            }
        }
        
        return val;
    }

    /**
     * setRating determines a rating factor from the total values.
     * 
     */
    private void setRating() {
        
        long val1 = 0;
        
        long val2 = 0;

        if (event1 != null) {
            
            val1 = event1.getTotalValue();            
        }

        if (event2 != null) {
            
            val2 = event2.getTotalValue();      
        }

        total = computeValue(val1, val2);

        boolean includePercent = (myMode == '/');
        
        scaling = Scaling.getAutoScaleOrder(total, includePercent);

        logger.info("setRating for this metric: " + getFormula(0));
        logger.info("setRating, L1 = " + val1 + ", L2 = " + val2 + ", total (ref) = " + total);
        logger.info("setRating done, is " + scaling);

    } // setRating

    /**
     * Set the scaling for this metric.
     * 
     * @param scaleVal is a string in array SCALE_ORDER
     */
    public void setScaling(String scaleVal) {
        
        scaling = Scaling.getCorrectScaling(scaleVal);

    }
 
    /**
     * Getter routine for the scaling string.
     * 
     * @return the current value of scaling for this metric
     */
    public String getScaling() {
        
        return scaling;
    }    
     
    /**
     * This routine increases the current scaling.
     * 
     */
    public void scaleUp() {
        
        scaling = Scaling.scaleUp(scaling);
    }
    
    /**
     * This routine descreases the current scaling to the next
     * lower level.
     * 
     */
    public void scaleDown() {
        
        scaling = Scaling.scaleDown(scaling);

    }
    /**
     * Getter routine for the name of the metric.
     * 
     * @return the name of this metric.
     */
    public String getName() {
        
        return myName;
    }
    
    /**
     * get the metric values by given counter values.
     * 
     * @param countedValues are all counted event values
     * @return the metric value given by this metric.
     */
    public double getValue(long[] countedValues) {

        long val1 = 0;
        
        long val2 = 0;

        if (event1 != null) {
            
            val1 = event1.getCounterValue(countedValues);           
        }

        if (event2 != null) {
            
            val2 = event2.getCounterValue(countedValues);
           
        }

        return computeValue(val1, val2);
    }

    /**
     * Get the metric value by given counter values in the scaled version.
     * E.g. it returns 16.3 instead of 0.163 if scaling is "%".
     * 
     * @param countedValues are all counted event values
     * @return the metric value given by this metric.
     */
    public double getScaledValue(long[] countedValues) {
        
        double val = getValue(countedValues);
        
        val = Scaling.getScaledValue(val, scaling);

        return val;
    }

    /**
     * get the relative value.
     * 
     * @param countedValues are all counted event values 
     * @return the relative value for the metric.
     */
    public double getRelativeValue(long[] countedValues) {

        if (scaling == null) {
            
            setRating();            
        }

        double value = getValue(countedValues);

        if (myMode == '/') { 
            
            // 2 * total is worse: 0.0 
            // 0.0 is best: 1.0

            return value / (2.0 * total);
        }

        return value / total;
    }

    /**
     * This function returns a string that represents the value of
     * the metric for a given set of event counter values.
     * 
     * @param countedValues is the array with values of all counted events
     * @return string with rounded and scaled value
     */
    public String getValueString(long[] countedValues) {
        
        // make sure that we have already a legend for the rating
        
        if (scaling == null) {
            
            setRating();            
        }

        double val = getValue(countedValues);

        return Scaling.getValueString(val, scaling);

    }

    /**
     * get a description of the metric used for help routines.
     * 
     * @return String containing the description of the metric.
     */
    public String getDescription() {

        String description = null;
        
        if (event1 == null) {
            
            description = UNDEFINED;
            
        } else if (event2 == null) {

            description = "~ " + event1.getName();
            
        } else {
            
            description = "~ " + event1.getName() + " " + myMode + " " + event2.getName();
        }

        return description;

    } // getDescription

    /**
     * get a full description of the metric used for output.
     * 
     * @return String containing the full description of the metric.
     */
    public String getFullDescription() {

        String description = null;
        
        if (event1 == null) {
            
            description = UNDEFINED;
            
        } else if (event2 == null) {

            description = rate1 + " " + event1.getName();
            
        } else {
            
            description =  myMode + " " + rate1 + " " + event1.getName();
            description += " " + rate2 + " " + event2.getName();
        }

        return description;

    } // getFullDescription

    
    /**
     * returns a formula for the metric (needed for PMS files)!
     * 
     * @param offset is an offset used for counter index
     * 
     * @return a formula expression as a string
     */
    public String getFormula(int offset) {

        String formule = null;
        
        if (event1 == null) {
           
            formule = UNDEFINED;
            
        } else {
            
            formule = rate1 + " " + (event1.getIndex() + offset);
            
            if (event2 != null) {
                formule = myMode + " " + formule;
                formule += " " + rate2 + " " + (event2.getIndex() + offset);
            }
        }
        
        return formule;
    }

    /**
     * {@inheritDoc}
     *
     * @see java.lang.Object#toString()
     */
    public String toString() {
        
        return myName;
    }

} // class CounterMetric

