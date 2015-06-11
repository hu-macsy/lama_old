/*
 * Scaling.java
 * 
 * Help class Scaling to make a good scaling of counter values (m, %, k, M, G).
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
 * Help class Scaling to make a good scaling of counter values (m, %, k, M, G).
 *
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */
public final class Scaling {

    /**
     * The value kiloVal specifies the amount of one kilo. An
     * alternative could be 1024.0 instead of 1000.0.
     * 
     */
    private static final double KILO_VAL = 1000.0;
    
    /**
     * The value milliVal specifies the amount of milli.
     */
    private static final double MILLI_VAL = 0.001;
    
    /**
     * The value percentageValue specifies how much 1 percent is.
     */
    private static final double PERCENTAGE_VAL = 0.01;
    
    /**
     * String used for no scaling.
     */
    private static final String NO_SCALING = "-";
    
    /**
     * String used for scaling with kilo.
     */
    private static final String KILO_SCALING = "K";
    
    /**
     * String used for scaling with Mega.
     */
    private static final String MEGA_SCALING = "M";
    
    /**
     * String used for scaling with Giga.
     */
    private static final String GIGA_SCALING = "G";
    
    /**
     * String used for scaling with milli.
     */
    private static final String MILLI_SCALING = "m";
    
    /**
     * String used for scaling with micro.
     */
    private static final String MICRO_SCALING = "u";
    
    /**
     * String used for scaling with percent, e.g. 5.14% instead of 0.0514.
     */
    private static final String PERCENT_SCALING = "%";
    
    /**
     * We order the scale strings to define scale up and scale down.
     */
    private static final String[] SCALE_ORDER = { MICRO_SCALING, MILLI_SCALING, PERCENT_SCALING, NO_SCALING, 
        KILO_SCALING, MEGA_SCALING, GIGA_SCALING };
    
    /**
     * Array will contain the scale values for the scale orders, e.g.
     * KILO_VAL for "K".
     */
    private static double[] scaleValues = null;
    
    /**
     * The array scaleRate will contain the inverse values of scaleValues and
     * are used for the scaling.
     */
    private static double[] scaleRates = null;
   
    /**
     * Logger for this class.
     */
    private static Logger logger = Logger.getLogger(Scaling.class);
    
    /**
     * If this string contains a value autoScaling will be disabled
     * and this scaling will be preferred always.
     */
    private static String globalScaling = null;

    
    /**
     * This value is used to fix the number of digits.
     */
    private static final double INIT_FIX_VAL = 100.0;

    /**
     * String constant used if auto scaling should be enabled.
     */
    public static final String AUTO_SCALING = "<auto>";    

    /**
     * Override the default constructor for this utility class.
     * 
     */
    private Scaling() {
        
    }    

    /**
     * This value is used to fix the number of digits.
     */
    private static double fixVal = INIT_FIX_VAL;
    
    /**
     * This routine returns an array for selections of global scaling.
     * 
     * @return a string array with the possible selections.
     */
    public static String[] getScaleOrders() {
        
        int len = SCALE_ORDER.length;
        
        String[] orders = new String[len + 1];
    
        orders[0] = AUTO_SCALING;
        
        for (int i = 0; i < len; i++) {
            
            orders[i + 1] = SCALE_ORDER[i]; 
        }
        
        return orders;
    }
    
    /**
     * This routine returns for a total value a scaling that will probably fit best for the
     * presentation.
     * 
     * @param total
     *            is the highest possible value that might appear
     * @param includePercent
     *            is a flag to enable also the "%" scaling
     * @return a String for scaling that fits best to the total value.
     */
    public static String getAutoScaleOrder(double total, boolean includePercent) {

        int len = SCALE_ORDER.length;

        // default value is the highest scaling

        String scaling = SCALE_ORDER[0];

        for (int i = 1; i < len; i++) {

            // if total value is lower than twice of scaling value, take it

            logger.debug("getAutoScaleOrder, total = " + total + " scaleValue = " + scaleValues[i]);

            if (total < 2 * scaleValues[i]) {

                // we take percent only if flag is enabled

                if ((scaling != Scaling.PERCENT_SCALING) || (includePercent)) {

                    break;
                }
            }

            // the scaling was too small, so we try the next one

            scaling = SCALE_ORDER[i];
        }

        return scaling;
    }
   

    /**
     * A help routine to make sure that we have a correct scaling.
     * 
     * @param scaleVal is the proposed scaling.
     * @return a correct scaling string.
     */
    public static String getCorrectScaling(String scaleVal) {
                
        String scaling = NO_SCALING;
        
        int pos = ArraySelection.findIndex(SCALE_ORDER, scaleVal);
        
        if (pos >= 0) {
            
            // this is a legal value
            
            scaling = scaleVal;
        
        } else if (scaleVal.equals("k")) {
            
            // we treat lower case k as upper class K
            
            scaling = KILO_SCALING;
            
        } else {
             
            logger.error("illegal scaling: " + scaleVal);
        }
        
        return scaling;
    }
    
    /**
     * This routine sets a global value for the scaling. If
     * this string is AUTO_SCALING, scaling arguments will 
     * be considered.
     * 
     * @param scaling is the new setting for global scaling.
     */
    public static void setGlobalScaling(String scaling) {
        
        if (scaling == AUTO_SCALING) {
            
            globalScaling = null;
            
        } else {
            
            globalScaling = getCorrectScaling(scaling);
        }
    }

    /**
     * This routine increases the current scaling.
     * 
     * @param lastScaling is the input scaling value.
     * @return a string with the next higher scaling.
     */
    public static String scaleUp(String lastScaling) {
        
        String scaling = lastScaling;
        
        int index = ArraySelection.findIndex(SCALE_ORDER, lastScaling);
        
        if ((index >= 0) && (index + 1 < SCALE_ORDER.length)) {
            
            scaling = SCALE_ORDER[index + 1];
        }

        return scaling;
    }
    
    /**
     * This routine decreases the current scaling.
     * 
     * @param lastScaling is the input scaling value.
     * @return a string with the next lower scaling.
     */
    public static String scaleDown(String lastScaling) {
        
        String scaling = lastScaling;
        
        int index = ArraySelection.findIndex(SCALE_ORDER, lastScaling);
                
        if (index != ArraySelection.NO_SELECTION) {
            
            scaling = SCALE_ORDER[index - 1];
        }

        return scaling;
    }
    
    /**
     * This routine returns for a double value a good string
     * represenation that observes the given scaling and the global
     * value for number of relevant digits.
     * 
     * @param theVal is the value that needs scaling
     * @param scaling is the scaling (will always be taken)
     * @return the scaled value as a string
     */
    private static double makeScaledValue(double theVal, String scaling) {

        double val = theVal;

        if (scaling.length() > 0) {
            
            int pos = ArraySelection.findIndex(SCALE_ORDER, scaling);
            
            if (pos >= 0) {
                
                val *= scaleRates[pos];      
            }            
        }

        return val;
    }

    /**
     * This routine returns for a double value the scaled value. Be
     * careful that the scaling becomes visible as suffix or in a legend.
     * 
     * @param theVal is the value that needs scaling
     * @param autoScaling is the scaling (globalScaling might be taken)
     * @return the scaled value as a double value
     */
    
    public static double getScaledValue(double theVal, String autoScaling) {
        
        if (globalScaling == null) {
            
            return makeScaledValue(theVal, autoScaling);
            
        } else {
            
            return makeScaledValue(theVal, globalScaling);
            
        }
    }

    /**
     * This routine sets the number of relevant digits.
     * 
     * @param noDigits is an int value for the number of relevant digits.
     */
    public static void setPrecision(int noDigits) {
        
        fixVal = 1.0;
        
        for (int i = 0; i < noDigits; i++) {
            
            fixVal = fixVal * 10.0;
        }

    }
    
    /**
     * This routine returns for a double value a good string
     * representation that observes the given scaling and the global
     * value for number of relevant digits.
     * 
     * @param val is the value that needs scaling
     * @param scaling is the scaling (will always be taken)
     * @return the scaled value as a string
     */
    
    private static String makeValueString(double val, String scaling) {
        
        double scaledVal = makeScaledValue(val, scaling);
        
        /* we make a transformation to get fixed value */

        long longVal = Math.round(scaledVal * fixVal);
        
        String valString;
        
        if (fixVal == 1.0) {
            
            valString = Long.toString(longVal);
            
        } else {
            
            scaledVal = longVal / fixVal;
            
            valString = Double.toString(scaledVal);
        }
        
        if (scaling.length() > 0) {
            
            valString += " " + scaling; 
            
        }

        return valString;

    }
    
    /**
     * This routine returns for a double value a good string
     * representation that observes the given scaling and the global
     * value for number of relevant digits.
     * 
     * @param val is the value that needs scaling
     * @param autoScaling is the scaling (will always be taken)
     * @return the scaled value as a string
     */
    
    public static String getValueString(double val, String autoScaling) {
        
        if (globalScaling == null) {
            
            return makeValueString(val, autoScaling); 
            
        } else {
            
            return makeValueString(val, globalScaling);
        }
 
    }
    
    /**
     * This routine returns for a double value a good string
     * representation that shows the percentage value. Global
     * scaling is not considered. Global value for number of
     * relevant digits will be taken.
     * 
     * @param val is the value that needs scaling
     * @return the scaled value in percentage as a string
     */
    
    public static String getPercentString(double val) {
        
        return makeValueString(val, PERCENT_SCALING);
    
    }
    
    static {
        
        int len = SCALE_ORDER.length;
        
        scaleRates  = new double[len];
        scaleValues = new double[len];
        
        for (int i = 0; i < len; i++) {
            
            char rate = SCALE_ORDER[i].charAt(0);
            
            double val  = 1.0;
            double rval = 1.0;
            
            switch (rate) {
            
            case 'G':
                
                val  *= KILO_VAL;
                rval /= KILO_VAL;
                
            case 'M':
                
                val  *= KILO_VAL;
                rval /= KILO_VAL;

            case 'k':
            case 'K':
                
                val  *= KILO_VAL;
                rval /= KILO_VAL;
                
                break;
                                
            case 'u':
                
                val  *= MILLI_VAL;
                rval /= MILLI_VAL;
                
            case 'm':
                
                val  *= MILLI_VAL;
                rval /= MILLI_VAL;
                
                break;
                
            case '%':
                
                val  *= PERCENTAGE_VAL;
                rval /= PERCENTAGE_VAL;
                break; 
                
            default:
                
                val  = 1.0;
                rval = 1.0;
                
            } // switch
        
            scaleValues[i] = val;
            scaleRates[i]  = rval;
        }
    }
        
} // class Scaling

