/*
 * CTProperties.java
 *
 * Properties for Calltree Performance Monitoring.
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

import org.apache.log4j.Logger;

import adaptor.General.ArraySelection;
import adaptor.General.CallNode;
import adaptor.General.CountedEvent;
import adaptor.General.CounterMetric;
import adaptor.General.PropertyMenu;
import adaptor.General.RegionDescriptor;
import adaptor.General.Scaling;

/**
 * CTProperites is a class containing all global values
 * needed for drawing of Calltree Performance Data.
 * <ul>
 * <li>
 * Depth (value in %) for showing relevant nodes and edges
 * <li>
 * Stress (stressing the relative value for colors)
 * <li>
 * Mode (can be absolute or relative)
 * <li>
 * Node Properties
 * <li>
 * Edge Properties
 * </ul>
 *
 * @author Thomas Brandes
 * @version $LastChangedRevision$
 */

public final class CTProperties
{

    /**
     * This is the color taken for edges.
     */
    private static final String EDGE_MAIN_COLOR = "black";

    /**
     * Enumeration type for all possible kind of modes.
     *
     * @version $LastChangedRevision$
     * @author Thomas Brandes
     */
    public static enum ModeKind { ABSOLUTE, RELATIVE, METRIC };

    /**
     * This is the title of the costs menu to select
     * how to display the costs of nodes and edges.
     */
    static final String COSTS_MENU = "Costs";

    /**
     * This is the title of the node menu to select properties
     * of the nodes (labels, color).
     */
    static final String NODE_MENU = "Node";

    /**
     * This is the title for the shape menu to select the shape
     * of nodes for the different kind of regions.
     */
    static final String SHAPE_MENU = "Shape";

    /**
     * This is the title for the edge menu to select
     * properties for the edges.
     */
    static final String EDGE_MENU = "Edge";

    /**
     * This is scale factor for percentage.
     */
    private static final double PERCENT = 100.0;

    /**
     * This array sorts the different menus in a order how
     * they should appear later.
     */
    private static final String[] MENU_HEADERS = { COSTS_MENU, NODE_MENU, SHAPE_MENU, EDGE_MENU };

    /**
     * String for INLCUSIVE costs (avoids misspelling).
     */
    private static final String INCLUSIVE = "inclusive";

    /**
     * String for EXCLUSIVE costs (avoids misspelling).
     */
    private static final String EXCLUSIVE = "exclusive";

    /**
     * String for SOURCE_CALLS (avoids misspelling).
     */
    private static final String SOURCE_CALLS = "Source Calls";

    /**
     * String for RUNTIME_CALLS (avoids misspelling).
     */
    private static final String RUNTIME_CALLS = "Runtime Calls";

    /**
     * String for none selection (avoids misspelling).
     */
    private static final String NONE_ITEM = "<none>";

    /**
     * String for none selection (avoids misspelling).
     */
    private static final String KEEP_ITEM = "<keep>";

    /**
     * String for none selection (avoids misspelling).
     */
    private static final String RANDOM = "random";

    /**
     * String for costs selection (avoids misspelling).
     */
    private static final String EVENT_COSTS = "Event Costs";

    /**
     * String for metric selection (avoids misspelling).
     */
    private static final String METRIC_COSTS = "Metric Costs";

    /**
     * String for cost selection (avoids misspelling).
     */
    private static final String COSTS = "costs";

    /**
     * String for carriage return in output file.
     */
    private static final String CR = "\\n";

    /**
     * This value is taken as default for node brightness.
     */
    private static final String BRIGHT_FULL = "1.0";

    /**
     * This is the value for high saturation in color.
     */
    private static final String SAT_HIGH = "1.0";

    /**
     * This is the value for middle saturation in color.
     */
    private static final String SAT_MIDDLE = "0.75";

    /**
     * This value is taken for low saturation in color.
     */
    private static final String SAT_LOW = "0.4";

    /**
     * Logger for this class.
     */
    private static Logger logger = Logger.getLogger( CTProperties.class );

    /**
    * These are the items for the Show Menu.
    */
    private static final String[] SHOW_ITEMS = { EVENT_COSTS, METRIC_COSTS, RUNTIME_CALLS, SOURCE_CALLS};

    /**
     * This is the property menu for all counted events. The selected event
     * is used to show the costs for nodes and edges.
     */
    private static PropertyMenu propEvent = null;

    /**
     * This is the property meun for all available counter metrics. The selected metric
     * can be used to label the nodes and edges with corresponding values.
     */
    private static PropertyMenu propMetric = null;

    /**
     * These are all possible values that can be selected for depth of presentation.
     */
    private static final double[] DEPTH_ITEMS = { 10.0, 5.0, 2.0, 1.0, 0.5, 0.2, 0.1, .01, 0.0 };

    /**
     * Tooltip for the depth property.
     */
    private static final String DEPTH_HELP = "Show only nodes with a certain percentage of costs";

    /**
     * This is the value selected for presentation depth.
     */
    private static PropertyMenu propDepth = new PropertyMenu( COSTS_MENU, "Depth", DEPTH_HELP, DEPTH_ITEMS, 1.0 );

    /**
     * These are all possible values that can be selected for the stress value.
     * As depth and stress are correlated, this array should have the same number of
     * entries that are greater or equal to zero.
     */
    private static final int[] STRESS_ITEMS = { -1, 0, 1, 2, 3, 4, 5, 6, 8, 10 };

    /**
     * Tooltip for the stress property.
     */
    private static final String STRESS_HELP = "Change the range of color values for costs";
    /**
     * The stress value changes the color scaling. The higher the value the more
     * smaller values become highlighted. The stress value will be set implicitly
     * with the depth. The default value for depth 1.0 is 4.
     */
    private static PropertyMenu propStress = new PropertyMenu( COSTS_MENU, "Stress", STRESS_HELP, STRESS_ITEMS, 4 );

    // Depth and Stress are correlated, Stress should have same number of entries

    // Mode for absolute or relative Counter values

    /**
     * Tooltip for the mode property.
     */
    private static final String MODE_HELP = "Choose between relative event costs, absolute event costs or metric costs";

    /**
     * Property menu for mode of cost representation.
     */
    private static PropertyMenu propMode = new PropertyMenu( COSTS_MENU, "Mode", MODE_HELP, ModeKind.values(), ModeKind.RELATIVE );

    /**
     * These are the possibilities of how a label for a node can be chosen.
     * <ul>
     * <li>
     * inclusive gives as label the inclusive costs for the selected counter or metric
     * <li>
     * exclusive gives as label the exclusive costs for the selected counter or metric
     * <li>
     * source calls will print as label the number of known source calls of the routine
     * <li>
     * runtime calls will print as label the number of calls for this node at runtime
     * <li>
     * distance will print the distance of the node to the selected node(s)
     * <li>
     * #subnodes will print as label the number of nodes that are grouped by this node
     * <li>
     * keep will keep the last selected label even if the counter is changed; by this
     * way it is possible to have labels for different counters
     * <li>
     * none will not print any label
     * </ul>
     */
    private static final String[] NODE_LABEL_ITEMS = { INCLUSIVE, EXCLUSIVE, SOURCE_CALLS, RUNTIME_CALLS,
                                  "distance", "#subnodes", KEEP_ITEM, NONE_ITEM
                                                     };

    /**
     * Tooltip for the node label1 property.
     */
    private static final String NODE_LABEL1_TIP = "Selection of info to be displayed as first label.";

    /**
     * Tooltip for the node label2 property.
     */
    private static final String NODE_LABEL2_HELP = "Selection of info to be displayed as second label.";

    /**
     * Tooltip for the node hue property.
     */
    private static final String NODE_HUE_TIP = "Choose how the hue of the node color will be determined.";

    /**
     * Tooltip for the node saturation property.
     */
    private static final String NODE_SAT_HELP = "Choose how the saturation of the node color will be determined.";

    /**
     * Tooltip for the node brightness property.
     */
    private static final String NODE_BRIGHTNESS_TIP = "Select the brightness of the color for the nodes.";

    /**
     * Property menu for node label 1.
     */
    private static PropertyMenu propNodeLabel1 = new PropertyMenu( NODE_MENU, "Label1", NODE_LABEL1_TIP, NODE_LABEL_ITEMS, INCLUSIVE );

    /**
     * Property menu for node label 2.
     */
    private static PropertyMenu propNodeLabel2 = new PropertyMenu( NODE_MENU, "Label2", NODE_LABEL2_HELP, NODE_LABEL_ITEMS, EXCLUSIVE );

    /**
     * This variable will contain the legend for node label 1. We need a private variable
     * as the first label might be kept.
     */
    private static String theNodeLegend1 = "";

    /**
     * This variable will contain the legend for node label 2. We need a private variable
     * as the second label might be kept.
     */
    private static String theNodeLegend2 = "";



    /**
     * These are the menu items for selection of hue for the nodes.
     */
    private static final String[] NODE_HUE_ITEMS = { INCLUSIVE, RANDOM, KEEP_ITEM };

    /**
     * This is the selected kind of how to compute the hue of the node.
     */
    private static PropertyMenu propNodeHue = new PropertyMenu( NODE_MENU, "Hue", NODE_HUE_TIP, NODE_HUE_ITEMS, INCLUSIVE );

    /**
     * These are all possible value of how the color saturation for a node is determined.
     */
    private static final String[] NODE_SAT_ITEMS = { EXCLUSIVE, SAT_HIGH, SAT_MIDDLE, SAT_LOW, KEEP_ITEM };

    /**
     * This is the property menu for node saturation.
     */

    private static PropertyMenu propNodeSaturation = new PropertyMenu( NODE_MENU, "Saturation", NODE_SAT_HELP, NODE_SAT_ITEMS, EXCLUSIVE );

    /**
     * These are the values offered in the brightness menu for nodes.
     */
    private static final double[] NODE_BRIGHTNESS_ITEMS = { 1.0, 0.9, 0.8, 0.7, 0.6 };

    /**
     * This ist the selected value for node brightness.
     */

    private static PropertyMenu propNodeBrightness = new PropertyMenu( NODE_MENU, "Brightness",
            NODE_BRIGHTNESS_TIP, NODE_BRIGHTNESS_ITEMS, 1.0 );


    /**
     * This value is the multiplicator before rouding to integer. It
     * determines the relevant digits for the fixed size.
     */
    private static double fixVal = 100.0;

    /**
     * This array contains possible values for number of relevant digits
     * to fix the number of digits.
     */
    private static final int[] FIX_VALUE_ITEMS = { 0, 1, 2, 3 };

    /**
     * Tool tip for the property menu to select the fix value.
     */
    private static final String NODE_FIX_TIP = "number of relevant digits for relative value";

    /**
     * Property menu to fix the node values.
     */
    private static PropertyMenu propNodeFixValue = new PropertyMenu( COSTS_MENU, "Fix", NODE_FIX_TIP, FIX_VALUE_ITEMS, 2 );

    private static PropertyMenu propScaling = new PropertyMenu( COSTS_MENU, "Scaling", "NoTip", Scaling.getScaleOrders(), Scaling.AUTO_SCALING );

    /**
     * These shapes are supported for nodes in the GRAPPA representation.
     */
    private static final String[] NODE_SHAPE_ITEMS = { "record", "Mrecord", "box", "ellipse", "hexagon", "octagon",
                                  "circle", "plaintext", "egg",
                                  "triangle", "diamond", "trapezium", "parallelogram", "doublecircle"
                                                     };

    /**
     * Each node or better region kind might have its own shape.
     */
    private static PropertyMenu [] propNodeShapes = null;

    /**
     * Tooltip for the edge label menu.
     */
    private static final String EDGE_LABEL_TIP = "Selection of info to be displayed at the edges.";

    /**
     * These are the menu items of how to compute the labels of edges.
     */
    private static final String[] EDGE_LABEL_ITEMS = { COSTS, SOURCE_CALLS, RUNTIME_CALLS };

    /**
     * This is the corresponding property menu for edge label.
     */
    private static PropertyMenu propEdgeLabel = new PropertyMenu( EDGE_MENU, "Label", EDGE_LABEL_TIP, EDGE_LABEL_ITEMS, COSTS );

    /**
     * This label is used as the legend for edges.
     */
    private static String edgeLegend = "";

    /**
     * These are the items offered in the Width menu for edges.
     */
    private static final int[] EDGE_WIDTH_ITEMS = { 1, 2, 3, 4, 5 };

    /**
     * This is the property menu for the edge width.
     */
    private static PropertyMenu propEdgeWidth = new PropertyMenu( EDGE_MENU, "Width", "", EDGE_WIDTH_ITEMS, 1 );

    /**
     * Array with possible items for coloring edges.
     */
    private static final String[] EDGE_COLOR_ITEMS = { COSTS, EDGE_MAIN_COLOR };

    /**
     * This is the slected item for the edge color.
     */
    private static PropertyMenu propEdgeColor = new PropertyMenu( EDGE_MENU, "Color", "", EDGE_COLOR_ITEMS, EDGE_MAIN_COLOR );

    /**
     * Private constructor hides the public constructor.
     *
     */
    private CTProperties()
    {

    }

    /**
     * This routine creates all property menus that are not available yet.
     * It is mandatory to call this routine to make sure that also such menus
     * will be created that depend on the input data (e.g. counters, metrics).
     *
     */
    public static void createPropertyMenus()
    {

        if ( propNodeShapes == null )
        {

            final String defaultShape = NODE_SHAPE_ITEMS[0];

            String[] kindItems = CallNode.getNodeKindStrings();

            int n = kindItems.length;

            propNodeShapes = new PropertyMenu[n + 1];

            String tip = "Selection of the default shape for nodes.";

            propNodeShapes[0] = new PropertyMenu( SHAPE_MENU, "default", tip, NODE_SHAPE_ITEMS, defaultShape );

            for ( int i = 0; i < n; i++ )
            {

                tip = "Selection of the shape for nodes that stand for a " + kindItems[i] + " region.";

                // all different region kinds get no selection so it will be the default

                propNodeShapes[i + 1] = new PropertyMenu( SHAPE_MENU, kindItems[i], tip, NODE_SHAPE_ITEMS );

            }

        }

        if ( propEvent == null )
        {

            final String help = "Selection of a counted event for the costs.";

            propEvent = new PropertyMenu( COSTS_MENU, "Event", help, CounterData.getCountedEvents() );

        }

        if ( propMetric == null )
        {

            final String help = "Selection of a counter metric to label nodes.";

            propMetric = new PropertyMenu( COSTS_MENU, "Metric", help, CounterData.getCounterMetrics() );
        }

    }
    /**
     * Sets default properties to display source call information.
     *
     */
    public static void setSourceCallDefaults()
    {

        propNodeBrightness.setItem( BRIGHT_FULL );

        propNodeLabel1.setItem( SOURCE_CALLS );
        propNodeLabel2.setItem( NONE_ITEM );

        propEdgeLabel.setItem( SOURCE_CALLS );

        propNodeHue.setItem( RANDOM );

        propNodeSaturation.setItem( SAT_MIDDLE );

    } // setSourceCallDefaults

    /**
     * Set default properties to display runtime call information.
     */

    public static void setRuntimeCallDefaults()
    {

        propNodeBrightness.setItem( BRIGHT_FULL );

        propNodeLabel1.setItem( RUNTIME_CALLS );
        propNodeLabel2.setItem( NONE_ITEM );

        propEdgeLabel.setItem( RUNTIME_CALLS );

        propNodeHue.setItem( RANDOM );

        propNodeSaturation.setItem( SAT_MIDDLE );


    } // setRuntimeCallDefaults

    /**
     * Set default properties to display costs for a counter.
     */

    public static void setEventCostDefaults()
    {

        // do not set defaults wihout a counter

        if ( propEvent.getSelectedObject() == null )
        {

            logger.error( "no counter selected, no costs" );
            return;
        }

        propMode.setItem( ModeKind.RELATIVE );

        propNodeBrightness.setItem( BRIGHT_FULL );

        propNodeLabel1.setItem( INCLUSIVE );
        propNodeLabel2.setItem( EXCLUSIVE );

        propEdgeLabel.setItem( COSTS );

        propNodeHue.setItem( INCLUSIVE );
        propNodeFixValue.setItem( 2 );

        propNodeSaturation.setItem( EXCLUSIVE );

    }

    /**
     * Set default properties to display costs for a counter.
     */

    public static void setMetricCostDefaults()
    {

        // do not set defaults wihout a counter

        if ( propMetric.getSelectedObject() == null )
        {

            logger.error( "no metric selected, no costs" );
            return;
        }

        propMode.setItem( ModeKind.METRIC );

        propNodeBrightness.setItem( BRIGHT_FULL );

        propNodeLabel1.setItem( INCLUSIVE );
        propNodeLabel2.setItem( EXCLUSIVE );

        propEdgeLabel.setItem( COSTS );

        propNodeHue.setItem( INCLUSIVE );
        propNodeFixValue.setItem( 2 );

        propNodeSaturation.setItem( EXCLUSIVE );

    }

    /**
     * Set default values for all properties that are used for CalltreePM.
     *
     */
    public static void setDefaults()
    {

        RegionDescriptor.RegionKind[] kindItems = RegionDescriptor.RegionKind.values();

        createPropertyMenus();

        logger.info( "set defaults for CT properties" );

        propNodeBrightness.setItem( BRIGHT_FULL );

        propNodeShapes[0].setItem( NODE_SHAPE_ITEMS[0] );

        for ( int i = 0; i < kindItems.length; i++ )
        {

            propNodeShapes[i + 1].setNone();
        }

        // take the latest event as the selected event

        propEvent.setLastItem();

        propMetric.setLastItem();

        if ( propEvent.getSelectedObject() != null )
        {

            setEventCostDefaults();

        }
        else
        {

            setRuntimeCallDefaults();

        }

    } // setDefaults

    /**
     * Default shape of a node.
     *
     * @return String for default node shape
     */

    public static String getNodeDefaultShape()
    {

        return propNodeShapes[0].getSelection();

    } // getNodeDefaultShape


    /**
     * The routine getNodeLegend gives the node legend for a given
     * node label item.
     *
     * @param oldLegend is the old legend string needed in case of 'keep' has been selected
     * @param labelItemSt is the kind of node label chosen
     * @param costName is the name of counter/metric used for inclusive/exclusive costs
     * @return the legend for the node as printable String
     */
    static String getNodeLegend( String oldLegend, String labelItemSt, String costName )
    {

        String newLegend = null;

        // inclusive (0), exclusive (1)

        if ( labelItemSt.equals( INCLUSIVE ) || labelItemSt.equals( EXCLUSIVE ) )
        {

            // inclusive or exclusive has been chosen

            newLegend = costName + ":" + " " + labelItemSt;

        }
        else if ( labelItemSt.equals( KEEP_ITEM ) )
        {

            newLegend = oldLegend;

        }
        else if ( labelItemSt.equals( NONE_ITEM ) )
        {

            newLegend = "";

        }
        else
        {

            // SOURCE_CALLS, RUNTIME_CALLS, DISTANCE, #SUBNODES

            newLegend = labelItemSt;
        }

        return newLegend;

    } // getNodeLegend

    /**
     * This return returns a legend for the current selected counter/metric.
     *
     * @return String with desciption
     */
    static String getCounterLegend()
    {

        if ( propMetric.getSelectedObject() != null )
        {

            // Counter Metric has been selected

            return propMetric.getSelection();

        }
        else if ( propEvent.getSelectedObject() != null )
        {

            return propEvent.getSelection();

        } // Counted Event has been selected

        return "";

    } // getCounterLegend

    /**
     * Get the Title for the Call Graph (based on the properties).
     *
     * @return String with the title
     */
    static String getTitle()
    {

        String title = "";
        String event = "";



        Object mode = propMode.getSelectedObject();

        if ( mode == ModeKind.METRIC )
        {

            CounterMetric selectedMetric = ( CounterMetric ) propMetric.getSelectedObject();

            if ( selectedMetric != null )
            {

                // Counter Metric has been selected

                event = selectedMetric.getName();
                title = event + ":" + " " + selectedMetric.getDescription();
                title += CR;

            }

        }
        else
        {

            CountedEvent  selectedEvent  = ( CountedEvent )  propEvent.getSelectedObject();

            if ( selectedEvent != null )
            {

                // Counter has been selected

                event = selectedEvent.getName();
                title = event + ": " + selectedEvent.getDescription();
                title += " (Depth=" + propDepth.getSelection() + ")";
                title += CR;

            } // Counted Event has been selected

        }

        title += "Node [label = ";

        theNodeLegend1 = getNodeLegend( theNodeLegend1, propNodeLabel1.getSelection(), event );

        logger.info( "NodeLegend1 = " + theNodeLegend1 + ", Item = " + propNodeLabel1.getSelection() );

        theNodeLegend2 = getNodeLegend( theNodeLegend2, propNodeLabel2.getSelection(), event );

        logger.info( "NodeLegend2 = " + theNodeLegend2 + ", Item = " + propNodeLabel2.getSelection() );

        if ( ( theNodeLegend1.length() == 0 ) || ( theNodeLegend2.length() == 0 ) )
        {

            title += theNodeLegend1 + theNodeLegend2;

        }
        else
        {

            title += theNodeLegend1 + " | " + theNodeLegend2;
        }

        if ( propEdgeLabel.isSelected( COSTS ) )
        {

            edgeLegend = event + ": " + COSTS;

        }
        else
        {

            edgeLegend = propEdgeLabel.getSelection();
        }

        title += "]" + CR + "Edge [label = " + edgeLegend + "]";

        return title;

    } // getTitle

    /**
     * This routine acts on the menu selection that has been chosen.
     *
     * @param cmd is the selected command
     * @param title is the menu from which command comes
     * @param header is the header of the menu from which the command comes
     * @param main is pointer back to the main window
     */
    static void setProperty( String cmd, String title, String header )
    {

        boolean redraw = false;

        logger.info( "setProperty, cmd = " + cmd + ", menu = " + title + ", header = " + header );

        PropertyMenu menu = PropertyMenu.findMenu( header, title );

        /* no longer available
         *

        if (menu == propShow) {

            // this menu is not really a property but a multiple setting

            setShowDefaults(cmd);

            redraw = true;           // redraw in any case

        } else

         */

        if ( menu != null )
        {

            // okay, we have found the menu and so we change it

            redraw = menu.changeItem( cmd );

            if ( menu == propNodeFixValue )
            {

                int digits = Integer.parseInt( propNodeFixValue.getSelection() );

                logger.info( "Labels will now have " + digits + " relevant digits" );

                Scaling.setPrecision( digits );

            }

            if ( menu == propScaling )
            {

                String globalScaling = propScaling.getSelection();

                logger.info( "Scaling is now: " + propScaling.getSelection() );

                Scaling.setGlobalScaling( globalScaling );
            }

            // we have to handle the depth menu as a special case

            if ( redraw && title.equals( propDepth.getTitle() ) )
            {

                // change of depth implies change of stress

                int index = ArraySelection.findIndex( propDepth.getItems(), cmd );

                propStress.setItem( STRESS_ITEMS[index + 1] + "" );
            }

        }

    } // setProperty

    public static String[] getShowItems()
    {

        return SHOW_ITEMS;
    }

    /**
     * This routine is used to set the default properties if an item of
     * the show menu has been selected.
     *
     * @param cmd is an element of getShowItems()
     */
    public static void setShowDefaults( String cmd )
    {

        if ( cmd.equals( METRIC_COSTS ) )
        {

            setMetricCostDefaults();

        }
        else if ( cmd.equals( EVENT_COSTS ) )
        {

            setEventCostDefaults();

        }
        else if ( cmd.equals( SOURCE_CALLS ) )
        {

            setSourceCallDefaults();

        }
        else
        {

            setRuntimeCallDefaults();
        }
    }

    /**
     * This routine stresses a value between 0.0 and 1.0. For
     * stress == 1, values remain unchanged. For stress > 1, values
     * become larger and closer to 1. E.g. 0.5 is 0.75 for stress =2,
     * 0.9375 for stress =3 and so on. This routine is used to highlight or stress
     * also smaller values.
     *
     * @param x is the input value in range from 0.0 to 1.0
     * @return output value even in the same range but closer to 1.0
     */
    private static double rateByStress( double x )
    {

        // at first we build the inverse

        double inv = 1.0 - x;

        // square the inverse n-times, where n is the stress

        int stress = Integer.parseInt( propStress.getSelection() );

        for ( int i = 0; i < stress; i++ )
        {

            inv = inv * inv;
        }

        // reinverse if stress was positive

        if ( stress >= 0 )
        {

            inv = 1.0 - inv;

        }

        return inv;

    }

    /**
     * This is the inverse routine to rateByStress. It is needed
     * to make a useful color legend where we need the threshold values
     * for the colors.
     *
     * @param x is the input value
     * @return y with rateByStress(y) == x
     */
    public static double invByStress( double x )
    {

        // at first we build the inverse

        double inv = 1.0 - x;

        // square the inverse n-times, where n is the stress

        int stress = Integer.parseInt( propStress.getSelection() );

        for ( int i = 0; i < stress; i++ )
        {

            inv = Math.sqrt( inv );
        }

        inv = 1.0 - inv;

        return inv;

    }

    /**
     * This function return returns the relative value for the current selected counter
     * or metric. It should return a value between 0.0 and 1.0.
     *
     * @param counterVals are the counted cost values
     * @return relative value compared to the total costs
     */
    private static double getRelativeValue( long[] counterVals )
    {

        CounterMetric selectedMetric = ( CounterMetric ) propMetric.getSelectedObject();
        CountedEvent  selectedEvent  = ( CountedEvent )  propEvent.getSelectedObject();

        if ( propMode.isSelected( ModeKind.METRIC ) && selectedMetric != null )
        {

            return selectedMetric.getRelativeValue( counterVals );

        }
        else if ( selectedEvent != null )
        {

            return selectedEvent.getRelativeValue( counterVals );

        }
        else
        {

            return PERCENT;
        }
    }

    /**
     * This routine computes the relative value with the stress factor.
     *
     * @param counterVals are the counter values
     * @return the relative value (between 0.0 and 1.0) and stressed
     */
    private static double getStressedValue( long[] counterVals )
    {


        double relValue    = getRelativeValue( counterVals );
        double stressValue = rateByStress( relValue );

        logger.debug( "relative value " + relValue + " (stressed " + stressValue + ")" );

        if ( relValue < 0.0 )
        {

            relValue = 0.0;
            stressValue = 0.0;

        }
        else if ( relValue > 1.0 )
        {

            relValue = 1.0;
            stressValue = 1.0;
        }

        return stressValue;

    } // getStressedValue

    /**
     * This routine returns the saturation color value for a node depending on the
     * selected kind.
     *
     * @param counterVals are the values of the exlusive counter
     * @param oldSaturation is the old saturation value
     * @return new saturation value
     */
    static double getNodeSaturation( long[] counterVals, double oldSaturation )
    {

        // default is to keep the old saturation value

        double newSaturation = oldSaturation;

        if ( propNodeSaturation.isSelected( EXCLUSIVE ) )
        {

            // "exclusive" counter value is taken

            final double minSat = 0.1;
            final double maxSat = 1.0;

            // value by Counter

            double x = getStressedValue( counterVals );

            // scale X linear from SatMin to SatMax

            newSaturation = minSat + x * ( maxSat - minSat );

        }
        else if ( !propNodeSaturation.isSelected( KEEP_ITEM ) )
        {

            // one of the value "1.0", "0.1", "0.7"

            double val = Double.parseDouble( propNodeSaturation.getSelection() );

            newSaturation = val;

        }

        return newSaturation;

    }

    /**
     * This routine computes the hue value with properties of the node.
     *
     * @param counterVals are counter values
     * @param name as name of the node used for random color
     * @param oldHue is old value taken for keeping selection
     * @return a new hue value in range 0.0 to 1.0
     */
    static double getNodeHue( long[] counterVals, String name, double oldHue )
    {

        double newHue = oldHue;

        if ( propNodeHue.isSelected( INCLUSIVE ) )
        {

            // Hue value is determined by costs and stress

            final double minHue = 0.3;
            final double maxHue = 0.0;

            // value by Counter

            double stressVal = CTProperties.getStressedValue( counterVals );

            // scale X from HueMin = 0.3 (green) to HueMax = 0.0 red

            newHue = minHue + stressVal * ( maxHue - minHue );

        }
        else if ( propNodeHue.isSelected( RANDOM ) )
        {

            // Hue value is random


            final int initRandomVal = 6191;
            final int reminder      = 13191;
            final int adder         = 131;

            int randomVal = initRandomVal;

            // int C = Name.hashCode (); not very convenient

            for ( int i = 0; i < name.length(); i++ )
            {

                int posVal = name.charAt( i );

                randomVal = ( posVal * randomVal + adder ) % reminder;

            }

            newHue = ( double ) randomVal / ( double ) reminder;

        }

        return newHue;

    }

    /**
     * Asking for the current value of brightness.
     *
     * @return value to be taken for node brightness color
     */
    public static double getNodeBrightness()
    {

        return Double.parseDouble( propNodeBrightness.getSelection() );

    }

    /**
     * Routine to check if the costs of the selected event are about the selected
     * depth.
     *
     * @param countersVals are the counterValues
     * @param isEdge for true an edge is assumed
     * @return true if costs are about threshold
     */
    static boolean aboutDepth( long[] countersVals, boolean isEdge )
    {

        final double nodeThreshold = PERCENT;
        final double edgeThreshold = 4 * PERCENT; // edges are better visible

        boolean about = true;

        CountedEvent  selectedEvent  = ( CountedEvent )  propEvent.getSelectedObject();

        if ( selectedEvent != null )
        {

            double depth = Double.parseDouble( propDepth.getSelection() );

            double x = selectedEvent.getRelativeValue( countersVals );

            if ( isEdge )
            {

                about = x * edgeThreshold >= depth;

            }
            else
            {

                about = x * nodeThreshold >= depth;
            }

        }

        return about;

    } // aboutDepth

    /**
     * This routine returns for counter values the selected counter or metric value.
     *
     * @param counterVals are the given costs
     * @return the corresponding value for selected counter or metric
     */
    static double getValue( long[] counterVals )
    {

        CounterMetric selectedMetric = ( CounterMetric ) propMetric.getSelectedObject();
        CountedEvent  selectedEvent  = ( CountedEvent )  propEvent.getSelectedObject();

        if ( selectedMetric != null )
        {

            return selectedMetric.getValue( counterVals );

        }
        else if ( selectedEvent != null )
        {

            return ( double ) selectedEvent.getCounterValue( counterVals );

        }
        else
        {

            return 0.0;
        }

    }

    /**
     * Returns for the selected counter or metric the value as a string.
     *
     * @param counterVals are the actual counter values
     * @return String with value for the current counter-metric selection
     */
    static String getLabel( long[] counterVals )
    {

        String label = "";


        CountedEvent  selectedEvent  = ( CountedEvent )  propEvent.getSelectedObject();

        if ( propMode.isSelected( ModeKind.METRIC ) )
        {

            CounterMetric selectedMetric = ( CounterMetric ) propMetric.getSelectedObject();

            if ( selectedMetric != null )
            {

                label = selectedMetric.getValueString( counterVals );

            }
            else
            {

                label = "~~";
            }

        }
        else if ( selectedEvent == null )
        {

            label = "--";

        }
        else if ( propMode.isSelected( ModeKind.RELATIVE ) )
        {

            // percentage value of relative costs

            // label = getPercentage(selectedEvent.getRelativeValue(counterVals));

            label = selectedEvent.getRelValueString( counterVals );

        }
        else
        {

            // absolute costs

            // label = "" + selectedEvent.getCounterValue(counterVals);

            label = selectedEvent.getAbsValueString( counterVals );
        }

        return label;

    } // getLabel

    /**
     * Make a label for the node.
     *
     * @param second should be set true for the second label
     * @param oldLabel is the old label (might be kept)
     * @param sourceCalls is the number of source calls
     * @param runCalls is the number of runtime calls
     * @param distance is the distance
     * @param noSubnodes is the number of subnodes
     * @param inclusiveVals are the inclusive costs
     * @param exclusiveVals are the exclusive costs
     * @return string used as label for the node in dot file
     */
    static String getNodeLabel( boolean second, String oldLabel, long sourceCalls, long runCalls, int distance, int noSubnodes,
                                long[] inclusiveVals, long[] exclusiveVals )
    {

        PropertyMenu propMenu = propNodeLabel1;

        if ( second )
        {

            propMenu = propNodeLabel2;
        }

        String label = oldLabel;

        if ( propMenu.isSelected( INCLUSIVE ) )
        {

            label = getLabel( inclusiveVals );

        }
        else if ( propMenu.isSelected( EXCLUSIVE ) )
        {

            label = getLabel( exclusiveVals );

        }
        else if ( propMenu.isSelected( SOURCE_CALLS ) )
        {

            label = "" + sourceCalls;

        }
        else if ( propMenu.isSelected( RUNTIME_CALLS ) )
        {

            label = "" + runCalls;

        }
        else if ( propMenu.isSelected( "distance" ) )
        {

            label = "" + distance;

        }
        else if ( propMenu.isSelected( "#subnodes" ) )
        {

            label = "" + noSubnodes;

        }
        else if ( propMenu.isSelected( NONE_ITEM ) )
        {

            label = "";

        }
        else if ( propMenu.isSelected( KEEP_ITEM ) )
        {

            label = oldLabel;

        }
        else
        {

            logger.error( "getNodeLabel: illegal item " + propMenu.getSelection() );
        }

        return label;

    } // getNodeLabel

    /**
     * Getter routine for the currently selected edge width.
     *
     * @return current edge width value
     */
    static int getEdgeWidth()
    {

        return Integer.parseInt( propEdgeWidth.getSelection() );
    }

    /**
     * This routine returns a label for the edge. The kind of label depends
     * on the current selection for edge label kind.
     *
     * @param sourceCalls is the number of source calls
     * @param runCalls is the number of runtime calls
     * @param counterVals are the cost values of the edge
     * @return String that woks as label for the edge
     */
    static String getEdgeLabel( long sourceCalls, long runCalls, long[] counterVals )
    {

        if ( propEdgeLabel.isSelected( COSTS ) )
        {

            return getLabel( counterVals );

        }
        else if ( propEdgeLabel.isSelected( SOURCE_CALLS ) )
        {

            return "" + sourceCalls;
        }

        return "" + runCalls;

    } // getEdgeLabel

    /**
     * Asking for the color of an edge. The counter values might
     * be needed if color depends on the costs.
     *
     * @param counterValues might be needed for the color value
     * @return String containing the color value.
     */
    public static String getEdgeColor( long[] counterValues )
    {

        if ( propEdgeColor.isSelected( EDGE_MAIN_COLOR ) )
        {

            return EDGE_MAIN_COLOR;
        }

        final double minHue = 0.3;
        final double maxHue = 0.0;

        // value by Counter

        double x = getStressedValue( counterValues );

        // scale X from HueMin = 0.3 (green) to HueMax = 0.0 red

        double hue = minHue + x * ( maxHue - minHue );

        return makeColorString( hue, 1.0, 1.0 );

    }

    /**
     * Truncates the double value to a certain number of relevant
     * digits. E.g. 3.14145 returns 3.141
     *
     * @param x is the double value
     * @return x with a fixed number of digits
     */
    private static String fixValue( double x )
    {

        // scaleVal determines the number of positions

        final long scaleVal = 1000;

        long xLong = ( long ) ( x * scaleVal );

        double xFix = ( ( double ) xLong ) / ( ( double ) scaleVal );

        return xFix + "";

    }

    /**
     * Make a printable color string by hue, saturation and brightness.
     *
     * @param hue is the hue value
     * @param sat is the saturation between 0.0 and 1.0
     * @param br is the brightness between 0.0 and 1.0
     * @return String for color that can be used in dot files
     */
    static String makeColorString( double hue, double sat, double br )
    {

        String color = fixValue( hue ) + "," + fixValue( sat ) + "," + fixValue( br );

        return "\"" + color + "\"";
    }

    /**
     * This routine returns the node label to be used for the graph files.
     *
     * @param name is the name of the node
     * @param values is a value string
     * @param kind is the kind of the node
     * @return the label as a string
     */
    static String makeNodeLabel( String name, String values, int kind )
    {

        String label = name;

        if ( values.length() > 0 )
        {

            String shape = getNodeShapeString( kind );

            if ( shape.equals( "record" ) || shape.equals( "Mrecord" ) )
            {

                // this is a record

                label = "{ " + name + " | { " + values + " } } ";

            }
            else
            {

                label = name + CR + values;
            }

        }

        return label;
    }

    /**
     * This routine returns the selected shape for a certain node kind.
     *
     * @param kind is the kind of the node
     * @return an element of the array NODE_SHAPE_ITEMS
     */
    private static String getNodeShapeString( int kind )
    {

        String shape;

        // be careful: propNodeShapes[0] is reserved for default

        if ( propNodeShapes[kind + 1].getSelectedObject() == null )
        {

            // take the default shape

            shape = propNodeShapes[0].getSelection();

        }
        else
        {

            shape = propNodeShapes[kind + 1].getSelection();
        }

        return shape;

    }

    /**
     * This routine returns the selected shape for a certain node kind.
     *
     * @param kind is the kind of the node (starting with 0)
     * @return an element of the array NODE_SHAPE_ITEMS or empty string for the default.
     */
    public static String getNodeShape( int kind )
    {

        String shape = getNodeShapeString( kind );

        if ( propNodeShapes[0].isSelected( shape ) )
        {

            // the selection is the default one, so we can return empty string

            shape = "";
        }

        return shape;

    }

    /**
     * Routine to access all headers that are defined in this class.
     *
     * @return array of headers for which menus are available.
     */
    public static String[] getAllHeaders()
    {

        return MENU_HEADERS;

    }

} // CTProperties
