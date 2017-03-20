/*
 * ChartProperties.java
 *
 * Properties of charts used for ShowPM and routines for the generation
 * of charts using these properties.
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

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;

import javax.swing.JFrame;

import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.CategoryAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.labels.PieToolTipGenerator;
import org.jfree.chart.labels.StandardCategoryToolTipGenerator;
import org.jfree.chart.labels.StandardPieSectionLabelGenerator;
import org.jfree.chart.labels.StandardPieToolTipGenerator;
import org.jfree.chart.plot.CategoryPlot;
import org.jfree.chart.plot.MultiplePiePlot;
import org.jfree.chart.plot.PiePlot;
import org.jfree.chart.plot.Plot;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.renderer.category.BarRenderer;
import org.jfree.chart.renderer.category.CategoryItemRenderer;
import org.jfree.chart.renderer.category.StackedBarRenderer;
import org.jfree.chart.renderer.category.StatisticalBarRenderer;
import org.jfree.data.category.CategoryDataset;
import org.jfree.ui.RefineryUtilities;
import org.jfree.util.TableOrder;

/**
 * ChartProperties is a class that manages global properties for the generation of
 * charts as used for ShowPM.
 *
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */
public final class ChartProperties
{

    /**
     * Menu items for the orientation of the charts.
     */
    public static final String[] ITEMS_ORIENTATION = {"HORIZONTAL", "VERTICAL" };

    /**
     * Menu items for labels on or off.
     */
    public static final String[] ITEMS_LABEL = {"Labels off", "Labels on" };

    /**
     * Preferred size for charts.
     */
    private static final Dimension MY_PREFERRED_SIZE = new Dimension( 640, 420 );

    /**
     * Currently selected orientation for plots.
     */
    private static PlotOrientation orientation = PlotOrientation.HORIZONTAL;

    /**
     * Currently selection for labels on or off in pie charts.
     */
    private static boolean showLabels = true;

    /**
     * This font is used for the title of charts. Currently it cannot be
     * changed.
     */
    private static Font titleFont = JFreeChart.DEFAULT_TITLE_FONT;

    /**
     * Dummy constructor to overwrite default constructor. An object of this
     * class cannot be generated.
     */
    private ChartProperties()
    {
    }

    /**
     * This routine is called to set a property for charts.
     *
     * @param cmd is the item for the property that should be set
     */
    static void setProperty( String cmd )
    {

        if ( cmd.equals( ITEMS_LABEL[1] ) )
        {

            showLabels = true; // Labels on
        }

        if ( cmd.equals( ITEMS_LABEL[0] ) )
        {

            showLabels = false; // Labels off
        }

        if ( cmd.equals( "HORIZONTAL" ) )
        {

            orientation = PlotOrientation.HORIZONTAL;
        }

        if ( cmd.equals( "VERTICAL" ) )
        {

            orientation = PlotOrientation.VERTICAL;
        }

    }

    /**
     * ShowChart is a routine to show the chart plot in a frame.
     *
     * @param frameTitle becomes the title of the frame
     * @param chartTitle becomes the title of the chart
     * @param plot is the plot of the chart to be displayed
     */
    private static void showChart( String frameTitle, String chartTitle, Plot plot )
    {

        final JFreeChart chart;

        // titleFont is a property of this class

        chart = new JFreeChart( chartTitle, titleFont, plot, true );

        // set the background color for the chart...

        chart.setBackgroundPaint( Color.lightGray );

        // add the chart to a panel...

        final ChartPanel chartPanel = new ChartPanel( chart );

        chartPanel.setPreferredSize( MY_PREFERRED_SIZE );

        JFrame chartFrame = new JFrame( frameTitle );
        chartFrame.setContentPane( chartPanel );
        chartFrame.pack();
        RefineryUtilities.centerFrameOnScreen( chartFrame );
        chartFrame.setVisible( true );

    } // ShowChart

    /**
     * This routine creates for a category data set a bar chart.
     *
     * @param title is the title of the chart
     * @param domainLabel label for domain
     * @param rangeLabel label for range
     * @param dataset input data set for which chart is created
     */
    static void makeBarChart( String title, String domainLabel, String rangeLabel, CategoryDataset dataset )
    {

        CategoryAxis categoryAxis = new CategoryAxis( domainLabel );
        ValueAxis valueAxis = new NumberAxis( rangeLabel );

        BarRenderer renderer = new BarRenderer();

        renderer.setBaseToolTipGenerator( new StandardCategoryToolTipGenerator() );

        CategoryPlot plot = new CategoryPlot( dataset, categoryAxis, valueAxis, renderer );

        plot.setOrientation( ChartProperties.orientation );

        showChart( "ShowPM: Bar Chart", title, plot );

    } // MakeBarChart

    /**
     * This routine creates for a category data set a stacked bar chart.
     *
     * @param title is the title of the chart
     * @param domainLabel label for domain
     * @param rangeLabel label for range
     * @param dataset input data set for which chart is created
     */

    static void makeStackedBarChart( String title, String domainLabel, String rangeLabel, CategoryDataset dataset )
    {

        CategoryAxis categoryAxis = new CategoryAxis( domainLabel );
        ValueAxis valueAxis = new NumberAxis( rangeLabel );

        StackedBarRenderer renderer = new StackedBarRenderer();

        // tooltips are added

        renderer.setBaseToolTipGenerator( new StandardCategoryToolTipGenerator() );

        CategoryPlot plot = new CategoryPlot( dataset, categoryAxis, valueAxis, renderer );

        plot.setOrientation( orientation );

        showChart( "ShowPM: Stacked Bar Chart", title, plot );

    } // MakeStackedBarChart

    /**
     * This routine creates for a category data set a statistical bar chart.
     *
     * @param title is the title of the chart
     * @param domainLabel label for domain
     * @param rangeLabel label for range
     * @param dataset input data set for which chart is created
     */

    static void makeStatBarChart( String title, String domainLabel, String rangeLabel, CategoryDataset dataset )
    {

        final CategoryAxis xAxis = new CategoryAxis( domainLabel );
        final ValueAxis yAxis = new NumberAxis( rangeLabel );

        // define the plot

        final CategoryItemRenderer renderer = new StatisticalBarRenderer();

        final CategoryPlot plot = new CategoryPlot( dataset, xAxis, yAxis, renderer );

        plot.setOrientation( ChartProperties.orientation );

        showChart( "ShowPM: Statistical Bar Chart", title, plot );

    } // MakeStatBarChart

    /**
     * This routine creates for a category data set a pie chart.
     *
     * @param title is the title of the chart
     * @param domainLabel label for domain
     * @param rangeLabel label for range
     * @param dataset input data set for which chart is created
     * @param order can be BY_ROW (chart for every counter) or BY_COLUMN (chart for every region)
     */

    static void makePieChart( final String title, final String domainLabel, final String rangeLabel, final CategoryDataset dataset,
                              final TableOrder order )
    {

        // BY_ROW    : chart for every counter
        // BY_COLUMN : chart for every region

        MultiplePiePlot plot = new MultiplePiePlot( dataset );

        plot.setDataExtractOrder( order );
        plot.setBackgroundPaint( null );
        plot.setOutlineStroke( null );

        PieToolTipGenerator tooltipGenerator = new StandardPieToolTipGenerator();

        PiePlot p = ( PiePlot ) plot.getPieChart().getPlot();
        p.setToolTipGenerator( tooltipGenerator );

        if ( showLabels )
        {

            final double gap = 0.30;

            p.setLabelGenerator( new StandardPieSectionLabelGenerator( "{0}" ) );
            p.setLabelFont( new Font( "SansSerif", Font.PLAIN, 8 ) );
            p.setInteriorGap( gap );

        }
        else
        {

            p.setLabelGenerator( null );

        }

        // Note: PiePlot has not a method setOrientation

        showChart( "ShowPM: Pie Chart", title, plot );

    } // MakePieChart

} // ChartProperties
