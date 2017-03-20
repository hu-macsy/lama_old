/*
 * ExtraPanel.java
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

import java.io.File;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.SwingConstants;

import org.apache.log4j.Logger;

/***************************************************************************************************
 * * ExtraPane: Pane used to specify data profiling, trace, RIF * *
 **************************************************************************************************/

/**
 * The class ExtraPanel defines a panel for all attributes that are not related to counter.
 *
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */

public class ExtraPanel extends JPanel implements ActionListener
{

    /**
     * String value needed when something is enabled.
     */
    private static final String ENABLED_STRING = "enabled";

    /**
     * String value needed when something is disabled.
     */
    private static final String DISABLED_STRING = "disabled";

    /**
     * Logger for this class.
     */
    private static Logger logger = Logger.getLogger( ExtraPanel.class );

    /**
     * Pointer to get infos about the PM configuration.
     */
    private PMConfiguration myPMConfiguration = null;

    /**
     * Command is used to display the current execution command, e.g.
     * hello_world -pm=pm_config.
     */
    private JLabel executeCommand = null;

    /**
     * Check box for profiling depth.
     */
    private JCheckBox myProfDepthBox;

    /**
     * Value button for the profiling depth.
     */
    private ValueButton depthButton;

    /**
     * Check box for data profiling.
     */
    private JCheckBox myDataProfCheckBox;

    /**
     * Value button for the rate of data profiling.
     */
    private ValueButton rateButton;

    /**
     * Value button for the size of data profiling.
     */
    private ValueButton sizeButton;

    /**
     * Check box for tracing.
     */
    private JCheckBox traceEnableBox;

    /**
     * Label for tracing.
     */
    private JLabel    traceLabel;

    /**
     * Combo box to select how regions are grouped for tracing.
     */
    private JComboBox traceGroup;

    /**
     * Check box if calltree data is collected.
     */
    private JCheckBox myCalltreeCheckBox;

    /**
     * Label for calltree flag.
     */
    private JLabel calltreeLabel;

    /**
     * Check box if calltree data for data arrays is collected.
     */
    private JCheckBox myCalldataCheckbox;

    /**
     * Label for corresponding check box.
     */
    private JLabel calldataLabel;

    /**
     * Check box for RIF filename. If enabled, an
     * RIF filename is specified.
     */
    private JCheckBox myRIFCheckBox;

    /**
     * Value button for the current RIF file.
     */
    private ValueButton myRIFNameButton;

    /**
     * Check box for output file name.
     */
    private JCheckBox myOutCheckBox;

    /**
     * Value button for the output file name.
     */
    private ValueButton myOutfileNameButton;

    /**
     * Constructor for this extra panel that
     * allows to define special value for the PM configuration.
     *
     * @param config is reference to the current configuration
     */
    public ExtraPanel( PMConfiguration config )
    {

        super( new BorderLayout( 3, 3 ) );

        final int defaultSamplingRate = 128;
        final int defaultSamplingSize = 50000;
        final int defaultDepth        = Integer.MAX_VALUE;

        myPMConfiguration = config;

        setBackground( Color.black );

        JPanel mainPanel = new JPanel();
        mainPanel.setLayout( new GridLayout( 7, 3, 2, 2 ) );

        String name = myPMConfiguration.getExecutable();
        executeCommand = new JLabel( name, SwingConstants.CENTER );
        executeCommand.setForeground( Color.yellow );
        executeCommand.setBackground( Color.black );

        // Row 1:  OutputFile   <name>    --

        myOutCheckBox = new JCheckBox( "Output File", null, false );
        myOutCheckBox.addActionListener( this );
        mainPanel.add( myOutCheckBox );

        myOutfileNameButton = new ValueButton( "file name", "<default>", this );
        mainPanel.add( myOutfileNameButton );

        mainPanel.add( new JLabel( "" ) );

        /* next row: profiling depth */

        myProfDepthBox = new JCheckBox( "Profiling Depth", null, false );
        myProfDepthBox.addActionListener( this );
        mainPanel.add( myProfDepthBox );

        depthButton = new ValueButton( "max profiling depth", new Integer( defaultDepth ), this );
        mainPanel.add( depthButton );

        mainPanel.add( new JLabel( "" ) );

        /* next row: data profiling */

        myDataProfCheckBox = new JCheckBox( "Data Profiling", null, false );
        myDataProfCheckBox.addActionListener( this );
        mainPanel.add( myDataProfCheckBox );

        rateButton = new ValueButton( "sampling rate", new Integer( defaultSamplingRate ), this );
        mainPanel.add( rateButton );

        sizeButton = new ValueButton( "sampling size", new Integer( defaultSamplingSize ), this );
        mainPanel.add( sizeButton );

        traceEnableBox = new JCheckBox( "Trace", null, false );
        traceEnableBox.addActionListener( this );
        mainPanel.add( traceEnableBox );

        traceLabel = new JLabel( DISABLED_STRING );
        mainPanel.add( traceLabel );

        traceGroup = new JComboBox();

        traceGroup.addItem( "GROUP=ALL" );
        traceGroup.addItem( "GROUP=FILE" );
        traceGroup.addItem( "GROUP=CLASS" );
        traceGroup.addItem( "GROUP=SINGLE" );
        traceGroup.setToolTipText( "select item which regions build a group" );
        traceGroup.addActionListener( this );

        traceGroup.setSelectedIndex( myPMConfiguration.getGroupIndex() );

        mainPanel.add( traceGroup );

        myCalltreeCheckBox = new JCheckBox( "Calltree", null, false );
        myCalltreeCheckBox.addActionListener( this );
        mainPanel.add( myCalltreeCheckBox );

        calltreeLabel = new JLabel( DISABLED_STRING );
        mainPanel.add( calltreeLabel );

        mainPanel.add( new JLabel( "" ) );

        // start Grid line for Calldata

        myCalldataCheckbox = new JCheckBox( "Calldata", null, false );
        myCalldataCheckbox.addActionListener( this );
        mainPanel.add( myCalldataCheckbox );

        calldataLabel = new JLabel( DISABLED_STRING );
        mainPanel.add( calldataLabel );

        mainPanel.add( new JLabel( "" ) );

        myRIFCheckBox = new JCheckBox( "Region Information File", null, false );
        myRIFCheckBox.addActionListener( this );
        mainPanel.add( myRIFCheckBox );

        myRIFNameButton = new ValueButton( "file name", new File( "RIF" ), this );

        mainPanel.add( myRIFNameButton );

        mainPanel.add( new JLabel( "" ) );

        showValues();

        add( "Center", mainPanel );
        add( "North", executeCommand );

    }

    /**
     * This routine returns the string for enabled or disabled by
     * a boolean flag.
     *
     * @param val decieds about enabled or disabled
     * @return string representation for enabled or disabled
     */
    private static String enabledText( boolean val )
    {

        if ( val )
        {
            return ENABLED_STRING;
        }
        else
        {
            return DISABLED_STRING;
        }
    }

    /**
     * This routine shows all the extra values in the panel.
     *
     */
    public void showValues()
    {

        boolean enabled;

        // profiling depth

        enabled = myPMConfiguration.getProfilingDepth() >= 0;

        if ( enabled )
        {

            int depth = myPMConfiguration.getProfilingDepth();
            depthButton .setValue( new Integer( depth ) );

        }

        depthButton.setEnabled( enabled );
        myProfDepthBox.setSelected( enabled );

        // dataProfiling

        enabled = myPMConfiguration.doDataProfling();

        myDataProfCheckBox.setSelected( enabled );

        if ( enabled )
        {

            int rate = myPMConfiguration.getSamplingRate();
            int size = myPMConfiguration.getSamplingSize();

            rateButton.setValue( new Integer( rate ) );
            sizeButton.setValue( new Integer( size ) );

        }

        rateButton.setEnabled( enabled );
        sizeButton.setEnabled( enabled );

        // calltree

        enabled = myPMConfiguration.doCalltree();
        myCalltreeCheckBox.setSelected( enabled );
        calltreeLabel.setText( enabledText( enabled ) );

        // calldata

        enabled = myPMConfiguration.doCalldata();
        myCalldataCheckbox.setSelected( enabled );
        calldataLabel.setText( enabledText( enabled ) );

        // trace

        enabled = myPMConfiguration.doTracing();
        traceEnableBox.setSelected( enabled );
        traceLabel.setText( enabledText( enabled ) );
        traceGroup.setSelectedIndex( myPMConfiguration.getGroupIndex() );

        // rifFile

        String rifFile = myPMConfiguration.getRIF();
        enabled = ( rifFile != null );

        myRIFCheckBox.setSelected( enabled );
        myRIFNameButton.setEnabled( enabled );

        if ( enabled )
        {
            myRIFNameButton.setValue( new File( rifFile ) );
        }

        // outFile

        String outFile = myPMConfiguration.getOutFile();
        myOutfileNameButton.setEnabled( outFile != null );

        // command for profiling

        String command = myPMConfiguration.getExecutable();
        String configFile = myPMConfiguration.getConfigFileName();

        if ( configFile != null )
        {
            if ( !configFile.equals( "" ) )
            {
                command += " -pm=" + configFile;
            }
        }

        executeCommand.setText( command );

    }

    /**
     * This routine is called when a checkBox has changed its status.
     * Its new status can be asked for by its methods.
     *
     * @param checkBox is the selected checkBox
     */
    private void actualizeCheckBox( JCheckBox checkBox )
    {

        boolean enabled = checkBox.isSelected();

        if ( checkBox == myDataProfCheckBox )
        {

            myPMConfiguration.setDataProfling( enabled );

        }
        else if ( checkBox == myCalltreeCheckBox )
        {

            myPMConfiguration.setCalltree( enabled );

        }
        else if ( checkBox == myCalldataCheckbox )
        {

            myPMConfiguration.setCalldata( enabled );

        }
        else if ( checkBox == traceEnableBox )
        {

            myPMConfiguration.setTracing( enabled );

        }
        else if ( checkBox == myProfDepthBox )
        {

            if ( enabled )
            {

                Integer depth = ( Integer ) depthButton.getValue();

                int val = depth.intValue();

                if ( val < 0 )
                {

                    val = 0;
                }

                myPMConfiguration.setProfilingDepth( val );

            }
            else
            {

                myPMConfiguration.setProfilingDepth( -1 );
            }


        }
        else if ( checkBox == myOutCheckBox )
        {

            if ( enabled )
            {

                myPMConfiguration.setOutFile( myOutfileNameButton.getText() );

            }
            else
            {

                myPMConfiguration.setOutFile( null );
            }

        }
        else if ( checkBox == myRIFCheckBox )
        {

            if ( enabled )
            {

                myPMConfiguration.setRIF( myRIFNameButton.getValue() + "" );

            }
            else
            {

                myPMConfiguration.setRIF( null );
            }
        }

        showValues();

    }

    /**
     * This routine is called when one of my value buttons has been
     * changed. The new value has to be set in the configuration
     * data.
     *
     * @param button is the button with a new value
     */
    private void actionValue( ValueButton button )
    {

        if ( button == depthButton )
        {

            Integer depth = ( Integer ) button.getValue();

            myPMConfiguration.setProfilingDepth( depth.intValue() );

            logger.info( "depth = " + depth.intValue() );

        }
        else if ( button == rateButton )
        {

            Integer val = ( Integer ) button.getValue();

            myPMConfiguration.setSamplingRate( val.intValue() );

            logger.info( "rateSampling = " + val.intValue() );

        }
        else if ( button == sizeButton )
        {

            Integer val = ( Integer ) button.getValue();

            myPMConfiguration.setSamplingSize( val.intValue() );

            logger.info( "sizeSampling = " + val.intValue() );

        }
        else if ( button == myRIFNameButton )
        {

            File f = ( File ) button.getValue();

            myPMConfiguration.setRIF( f.getName() );

            logger.info( "RIF file = " + f.getName() );

        }
        else if ( button == myOutfileNameButton )
        {

            String s = ( String ) button.getValue();

            myPMConfiguration.setOutFile( s );

            logger.info( "Out file = " + s );

        }


    }

    /**
     * {@inheritDoc}
     *
     * @see java.awt.event.ActionListener#actionPerformed(java.awt.event.ActionEvent)
     */
    public void actionPerformed( ActionEvent e )
    {

        Object obj = e.getSource();

        // JCheckBox is the only kind of items we support

        if ( obj instanceof JCheckBox )
        {

            actualizeCheckBox( ( JCheckBox ) obj );

        }
        else if ( obj instanceof ValueButton )
        {

            ValueButton button = ( ValueButton ) obj;

            logger.info( "Value Button has changed, name = "
                         + button.getValue() );

            actionValue( button );

        }
        else if ( obj == traceGroup )
        {

            int index = traceGroup.getSelectedIndex();

            logger.info( "trace group is now " + index );

            myPMConfiguration.setGroupIndex( index );
        }

    }

}
