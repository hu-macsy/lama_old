/*
 * CallDescriptor.java
 *
 * Descriptor for a call in the RIF file
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

import java.io.File;
import java.io.FileInputStream;
import java.io.InputStream;

import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;

import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JPanel;

import org.apache.log4j.Logger;

import org.jfree.ui.ExtensionFileFilter;

import adaptor.General.OptionsAdaptor;
import adaptor.General.RIFData;
import adaptor.General.ShowInfo;
import att.grappa.Graph;
import att.grappa.Parser;
import att.grappa.Subgraph;

/**
 * CalltreePM: class reads callgraph files and prints output.
 *
 * TODO: legend for the colors
 *
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */
public class CalltreePM implements CTInterface
{

    /**
     * Logger for this class.
     */
    private static Logger logger = Logger.getLogger( CalltreePM.class );

    /**
     * Array of file suffixes for which we have export support here.
     */
    private static final String[] EXPORT_SUFFIXES = { "png", "jpg", "fig", "ps", RIFData.SUFFIX, CTDataFile.PM_SUFFIX, CTDataFile.SUFFIX, "Bar Chart", "Pie Chart" };

    /**
     * Default suffix for calltree files generated at runtime.
     */
    private static final String SUFFIX = ".ct";

    /**
     * Calltree data for each file.
     */
    private static CTData[] myCTFiles = null;

    /**
     * noGroups is a constant for the number of groups supported by NodeSelectionFrame.
     */
    private static final int NO_GROUPS = NodeSelectionPanel.getGroupItems().length;

    /**
     * Calltree data for each group.
     */
    private static CTData[] myCTGroups = new CTData[NO_GROUPS];

    /**
     * This is the currently selected CTData.
     */
    private static CTData myCT = null; // currently selected CTData

    /**
     * Calltree data of regions.
     */
    private static CTData myCTRegion = null;

    /**
     * Calltree data of current group (can be null).
     */
    private static CTData myCTGroup = null;

    /**
     * GraphFrame for the display of the Call Graph.
     */
    private static GraphFrame myGraphFrame = null;

    /**
     * Additional panel to be used with calltree and its menu.
     */
    private static CallEnvironmentPanel myCallEnvironmentPanel = null;

    /**
     * Additional frame to be used with calltree and its menu.
     */
    private static NodeSelectionPanel myNodeSelectionPanel = null;

    /**
     * Additional frame to be used with calltree and its menu.
     */
    private static FileSelectionPanel myFileSelectionPanel = null;

    /**
     *
     */
    private int doGroup = 0;

    /**
     * {@inheritDoc}
     *
     * @see adaptor.Calltree.CTInterface#getExportStrings()
     */
    public String[] getExportStrings()
    {

        return EXPORT_SUFFIXES;

    } // getExportStrings

    /**
     * {@inheritDoc}
     *
     * @see adaptor.Calltree.CTInterface#showGraphFile(java.lang.String)
     */
    public void showGraphFile( String outputFileName )
    {

        // result of layout is in output file output_name

        Parser program = null;

        try
        {

            InputStream input = new FileInputStream( outputFileName );
            program = new Parser( input, System.err );
            program.parse();

        }
        catch ( Exception ex )
        {

            final int exitError = -1;

            logger.error( ex.getMessage(), ex );
            System.exit( exitError );

        } // try - catch

        // now we have parsed the graph and get it

        logger.info( "showGraphFile: output Graph parsed succesfully" );

        Graph g = program.getGraph();

        // now we show the Graph in the CT Frame

        showGraph( g );

    } // showGraphFile

    /**
     * {@inheritDoc}
     *
     * @see adaptor.Calltree.CTInterface#redraw()
     */
    public void redraw()
    {

        // runs as an own thread

        String infile = "graph_in.dot";

        logger.info( "write the current graph in file " + infile );

        myCT.printGraph( infile );

        logger.info( "written the graph, now make layout thread" );

        LayoutThread localThread = new LayoutThread( infile, this );

        localThread.start();

        logger.info( "layout will be printed when ready, but enable also kill button" );

        myGraphFrame.enableKill( localThread );

    } // newLayout

    /**
     * {@inheritDoc}
     *
     * @see adaptor.Calltree.CTInterface#recalc()
     */
    public void recalc()
    {

        // get the node selection

        boolean onlyRoot = myNodeSelectionPanel.getRootFlag();
        boolean onlyLeaf = myNodeSelectionPanel.getLeafFlag();

        int selectedGroup = myNodeSelectionPanel.getSelectedGroup();


        // Step 1: select the nodes according to the group

        setNodeSelection( selectedGroup, onlyRoot, onlyLeaf );

        String pattern = myNodeSelectionPanel.getPattern();

        // Step 2: apply the search pattern

        if ( pattern != null )
        {

            searchNodes( pattern );

        }

        // Step 3: select also the environment

        if ( myNodeSelectionPanel.getNeighborsEnabled() )
        {

            addNeighbors();

        }

        redraw();

    } // redraw

    /**
     * This routine will set the corresponding selection of nodes for the graph frame.
     *
     * @param groupSelection specifies how to group the nodes, 0 is no grouping
     * @param rootFlag if true only root nodes are selected
     * @param leafFlag if true only leaf nodes are selected
     */

    private void setNodeSelection( int groupSelection, boolean rootFlag, boolean leafFlag )
    {

        // make the correct selection for the group

        logger.info( "setNodeSelection (groupSelection=" + groupSelection + ")" );

        if ( groupSelection == 0 )
        {

            if ( myFileSelectionPanel != null )
            {

                boolean[] filesEnabled = myFileSelectionPanel.getSelectedFiles();

                logger.info( filesEnabled.length + " files, inherit costs from the enabled" );

                myCTRegion.resetCosts();

                for ( int i = 0; i < filesEnabled.length; i++ )
                {

                    if ( filesEnabled[i] )
                    {

                        logger.info( "inherit costs from file " + i );

                        myCTRegion.inheritCosts( myCTFiles[i] );

                    }
                }
            }

            myCT = myCTRegion;
            myCTGroup = null;

        }
        else
        {

            if ( myCTGroups[groupSelection] == null )
            {

                // Calltree for this group not available yet

                myCTGroups[groupSelection] = new CTDataGroup( myCTRegion, groupSelection );

            }
            else
            {

                myCTGroups[groupSelection].resetView();
            }

            myCTGroup = myCTGroups[groupSelection];
            myCT = myCTGroup;

        } // GroupSelection

        // set for the current Calltree the parent correctly

        myCT.setParent();

        // make sure that the nodes of the
        if ( rootFlag || leafFlag )
        {

            myCT.markRootOrLeaf( rootFlag, leafFlag );

        }
        else
        {

            myCT.resetSubGraph();
        }

    }

    /**
     * {@inheritDoc}
     *
     * @see adaptor.Calltree.CTInterface#showCallEnvironment()
     */
    public void showCallEnvironment()
    {

        Subgraph subg = myGraphFrame.getSubgraph();

        CTNode[] selectedNodes = MySubgraph.getSelectedNodes( subg );

        if ( selectedNodes == null )
        {

            JOptionPane.showMessageDialog( null, "No nodes are selected", "showCallEnvironment", JOptionPane.ERROR_MESSAGE );

            return;
        }

        MySubgraph.fixSelectedNodes( selectedNodes );

        myCT.buildSubGraph( selectedNodes, myCallEnvironmentPanel.getCalledDepth(), myCallEnvironmentPanel.getCallingDepth() );

        // be careful: redraw would also make a new selection that we do not want here.

        redraw();

    } // showCallEnvironment

    /**
     * This routine is used to restrict the selected nodes by a search pattern.
     *
     * @param pattern is the pattern to be applied for the nodes
     */
    private void searchNodes( String pattern )
    {

        logger.info( "search nodes by this pattern: " + pattern );

        myCT.searchNodes( pattern );

    }

    /**
     * {@inheritDoc}
     *
     * @see adaptor.Calltree.CTInterface#searchNodes(java.lang.String, boolean)
     */
    public void addNeighbors()
    {

        int callingDepth = myCallEnvironmentPanel.getCallingDepth();
        int calledDepth = myCallEnvironmentPanel.getCalledDepth();

        logger.info( "addNeightbors, callingDepth = " + callingDepth + ", calledDepth = " + calledDepth );

        myCT.markEnvironment( calledDepth, callingDepth );

    }


    /**
     * This routine is used to choose a file with a given extension. It returns null if the
     * operation has been canceled otherwise the name of the chosen file.
     *
     * @param description
     *            is an additional description for the selection (hint).
     * @param extension
     *            specifies the extension type for which a file should be chosen.
     * @return the chosen filename (null if canceled).
     */
    private String getFileName( String description, String extension )
    {

        String currentDir = System.getProperty( "user.dir" );

        JFileChooser chooser = new JFileChooser( currentDir );

        chooser.setDialogTitle( "Choose export file" );

        ExtensionFileFilter filter = new ExtensionFileFilter( description, extension );

        chooser.setFileFilter( filter );

        int result = chooser.showSaveDialog( myGraphFrame );

        if ( result == JFileChooser.CANCEL_OPTION )
        {

            return null;
        }

        File selectedFile = chooser.getSelectedFile();

        // make sure that that really a file has been selected

        if ( selectedFile == null )
        {

            return null;
        }

        return selectedFile.getAbsolutePath();

    }

    /**
     * {@inheritDoc}
     *
     * @see adaptor.Calltree.CTInterface#export(java.lang.String)
     */
    public void export( String suffix )
    {

        logger.info( "Export current graph, suffix = " + suffix );

        if ( suffix.equals( CTDataFile.PM_SUFFIX ) )
        {

            // Attention: we cannot export grouped nodes

            if ( myCTRegion instanceof CTDataFile )
            {

                String help = "Performance Monitoring Statistic (*." + CTDataFile.PM_SUFFIX + ")";

                String filename = getFileName( help, CTDataFile.PM_SUFFIX );

                if ( filename != null )
                {

                    ( ( CTDataFile ) myCTRegion ).exportPMS( filename );
                }

            }
            else
            {

                JOptionPane.showMessageDialog( null, "No Calltree performance file", "export PMS file", JOptionPane.ERROR_MESSAGE );

            }

        }
        else if ( suffix.equals( CTDataFile.SUFFIX ) )
        {

            // Attention: we cannot export grouped nodes

            if ( myCTRegion instanceof CTDataFile )
            {

                String filename = getFileName( "CallTree Data File (*.ct)", "ct" );

                if ( filename != null )
                {

                    ( ( CTDataFile ) myCTRegion ).output( filename );
                }

            }
            else
            {

                JOptionPane.showMessageDialog( null, "No Calltree performance file", "export CT file", JOptionPane.ERROR_MESSAGE );

            }

        }
        else if ( suffix.equals( RIFData.SUFFIX ) )
        {

            // Attention: we cannot export grouped nodes ?

            String filename = getFileName( "Region Information File (*." + RIFData.SUFFIX + ")", RIFData.SUFFIX );

            if ( filename == null )
            {

                logger.info( "Export RIF file has been canceled" );

            }
            else
            {

                myCTRegion.exportRIF( filename );

                logger.info( "have exported RIF file " + filename );
            }

        }
        else if ( suffix.equals( "Bar Chart" ) )
        {

            myCTRegion.makeChart( true );

        }
        else if ( suffix.equals( "Pie Chart" ) )
        {

            myCTRegion.makeChart( false );

        }
        else
        {

            // generates the graph with the current selected properties

            String infile = "graph_out.dot";

            String outfile = getFileName( "Output file for graph (*." + suffix + ")", suffix );

            if ( outfile == null )
            {

                logger.info( "Export of graph file has been canceled" );

            }
            else
            {

                myCT.printGraph( infile );

                LayoutThread thread = new LayoutThread( infile, outfile, suffix );

                thread.start(); // start the layout thread

                logger.info( "Export thread has been started" );

                // export Thread cannot be killed at this time:w

            }

        }

    } // export

    /**
     * {@inheritDoc}
     *
     * @see adaptor.Calltree.CTInterface#resetUnfolding()
     */
    public void resetUnfolding()
    {

        if ( myCT.resetUnfolding() )
        {

            redraw();
        }

    } // resetUnfolding

    /**
     * {@inheritDoc}
     *
     * @see adaptor.Calltree.CTInterface#showInfo(java.lang.Object)
     */
    public void showInfo( Object infoObject )
    {

        if ( infoObject instanceof CTNode )
        {

            CTNode node = ( CTNode ) infoObject;

            ShowInfo.showInfoText( node.infoNode() );

        }
        else
        {

            // for all other objects we show full info about current CT data

            ShowInfo.showInfoText( myCT.getInfo() );

        }

    } // showInfo

    /**
     * Displays a call graph (internal data structure of a graph with given layout).
     *
     * @param aGraph
     *            is the call graph to be displayed.
     */
    private void showGraph( Graph aGraph )
    {

        myGraphFrame.disableKillButton();

        myGraphFrame.showGraph( aGraph, this );

        logger.info( "showGraph of thread has been done" );

    }

    /**
     * This routine is used to get the node in the CT data by the identification of the node.
     *
     * @param printIdent
     *            is the internal identification of the node
     * @return the CTNode with the given identification (can be null)
     */
    public static CTNode getNode( String printIdent )
    {

        // our first try is in current calltree graph

        CTNode node = myCT.getNode( printIdent );

        if ( node == null )
        {

            // our second try is on the full calltree graph

            if ( myCTRegion != null )
            {
                if ( myCTRegion != myCT )
                {
                    node = myCTRegion.getNode( printIdent );
                }
            }
        }

        if ( node == null )
        {

            // our last try is on the group calltree graph

            if ( myCTGroup != null )
            {
                if ( myCTGroup != myCT )
                {
                    node = myCTGroup.getNode( printIdent );
                }
            }
        }

        return node;

    } // getNode

    /**
     * Checks whether we have a valid filename.
     *
     * @param arg
     *            is the name of the file to be opened.
     * @return a valid filename that can be opened or null.
     */
    String validFilename( String arg )
    {

        // check for a valid file

        String validName = arg;

        File testFile = new File( validName );

        if ( !testFile.exists() )
        {

            // so we give it a try on filename with Suffix

            validName = arg + SUFFIX;

            testFile = new File( validName );

            if ( !testFile.exists() )
            {

                // did not work, so we set it to null

                validName = null;
            }
        }

        return validName;
    }

    /**
     * This routine is called after the evaluation of the arguments.
     *
     * @param theOptions
     *            are the evaluated options
     * @param filenames
     *            are the input filenames
     */
    void doit( OptionsAdaptor theOptions, String[] filenames )
    {

        // -root will show only root nodes

        boolean onlyRoot = theOptions.isOptionSet( "-root" );

        // -leaf will show only leaf nodes

        boolean onlyLeaf = theOptions.isOptionSet( "-leaf" );

        boolean useFictiveCounter = theOptions.isOptionSet( "-fictive" );

        boolean isRIF = theOptions.isOptionSet( "-rif" );

        String groupVal = theOptions.getOptionVal( "-group" );

        doGroup = 0;

        String[] groups = NodeSelectionPanel.getGroupItems();

        if ( groupVal != null )
        {

            for ( int i = 0; i < groups.length; i++ )
            {

                if ( groupVal.equals( groups[i] ) )
                {

                    doGroup = i;
                }
            }
        }

        if ( filenames.length < 1 )
        {

            logger.error( "no input files specified" );
            theOptions.printUsage();
            System.exit( 0 );
        }

        if ( isRIF )
        {

            if ( filenames.length > 1 )
            {

                logger.error( "multiple RIF files not supported yet" );
                theOptions.printUsage();
            }

            String nameRIF = filenames[0];

            logger.info( "construct RIF data from file " + nameRIF );

            RIFData theRIFData = new RIFData( nameRIF );

            myCTRegion = new CTDataRIF( theRIFData, useFictiveCounter );

            CTProperties.setDefaults();

            // without fictive costs we only have source calls

            if ( !useFictiveCounter )
            {

                CTProperties.setSourceCallDefaults();
            }

        }
        else
        {

            int noFiles = filenames.length;

            myCTFiles = new CTData[noFiles];

            for ( int i = 0; i < noFiles; i++ )
            {

                logger.info( "make CT data from file (" + ( i + 1 ) + " of " + noFiles + ") " + filenames[i] );

                CTDataFile ctDataFile = new CTDataFile( filenames[i] );

                // resolve hexadeciaml addresses

                ctDataFile.addr2line();

                myCTFiles[i] = ctDataFile;

            }

            if ( noFiles > 1 )
            {

                // make a file selection menu (will be used later)

                myCTRegion = new CTDataUnion( myCTFiles );

            }
            else
            {

                myCTRegion = myCTFiles[0];
            }

            CTProperties.setDefaults();

        }

        String outFile = theOptions.getOptionVal( "-o" );

        if ( outFile != null )
        {

            logger.info( "write compressed / resolved data out to file " + outFile );

            if ( myCTFiles == null )
            {

                logger.error( "sorry, no CT data available for output" );

            }
            else
            {

                for ( int i = 0; i < myCTFiles.length; i++ )
                {

                    CTData ctData = myCTFiles[i];

                    if ( ctData instanceof CTDataFile )
                    {

                        CTDataFile ctDataFile = ( CTDataFile ) ctData;

                        if ( myCTFiles.length > 1 )
                        {

                            ctDataFile.output( outFile + "." + ( i + 1 ) );

                        }
                        else
                        {

                            ctDataFile.output( outFile );
                        }

                    }
                    else
                    {

                        logger.error( "no CT Data file" );
                    }

                }

            }

            System.exit( 0 );
        }

        logger.info( "data has been read in, now set up the GUI" );

        myCT = myCTRegion;

        myCallEnvironmentPanel = new CallEnvironmentPanel();

        myNodeSelectionPanel = new NodeSelectionPanel( this, doGroup, onlyRoot, onlyLeaf );

        if ( filenames.length > 1 )
        {

            myFileSelectionPanel = new FileSelectionPanel( filenames );
        }

        myGraphFrame = new GraphFrame( new CTMenuBar( this ) );

        // now we create an array of selection panels

        JPanel[] selectionPanels = { myNodeSelectionPanel, myCallEnvironmentPanel, myFileSelectionPanel };

        myGraphFrame.addTab( "Properties", new CTPropertiesPanel( this, selectionPanels ) );

        logger.info( "all tabbed panels added" );

        // now make the layout of the first graph and show it

        setNodeSelection( doGroup, onlyRoot, onlyLeaf );

        recalc();

        myGraphFrame.setDefaultCloseOperation( JFrame.EXIT_ON_CLOSE );

        myGraphFrame.setVisible( true );

    }

    /**
     * The main routine of CalltreePM is used for running it as an application. Reads the command
     * line arguments. Reads the inputfiles and shows the Mainframe with the call tree graph
     *
     * @param arguments
     *            [-rif] [-group=xxx] [-root] [-leaf] filename ...
     */
    public static void main( String[] arguments )
    {

        final int errorCode = -1;

        // Main program should start with an invocation of OptionsAdaptor

        String help = "Call Graph with Performance Monitoring Data";

        OptionsAdaptor theOptions = new OptionsAdaptor( "CalltreePM", help );

        // Now we defined the options that can be used

        theOptions.addFlagOption( "-rif", "inputfile is RIF file" );
        theOptions.addFlagOption( "-root", "display only root nodes" );
        theOptions.addFlagOption( "-leaf", "display only leaf nodes" );
        theOptions.addFlagOption( "-fictive", "add fictive performance event" );
        theOptions.addValueOption( "-group", "build group of nodes", NodeSelectionPanel.getGroupItems() );
        theOptions.addValueOption( "-o", "output file " );

        String[] inputFiles = null;

        try
        {

            inputFiles = theOptions.parseArguments( arguments );

        }
        catch ( IllegalArgumentException e )
        {

            logger.error( e.getMessage() );
            theOptions.printUsage();
            System.exit( errorCode );

        }

        logger.info( "main program starts" );

        CalltreePM mainCT = new CalltreePM();

        mainCT.doit( theOptions, inputFiles );

        logger.info( "main program terminates, GUI remains" );

    }

} // class CalltreePM
