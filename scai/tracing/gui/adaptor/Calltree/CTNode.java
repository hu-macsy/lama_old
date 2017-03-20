/*
 * CTNode.java
 *
 * A class used for representing nodes in the Calltree.
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
import java.util.List;
import java.util.Vector;

import org.apache.log4j.Logger;

import adaptor.General.CallNode;
import adaptor.General.FileTable;
import adaptor.General.RegionDescriptor;
import adaptor.General.Scaling;

/**
 * CTNode contains data and methods for nodes in the call tree.
 *
 * @version $LastChangedRevision$
 * @author Dr. Thomas Brandes
 */

public class CTNode
{

    /**
     * The names, labels and colors of nodes must be included apostrophs.
     * This special string is used for this to have a unique representation.
     */
    private static final String APO = "\"";

    /**
     * Logger for this class.
     */
    private static Logger logger = Logger.getLogger( CTNode.class );

    /**
     * myCallNode is taken as entity (actually we could have extended
     * this class).
     */

    private CallNode myCallNode = null;

    /**
     * myIdent will be a unique identification of a node in one CT data set as the name of node or
     * region is not unique.
     */
    private String myIdent;

    /**
     * Value for how often this node is called in source code.
     */
    private long noSourceCalls = 0;

    /**
     * Value for how often this node has been called at runtime.
     */
    private long noRuntimeCalls = 0;

    /**
     * Exclusive costs for this node, value for each counter.
     */
    private long[] myExclusiveCounter;

    /**
     * Inclusive costs for this node, value for each counter.
     */
    private long[] myInclusiveCounter;

    /**
     * String for the first label of this node.
     */
    private String myLabel1 = "";

    /**
     * String for the second label of this node.
     */
    private String myLabel2 = "";

    // Color values if color is fixed

    /**
     * Hue value for the color of this node.
     */
    private double myHue = 0.0;

    /**
     * Saturation value for the color of this node.
     */
    private double mySaturation = 1.0;

    /**
     * Needed for algorithms on this node.
     */
    private boolean myVisited  = false;

    /**
     * Mark flag is used to hide or unselected nodes. If a
     * node is not marked it will not be displayed.
     */
    private boolean myMarked   = true;

    /**
     * Distance value to nodes that have the distance 0.
     */
    private int     myDistance = 0;

    // Variable needed for unfolding

    /**
     * If the flag is true, this group node will be unfolded.
     */
    private boolean myDoUnfold = false;

    // Nodes can build a Hierarchy

    /**
     * This is the group node to which I belong. It can be null
     * if this node is not grouped at all.
     */
    private CTNode myParent   = null;

    /**
     * Vector of all nodes that are grouped to this node.
     */
    private List mySubNodes = null;

    /**
     * This constructor creates a Calltree node from a call node.
     *
     * @param theCallNode is the call node
     */
    public CTNode( CallNode theCallNode )
    {

        int numberCounters = CounterData.numberEvents();

        // construct a node and init all counter values with 0

        myExclusiveCounter = new long[numberCounters];
        myInclusiveCounter = new long[numberCounters];

        for ( int i = 0; i < numberCounters; i++ )
        {

            myExclusiveCounter[i] = 0;
            myInclusiveCounter[i] = 0;

        } // counter intialization

        this.myCallNode = theCallNode;

    } // CTNode

    /**
     * Getter routine for the number of runtime calls of this node.
     * @return Returns the myNoRuntimeCalls.
     */
    public long getRuntimeCalls()
    {

        return noRuntimeCalls;
    }

    /**
     * Setter routine for the number of runtime calls.
     * @param myNoRuntimeCalls The myNoRuntimeCalls to set.
     */

    public void setRuntimeCalls( long myNoRuntimeCalls )
    {

        this.noRuntimeCalls = myNoRuntimeCalls;
    }

    /**
     * Setter routine for the unique id of the node. It is
     * set by the CT data.
     *
     * @param ident is the unique identification.
     */
    public void setIdent( String ident )
    {

        myIdent = ident;
    }

    /**
     * Asks whether this node is a source root node (main program).
     * It is assumed to be a source root if there are no source
     * calls for this node.
     *
     * @return true if this node is never called in the source code
     */
    public boolean isSourceRoot()
    {

        return ( noSourceCalls == 0 );
    }

    /**
     * Getter routine for the exclusive costs of a node.
     *
     * @return array of all exclusive costs
     */
    public long[] getExclusiveCosts()
    {

        return myExclusiveCounter;
    }

    /**
     * Getter routine for the inclusive costs of a node.
     *
     * @return array of all inclusive costs
     */
    public long[] getInclusiveCosts()
    {

        return myInclusiveCounter;
    }

    /**
     * Setter routine for the exclusive costs.
     *
     * @param values are the exclusvie costs for this node
     */
    public void setExclusiveCosts( long[] values )
    {

        for ( int i = 0; i < myExclusiveCounter.length; i++ )
        {

            myExclusiveCounter[i] = values[i];
        }

    }

    /**
     * Setter routine for the inclusive costs.
     *
     * @param values are the inclusvie costs for this node
     */

    public void setInclusiveCosts( long[] values )
    {

        for ( int i = 0; i < myInclusiveCounter.length; i++ )
        {

            myInclusiveCounter[i] = values[i];
        }
    }

    /**
     * This routine can be used to add inclusive costs for a node.
     *
     * @param values are the counter values for each event
     */
    public void addInclusiveCosts( long[] values )
    {

        for ( int i = 0; i < myInclusiveCounter.length; i++ )
        {

            myInclusiveCounter[i] += values[i];
        }
    }

    /**
     * This routine can be used to add exclusive costs for a node.
     *
     * @param values are the counter values for each event
     */

    public void addExclusiveCosts( long[] values )
    {

        for ( int i = 0; i < myExclusiveCounter.length; i++ )
        {

            myExclusiveCounter[i] += values[i];
        }
    }

    /**
     * Reset all runtime costs (exclusive and inclusive costs,
     * number of runtime calls).
     *
     */
    public void resetCosts()
    {

        for ( int i = 0; i < myInclusiveCounter.length; i++ )
        {

            myInclusiveCounter[i] = 0;
            myExclusiveCounter[i] = 0;

        }

        noRuntimeCalls = 0;

    }

    /**
     * Adding a certain number of runtime calls for this node.
     *
     * @param calls is the number of calls to add
     */
    public void addNoRunCalls( long calls )
    {

        noRuntimeCalls += calls;

    }


    /**
     * Adding a certain number of source calls for this node.
     *
     * @param calls is the number of calls to add
     */

    public void addNoSourceCalls( long calls )
    {

        noSourceCalls += calls;

    }

    /**
     * There is an additional source call of this node.
     *
     */
    public void incSourceCall()
    {

        noSourceCalls++;
    }

    /**
     * This routine adds all the number of calls and costs from the input node N
     * to the own costs.
     *
     * @param node is used for adding costs and calls
     */
    public void addCounters( CTNode node )
    {

        for ( int i = 0; i < myExclusiveCounter.length; i++ )
        {

            myExclusiveCounter[i] = node.myExclusiveCounter[i];
            myInclusiveCounter[i] = node.myInclusiveCounter[i];

        }

        noSourceCalls += node.noSourceCalls;
        noRuntimeCalls += node.noRuntimeCalls;

    } // addCounters

    /**
     * Asking for the visibility of this node.
     *
     * @param inclusiveFlag if true inclusive costs must be about threshold, otherwise
     *        the exclusive costs
     * @return true if node is visible under the given conditions
     */
    public boolean isVisible( boolean inclusiveFlag )
    {

        boolean visible = false;

        if ( myMarked )
        {

            if ( inclusiveFlag )
            {

                visible = CTProperties.aboutDepth( myInclusiveCounter, false );

            }
            else
            {

                visible = CTProperties.aboutDepth( myExclusiveCounter, false );
            }
        }

        return visible;

    } // isVisible

    /**
     * Routine to ask for the visibility of this node. It must be marked
     * and the inclusive costs must be about the depth threshold.
     *
     * @return true if node is visible
     */
    public boolean isVisible()
    {

        return isVisible( true );

    }

    /**
     * Getter routine for set unfold flag.
     *
     * @return the flag that indicates unfolding
     */
    public boolean isUnfolded()
    {

        return myDoUnfold;

    }

    /**
     * Setter routine for the unfold flag.
     *
     * @param flag is vlaue for unfold flag
     */
    public void setUnfolded( boolean flag )
    {

        myDoUnfold = flag;
    }


    /**
     * Setter routine for the mark flag.
     *
     * @param mark is value for the mark flag
     */
    public void setMarked( boolean mark )
    {

        myMarked = mark;
    }

    /**
     * Getter routine for the mark flag.
     *
     * @return the mark flag
     */
    public boolean isMarked()
    {

        return myMarked;
    }

    /**
     * Setter routine for the visit flag.
     *
     * @param flag is the value for the visit flag
     */
    public void setVisited( boolean flag )
    {

        myVisited = flag;
    }

    /**
     * Getter routine for the visit flag.
     *
     * @return true if node has already been visited
     */
    public boolean isVisited()
    {

        return myVisited;
    }

    /**
     * Reset visit flag, mark flag and set distance to 0.
     *
     */
    public void reset()
    {

        myVisited = false;
        myMarked  = false;
        myDistance = 0;
    }

    /**
     * Getter routine for the distance.
     *
     * @return distance value of this node
     */
    public int getDistance()
    {

        return myDistance;
    }

    /**
     * Setter routine for the parent of the node.
     *
     * @param parent is the group node
     */
    public void setParent( CTNode parent )
    {

        myParent = parent;
    }

    /**
     * This routine is used to set the labels of this node. The labels
     * depend on the CT properties.
     *
     */
    public void setLabel()
    {

        myLabel1 = CTProperties.getNodeLabel( false, myLabel1, noSourceCalls, noRuntimeCalls, myDistance,
                                              getNoSubNodes(), myInclusiveCounter, myExclusiveCounter );

        myLabel2 = CTProperties.getNodeLabel( true, myLabel2, noSourceCalls, noRuntimeCalls, myDistance,
                                              getNoSubNodes(), myInclusiveCounter, myExclusiveCounter );

    } // setLabel

    /**
     * This routine gets the exclusive value for the current selected counter
     * or metric.
     *
     * @return value of exclusive counter/metric
     */
    public double getExclusiveValue()
    {

        // return the exclusive value of the current selected counter/metric

        return CTProperties.getValue( myExclusiveCounter );

    }

    /**
     * This routine returns the label that should be used for this node.
     *
     * @return the concatenation of the two labels for this node
     */
    public String getLabel()
    {

        String fullLabel = myLabel1;

        if ( myLabel2.length() > 0 )
        {

            if ( fullLabel.length() > 0 )
            {

                fullLabel += " | ";
            }

            fullLabel += myLabel2;
        }

        // label string itself depends on shape of the node

        return CTProperties.makeNodeLabel( getName(), fullLabel, myCallNode.getNodeKind() );

    } // getLabel

    /**
     * This routine calculates the color for this node and sets is hue
     * and saturation value. These values can be fixed for a certain time.
     *
     */
    public void setColor()
    {

        // this routine is called for all nodes

        myHue = CTProperties.getNodeHue( myInclusiveCounter, getName(), myHue );

        mySaturation = CTProperties.getNodeSaturation( myExclusiveCounter, mySaturation );

    } // setColor

    /***********************************************************************************************
     * * String getColor () * * - combine values Hue (local) / Saturation (local) * - and global
     * Brightness * *
     **********************************************************************************************/

    /**
     * This routine returns a string that stands for the current color of this node.
     *
     * @return String for the color
     */
    public String getColor()
    {

        double theBrightness = CTProperties.getNodeBrightness();

        return CTProperties.makeColorString( myHue, mySaturation, theBrightness );
    }

    /**
     * This routine return the name of the node used for presentation. The name
     * must not be unique among all existing nodes.
     *
     * @return String containing the name of this node
     */
    public String getName()
    {

        return myCallNode.getName();

    }

    /**
     * Getter routine for the unique representation of the node.
     *
     * @return unique id as a string
     */
    public String getPrintIdent()
    {

        // getIdent returns the unique id

        return myIdent;

    }

    /***********************************************************************************************
     * * String getGroupName (int kind) * * - get the group name used for grouping of regions * - 1 :
     * group by class_name "class" * - 2 : group by file_name "file" * - 3 : group by dir_name
     * (level1) "dir-1" * - 4 : group by dir_name (level2) "dir-2" * *
     **********************************************************************************************/


    /**
     * This routine returns the name of the group for a node.
     *
     * <ul>
     * <li>
     * 1 stands for class name
     * <li>
     * 2 group by name of file
     * <li>
     * 3 group by last item of the path
     * <li>
     * 4 group by seond last item of the path
     * </ul>

     *
     * @param kind is the kind of group to build
     * @return the group name
     */
    public String getGroupName( int kind )
    {

        String groupName = null;

        RegionDescriptor dsp = myCallNode.getRegion();

        if ( dsp == null )
        {

            // we take internal, external, dummy  as group name

            groupName = myCallNode.getNodeKindString();

        }
        else if ( kind == 1 )
        {

            groupName = dsp.getClassName();

        }
        else if ( kind == 2 )
        {

            groupName = dsp.getFile().getShortFileName();

        }
        else if ( kind == 3 )
        {

            groupName = dsp.getFile().getFilePath( 1, 1 );

        }
        else
        {

            groupName = dsp.getFile().getFilePath( 2, 2 );
        }

        return groupName;
    }

    /**
     * This routine returns the selected shape for this node kind.
     *
     * @return string containing the shape
     */
    public String getShape()
    {

        return CTProperties.getNodeShape( myCallNode.getNodeKind() );

        // shape of the node depends on the kind and user settings

    }

    /***********************************************************************************************
     * * Methods for SubNodes * *
     **********************************************************************************************/

    /**
     * This routine adds to a group node a node that is in its group.
     *
     * @param node is a node for this group node
     */
    public void addSubNode( CTNode node )
    {

        if ( mySubNodes == null )
        {

            mySubNodes = new Vector();

        }

        mySubNodes.add( node );

    } // add SubNode

    /**
     * This routine returns for a group node the number of nodes
     * that are grouped by this node.
     *
     * @return the number of nodes grouped by this node
     */
    public int getNoSubNodes()
    {

        int n = 0;

        if ( mySubNodes != null )
        {

            n = mySubNodes.size();

        }

        return n;

    }

    /**
     * Access to the nodes of this group node.
     *
     * @param i is the index of subnode
     * @return the i-th subnode
     */
    public CTNode getSubNode( int i )
    {

        return ( CTNode ) mySubNodes.get( i );

    } // getSubNode

    /**
     * Unfold this node.
     *
     */
    public void unfold()
    {

        // do not unfold this node if there are no subnodes

        if ( mySubNodes == null )
        {

            logger.error( "Node " + getName() + " cannot be unfolded" );
            return;
        }

        myDoUnfold = true;
    }

    /**
     * Fold back this node by setting the unfold flag to false.
     * This will also be done for the parent if available.
     *
     * @return true if the node or parent was unfolded before
     */
    public boolean fold()
    {

        boolean done = myDoUnfold;

        myDoUnfold = false;

        if ( myParent != null )
        {

            done = done || myParent.fold();
        }

        return done;
    }

    /**
     * Getter routine for the parent of a CT node.
     *
     * @return the parent of this CT node (maybe null if not grouped)
     */
    public CTNode getParent()
    {

        return myParent;
    }

    /**
     * Get the group node to which this node belongs
     * if this group node is not unfolded.
     *
     * @return the parent CT node
     */
    public CTNode getFoldedParent()
    {

        // by default: parent is the node itself

        CTNode parent = this;

        if ( myParent != null )
        {

            // okay, I have a parent

            if ( !myParent.myDoUnfold )
            {

                // and it is really folded

                parent = myParent;
            }
        }

        return parent;
    }

    // try to find existing Edge

    /**
     * This routine returns for this node the node that should be printed.
     * For level = 0 its me, for any higher level it is the corresponding
     * group node to which the node belongs.
     *
     * @param level is the number of groups that are used.
     * @return the corresponding group node
     */
    public CTNode getPrintedNode( int level )
    {

        // in most cases its myself who is printed

        CTNode printedNode = this;

        // may be its not me to be printed
        // if we want to be at a higher level

        if ( level > 0 )
        {

            // okay, we want to look at a higher level

            if ( myParent != null )
            {

                // so I am really grouped in a set of nodes

                if ( !myParent.myDoUnfold )
                {

                    // and my parent node is not unfoled

                    printedNode = myParent.getPrintedNode( level - 1 );


                }

            }

        }

        return printedNode;

    } // getPrintedNode

    /**
     * Writing this node in a buffered write for the dot file.
     *
     * @param buff is the writer
     * @throws IOException in case of any IO problem
     */
    public void writeNode( BufferedWriter buff ) throws IOException
    {

        setColor(); // this is always done
        setLabel(); // this is always done

        if ( !isVisible() )
        {

            return;
        }

        String name = APO + getPrintIdent() + APO;

        if ( myDoUnfold )
        {

            // so we will print all subNodes
            // give it a try to a subgraph

            name = APO + "cluster_" + getName() + APO;

            buff.write( " subgraph " + name + " {" );
            buff.newLine();
            buff.write( "   graph [label=" + APO + getName() + APO + "]; " );
            buff.newLine();

            for ( int i = 0; i < getNoSubNodes(); i++ )
            {

                getSubNode( i ).writeNode( buff );

            }

            buff.write( "}" );
            buff.newLine();

        }
        else
        {

            buff.write( "      " + name );

            String label = APO + getLabel() + APO;

            buff.write( " [label=" + label );

            String color = getColor();

            buff.write( ", color=" + color );

            String shape = getShape();

            if ( shape.length() > 0 )
            {

                buff.write( ", shape=" + shape );
            }

            buff.write( "];" );

            buff.newLine();

        }

    } // writeNode

    /**
     * This routine writes the counter values into a buffered writer.
     *
     * @param buff is the buffer to which we write the counter values
     * @throws IOException in case of IO exception
     */
    public void writeCounterVals( BufferedWriter buff ) throws IOException
    {

        int noCounters = myInclusiveCounter.length;

        for ( int i = 0; i < noCounters; i++ )
        {

            buff.write( myInclusiveCounter[i] + " " );
        }


        for ( int i = 0; i < noCounters; i++ )
        {

            buff.write( "" + myExclusiveCounter[i] );

            if ( i + 1 < noCounters )
            {

                buff.write( " " );
            }
        }

        buff.newLine();

    } // writeCounterVals

    /**
     * Action for this node if it should be shown.
     *
     */
    public void showNode()
    {

        logger.info( "show CT Node " + getName() );

        RegionDescriptor myRegion = myCallNode.getRegion();

        if ( myRegion == null )
        {

            return;
        }

        logger.info( "CT node stands for region " + myRegion.getName() );

        // okay, is a node with a region so we can open it

        FileTable.showFile( myRegion.getFile(), myRegion.getFirstLine(), myRegion.getLastLine() );

    }

    /**
     * This routine returns an info string for this node.
     *
     * @return String containing relevant node information
     */
    public String infoNode()
    {

        final String ident = "  ";

        StringBuffer info = new StringBuffer();

        info.append( "Info of CT node <" + getName() );
        info.append( ":" + myCallNode.getNodeKindString() + ">" );
        info.append( "\n" );

        info.append( "  Identification = " + myIdent + "\n" );

        RegionDescriptor myRegion = myCallNode.getRegion();

        if ( myRegion != null )
        {

            myRegion.appendInfo( info, ident );

        }

        // now we print all the counter information available

        int numberCounters = CounterData.numberEvents();

        long value;

        double dval;

        String pVal;

        for ( int indexCounter = 0; indexCounter < numberCounters; indexCounter++ )
        {

            info.append( ident );
            info.append( CounterData.getEventName( indexCounter ) );
            value = myInclusiveCounter[indexCounter];
            dval  = CounterData.getRelativeValue( indexCounter, value );
            pVal  = Scaling.getValueString( dval, "m" );
            info.append( ": " + value );
            info.append( " (incl, " + pVal + ")" );
            value = myExclusiveCounter[indexCounter];
            dval  = CounterData.getRelativeValue( indexCounter, value );
            pVal  = Scaling.getValueString( dval, "m" );
            info.append( ", " + value );
            info.append( " (excl, " + pVal + ")\n" );
        }

        numberCounters = CounterData.numberMetrics();

        for ( int counter = 0; counter < numberCounters; counter++ )
        {

            info.append( ident );
            info.append( CounterData.getMetricName( counter ) );
            pVal = CounterData.getMetricValueString( counter, myInclusiveCounter );
            info.append( ": " + pVal + " (incl), " );
            pVal = CounterData.getMetricValueString( counter, myExclusiveCounter );
            info.append( pVal + " (excl)\n" );

        }

        if ( mySubNodes != null )
        {

            info.append( "has " + getNoSubNodes() + " subnodes\n" );

            for ( int i = 0; i < getNoSubNodes(); i++ )
            {

                info.append( getSubNode( i ).infoNode() );
            }
        }

        return info.toString();

    } // infoNode

    /**
     * Get the region for which a CT node can stand.
     *
     * @return Returns the region of this node.
     */
    public RegionDescriptor getRegion()
    {

        return myCallNode.getRegion();
    }

    /**
     * This routine is used to update the distance value of this node. If this node
     * has not been visited yet, the distance value is taken. Otherwise the new
     * value is only taken, if the absolute distance becomes shorter.
     *
     * @param distance is a known distance value for this node
     */
    public void visitWithDistance( int distance )
    {

        if ( myVisited )
        {

            // I have already a distance so we check for a better value

            if ( ( distance > 0 ) && ( distance < myDistance ) )
            {

                myDistance = distance;
            }

            if ( ( distance < 0 ) && ( distance > myDistance ) )
            {

                myDistance = distance;
            }

        }
        else
        {

            // this is my first value for distance

            myVisited = true;
            myDistance = distance;

        }

    }

    /**
     * This routine marks this node is the distance of the node
     * is in a certain range. -backward <= distance <= forward)
     *
     * @param forward is the maximal distance forward
     * @param backward is the maximal distance backward
     */
    public void markDistance( int forward, int backward )
    {

        // This node has only a valid distance if it has been visited.

        if ( myVisited )
        {

            if ( ( myDistance >= 0 ) && ( myDistance <= forward ) )
            {

                myMarked = true;
            }

            if ( ( myDistance <= 0 ) && ( myDistance >= ( -backward ) ) )
            {

                myMarked = true;
            }
        }
    }

    /**
     * Getter routine for the call node.
     *
     * @return the call node for which the CT node stands
     */
    public CallNode getCallNode()
    {

        return myCallNode;
    }


} // CTNode
