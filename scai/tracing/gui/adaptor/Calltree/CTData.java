/*
 * CTData.java
 *
 * Data of a call tree
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
import java.io.FileWriter;
import java.io.IOException;
import java.util.regex.Pattern;

import java.awt.Color;
import java.awt.Font;

import javax.swing.JFrame;

import org.apache.log4j.Logger;

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
import org.jfree.data.category.DefaultCategoryDataset;
import org.jfree.ui.RefineryUtilities;
import org.jfree.util.TableOrder;

import adaptor.General.ArraySelection;
import adaptor.General.CallNode;
import adaptor.General.FileDescriptor;
import adaptor.General.FileTable;
import adaptor.General.RIFData;
import adaptor.General.RegionDescriptor;

/**
 * The class CTData is a class containing all the data to generate a call tree.
 * CTData contains the nodes and the edges.
 *
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */
public class CTData
{

    /**
     * This string is used to seperate the two items of the print id.
     */
    private static final String ID_SEPERATOR = "_";

    /**
     * initial size of the array for the CT nodes.
     */
    private static final int MAX_NODES = 30;

    /**
     * initial size of the array for the CT edges.
     */
    private static final int MAX_EDGES = 30;

    /**
     * Size of the hash table for nodes and edges.
     */
    private static final int HASH_TABLE_SIZE = 311;

    /**
     * Logger for this class.
     */
    private static Logger logger = Logger.getLogger( CTData.class );

    /**
     * Identification to distinguish different CTData sets.
     */
    protected int myIdentification = 0;

    /**
     * CTNode contains all nodes of the calltree. This array is likely
     * to have gaps for certain entries. The position in the node
     * array corresponds to the node identification.
     */
    private CTNode[] myNodes = new CTNode[MAX_NODES];

    /**
     * This is the number of edges currently in use.
     */
    private int noEdges = 0;

    /**
     * This is the array of all edges of this CT data.
     */
    private CTEdge[] myEdges = new CTEdge[MAX_EDGES];

    /**
     * Nodes can be accessed via a hash table to have very fast
     * access to a node by the given name.
     */
    private MyHashTable myCTNodeHashTable = new MyHashTable( HASH_TABLE_SIZE );

    /**
     * Edges can be accessed via a hash table to have very fast access
     * to an edge by the names of source and target node.
     */
    private MyHashTable myCTEdgeHashTable = new MyHashTable( HASH_TABLE_SIZE );

    /**
     * The printed id of a node is the combination of the CT data identification
     * and of the node id.
     *
     * @param nodeId is the unique node id of this CT data
     * @return a unique id among all nodes
     */
    private String makePrintIdent( int nodeId )
    {

        return myIdentification + ID_SEPERATOR + nodeId;
    }

    /**
     * This routine returns for a the unique calltree node
     * identification (printIdent) the corresponding CT node.
     *
     * @param printIdent is the unique id of the node
     * @return full CT node with the given name (may be null)
     */
    public CTNode getNode( String printIdent )
    {

        String[] items = printIdent.split( ID_SEPERATOR );

        if ( items.length < 2 )
        {

            return null;
        }

        int pos = Integer.parseInt( items[0] );

        if ( pos != myIdentification )
        {

            return null;
        }

        pos = Integer.parseInt( items[1] );

        return myNodes[pos];

    }

    /**
     * This routine returns the CT Node belonging to a given
     * region by its descriptor. A very fast access is assumed.
     * This routine does not create a new CT node.
     *
     * @param dsp is the region descriptor for which the CT node is searched
     * @return the CT node for the region
     */
    public CTNode getNode( RegionDescriptor dsp )
    {

        CTNode node = null;

        // we have a quick access to the node via the region id

        int id = dsp.getRegionId();

        if ( ( id > 0 ) && ( id <= myNodes.length ) )
        {

            node = myNodes[id - 1];
        }

        return node;
    }

    /**
     * CT nodes of CT data can be accessed by the position in the
     * array that is also the node identification.
     *
     * @param pos is the position
     * @return is the CT node at the given position (can be null)
     */
    public CTNode getNode( int pos )
    {

        if ( ( 0 <= pos ) && ( pos < myNodes.length ) )
        {

            return myNodes[pos];
        }

        // illegal position >= pos can be possible

        return null;
    }

    /**
     * This routine finds for a given node the index position.
     *
     * @param node is the Calltree node
     * @return the index position (is less zero if not found)
     */
    public int getNodeIndex( CTNode node )
    {

        // we have the node index in the print ident

        String[] items = node.getPrintIdent().split( ID_SEPERATOR );

        if ( items.length < 2 )
        {

            return ArraySelection.NO_SELECTION;
        }

        return Integer.parseInt( items[1] );

    }

    /**
     * The routine findNode tries to find a calltree node whose
     * call node equals to the given input node.
     *
     * @param node is a call node that is searched in the CT data
     * @return a CT node where the call node equals the input node (null otherwise).
     */
    public CTNode findNode( CallNode node )
    {

        CTNode cnode = null;

        RegionDescriptor dsp = node.getRegion();

        if ( dsp != null )
        {

            // with a region we have the region id that gives
            // direct access to the CT node

            cnode = getNode( dsp );

            if ( cnode == null )
            {

                // this can happen when we union a data set to an empty set

                dsp  = null;

            }
            else if ( !cnode.getCallNode().equals( node ) )
            {

                // be careful: it might happen that region identifiers are different

                logger.warn( "descriptors of nodes do not match" );

                // descriptor of the input node was nothing worth

                dsp   = null;
                cnode = null;
            }

        }

        if ( dsp == null )
        {

            // no descriptor, so we find the node by the hash routine

            int pos = myCTNodeHashTable.queryEntry( node );

            if ( pos >= 0 )
            {

                cnode = myNodes[pos];
            }

        }

        return cnode;

    }


    /**
     * Routine returns the maximal node identification.
     *
     * @return maximal node identification in use
     */
    public int getNumberNodes()
    {

        return myNodes.length;
    }

    /**
     * Asking for the number of edges in the CT data.
     *
     * @return the exact number of edges of the CT data
     */
    public int getNumberEdges()
    {

        return noEdges;
    }

    /**
     * Edges of the CT data can be accessed by the position.
     *
     * @param pos with 0 <= pos < getNumberEdges()
     * @return the edge at position pos
     */
    public CTEdge getEdge( int pos )
    {

        return myEdges[pos];
    }
    /**
     * This node returns the edge between two nodes. If the
     * edge is not available, a new one will be created.
     *
     * @param node1 is the source node
     * @param node2 is the target node
     * @return an edge between the nodes belonging to the CT data
     */
    public CTEdge getEdge( CTNode node1, CTNode node2 )
    {

        CTEdge edge;

        String name = node1.getName() + "-" + node2.getName();

        int pos = myCTEdgeHashTable.queryEntry( name );

        if ( pos < 0 )
        {

            // we have to add a new edge

            edge = new CTEdge( node1, node2 );

            if ( noEdges >= myEdges.length )
            {

                extendEdges();

            }

            myCTEdgeHashTable.addEntry( name, noEdges );

            myEdges[noEdges++] = edge;

            return edge;

        } // getEdge

        // pos >= 0 -> Edge was already defined

        return myEdges[pos];

    } // getEdge

    /**
     * Dynamic extension of the array of nodes.
     *
     * @param nr is the minimal size that myNodes must have
     */
    private void extendNodes( int nr )
    {

        int sizeNodes = myNodes.length;

        while ( nr >= sizeNodes )
        {

            // double number of Nodes

            CTNode[] newNodes = new CTNode[2 * sizeNodes];

            for ( int i = 0; i < sizeNodes; i++ )
            {
                newNodes[i] = myNodes[i];
            }

            sizeNodes = 2 * sizeNodes;
            myNodes = newNodes;

        } // while not enough size

    } // extend_nodes


    /**
     * This routine adds a new node at a given position.
     *
     * @param pos is the position where to add
     * @param node is the CT node to add
     */
    public void addCTNode( int pos, CTNode node )
    {

        node.setIdent( makePrintIdent( pos ) );
        extendNodes( pos );
        myNodes[pos] = node;
        myCTNodeHashTable.addEntry( node.getCallNode(), pos );
        logger.info( "added CTNode " + node.getName() + " in hash table at pos " + pos );

    }

    /**
     * This routine doubles the size of the arrays with edges.
     *
     */
    private void extendEdges()
    {

        int sizeEdges = myEdges.length;

        logger.info( "extend array of edges, double current size of " + sizeEdges );

        // double number of Edges

        CTEdge[] newEdges = new CTEdge[2 * sizeEdges];

        for ( int i = 0; i < sizeEdges; i++ )
        {

            newEdges[i] = myEdges[i];
            newEdges[i + sizeEdges] = null;
        }

        myEdges = newEdges;

    }

    /**
     * This routine returns all CT nodes n with n.getName()
     * equal to a given input name.
     *
     * @param name is the input name for which we search nodes
     * @return an array with all matching CT nodes (can be null)
     */
    CTNode[] findNodes( String name )
    {

        // we make a dummy node that equals with all other call nodes of the same name

        CallNode dummy = new CallNode( CallNode.ANY_NODE, name );

        int[] pos = myCTNodeHashTable.queryEntries( dummy );

        if ( pos == null )
        {

            // make an efficient solution for this case

            return null;

        }

        if ( pos.length == 0 )
        {

            return null;

        }

        // the CT nodes we get by their positions

        CTNode[] nodes = new CTNode[pos.length];

        for ( int i = 0; i < nodes.length; i++ )
        {

            nodes[i] = myNodes[pos[i]];
        }

        return nodes;
    }

    /**
     * This routine is used to reset for all nodes certain attributes.
     * The node is not marked, not visited and the distance is set to zero.
     */
    void resetAllNodes()
    {

        for ( int i = 0; i < myNodes.length; i++ )
        {

            CTNode node = myNodes[i];

            if ( node != null )
            {

                node.reset();

            }

        } // for all nodes

    } // setAllNodes

    /**
     * Update the distance for two nodes between an edge
     * exists.
     *
     * @param node1 is the source node
     * @param node2 is the target node
     */
    private void updateDistance( CTNode node1, CTNode node2 )
    {

        if ( node1.isVisited() && node1.getDistance() >= 0 )
        {

            node2.visitWithDistance( node1.getDistance() + 1 );

        }

        if ( node2.isVisited() && node2.getDistance() <= 0 )
        {

            node1.visitWithDistance( node2.getDistance() - 1 );

        }

    }

    /**
     * This routine is used to mark all nodes in a certain
     * environment of a set of nodes. It starts with all nodes
     * that have the distance 0 and updates then the distance
     * for neighbored nodes.
     *
     * @param sourceDist is the distance from source
     * @param targetDist is the maximal distance to target
     */
    void markEnvironment( int sourceDist, int targetDist )
    {

        int depth = targetDist;

        if ( sourceDist > depth )
        {

            depth = sourceDist;
        }

        // Step 1: set visited = true && distance = 0 for all marked nodes

        for ( int i = 0; i < myNodes.length; i++ )
        {

            CTNode node = myNodes[i];

            if ( node != null )
            {

                if ( node.isMarked() )
                {

                    node.visitWithDistance( 0 );
                }
            }

        } // for all nodes

        for ( int d = 0; d < depth; d++ )
        {

            // loop at least depth times

            for ( int i = 0; i < noEdges; i++ )
            {


                CTEdge edge = myEdges[i];

                CTNode node1 = edge.getSource();
                CTNode node2 = edge.getTarget();

                updateDistance( node1, node2 );

            }
        }

        // now we mark all nodes

        for ( int i = 0; i < myNodes.length; i++ )
        {

            CTNode node = myNodes[i];

            if ( node != null )
            {

                myNodes[i].markDistance( targetDist, sourceDist );
            }

        } // for all nodes

    } // markEnvironment

    /**
     * This routine restricts marked nodes that satisfy a search pattern.
     *
     * @param pattern is the search pattern for the nodes
     */
    void searchNodes( String pattern )
    {

        Pattern searchPattern = Pattern.compile( pattern, Pattern.CASE_INSENSITIVE );

        for ( int i = 0; i < myNodes.length; i++ )
        {

            CTNode node = myNodes[i];

            if ( node != null )
            {

                if ( node.isMarked() )
                {

                    String name = node.getName();

                    node.setMarked( searchPattern.matcher( name ).matches() );

                }
            }
        }

    } // searchNodes

    /**
     * This routine marks the selected nodes and all other nodes in a
     * certain neighborhood (calling and called depth).
     *
     * @param selNodes are the nodes that haven been selected
     * @param sourceDist is the positive distance from selected node
     * @param targetDist is the negative distance from selected node
     */
    void buildSubGraph( CTNode[] selNodes, int sourceDist, int targetDist )
    {

        resetAllNodes();

        // If no nodes are specified we will mark none

        if ( selNodes != null )
        {

            for ( int i = 0; i < selNodes.length; i++ )
            {

                selNodes[i].setMarked( true );
                selNodes[i].setVisited( true );

            }
        }

        markEnvironment( sourceDist, targetDist );

    } // buildSubGraph

    /**
     * Routine marks the root or leaf nodes of the CTData. A root is a node
     * with no incoming edge, a leaf is a node without outgoing edge.
     *
     * @param rootFlag if set true will mark only root nodes
     * @param leafFlag if set true will mark only leaf nodes
     */
    public void markRootOrLeaf( boolean rootFlag, boolean leafFlag )
    {

        // mark all nodes

        for ( int i = 0; i < myNodes.length; i++ )
        {

            CTNode node = myNodes[i];

            if ( node != null )
            {

                node.setMarked( true );
            }

        }

        // unmark all nodes that are source/target of an edge

        for ( int i = 0; i < noEdges; i++ )
        {

            CTEdge edge = myEdges[i];

            // unmark source node

            if ( leafFlag )
            {

                edge.getSource().setMarked( false );
            }

            if ( rootFlag )
            {

                edge.getTarget().setMarked( false );
            }

        } // for all myEdges

        return;

    } // mark RootOrLeaf

    /**
     * update the number of calls for a non-visited node
     * updates also the number of calls for all incoming edges
     * and sets now correct values for inclusive/exclusive/call costs.
     *
     * @param node is the node that is now visited
     */
    void calcFictiveRuntimeCalls( CTNode node )
    {

        if ( node.isVisited() )
        {

            return;
        }


        node.setVisited( true );
        node.setMarked( true ); // marked as long as on call stack

        node.setRuntimeCalls( 0 );

        // look for incoming edges (calls of this node N)
        // important: all edges will be considered only once

        for ( int i = 0; i < noEdges; i++ )
        {

            CTEdge edge = myEdges[i];

            // only incoming edges

            if ( edge.getTarget() != node )
            {

                continue;
            }


            // avoid recursive calls

            if ( edge.getSource().isMarked() )
            {

                continue;
            }


            // make sure that we know the number of runtime calls
            // for the source node of this edge

            calcFictiveRuntimeCalls( edge.getSource() );

            // node N is then called as follows

            long myCalls = edge.getSource().getRuntimeCalls() * edge.getSourceCalls();

            edge.addNoRunCalls( myCalls );

            node.setRuntimeCalls( node.getRuntimeCalls() + myCalls );

        } // for all edges

        if ( node.getRuntimeCalls() == 0 )
        {

            node.setRuntimeCalls( 1 );

        }

        node.setMarked( false ); // marked as long as on call stack

    } // calcFictiveRuntimeCalls

    /**
     * Calculate fictive number of runtime calls.
     * e.g.: A calls B (NoSourceCalls = 3) becomes (NoCalls = 6)
     *       if A itself is called two times
     */
    void calcFictiveRuntimeCalls()
    {

        logger.info( "calc runtime calls for all nodes/edges" );

        // all nodes are marked as not visited yet

        resetAllNodes();

        // actualize number of runtime calls for all nodes
        // (bottom up traversing)

        for ( int i = 0; i < myNodes.length; i++ )
        {

            CTNode node = myNodes[i];

            if ( node != null )
            {

                calcFictiveRuntimeCalls( node );
            }

        }

    }

    /**
     * routine that computes the inclusive costs for a node
     * respects the cycles in the callgraph.
     *
     * @param node is the node that is now visited
     */
    void calcInclCosts( CTNode node )
    {

        if ( node.isVisited() )
        {
            return;
        }


        node.setVisited( true );
        node.setMarked( true ); // marked as long as on call stack

        node.setRuntimeCalls( 0 );

        node.setInclusiveCosts( node.getExclusiveCosts() );

        // look for incoming edges (calls of this node N)
        // important: all edges will be considered only once

        for ( int i = 0; i < noEdges; i++ )
        {

            CTEdge edge = myEdges[i];

            if ( edge.getSource() == node )
            {

                if ( edge.getTarget().isMarked() )
                {

                    continue;
                }

                calcInclCosts( edge.getTarget() );

                node.addInclusiveCosts( edge.getCallCosts() );

            }

        } // for all incoming edges

        node.setMarked( false ); // marked as long as on call stack

    } // calcInclCosts

    /**
     * This routine calculates for all CT nodes the inclusive costs.
     *
     */
    void calcInclCosts()
    {

        logger.info( "calculate inclusive costs" );

        // all nodes are marked as not visited yet

        resetAllNodes();

        // actualize number of calls for all nodes
        // (top down traversing)

        for ( int i = 0; i < myNodes.length; i++ )
        {

            CTNode node = myNodes[i];

            if ( node != null )
            {

                calcInclCosts( node );
            }


        }

    } // calcInclCosts

    /**
     * calculate fictive runtime costs for a node and its
     * outgoing edges.
     * This routine is called recursively (depth first search).
     *
     * @param node is the calltreeNode currently visited
     */
    void calcFictiveCosts( CTNode node )
    {

        if ( node.isVisited() )
        {

            return;
        }

        node.setVisited( true );
        node.setMarked( true ); // needed to avoid call cycles

        long costs = node.getRuntimeCalls();

        // set inclusive and exclusive costs, one counter only

        long [] allCosts = { costs };

        node.setInclusiveCosts( allCosts );
        node.setExclusiveCosts( allCosts );

        for ( int i = 0; i < noEdges; i++ )
        {

            CTEdge edge = myEdges[i];

            // consider only outgoing edges from this node

            if ( edge.getSource() != node )
            {

                continue;
            }

            // if E.Source is on stack we avoid recursion

            CTNode calledNode = edge.getTarget();

            if ( calledNode.isMarked() )
            {

                continue;
            }


            calcFictiveCosts( calledNode );

            // get the inclusive costs of the called routine for one call

            long[] inclVals = calledNode.getInclusiveCosts();

            long callCosts = inclVals[0] / calledNode.getRuntimeCalls();

            callCosts *= edge.getRuntimeCalls();

            long[] edgeCallCosts = edge.getCallCosts();

            edgeCallCosts[0] = callCosts;

            long[] allCallCosts = { callCosts };

            node.addInclusiveCosts( allCallCosts ) ;

        } // for all outgoing edges

        node.setMarked( false );

    } // calcFictiveCosts

    /**
     * TODO: brandes: Enter comment!
     *
     */
    void calcFictiveCosts()
    {

        logger.info( "calcFictiveCosts: started" );

        // all nodes are marked as not visited yet

        resetAllNodes();

        for ( int i = 0; i < myNodes.length; i++ )
        {

            CTNode node = myNodes[i];

            if ( node != null )
            {

                // we only visit the root nodes

                if ( node.isSourceRoot() )
                {

                    calcFictiveCosts( node );
                }

            }

        } // for all nodes

        logger.info( "calcFictiveCosts: ready" );

    } // calcFictiveCosts

    /**
     * TODO: brandes: Enter comment!
     *
     * @return array with all total values
     */
    long[] calcTotalCosts()
    {

        int noCounters = CounterData.numberEvents();

        long[] totalCosts = new long[noCounters];

        for ( int i = 0; i < myNodes.length; i++ )
        {

            CTNode node = myNodes[i];

            if ( node != null )
            {

                long [] costs = node.getExclusiveCosts();

                for ( int j = 0; j < noCounters; j++ )
                {

                    totalCosts[j] += costs[j];
                }

            }

        } // for all nodes

        return totalCosts;

    } // calcTotalCosts

    /**
     * This routine resets selected nodes and distances.
     *
     */
    void resetSubGraph()
    {

        int i;

        int numberNodes = myNodes.length;

        for ( i = 0; i < numberNodes; i++ )
        {

            CTNode node = myNodes[i];

            if ( node != null )
            {

                node.setMarked( true );
                node.setVisited( false );
            }

        } // for all nodes

    } // resetSubGraph

    /**
     * This routine sets for all nodes the parent.
     *
     */
    void setParent()
    {

        int noNodes = myNodes.length;

        for ( int i = 0; i < noNodes; i++ )
        {

            CTNode node = myNodes[i];

            if ( node != null )
            {

                int noSubNodes = node.getNoSubNodes();

                if ( noSubNodes == 0 )
                {

                    node.setParent( null );

                }
                else
                {

                    for ( int k = 0; k < noSubNodes; k++ )
                    {

                        CTNode subNode = node.getSubNode( k );
                        subNode.setParent( node );
                    }
                }

            }

        } // for all nodes

    } // setParent

    /**
     * This routine resets unfolded nodes.
     *
     * @return true if there was anz unfolded node.
     */
    boolean resetUnfolding()
    {

        int i;

        boolean done = false;

        int counter = 0;

        int noNodes = myNodes.length;

        for ( i = 0; i < noNodes; i++ )
        {

            CTNode node = myNodes[i];

            if ( node != null )
            {

                if ( node.isUnfolded() )
                {

                    done = true;
                    counter += 1;
                    node.setUnfolded( false );
                }
            }

        } // for all nodes

        logger.info( counter + " nodes of CT were unfolded" );

        return done;

    } // resetUnfolding

    /**
     * reset unfold, hide attributes.
     */
    void resetView()
    {

        int numberNodes = myNodes.length;

        for ( int i = 0; i < numberNodes; i++ )
        {

            CTNode node = myNodes[i];

            if ( node != null )
            {

                node.setUnfolded( false );
                node.setMarked( true );
            }

        } // for all nodes

    } // resetUnfolding

    /**
     * This routine prints the calltree data into a dot file.
     *
     * @param filename is the name of the output file
     */
    void printGraph( String filename )
    {

        int i;

        try
        {

            FileWriter file = new FileWriter( filename );
            BufferedWriter buff = new BufferedWriter( file );

            buff.write( "digraph \"Test Graph\" {" );
            buff.newLine();

            buff.write( " node [ shape = " + CTProperties.getNodeDefaultShape() + ", style=filled ];" );

            buff.write( " edge [ style = \"setlinewidth(" + CTProperties.getEdgeWidth() + ".0)\" ];" );
            buff.newLine();

            // bgcolor=cyan has no importance for GRAPPA
            // fontsize=20 applies to Title of the drawing

            buff.write( " graph [fontsize=20,label=\"" + CTProperties.getTitle() + "\"]; " );
            buff.newLine();

            for ( i = 0; i < myNodes.length; i++ )
            {

                CTNode node = myNodes[i];

                if ( node != null )
                {

                    node.writeNode( buff );
                }


            } // for all Regions

            for ( i = 0; i < noEdges; i++ )
            {

                CTEdge edge = myEdges[i];

                edge.writeEdge( buff, 0 );

            }

            buff.write( "}" );
            buff.newLine();
            buff.close();

        }
        catch ( IOException e )
        {

            logger.error( "IO Exception  -- " + e.toString() );
        }

        logger.info( ">>>>> output done" );

    }

    /**
     * Routine delivers information about this Calltree data.
     *
     * @return some Info about this CTData
     */
    String getInfo()
    {

        int i;
        int noNodes = 0;
        int visibleNodes = 0;
        int visibleEdges = 0;

        for ( i = 0; i < myNodes.length; i++ )
        {

            CTNode node = myNodes[i];

            if ( node != null )
            {

                noNodes++;

                if ( node.isVisible() )
                {

                    visibleNodes++;
                }
            }

        } // for all Nodes

        for ( i = 0; i < noEdges; i++ )
        {

            CTEdge edge = myEdges[i];

            if ( edge.isVisible() )
            {

                visibleEdges++;
            }
        }

        String info = "CT has " + noNodes + " Nodes (" + "visible " + visibleNodes + ")\n";

        info += "   and " + noEdges + " Edges (" + "visible " + visibleEdges + ")\n";

        return info;

    } // getInfo


    /**
     * Simple help function to get the number of region nodes.
     *
     * @return number of region nodes
     */
    public int countNodeRegions()
    {

        int counterRegions = 0;

        for ( int i = 0; i < myNodes.length; i++ )
        {

            CTNode node = myNodes[i];

            if ( node != null )
            {

                if ( node.getRegion() != null )
                {

                    counterRegions++;
                }
            }
        }

        return counterRegions;
    }


    /**
     * Export the calltree data as a RIF file.
     *
     * @param filename
     *            is the name of the output file.
     */
    public void exportRIF( String filename )
    {

        // create RIF data from this CT data

        int noRegions = countNodeRegions();

        logger.info( "exportRIF: " + noRegions + " regions" );

        RegionDescriptor[] myRegions = new RegionDescriptor[noRegions];

        int index = 0;

        for ( int i = 0; i < myNodes.length; i++ )
        {

            CTNode node = myNodes[i];

            if ( node != null )
            {

                RegionDescriptor region = node.getRegion();

                if ( region != null )
                {

                    myRegions[index++] = region;

                    // set all flags to overwrite flags from previous export

                    region.setProfEnabled( node.isVisible() );
                    region.setNestProfiling( 0 ); // do not enter nested regions

                }
            }

        } // for all Regions

        // By default profiling of called regions is disabled. But
        // we enable it if there is at least one edge that goes out
        // from the region node

        for ( int i = 0; i < noEdges; i++ )
        {

            CTEdge edge = myEdges[i];

            // we consider only visible edges

            CTNode source = edge.getSource();
            CTNode target = edge.getTarget();

            if ( target.isVisible() && source.isVisible() )
            {

                RegionDescriptor region = source.getRegion();

                if ( region != null )
                {

                    region.setNestProfiling( -1 );
                }
            }

        }

        // ToDo: generate new file data

        FileTable newFileData = new FileTable();

        for ( int i = 0; i < myRegions.length; i++ )
        {

            FileDescriptor f = myRegions[i].getFile();

            newFileData.defineFile( f.getFileId(), f );

        }

        RIFData myRIFData = new RIFData( newFileData.getFiles(), myRegions );

        myRIFData.writeRIF( filename );

    }

    /**
     * This routine opens a new chart frame with the exlusive values of the selected event.
     *
     * @param barFlag
     *            if true makes a bar chart otherwise a pie chart
     */
    void makeChart( boolean barFlag )
    {

        final DefaultCategoryDataset dataset = new DefaultCategoryDataset();

        // String CounterName = CTProperties.getSelectedCounter ();

        String nameCounter = CTProperties.getCounterLegend();

        for ( int i = 0; i < myNodes.length; i++ )
        {

            CTNode node = myNodes[i];

            if ( node != null )
            {

                if ( node.isVisible( false ) )
                {

                    dataset.addValue( node.getExclusiveValue(), nameCounter, node.getName() );
                }

            }
        } // for all nodes

        Plot drawPlot;

        if ( barFlag )
        {

            CategoryAxis categoryAxis = new CategoryAxis( "" );
            ValueAxis valueAxis = new NumberAxis( "" );

            BarRenderer renderer = new BarRenderer();

            renderer.setBaseToolTipGenerator( new StandardCategoryToolTipGenerator() );

            CategoryPlot plot = new CategoryPlot( dataset, categoryAxis, valueAxis, renderer );

            plot.setOrientation( PlotOrientation.HORIZONTAL );

            drawPlot = plot;

        }
        else
        {

            final double interiorGap = 0.30;

            MultiplePiePlot plot = new MultiplePiePlot( dataset );

            plot.setDataExtractOrder( TableOrder.BY_ROW );
            plot.setBackgroundPaint( null );
            plot.setOutlineStroke( null );

            PieToolTipGenerator tooltipGenerator = new StandardPieToolTipGenerator();

            PiePlot p = ( PiePlot ) plot.getPieChart().getPlot();
            p.setToolTipGenerator( tooltipGenerator );

            p.setLabelGenerator( new StandardPieSectionLabelGenerator( "{0}" ) );
            p.setLabelFont( new Font( "SansSerif", Font.PLAIN, 8 ) );
            p.setInteriorGap( interiorGap );

            // use the following if no labels are wanted
            // p.setLabelGenerator (null);

            drawPlot = plot;
        }

        final JFreeChart chart = new JFreeChart( "", JFreeChart.DEFAULT_TITLE_FONT, drawPlot, true );

        // set the background color for the chart...

        chart.setBackgroundPaint( Color.lightGray );

        // add the chart to a panel...

        final ChartPanel chartPanel = new ChartPanel( chart );

        final int preferredWidth = 640;
        final int preferredHeight = 420;

        chartPanel.setPreferredSize( new java.awt.Dimension( preferredWidth, preferredHeight ) );

        JFrame chartFrame = new JFrame( "CalltreePM Chart" );

        chartFrame.setContentPane( chartPanel );
        chartFrame.pack();

        RefineryUtilities.centerFrameOnScreen( chartFrame );

        chartFrame.setVisible( true );

    }

    /**
     * This routine takes the costs from given data set to set
     * the costs of this CT data. For the inheritance the name
     * of the node is taken.
     *
     * @param inCT is the data from which we inherit costs
     */
    public void inheritCosts( CTData inCT )
    {

        int noNodes = inCT.myNodes.length;

        for ( int i = 0; i < noNodes; i++ )
        {

            CTNode ctNode = inCT.getNode( i );

            if ( ctNode == null )
            {

                continue;

            }

            // get the counterpart of ctNode in this CT data

            CTNode node1 = findNode( ctNode.getCallNode() );

            // inherit the self costs from the node

            if ( node1 != null )
            {

                long [] costs = ctNode.getExclusiveCosts();

                node1.addInclusiveCosts( costs );
                node1.addExclusiveCosts( costs );

            }

        } // for all nodes

        int noInheritedEdges = inCT.noEdges;

        for ( int i = 0; i < noInheritedEdges; i++ )
        {

            CTEdge edge = inCT.myEdges[i];

            CTNode node1 = edge.getSource();
            CTNode node2 = edge.getTarget();

            // find the counterparts of N1 and N2

            CTNode parent1 = findNode( node1.getCallNode() );
            CTNode parent2 = findNode( node2.getCallNode() );

            // find/insert the Edge for the supernodes

            CTEdge parentEdge = getEdge( parent1, parent2 );

            // add the costs / calls

            parentEdge.inheritCosts( edge );

        } // for all Edges

    } // inheritCosts

    /**
     * The routine resetCosts resets the cost counter values
     * for all nodes and edges.
     */
    void resetCosts()
    {

        int noNodes = myNodes.length;

        for ( int i = 0; i < noNodes; i++ )
        {

            CTNode node = myNodes[i];

            if ( node != null )
            {

                node.resetCosts();
            }

        } // for all nodes

        for ( int i = 0; i < noEdges; i++ )
        {

            CTEdge edge = myEdges[i];

            edge.resetCosts();

        } // for all Edges

    } // resetCosts

} // class CTData
