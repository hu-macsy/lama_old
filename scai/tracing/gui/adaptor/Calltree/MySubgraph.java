/*
 * MySubgraph.java
 *
 * Frame that displays the call graph.
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

import java.util.List;
import java.util.Vector;

import org.apache.log4j.Logger;

import att.grappa.Element;
import att.grappa.Grappa;
import att.grappa.GrappaConstants;
import att.grappa.Node;
import att.grappa.Subgraph;

/**
 * The class MySubgraph extends the class Subgraph with new methods for
 * selections.
 *
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */
public final class MySubgraph
{

    /**
     * This is an array with the names of all nodes of the
     * last selection.
     *
     */
    private static String[] lastSelection = null;

    /**
     * Logger for this class.
     */
    private static Logger logger = Logger.getLogger( MySubgraph.class );

    /**
     * Overwriting the default constructor.
     *
     */
    private MySubgraph()
    {

    }
    /**
     * This routine returns from a Subgraph the selected node and
     * returns the corresponding CT node.
     *
     * @param graph is the graph where one node should be selected
     * @return the CT node or null if no or multiple nodes are selected
     */
    static CTNode getSelectedNode( Subgraph graph )
    {

        if ( graph.currentSelection instanceof Element )
        {

            Element selElem = ( Element ) graph.currentSelection;

            logger.info( "selected element = " + selElem );

            if ( selElem.getType() == GrappaConstants.NODE )
            {

                Node node = ( Node ) selElem;
                logger.info( "NODE " + node.getName() + " clicked" );
                CTNode nodeCT = CalltreePM.getNode( node.getName() );

                return nodeCT;

            }

        }

        return null;

    } // getSelectedNode

    /**
     * This routine returns the selected nodes in a graph
     * presentation but the corresponding ones of the CT data.
     *
     * @param graph is the layout graph of GRAPPA
     * @return an array of CT nodes that are selected
     */
    static CTNode[] getSelectedNodes( Subgraph graph )
    {

        CTNode[] selectedCTNodes = null;

        if ( graph.currentSelection instanceof Element )
        {

            CTNode nodeCT = getSelectedNode( graph );

            if ( nodeCT != null )
            {

                selectedCTNodes = new CTNode[1];

                selectedCTNodes[0] = nodeCT;

                return selectedCTNodes;
            }

        }
        else if ( graph.currentSelection instanceof Vector )
        {

            List vec = ( Vector ) graph.currentSelection;

            List vecCT = new Vector();

            for ( int i = 0; i < vec.size(); i++ )
            {

                Element elem = ( Element ) vec.get( i );

                if ( elem.getType() == GrappaConstants.NODE )
                {

                    Node node = ( Node ) elem;

                    CTNode nodeCT = CalltreePM.getNode( node.getName() );

                    if ( nodeCT != null )
                    {

                        vecCT.add( nodeCT );
                    }

                }

            }

            int size = vecCT.size();

            selectedCTNodes = new CTNode[size];

            for ( int i = 0; i < size; i++ )
            {

                selectedCTNodes[i] = ( CTNode ) vecCT.get( i );

            }

        }

        return selectedCTNodes;

    } // getSelectedNodes

    /**
     * This routine takes the given CT nodes and makes them to
     * the last selected nodes.
     *
     * @param selectedNodes is the array with selected nodes
     */
    public static void fixSelectedNodes( CTNode[] selectedNodes )
    {

        int noNodes = selectedNodes.length;

        lastSelection = new String[noNodes];

        for ( int i = 0; i < noNodes; i++ )
        {

            lastSelection[i] = selectedNodes[i].getPrintIdent();

        }

    } // fixSelectedNodes

    /**
     * This routine sets in a subgraph the last selected nodes
     * to the selected nodes of this graph.
     *
     * @param theGraph is the graph where we set the latest selection
     */
    public static void setLastSelection( Subgraph theGraph )
    {

        if ( lastSelection == null )
        {

            logger.info( "setLastSelection: not done, was null" );
            return;

        }

        if ( theGraph.currentSelection != null )
        {

            logger.info( "setLastSelection: not done, G has already selection" );
            return;

        }

        int numberNodes = lastSelection.length;

        if ( numberNodes == 0 )
        {

            logger.info( "setLastSelection: not done, 0 elements" );
            return;

        }

        List selectionVector = new Vector();

        for ( int i = 0; i < numberNodes; i++ )
        {

            Node graphNode = theGraph.findNodeByName( lastSelection[i] );

            if ( graphNode == null )
            {

                logger.info( "setLastSelection: node " + lastSelection[i] + " not found" );

            }
            else
            {

                graphNode.highlight |= Grappa.SELECTION_MASK;

                selectionVector.add( graphNode );

            }

        } //

        logger.info( "setLastSelection: done, " + selectionVector.size() + " elements selected" );

        theGraph.currentSelection = selectionVector;

    } // setLastSelection

}
