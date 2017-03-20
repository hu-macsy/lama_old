/*
 * CTDataGroup.java
 *
 * Calltree build by grouping nodes and edges of CTData.
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

import adaptor.General.CallNode;

/**
 * CTDataGroup builds CTData by grouping the nodes.
 *
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */
public class CTDataGroup extends CTData
{

    /**
     * Logger for this class.
     */
    private static Logger logger = Logger.getLogger( CTDataGroup.class );

    /**
     * A CTData group is created by given CTData.
     * Different possibilities of grouping are supported by the
     * value group_kind.
     *
     * @param theData is the underlying calltree data
     * @param groupKind is the kind of group to build
     */
    CTDataGroup( CTData theData, int groupKind )
    {

        myIdentification = groupKind;

        logger.info( "collapsing CTData, group kind = " + groupKind );

        // we start with the default sizes

        int i;

        int noGroupNodes = 0;

        int noNodes = theData.getNumberNodes();

        logger.info( "Original CT has " + noNodes + " nodes" );

        for ( i = 0; i < noNodes; i++ )
        {

            CTNode node = theData.getNode( i );

            if ( node == null )
            {

                continue;
            }

            // the new node gets the name of the file

            String groupName = node.getGroupName( groupKind );

            CallNode cnode = new CallNode( CallNode.GROUP_NODE, groupName );

            CTNode groupNode = findNode( cnode );

            if ( groupNode == null )
            {

                groupNode = new CTNode( cnode );
                addCTNode( noGroupNodes, groupNode );
                noGroupNodes++;
            }

            // add the node to this CT Data

            groupNode.addSubNode( node );

            long [] costs = node.getExclusiveCosts();

            groupNode.addInclusiveCosts( costs );
            groupNode.addExclusiveCosts( costs );

            // we set the parent of the node

            node.setParent( groupNode );

        }

        logger.info( "collapsed graph has " + noGroupNodes + " nodes" );

        int nEdges = theData.getNumberEdges();

        logger.info( "Original CT has " + nEdges + " edges" );

        for ( i = 0; i < nEdges; i++ )
        {

            CTEdge edge = theData.getEdge( i );

            CTNode sourceNode = edge.getSource();
            CTNode targetNode = edge.getTarget();

            // find the supernodes of N1 and N2

            CTNode groupSourceNode = sourceNode.getParent();
            CTNode groupTargetNode = targetNode.getParent();

            // find/insert the Edge for the supernodes

            CTEdge groupEdge = getEdge( groupSourceNode, groupTargetNode );

            groupEdge.addSubEdge( edge );

            // add the costs / calls

            groupEdge.inheritCosts( edge );

        }

        logger.info( "the new graph has " + nEdges + " edges" );

    }

} // class CTDataGroup
