/*
 * CTDataUnion.java
 *
 * Union of CTData sets to a single data set
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

/**
 * CTDataUnion is CTData build by the union of many CTData sets.
 *
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */
public class CTDataUnion extends CTData
{

    /**
     * Logger for this class.
     */
    private static Logger logger = Logger.getLogger( CTDataUnion.class );

    /**
     * constructor that builds a new CT data set that is union of multiple CT data sets.
     *
     * @param dataSets is the array with all CT data sets
     */
    public CTDataUnion( CTData[] dataSets )
    {

        int nCT = dataSets.length;

        logger.info( "build union of CTData of " + nCT + " CTData sets" );

        int noAllNodes = 0;

        for ( int k = 0; k < nCT; k++ )
        {

            CTData data = dataSets[k];

            int noNodes = data.getNumberNodes();

            for ( int i = 0; i < noNodes; i++ )
            {

                CTNode node = data.getNode( i );

                if ( node == null )
                {

                    continue;
                }

                CTNode partner = findNode( node.getCallNode() );

                if ( partner == null )
                {

                    partner = new CTNode( node.getCallNode() );
                    addCTNode( noAllNodes, partner );
                    noAllNodes++;
                }

                // add the costs  to this CT Data

                partner.addCounters( node );

            } // for all nodes

            int noEdges = data.getNumberEdges();

            for ( int i = 0; i < noEdges; i++ )
            {

                CTEdge edge = data.getEdge( i );

                CTNode node1 = edge.getSource();
                CTNode node2 = edge.getTarget();

                // find the counterparts of N1 and N2

                CTNode partner1 = findNode( node1.getCallNode() );
                CTNode partner2 = findNode( node2.getCallNode() );

                // find/insert the Edge for the supernodes

                CTEdge partnerEdge = getEdge( partner1, partner2 );

                // add the costs / calls

                partnerEdge.addCosts( edge );

            } // for all Edges

        } // for all CT Files

        logger.info( "built union of CTData ready" );
        logger.info( "CT has " + noAllNodes + " nodes + " + getNumberEdges() + " Edges" );
    }

} // class CTDataUnion
