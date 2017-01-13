/*
 * CTDataRIF.java
 *
 * Calltree build from RIF file.
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

import adaptor.General.CallDescriptor;
import adaptor.General.CallNode;
import adaptor.General.RIFData;
import adaptor.General.RegionDescriptor;
import adaptor.General.RegionTable;

/**
 * CTData build from RIFData.
 *
 * @version $LastChangedRevision$
 * @author Dr. Thomas Brandes
 */

public class CTDataRIF extends CTData
{

    /**
     * Logger for this class.
     */
    private static Logger logger = Logger.getLogger( CTDataRIF.class );

    /**
     * Constructor for Calltree data from the RIF data (read from RIF file).
     *
     * @param theRIFData is the data of the RIF
     * @param fictiveCosts if true calculates fictive runtime calls and costs
     */
    public CTDataRIF( RIFData theRIFData, boolean fictiveCosts )
    {

        int i;

        RegionTable myRegionData = theRIFData.getRegionData();

        // we keep a counter for the nodes as we might create new ones

        int noNodes = myRegionData.noRegions();

        logger.info( "CTData from RIFData: " + noNodes + " Regions" );

        // define the fictive counter to make sure that arrays for
        // the counter values will already be allocated properly

        if ( fictiveCosts )
        {

            CounterData.defineEvent( 0, "FICTIVE", "fictive calculated costs" );
        }

        // make nodes from all regions

        for ( i = 0; i < noNodes; i++ )
        {

            RegionDescriptor region = myRegionData.getRegion( i );

            CTNode node = new CTNode( new CallNode( region ) );

            // check for duplicate name before adding the CT node

            CTNode[] otherNodes = findNodes( node.getName() );

            if ( otherNodes != null )
            {

                // so we have another node with the same name

                logger.error( ">>>> Name conflict for : " + region.getDescription() );

                // Note: we know that all other nodes are also regions

                for ( int j = 0; j < otherNodes.length; j++ )
                {

                    logger.error( "     already defined   : " + otherNodes[j].getRegion().getDescription() );
                }

            }

            addCTNode( i, node );

        } // for all Regions

        // number of edges can be less than calls for multiple calls
        // e.g. MAIN calls SUB in two different lines

        CallDescriptor[] calls = theRIFData.getCalls();

        final int progressPrintStep = 1000;

        boolean hasIndirectCalls = false;

        for ( i = 0; i < calls.length; i++ )
        {

            if ( ( i % progressPrintStep ) == 0 )
            {

                logger.info( "now at call " + i );
            }

            CallDescriptor call = calls[i];

            CallNode indCallee = call.getIndCallee();

            if ( indCallee != null )
            {

                hasIndirectCalls = true;

            }
            else
            {

                noNodes = addRIFCall( noNodes, call );

            }

        } // for all static calls in the RIF file

        // we have another loop for the indirect calls

        if ( hasIndirectCalls )
        {

            for ( i = 0; i < calls.length; i++ )
            {

                CallDescriptor call = calls[i];

                CallNode indCallee = call.getIndCallee();

                if ( indCallee != null )
                {

                    noNodes = addIndRIFCall( noNodes, call );

                }
            }
        }

        // calculate fictive costs for a single call

        if ( fictiveCosts )
        {

            // calculate fictive number of runtime calls

            calcFictiveRuntimeCalls();
            calcFictiveCosts();

            // calcFictiveCosts ();

            long[] totalValues = calcTotalCosts();

            CounterData.setTotals( totalValues );

        } // ficitiveCosts have been enabled

    } // constructor CTData from RIFData

    /**
     * Additional CT data (edge and maybe nodes) is generated from
     * a call descriptor given in the RIF file.
     *
     * @param noNodes is the number of CT nodes already in use.
     * @param call is the call descriptor to add
     * @return the number of CT nodes of the CT data
     */
    private int addRIFCall( int noNodes, CallDescriptor call )
    {

        int n = noNodes;

        // It is sure that we have already a CT node for every region

        CTNode node1 = getNode( call.getCallingRegion() );

        if ( node1 == null )
        {

            throw new RuntimeException( "region node not found" );

        }

        CallNode callee = call.getCallee();

        CTNode node2 = getCTNode( callee );

        if ( node2 == null )
        {

            node2 = new CTNode( callee );

            addCTNode( noNodes, node2 );

            n++;
        }

        CTEdge edge = getEdge( node1, node2 );

        edge.incSourceCall();

        return n;
    }

    /**
     * Additional CT data (edge and maybe nodes) is generated from
     * an indirect call descriptor given in the RIF file.
     *
     * @param noNodes is the number of CT nodes already in use.
     * @param call is the indirect call descriptor to add
     * @return the number of CT nodes of the CT data
     */
    private int addIndRIFCall( int noNodes, CallDescriptor call )
    {

        int n = noNodes;

        // It is sure that we have already a CT node for every region

        CTNode node1 = getNode( call.getCallingRegion() );

        if ( node1 == null )
        {

            throw new RuntimeException( "region node not found" );

        }

        CallNode callee = call.getCallee();

        CTNode node2 = getCTNode( callee );

        if ( node2 == null )
        {

            node2 = new CTNode( callee );

            addCTNode( noNodes, node2 );

            n++;
        }

        // For certain: node1 calls node3 but we have now
        // to find dummy in node3 where node2 will be the actual

        CallNode indCallee = call.getIndCallee();

        CTNode node3 = getCTNode( indCallee );

        if ( node3 == null )
        {

            node3 = new CTNode( indCallee );

            addCTNode( noNodes, node3 );

            n++;
        }

        CTEdge edge = getEdge( node3, node2 );

        edge.incSourceCall();

        return n;
    }

    /**
     * Efficient routine to get the CT node for a call node.
     *
     * @param callee is the call node
     * @return the CT node for the given call node (null if not found)
     */
    private CTNode getCTNode( CallNode callee )
    {

        CTNode node = findNode( callee );

        // in case of external nodes we try to find a matching region

        if ( node == null )
        {

            // we try to find a region node with the name

            CTNode[] callees = findNodes( callee.getName() );

            if ( callees != null )
            {

                if ( callees.length == 1 )
                {

                    node = callees[0];

                }
                else
                {

                    // give a warning for mutliple possibilities

                    logger.warn( "found " + callees.length + " matches for " + callee.getDescription() );
                }

            }

        }

        return node;

    }

} // class CTDataRIF
