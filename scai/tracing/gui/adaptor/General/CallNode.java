/*
 * CallNode.java
 *
 * This class is needed for nodes in call graphs as these
 * nodes are not always standing just for a region but also
 * for intrinsic, external and dummy rotuines.
 *
 * Members of this class are used within CallDescriptors
 * and for Calltree nodes.
 *
 * Created: 18.05.2006 brandes <email>
 * Changed:
 *
 * $Id$
 *
 * Copyright (C) 2006 Fraunhofer SCAI, Germany
 *
 * All rights reserved
 *
 * http://www.rcenviroment.de/
 */

package adaptor.General;

import org.apache.log4j.Logger;

import adaptor.Calltree.CalltreePM;

/**
 * This class is needed for nodes in call graphs as these
 * nodes are not always standing just for a region but also
 * for intrinsic, external and dummy routines.
 *
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */
public class CallNode
{

    /**
     * Constant for any kind of call node.
     */
    public static final int ANY_NODE = -1;

    /**
     * This kind specifieds a call node that stands for a region.
     */
    public static final int REGION_NODE = 0;

    /**
     * This kind specifieds a call node that stands for an intrinsic node.
     */
    public static final int INTRINSIC_NODE = 1;

    /**
     * This kind specifieds a call node that stands for an external node.
     */
    public static final int EXTERNAL_NODE = 2;

    /**
     * This kind specifieds a call node that stands for a dummy node.
     */
    public static final int DUMMY_NODE = 3;

    /**
     * This kind specifieds a call node that stands for a group node.
     */
    public static final int GROUP_NODE = 4;

    /**
     * These are the additional kind of nodes if the node is not
     * a region at all.
     */
    private static final String[] CALL_NODE_KIND = { "region", "intrinsic", "external", "dummy", "grouped" };

    /**
      * These will be CALL_NODE_KIND + REGION_KIND.
      */
    private static String[] theKindItems = null;

    /**
     * Logger for this class.
     */
    private static Logger logger = Logger.getLogger( CalltreePM.class );

    /**
     * This value gives more information about the node
     * that might be specified only by the name.
     *
     * kind takes one of the values ANY_NODE, ..., REGION_NODE
     */
    private int myKind = CallNode.ANY_NODE;

    /**
     * explicit name of the call node.
     */
    private String myName = null;

    /**
     * Only for myKind = 0: we have directly the region descriptor.
     */
    private RegionDescriptor myRegion = null;

    /**
     * Constructor for a CallNode by a given region.
     *
     * @param region is the region descriptor for which CallNode is created.
     */
    public CallNode( RegionDescriptor region )
    {

        myKind = REGION_NODE;
        myRegion = region;

    }

    /**
     * Constructor for a CallNode with name and kind only (no
     * region available).
     *
     * @param kind is the kind of the called routine
     * @param name is the name of the called routine
     */
    public CallNode( int kind, String name )
    {

        myKind = kind;

        if ( myKind == REGION_NODE )
        {

            throw new IllegalArgumentException( "CallNode for region not with this constructor" );

        }
        else
        {

            myName = name;

        }

    }

    /**
     * This routine returns the name of the CallNode as string but
     * for regions it returns the region id as a string.
     *
     * @return Representation for a CallNode as a string
     */
    public String toString()
    {

        if ( myKind == REGION_NODE )
        {

            return Integer.toString( myRegion.getRegionId() );

        }
        else
        {

            return myName;
        }
    }

    /**
     * This routine return just the name of the called node
     * (for regions it is the region name).
     *
     * @return the name of this call node as ident string
     */
    public String getName()
    {

        if ( myKind == REGION_NODE )
        {

            return myRegion.getName();

        }
        else
        {

            return myName;
        }
    }

    /**
     * This routine returns a detailled description of a call node.
     *
     * @return description of a call node as a String
     */
    public String getDescription()
    {

        return getName() + ":" + getNodeKindString();

    }

    /**
     * Getter routine for the call node kind. It returns 0
     * for a region.
     *
     * @return the kind of the call node as integer
     */
    public int getKind()
    {

        return myKind;
    }

    /**
     * This routine returns an integer value for the kind of a node.
     * The kind of a node will decide about its shape. The name of the
     * kind can be asked for with getNodeKindString.
     *
     * @return the kind of this node
     */
    public int getNodeKind()
    {

        int kind = 0;

        if ( myKind == REGION_NODE )
        {

            kind = myRegion.getKind();

        }
        else
        {

            kind = RegionDescriptor.RegionKind.values().length + ( myKind - 1 );
        }

        return kind;

    } // getNodeKind

    /**
     * This routine returns the call node kind.
     *
     * @return the call node kind as a string.
     */
    public String getNodeKindString()
    {

        int kind = getNodeKind();

        if ( theKindItems == null )
        {

            setCallNodeItems();

        }

        if ( kind < 0 )
        {

            return "ANY";

        }
        else
        {

            return theKindItems[kind];
        }

    }

    /**
     * Initialization routine to set up theKindItems.
     *
     */
    private static void setCallNodeItems()
    {

        int numberRegionKinds = RegionDescriptor.RegionKind.values().length;

        // for the other kinds we do not count the region nodes

        int numberOtherKinds = CALL_NODE_KIND.length - 1;

        String [] items = new String[numberRegionKinds + numberOtherKinds];

        for ( int i = 0; i < numberRegionKinds; i++ )
        {

            items[i] = RegionDescriptor.RegionKind.values()[i].toString();
        }

        for ( int i = 0; i < numberOtherKinds; i++ )
        {

            items[numberRegionKinds + i] = CALL_NODE_KIND[i + 1];
        }

        theKindItems = items;
    }

    /**
     * This routine returns all possible kinds of call nodes.
     *
     * @return an array with strings for all possible node kinds.
     */
    public static String[] getNodeKindStrings()
    {

        if ( theKindItems == null )
        {

            setCallNodeItems();
        }

        return theKindItems;
    }

    /**
     * Getter routine for the region of a call node. It can be
     * null if no region is available.
     *
     * @return the region descriptor of the region for which the call node stands
     */
    public RegionDescriptor getRegion()
    {

        return myRegion;
    }

    /**
     * This routine overrides the hashCode routine of Object.
     *
     * {@inheritDoc}
     *
     * @see java.lang.Object#hashCode()
     */
    public int hashCode()
    {

        logger.debug( "hash code of " + getName() + " = " + getName().hashCode() );

        return getName().hashCode();

    }

    /**
     * {@inheritDoc}
     *
     * @see java.lang.Object#equals(java.lang.Object)
     */
    public boolean equals( Object x )
    {

        boolean equal = false;

        if ( x instanceof CallNode )
        {

            // okay, we can compare the callnodes now

            CallNode otherNode = ( CallNode ) x;

            equal = getName().equals( otherNode.getName() );

            // the same name is not sufficient, both nodes must have the same kind
            // where kind == -1 matches with everything

            if ( ( myKind != ANY_NODE ) && ( otherNode.myKind != ANY_NODE ) )
            {

                equal = equal && ( myKind == otherNode.myKind );
            }
        }

        return equal;

    }
}
