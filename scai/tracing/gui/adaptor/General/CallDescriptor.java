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

package adaptor.General;

import java.io.IOException;

import org.apache.log4j.Logger;

/**
 * class that desribes a call of a routine in a region.
 *
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */
public class CallDescriptor
{

    /**
     * Logger for this class.
     */
    private static Logger logger = Logger.getLogger( CallDescriptor.class );

    /**
     * region id of the calling region.
     */
    private RegionDescriptor myCallingRegion;

    /**
     * File id of file in which call appears. Even it is usually
     * the same file as the calling region, it could be some
     * included file.
     */
    private int myFileId = 0;

    /**
     * This is the line of the call.
     */
    private int myLine = 0;

    /**
     * This is routine we call.
     */
    private CallNode myCallee = null;

    /**
     * This might be the routine in case of an indirect call.
     * SUB1: call SUB2 (SUB3) : callee: SUB3, ind callee: SUB2
     */
    private CallNode myIndCallee = null;

    /**
     * Create a call descriptor from input line in RIF file. In the
     * call descriptor the region descriptors are set and not the region
     * identifiers.
     *
     * @param inputLine is line of RIF file
     * @param regions is an array of region descriptors where the region
     *        with regionId i is found at region[i-1]
     * @throws IOException in case of illegal line
     */
    public CallDescriptor( String inputLine, RegionTable regions ) throws IOException
    {

        try
        {

            logger.debug( "make call descriptor from this line: " + inputLine );

            int regionId = Utilities.readInt( inputLine, "region" );

            myCallingRegion = regions.getRegion( regionId - 1 );
            myFileId = Utilities.readInt( inputLine, "file" );
            myLine = Utilities.readInt( inputLine, "line" );

            int    kind = Utilities.readInt( inputLine, "kind", 2 );
            String name = Utilities.readString( inputLine, "name" );

            myCallee = makeCallNode( kind, name, regions );

            // now check for an indirect call
            // ANY_NODE as default can be used to see that there is no indirect callee

            kind = Utilities.readInt( inputLine, "bkind", CallNode.ANY_NODE );

            if ( kind != CallNode.ANY_NODE )
            {

                name = Utilities.readString( inputLine, "bname" );
                myIndCallee = makeCallNode( kind, name, regions );

            }

        }
        catch ( IOException e )
        {

            logger.error( e.getMessage(), e );

            throw e;

        }

    } // CallDescriptor

    /**
     * This routine generates a call node from the input data and region array.
     *
     * @param kind is the kind of the callee node
     * @param name is the name of the callee (is integer)
     * @param regions is an array with all region descriptors
     * @return a new call node generated from the input data
     */
    private CallNode makeCallNode( int kind, String name, RegionTable regions )
    {

        CallNode node;

        if ( kind == CallNode.REGION_NODE )
        {

            int regionId = Integer.parseInt( name );

            RegionDescriptor dsp = regions.getRegion( regionId - 1 );

            if ( dsp.getRegionId() != regionId )
            {

                throw new RuntimeException( "illegal region ids" );
            }

            node = new CallNode( regions.getRegion( regionId - 1 ) );

        }
        else
        {

            /* callee is not known as a region */

            node = new CallNode( kind, name );
        }

        return node;
    }

    /**
     * the region id of the called region (0 if region is not known).
     * @return callee as a CallNode
     */
    public CallNode getCallee()
    {

        return myCallee;
    }

    /**
     * the region id of the called region (0 if region is not known).
     * @return region identifiction of the called region
     */
    public CallNode getIndCallee()
    {

        return myIndCallee;
    }

    /**
     * the region id of the calling region.
     * @return id of the calling region.
     */
    public RegionDescriptor getCallingRegion()
    {

        return myCallingRegion;
    }

    /**
     * get the id of the file where the call is placed.
     * @return myFileId.
     */
    public int getFileId()
    {
        return myFileId;
    }

    /**
     * get the line in the file where the call is placed.
     * @return myLine.
     */
    public int getLine()
    {
        return myLine;
    }


    /**
     * This routine makes a string from the call node to be used
     * when writing call descriptor in the RIF file.
     *
     * @param prefix allows to use this routine for callee and indirect callee
     * @param node is the call node to write
     * @return a string that represents the call node
     */
    private String makeString( String prefix, CallNode node )
    {

        String callee = " " + prefix + "name=";

        callee += node.toString();

        callee += " " + prefix + "kind=" + node.getKind();

        return callee;
    }

    /**
     * make a string representation of this CallDescriptor.
     *
     * @return string with call description for RIF file
     */
    public String makeCallString()
    {

        // the call string starts with the region id of the calling region

        String callLine = "region=" + myCallingRegion.getRegionId();

        // then we put file number and line number

        callLine += " file=" + myFileId + " line=" + myLine;

        // then we add the callee

        callLine += makeString( "", myCallee );

        if ( myIndCallee != null )
        {

            callLine += makeString( "b", myIndCallee );
        }

        return callLine;

    } // makeCallString

} // class CallDescriptor
