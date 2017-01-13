/*
 * RegionDescriptor.java
 *
 * RegionDescriptor is a descriptor for a region
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

import org.apache.log4j.Logger;

/**
 * RegionDescriptor contains attributes and methods for a source code region.
 *
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */

public class RegionDescriptor
{

    /**
     * name of the attributes belong to a region (needed for RegionTable).
     */
    public static final String[] REGION_ENTRIES = { "Id", "RegionName", "ClassName",
                                 "FileId", "Lines", "Kind", "Enable", "Nested", "Data"
                                                  };

    /**
     * Enumeration type for all possible kind of regions.
     *
     * @version $LastChangedRevision$
     * @author Thomas Brandes
     */
    public static enum RegionKind {PROGRAM, SUBROUTINE, FUNCTION, USER, PARALLEL, LOOP, CALL, IO, ENTRY, DATA, CYGFUNC, CYGUSER};

    /**
     * This string is used to separate first and last line in String representation.
     * E.g. it could be " ", ":" or "-".
     */
    public static final String LINE_SEPARATOR = ":";

    /**
     * Logger for this class. Here it is very helpful to give
     * error messages about wrong constructions of region descriptors.
     */
    private static Logger logger = Logger.getLogger( PropertyMenu.class );

    /**
     * Unique identification of a region (integer value > 0). This
     * identification is not necessarily an index for the region table.
     */
    private int myRegionId;

    /**
     * kind of the region, 0 <= myRegionKind < regionKindTable.length.
     */
    private RegionKind myRegionKind;

    /**
     * This is the class name for the region (allows grouping of regions).
     */
    private String myClassName;

    /**
     * This the name of this region.
     */
    private String myRegionName;

    /**
     * This is the file identification and identifies the file where the
     * source code location of this region is.
     */
    private FileDescriptor myFile;

    /**
     * This is the line number where the region in the file starts.
     */
    private int myLineStart;

    /**
     * This is the line number where the region in the file stops.
     */
    private int myLineStop;

    /**
     * This flag indicates whether profiling for this region is enabled.
     */
    private boolean enableThisProfiling;

    /**
     * This flag indicates that profiling for nested regions is enabled.
     */
    private int nestedProfiling;

    /**
     * This flag is set true if sampling of data is enabled for this region.
     */
    private boolean enableDataSampling;

    /**
     * This constructor make a region from region and file id and other items.
     *
     * @param regionId is integer value for unique identification
     * @param regionName is name of the region
     * @param className is name of the class or module
     * @param regionKind is integer value for region kind
     * @param file is the source/binary file descriptor where this region is defined
     * @param startLine is integer value for start line
     * @param stopLine is integer value for stop line
     */
    public RegionDescriptor( int regionId, String regionName, String className,
                             int regionKind, FileDescriptor file, int startLine, int stopLine )
    {

        myRegionId   = regionId;
        myRegionName = regionName;
        myClassName  = className;
        myRegionKind = RegionKind.values()[regionKind];
        myFile       = file;  // has already been defined before
        myLineStart  = startLine;
        myLineStop   = stopLine;

        enableThisProfiling = false;  // no profiling
        nestedProfiling     = -1;     // do not restrict nested profiling
        enableDataSampling  = false;  // no data sampling

    }

    /**
     * Getter for the file identification.
     *
     * @return the file identification of the region
     */
    public FileDescriptor getFile()
    {

        return myFile;
    }

    /**
     * Getter for the region identification.
     *
     * @return the region identification of the region
     */
    public int getRegionId()
    {

        return myRegionId;
    }

    /**
     * This routine returns an attribute of the region.
     *
     * @param k (0 <= k < regionEntries.length)
     * @return an Object containing the k-th attribute
     */
    public Object getRegionObject( int k )
    {

        Object result;

        switch ( k )
        {

            case 0:
                // must fit with regionEntries[0] = "Id"
                result = new Integer( myRegionId );
                break;

            case 1:
                // must fit with regionEntries[0] = "RegionName"
                result = getName();
                break;

            case 2:
                // must fit with regionEntries[0] = "ClassName"
                result = myClassName;
                break;

            case 3:
                // must fit with regionEntries[0] = "FileId"
                result = new Integer( myFile.getFileId() );
                break;

            case 4:
                // must fit with regionEntries[0] = "Lines"
                result = myLineStart + LINE_SEPARATOR + myLineStop;
                break;

            case 5:
                // must fit with regionEntries[0] = "Kind"

                result = myRegionKind.toString();
                break;

            case 6:

                // must fit with regionEntries[0] = "Enable"
                if ( enableThisProfiling )
                {
                    result = Boolean.TRUE;
                }
                else
                {
                    result = Boolean.FALSE;
                }

                break;

            case 7:
                // must fit with regionEntries[0] = "Nested"
                result = new Integer( nestedProfiling );
                break;

            case 8:

                // must fit with regionEntries[0] = "Data"
                if ( enableDataSampling )
                {
                    result = Boolean.TRUE;
                }
                else
                {
                    result = Boolean.FALSE;
                }

                break;

            default:
                result = "";

        } // switch

        return result;

    } // method getRegionObject

    /**
     * Routine to set a certain attribute of the region at a given
     * position.
     *
     * @param value is the value to set
     * @param pos is the position where to set
     */
    public void setRegionObject( Object value, int pos )
    {

        boolean val;

        switch ( pos )
        {

            case 2:
                myClassName = ( String ) value;
                break;
            case 6:
                val = ( ( Boolean ) value ).booleanValue();
                enableThisProfiling = val;

                if ( !enableThisProfiling )
                {
                    enableDataSampling = false;
                }

                break;
            case 7:
                nestedProfiling = ( ( Integer ) value ).intValue();
                break;
            case 8:
                enableDataSampling = ( ( Boolean ) value ).booleanValue();
                break;
            default:
                // should not happen
                logger.fatal( "setRegionObject: should not happen" );

        } // switch

    } // setRegionObject

    /**
     * Check whether the k-th attribute of a region can be modified. At this
     * time we can only change the class name or the flags to enable profiling.
     *
     * @param k is the number of the attribute
     * @return true if attribute can be changed
     */
    public boolean isEditable( int k )
    {

        boolean is;

        is = ( k == 2 ) // "ClassName"
             || ( k == 6 ) // "Enable"
             || ( k == 7 ) // "Nested"
             || ( ( k == 8 ) && enableThisProfiling ); // "Data" only if enabled

        return is;

    } // method isEditable

    /**
     * make a one line description of the region for output info.
     *
     * @return String contains description
     */
    public String getDescription()
    {

        return "region=" + myRegionId + " file=" + myFile.getLongFileName()
               + " lines=" + myLineStart + LINE_SEPARATOR + myLineStop
               + " class=" + myClassName + " name=" + myRegionName + " kind=" + myRegionKind;

    } // getDescription

    /**
     * make a string of the region as used in the RIF file.
     * RegionDescriptor (R.makeRegionString()) is same region a R
     *
     * @return one line string for the region
     */
    public String makeRegionString()
    {

        int enable;

        if ( !enableThisProfiling )
        {
            enable = 0;
        }
        else if ( !enableDataSampling )
        {
            enable = 1;
        }
        else
        {
            enable = 2;
        }

        return    "region=" + myRegionId
                  + " file=" + myFile.getFileId()
                  + " lines=" + myLineStart + LINE_SEPARATOR + myLineStop
                  + " class=" + myClassName
                  + " name="  + myRegionName
                  + " kind=" + myRegionKind.ordinal()
                  + " enable=" + enable
                  + " depth=" + nestedProfiling;

    } // makeRegionString

    /**
     * make a string containing the main attributes (not the the profiling flags). This
     * representation is used in PMS files.
     *
     * @return String with attribute values separated by a space
     */
    public String getAttributeString()
    {

        String attribute = myRegionId + " " + myRegionName + " " + myClassName + " " + myRegionKind.ordinal();

        attribute += " " + myFile.getLongFileName() + " " + myLineStart + " " + myLineStop;

        return attribute;

    } // getAttributeString

    /**
     * Gets the kind of this region as a string.
     *
     * @return the kind of the region as string
     */
    public String getKindStringOld()
    {

        return myRegionKind.toString();

    }

    /**
     * Gets the kind of this region as an integer value.
     *
     * @return kind of this region
     */
    public int getKind()
    {

        return myRegionKind.ordinal();

    }

    /**
     * Gets the kind of this region as an enum value.
     *
     * @return kind of this region
     */
    public RegionKind getRegionKind()
    {

        return myRegionKind;

    }

    /**
     * Getter routine for the name of the region. In case
     * of no name we use the region id as a hexadecimal string.
     *
     * @return name of the region.
     */
    public String getName()
    {

        if ( myRegionName.startsWith( "?" ) )
        {

            return "0x" + Integer.toHexString( myRegionId );
        }

        return myRegionName;
    }

    /**
     * Setter routine for the name of the region. In case of
     * C++ name might contain also param and template arguments
     *
     * @param name is the new name for the region.
     */
    public void setName( String name )
    {

        int posArgList = name.indexOf( '(' );

        if ( posArgList > 0 )
        {

            // we do not take the argument list in the name

            name = name.substring( 0, posArgList );

        }

        int posSeparation = name.lastIndexOf( "::" );

        if ( posSeparation > 0 )
        {

            myClassName = name.substring( 0, posSeparation );
            myRegionName = name.substring( posSeparation + 2 );

        }
        else
        {

            myClassName  = "";
            myRegionName = name;

        }
    }

    /**
     * Getter routine for the class of the region.
     *
     * @return the class name of the region
     */
    public String getClassName()
    {

        return myClassName;

    }

    /**
     * Getter routine for first line of the region.
     *
     * @return first line of region in source file
     */
    public int getFirstLine()
    {

        return myLineStart;
    }

    /**
     * Getter routine for the last line of the region.
     *
     * @return last line of region in source file
     */
    public int getLastLine()
    {

        return myLineStop;
    }

    /**
     * This routine update the info about file location for a
     * descriptor. It is used when a region id has been identified
     * by the routine addr2line.
     *
     * @param file is the new descriptor for the file.
     * @param lineStart is the line number where the region starts.
     * @param lineStop is the line number where the region ends.
     */
    public void setFileInfo( FileDescriptor file, int lineStart, int lineStop )
    {

        myRegionId = myLineStart;
        myLineStart = lineStart;
        myLineStop = lineStop;
        myFile     = file;
    }

    /**
     * Setter routine for the profiling flag.
     *
     * @param flag is the new value for profiling flag.
     */
    public void setProfEnabled( boolean flag )
    {

        enableThisProfiling = flag;

    }

    /**
     * Setter routine for the nest flag.
     *
     * @param depth is the new value for the nest depth.
     */
    public void setNestProfiling( int depth )
    {

        nestedProfiling = depth;

    }

    /**
     * Setter routine for the data sampling flag.
     *
     * @param flag is the new value for the data sampling.
     */
    public void setDataEnabled( boolean flag )
    {

        enableDataSampling = flag;
    }

    /**
     * This routine appends info about the region at a SringBuffer.
     *
     *
     * @param buffer is the StringBuffer where info will be appended.
     * @param indentString is a indent String for each new line.
     */
    public void appendInfo( StringBuffer buffer, String indentString )
    {

        buffer.append( indentString );
        buffer.append( "Region: " );
        buffer.append( getName() );
        buffer.append( "\n" );
        buffer.append( indentString );
        buffer.append( "Class: " );
        buffer.append( getClassName() );
        buffer.append( "\n" );
        buffer.append( indentString );
        FileDescriptor file = getFile();
        buffer.append( "File: " );
        buffer.append( file.getLongFileName() );
        buffer.append( ", lines = " );
        buffer.append( getFirstLine() + ":" + getLastLine() + "\n" );
        buffer.append( indentString );
        buffer.append( "Kind: " + getRegionKind() + "\n" );

    }

} // class RegionDescriptor

