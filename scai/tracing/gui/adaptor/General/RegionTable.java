/*
 * RegionTable.java
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

package adaptor.General;

import org.apache.log4j.Logger;


/**
 * An object of the class RegionTable is an indexed list of the RegionDescriptors.
 *
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */
public class RegionTable {

    /**
     * This is the initial number of regions.
     */
    private static final int INIT_NO_REGIONS = 10;

    /**
     * Logger for this class.
     */
    private static Logger logger = Logger.getLogger(RegionTable.class);
    
    /**
     * This is the array containing all regions. The region id is
     * used as an index for this array. The array itself can have holes
     * in it.
     */
    private RegionDescriptor[] theRegions = new RegionDescriptor[INIT_NO_REGIONS];

    /**
     * This constructor creates a region table with at least noRegions entries.
     * But the table might be extended.
     * 
     * @param noRegions is the number of entries that are expected.
     */
    public RegionTable(int noRegions) {
        
        this.theRegions = new RegionDescriptor[noRegions];
        
        for (int i = 0; i < noRegions; i++) {
            
            this.theRegions[i] = null;
        }
    }

    /**
     * This constructor creates a region table with a default number of entries. 
     */
    public RegionTable() {
    
        this(INIT_NO_REGIONS);
    }
    
    /**
     * This routine takes an already existing array of region descriptors.
     * 
     * @param inRegions is the array of known region descriptors.
     */
    public RegionTable(RegionDescriptor[] inRegions) {

        theRegions = inRegions;

    } // setRegions 

    /**
     * This routine reallocates the region data in such a way that an entry
     * for a certain region id can be made.
     * 
     * @param nr is the region id for which an entry in the array is needed
     */
    private void extendRegions(int nr) {

        int noRegions = theRegions.length;
        
        while (nr >= noRegions) {

            // logger.error("get_region: maximal number of regions reached");
            // logger.info("double region table to " + 2*MaxRegions  + " entries");

            // double number of Regions

            RegionDescriptor[] newRegions = new RegionDescriptor[2 * noRegions];

            for (int i = 0; i < noRegions; i++) {
                
                newRegions[i]             = theRegions[i];
                newRegions[i + noRegions] = null;
            }

            noRegions = 2 * noRegions;
            theRegions = newRegions;
        }

    } // extend_regions

    /**
     * This routine returns the region descriptor by the index.
     * 
     * @param index is the index of the region.
     * @return the descriptor of the region (might be null)
     */
    public RegionDescriptor getRegion(int index) {

        RegionDescriptor region = null;
        
        if ((index >= 0) || (index < theRegions.length)) {
            
            region = theRegions[index];
            
        }

        return region;
    }

    /**
     * This routine sets a region descriptor at a given index position.
     * 
     * @param index is the index for the region descriptor in the table
     * @param region is the region to be defined
     */
    public void setRegion(int index, RegionDescriptor region) {
        
        // make sure that the region array is big enough
        
        extendRegions(index);

        theRegions[index] = region;

    }

    /**
     * This routine defines a region descriptor at a given index position.
     * 
     * @param index is the index for the region descriptor in the table
     * @param region is the region to be defined
     * @return the new region descriptor or a matching existing one.
     */
    public RegionDescriptor defineRegion(int index, RegionDescriptor region) {
        
        // make sure that the region array is big enough
        
        extendRegions(index);

        if (theRegions[index] == null) {
            
            theRegions[index] = region;
            
            return region;
            
        }
        
        boolean match = (theRegions[index].getRegionId() == region.getRegionId());
        
        if (!match) {

            logger.error("defineRegion at index " + index + ": region with id " + 
                         region.getRegionId() + " overrides region with id " +
                         theRegions[index].getRegionId());
            
            System.exit(0);            
                        
        }

        return theRegions[index];
                  
    }

    /**
     * This routine return the highest index for which a region is available.
     * 
     * @return the highest index with a valid region.
     */
    public int noRegions() {
        
        int n = 0;
        
        if (theRegions != null) {
            
            n = theRegions.length;
            
            // there might be unused regions at the end of the table
            
            for (int k = n - 1; k >= 0; k--) {
                
                if (theRegions[k] == null) {
                    
                    n = k;
                    
                } else {
                    
                    break;
                }
            }
        }
        
        return n;
    }
}
