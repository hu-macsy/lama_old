/*
 * EditRIFInterface.java
 *
 * Interface that has to be implemented by the RIF Editor.
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

package adaptor.EditRIF;

import adaptor.General.FileDescriptor;
import adaptor.General.RegionDescriptor;

/**
 * The interface for EditRIF specifies the routines that can be called for it. The routines will be
 * called especially from the menu bar.
 * 
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */

public interface EditRIFInterface {
    
    /**
     * Quits the RIF editor but makes also sure that a modified RIF file is saved.
     */
    void quitRIF();
    
    /**
     * This routine is used to get access to the file by the unique file identification.
     * 
     * @param fileId
     *            is the unique file identifictation
     * @return the file descriptor (may be null for a wrong file id)
     */
    FileDescriptor getFileDescriptor(int fileId);
    
    /**
     * This routine saves the current RIF file if one is open.
     * 
     */
    void save();
    
    /**
     * This routine saves the current RIF file and allows to give it a new name.
     * 
     */
    void saveAs();
    
    /**
     * This routine marks the current RIF data as having been modified. It guarantees that modified
     * RIF data will be written back.
     * 
     */
    void setModified();
    
    /**
     * This routine returns all region descriptors in an array. The unique region id can be used to
     * access a region descriptor.
     * @param index refers to the index in the region table 
     * @return the region descriptor at the given index position
     */
    RegionDescriptor getRegion(int index);
    
    /**
     * This routine returns the last index of a region descriptor in the region table.
     * 
     * @return the number of regions in the region table.
     */
    int noRegions();
    
}