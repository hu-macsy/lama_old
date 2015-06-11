/*
 * CTInterface.java
 * 
 * This file contains the interface of CalltreePM.
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

/**
 * The interface that must be implemented by the CalltreePM class that all methods of the other
 * classes work correctly.
 *
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */
public interface CTInterface {
    /**
     * The routine showGraph takes a file with a graph and displays it in the frame.
     * 
     * @param filename is the name of the file that contains the graph.
     */
    void showGraphFile(String filename);
         
    /**
     * showCallEnvironment prints the selected nodes and their calling environment. Depth for
     * calling and called routines is taken from the call environment frame.
     */
    void showCallEnvironment();
    
    /**
     * This routine will recalculate a calltree by the available properties.
     * Flags of nodes (mark, unfold) will be recalculated.
     */
    void recalc();
        
    /**
     * This routine will draw a new calltree with the new properties.
     */
    void redraw();
    
    /**
     * This routine returns the possible suffixes that
     * are supported to export data.
     * 
     * @return array of suffix strings
     */
    String[] getExportStrings();
    
    /**
     * Interface routine that has to be implemented
     */
    /**
     * This routine is called to export the current CT data 
     * to a file with this suffix. The suffix decides about the
     * format of the output file.
     * 
     * @param suffix must be a supported suffix
     */
    void export(String suffix);
    

    /**
     * This routine folds back all unfolded nodes.
     * 
     */
    void resetUnfolding();
    

    /**
     * This routine is used to show info on a certain object.
     * 
     * @param node is the node or list of nodes for which info will be shown
     */
    void showInfo(Object node);
    
} // interface CTInterface
