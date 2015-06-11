/*
 * CTEdge.java
 * 
 * A class used for representing edges in the Calltree.
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

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Vector;

import org.apache.log4j.Logger;

/**
 * A class used for representing edges in the Calltree (one node/routine calls another node/routine).
 *
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */
public class CTEdge {
    
    /**
     * Constant string used for apostrophing names. This is necessary
     * to allow special characters in names and labels. 
     */
    private static final String APO = "\""; 
    
    /**
     * Logger for this class.
     */
    private static Logger logger = Logger.getLogger(CTEdge.class);
    
    /**
     * These are the values for the costs of the corresponding calls.
     */
    private long[] myCallCounterValues;
    
    /**
     * The number of calls during the execution of the program.
     */
    private long myNoRuntimeCalls = 0;
    
    /**
     * The number of syntactical calls in source code.
     */
    private long myNoSourceCalls = 0;
    
    /**
     * mySource is the calling routine and source of the edge.
     */
    private CTNode mySource;
    
    /**
     * myTarget is the called routine and target of the edge.
     */
    private CTNode myTarget;
    
    /**
     * An edge can be a collection of subedges. In this case this
     * vector contains all the subedges.
     */
    private List<CTEdge> mySubEdges = null;
    
    /**
     * An edge is created by the specification of the source and the target.
     * 
     * @param theSource is the source node
     * @param theTarget is the target node
     */
    public CTEdge(CTNode theSource, CTNode theTarget) {
        
        int noCounters = CounterData.numberEvents();
        
        mySource = theSource;
        myTarget = theTarget;
        
        myCallCounterValues = new long[noCounters];
        
        for (int i = 0; i < noCounters; i++) {
            
            myCallCounterValues[i] = 0;
            
        } // counter intialization
        
    } // constrcutor CTEdge
    
    /**
     * Getter routine for the source of the edge.
     * 
     * @return the source node
     */
    CTNode getSource() {
        
        return mySource;
    }
    
    /**
     * Getter routine for the target of the edge.
     * 
     * @return the target node
     */
    CTNode getTarget() {
        
        return myTarget;
    }
    
    /**
     * This routine resets all runtime / call costs to zero.
     * 
     */
    void resetCosts() {
        
        for (int i = 0; i < myCallCounterValues.length; i++) {
            
            myCallCounterValues[i] = 0;            
        }
                    
        myNoRuntimeCalls = 0;
        
    } // resetCosts
    
    /**
     * Getter routine for the counter values of the call costs.
     * 
     * @return array with long value of all call costs
     */
    public long[] getCallCosts() {
        
        return myCallCounterValues;
        
    }
    
    /**
     * Getter routine for the number of runtime calls for this edge.
     * 
     * @return the number of calls at runtime
     */
    public long getRuntimeCalls() {
        
        return myNoRuntimeCalls;
    }
    
    /**
     * Getter routine for the number of source calls for this edge.
     * 
     * @return the number of known source calls
     */
    public long getSourceCalls() {
        
        return myNoSourceCalls;
    }
    
    /**
     * add costs to the call costs of the edge and to the inclusive costs of the calling node
     * (Source).
     * 
     * @param costs is the array containing cost values for each counted event
     */
    void addCallCosts(long[] costs) {
        
        // we have to add the values of the call to the inclusive costs
        
        for (int i = 0; i < costs.length; i++) {
            
            myCallCounterValues[i] += costs[i];
        }
        
        // add inclusive costs in case of non self-edge
        
        if (mySource != myTarget) {
            
            mySource.addInclusiveCosts(costs);
            
        }
        
    }
    
    /**
     * Increase the number of source code calls for this edge and for the calling node.
     */
    void incSourceCall() {
        
        myTarget.incSourceCall();
        
        myNoSourceCalls += 1;
        
    }
    
    /**
     * Add the number of runtime calls to an edge increases also the number of run calls for the
     * target.
     * 
     * @param calls is the number of runtime calls to be added
     */
    void addNoRunCalls(long calls) {
        
        myTarget.addNoRunCalls(calls);
        myNoRuntimeCalls += calls;
        
    }
    
    /**
     * Add the number of source calls to an edge increases also the number of source calls for the
     * target.
     * 
     * @param calls is the number of source calls to be added
     */
    void addNoSourceCalls(long calls) {
        
        myTarget.addNoSourceCalls(calls);
        myNoSourceCalls += calls;
        
    }
    
    /**
     * Add the costs of the edge E to this edge. This includes the number of source and runtime
     * calls
     * 
     * @param edge is the edge from which costs are just added
     */
    void addCosts(CTEdge edge) {
        
        myNoRuntimeCalls += edge.myNoRuntimeCalls;
        myNoSourceCalls += edge.myNoSourceCalls;
        
        for (int i = 0; i < myCallCounterValues.length; i++) {
            
            myCallCounterValues[i] += edge.myCallCounterValues[i];
        }
        
    }
    
    /**
     * Inherit the calls and costs of the edge E to this edge. Updates also the inclusive costs of
     * the source and the number of calls for the target node (if source and target are different)
     * 
     * @param edge is the edge of which the costs and calls are inherited
     */
    void inheritCosts(CTEdge edge) {
        
        // call costs are added to this edge in any case
        
        addCallCosts(edge.myCallCounterValues);
        
        if (mySource == myTarget)  {
            
            return;
        }
        
        addNoRunCalls(edge.getRuntimeCalls());
        addNoSourceCalls(edge.getSourceCalls());
        
    } // inheritCosts
       
    /**
     * Returns true if edge is visible in the current selection. An edge is only
     * visible if source and target are visible.
     * 
     * @return true if node is visible under the current conditions.
     */
    boolean isVisible() {
        
        // an edge is visible if source and target are visible
        
        boolean visible = mySource.isVisible() && myTarget.isVisible();
        
        if (visible) {

            visible = CTProperties.aboutDepth(myCallCounterValues, true);
        
        }

        return visible;
        
    }
    
    /**
     * Asking for the label of the edge. Its value depends on the
     * settings in CTProperites.
     * 
     * @return the label used for this edge in the dot file
     */
    String getCostLabel() {
        
        // class CTProperties has all the info about how
        // to make a label from the values of the Edge
        
        return CTProperties.getEdgeLabel(myNoSourceCalls, myNoRuntimeCalls, myCallCounterValues);
        
    } // getLabel
    
    /**
     * Asking for the color of this edge. The color depends on the 
     * settings in CTProperties.
     * 
     * @return the color to be used for this edge.
     */
    String getColor() {
        
        return CTProperties.getEdgeColor(myCallCounterValues);
            
    }
    
    /**
     * This routine makes an edge to a subedge of this edge.
     * 
     * @param subEdge is the edge that becomes of subedge
     */
    void addSubEdge(CTEdge subEdge) {
        
        if (mySubEdges == null) {

            mySubEdges = new Vector<CTEdge>();
        }
        
        mySubEdges.add(subEdge);
        
    } // add SubEdge
    
    /**
     * This routine returns the current number of subedges.
     * 
     * @return actual number of subedges
     */
    int getNoSubEdges() {

        int size = 0;
        
        if (mySubEdges != null) {

            size = mySubEdges.size();
        }

        return size;
        
    }
    
    /**
     * This routine allows to access the single subedges.
     * 
     * @param i with 0 <= i < getNoSubEdges()
     * @return the i-th subedge
     */
    CTEdge getSubEdge(int i) {
        
        return (CTEdge) mySubEdges.get(i);
        
    } // getSubEdge
    
    /**
     * This is one of the more complicated routines. It returns the edge or subedges 
     * dependenging whether source or target node are folded or not.
     * 
     * @return list of all subedges or folded ones
     */
    List<CTEdge> getFoldedEdges() {
        
        int noSubEdges = mySubEdges.size();
        
        // logger.info("getFoldedEdges: # subEdges = " + noSubEdges);
        // logger.info(mySource.getName() + " => " + myTarget.getName());
        
        List<CTEdge> foldedEdges = new Vector<CTEdge>();
        
        for (int i = 0; i < noSubEdges; i++) {
            
            CTEdge subEdge = (CTEdge) mySubEdges.get(i);
            
            CTNode source = subEdge.mySource;
            CTNode target = subEdge.myTarget;
            
            source = source.getFoldedParent();
            target = target.getFoldedParent();
            
            // try to find existing Edge
            
            // logger.debug("try to find edge " + source.getName() + " => " + target.getName());
            
            boolean found = false;
            
            for (int k = 0; k < foldedEdges.size(); k++) {
                
                CTEdge edge = (CTEdge) foldedEdges.get(k);
                
                if ((edge.mySource == source) && (edge.myTarget == target)) {
                    
                    // okay we have found an existing node
                    
                    edge.addCosts(subEdge);
                    
                    found = true;
                }
                
            } // for loop over existing Edges
            
            if (!found) {
                
                // so we have to add a new Edge
                
                CTEdge edge = new CTEdge(source, target);
                
                edge.addCosts(subEdge);
                
                foldedEdges.add(edge);
                
            }
            
        } // for all subEdges
        
        // logger.debug("getFoldedEdges done: # folded subEdges = " + foldedEdges.size());
        
        return foldedEdges;
        
    } // getFoldedEdges
    
    /**
     * This routine writes the edge data in a file.
     * 
     * @param buff is the ouput file writer
     * @param level specifies level to which unfolding is done
     * @throws IOException in case of problems
     */
    void writeEdge(BufferedWriter buff, int level) throws IOException {
        
        if (!isVisible()) {
            
            return;
        }

        
        String label = getCostLabel();
        String color = getColor();
        
        // Source and Target of a (Sub)Edge is a litte bit more complicated
        
        CTNode source = mySource.getPrintedNode(level);
        CTNode target = myTarget.getPrintedNode(level);
        
        String ident1 = source.getPrintIdent();
        String ident2 = target.getPrintIdent();
        
        // Attention: GRAPPA cannot deal with self-looping edges
        
        boolean doSubEdges = false;
        
        doSubEdges = (getNoSubEdges() > 0);
        
        if (!source.isUnfolded()) {
            
            if (!target.isUnfolded()) {
                
                doSubEdges = false;
            }
        }
        
        if (doSubEdges) {
            
            List<CTEdge> vecFoldedEdges = getFoldedEdges();
            
            for (int i = 0; i < vecFoldedEdges.size(); i++) {
                
                CTEdge edge = (CTEdge) vecFoldedEdges.get(i);
                
                edge.writeEdge(buff, level + 1);
                
            }
            
        } else if (source != target) {
            
            buff.write("    "
                       + APO + ident1 + APO 
                       + " -> " + APO + ident2 + APO
                       + " [label=" + APO + label + APO
                       + ", color=" + color + "];");
            buff.newLine();
            
        }
        
    } // writeEdge
    
} // CTEdge 

