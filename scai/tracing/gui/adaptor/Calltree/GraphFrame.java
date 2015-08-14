/*
 * GraphFrame.java
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

package adaptor.Calltree;

import java.awt.Component;
import java.awt.Dimension;
import java.awt.Point;

import javax.swing.JFrame;
import javax.swing.JMenuBar;
import javax.swing.JScrollPane;
import javax.swing.JTabbedPane;

import att.grappa.Graph;
import att.grappa.GrappaPanel;
import att.grappa.Subgraph;

/**
 * The class GraphFrame is a frame used to display the call graph generated for a calltree
 * data set.
 *
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */
public class GraphFrame extends JFrame {

    /**
     * The following dimension(width,height) is taken for the size
     * of the graph display frame.
     * 
     * (1000, 800) is good dimension for 1240 x 1028 screeen
     * (800, 600) is good dimension for 1024 x 800 screen
     */
    private static final Dimension SIZE_DIM = new Dimension(1000, 800);
    
    /**
     * The following point is taken as the left upper corner for
     * the graph display frame.
     */
    private static final Point ORIGIN_POINT = new Point(120, 100);
    
    /**
     * The panel to show the call graph via Grappa.
     */
    private GrappaPanel myGrappaPanel = null;
    
    /**
     * Enables the graph to be scrolled.
     */
    private JScrollPane myScrollPane = null;

    /**
     * Pane for the different tabs.
     */
    private JTabbedPane myTabbedPane = null;

    /**
     * A panel that allows to define additional buttons for the Calltree.
     */
    private CTPanel myCalltreePanel = null;

    /**
     * Constructor to genreate a Graph frame. It has no related call tree
     * at the beginning.
     * 
     * @param graphMenuBar must be an instance of the menu related to this frame.
     */
    GraphFrame(JMenuBar graphMenuBar) {

        super("GraphFrame");

        this.setSize(SIZE_DIM);
        this.setLocation(ORIGIN_POINT);

        myScrollPane = new JScrollPane();
        
        myTabbedPane = new JTabbedPane(JTabbedPane.TOP);
        
        myTabbedPane.addTab("Calltree", myScrollPane);
        
        myCalltreePanel = new CTPanel();

        getContentPane().add("North", graphMenuBar);
        getContentPane().add("Center", myTabbedPane);
        getContentPane().add("South", myCalltreePanel);

        setVisible(true);

    } // constructor GraphFrame

    /**
     * This routine adds a new selection panel as a tabbed pane.
     * 
     * @param title is the name of the new tab
     * @param newComponent is the panel that will be added
     */
    public void addTab(String title, Component newComponent) {
        
        myTabbedPane.add(title, newComponent);
        
    }
    
    /**
     * This routine is used to display the call graph in this frame. The CT interface
     * allows to call routines for certain actions.
     * 
     * @param graph is the graph to be displayed
     * @param theMain implements the CTInterface
     */
    void showGraph(Graph graph, CTInterface theMain) {

        myGrappaPanel = new GrappaPanel(graph);

        // gp.addGrappaListener(new GrappaAdapter());
        // attention: we take our own listener 

        myGrappaPanel.addGrappaListener(new MyGrappaListener(theMain));

        myGrappaPanel.setScaleToFit(false);

        MySubgraph.setLastSelection(myGrappaPanel.getSubgraph());

        myScrollPane.setViewportView(myGrappaPanel);
        
        myTabbedPane.setSelectedComponent(myScrollPane);

        repaint();

    } // showGraph

    /**
     * Getter routine for the current selected subgraph.
     * 
     * @return the selected Subgraph
     */
    public Subgraph getSubgraph() {
        
        if (myGrappaPanel == null) {
            return null;
        } else {
            return myGrappaPanel.getSubgraph();
        }
    }
    
    /**
     * Enables a kill button to kill thread T.
     * 
     * @param theThread is the thread that can be killed by the button
     */
    void enableKill(LayoutThread theThread) {

        myCalltreePanel.enableKillButton(theThread);
    }

    /**
     * Disable the kill button. Should be called when the thread has already
     * finised. 
     * 
     */
    void disableKillButton() {

        myCalltreePanel.disableKillButton();
    }

} // class GraphFrame
