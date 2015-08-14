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

package adaptor.EditPM;

import java.awt.GridLayout;

import javax.swing.JFrame;
import javax.swing.JPanel;

/**
 * The class EventFrame is a frame that contains all event buttons. It 
 * will be constructed for an array of performance events, but new 
 * events can be added later.
 *
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */

public class EventFrame extends JFrame {
    
    /**
     * This panel contains all the buttons for the performance events.
     */
    private JPanel myEventPanel; // panel contains button for every event
    
    /**
     * EditPM is needed to define target actions for new event buttons. 
     */
    private EditPM myEditPM;

    /**
     * Constructor for the event frame.
     * 
     * @param editor is pointer back to call routines of the PM edtior
     * @param events is an array for which event buttons will be available
     */
    EventFrame(EditPM editor, PerformanceEvent[] events) {
        
        super("Select Event for New Counter");
        
        myEditPM = editor;
        
        myEventPanel = new JPanel();
        
        final int noRows = 20;
        
        myEventPanel.setLayout(new GridLayout(noRows, 0, 3, 3));
        
        int i;
        
        for (i = 0; i < events.length; i++) {
            
            PerformanceEvent event = events[i];
            
            if (event != null) {
                
                EventButton button = new EventButton(editor, event);
                
                myEventPanel.add(button);                
            }
            
        }
        
        setContentPane(myEventPanel);
        
    } // constructor EventFrame
    
    /**
     * This routine must be called after having added a new event to
     * have also a button for this new event.
     * 
     * @param newPerformanceEvent is the new defined event
     */
    public void addEvent(PerformanceEvent newPerformanceEvent) {
        
        EventButton button = new EventButton(myEditPM, newPerformanceEvent);
        
        myEventPanel.add(button);
        
        myEventPanel.revalidate();
        
        pack();
        
    }
    
} // class EventFrame
