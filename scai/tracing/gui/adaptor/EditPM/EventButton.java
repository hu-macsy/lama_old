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

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;

/**
 * The class EventButton impelements a button assigned to each performance event.
 *
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */

public class EventButton extends JButton implements ActionListener
{

    /**
     * The performance event for which this button stands.
     */
    private PerformanceEvent myPerformanceEvent;

    /**
     * Pointer back to the PM editor.
     */
    private EditPM myEditPM;

    /**
     * Constructor for an event button.
     *
     * @param config is pointer for callback into the editor
     * @param event is the performance event for this button
     */
    public EventButton( EditPM config, PerformanceEvent event )
    {

        // name of the button is the name of the performance event

        super( event.getName() );

        // set the interface that we can call later method when button is selected

        myEditPM = config;

        setBackground( event.getColor() );

        myPerformanceEvent = event;

        setToolTipText( event.getDescription() );

        addActionListener( this );

    } // constructior EventButton

    /**
     * {@inheritDoc}
     *
     * @see java.awt.event.ActionListener#actionPerformed(java.awt.event.ActionEvent)
     */
    public void actionPerformed( ActionEvent e )
    {

        EventButton actionButton;

        actionButton = ( EventButton ) e.getSource();

        // we let the Editor decide what to do

        myEditPM.actionEvent( actionButton.myPerformanceEvent );

    }

} // class EventButton
