/*
 * DerivedFrame.java
 *
 * File contains class DerivedFrame.
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
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.SwingConstants;

/**
 * The class DerivedFrame provides a frame to define derived events.
 * 
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */

class DerivedFrame extends JFrame implements ActionListener {

    /**
     * Text field that contains the name of the new derived event.
     */
    private JTextField name;

    /**
     * Combo box used to select the performance event that will be derived.
     */
    private JComboBox event;

    /**
     * Menu bar for the deriving operation (wall time or cpu time).
     */
    private JMenuBar opMenuBar;
        
    /**
     * This is the label that contains the derivation operation.
     */
    private JLabel opLabel;

    /**
     * Button that closes this frame.
     */
    private JButton quitButton;
    
    /**
     * Button that defines a new derived event with the current settings.
     */
    private JButton defineButton;

    /**
     * Pointer to the current PM configuration.
     */
    private PMConfiguration myPMConfiguration;

    /**
     * Constructor for the frame to define derived events.
     * 
     * @param config is the pointer to the PM configuration
     */
    DerivedFrame(PMConfiguration config) {
    
        super("Derived Event");
    
        myPMConfiguration = config;
    
        JPanel pane = new JPanel();
    
        pane.setLayout(new GridLayout(4, 2, 10, 10));
        pane.setBackground(DerivedEvent.getEventColor());
    
        JLabel label;
    
        final int textFieldSize = 20;
        
        name = new JTextField(textFieldSize);
    
        quitButton = new JButton("Quit");
        defineButton = new JButton("Define");
    
        String[] events = myPMConfiguration.getBasicEventStrings();
    
        // make a Combo Box of all basic events
    
        event = new JComboBox();
        event.setBackground(BasicEvent.getEventColor());
        int i;
        for (i = 0; i < events.length; i++) {
            event.addItem(events[i]);
        }
    
        defineButton.addActionListener(this);
        quitButton.addActionListener(this);
    
        JMenu opMenu = new JMenu("Op");
    
        JMenuItem opMenuItem;
    
        opMenuItem = new JMenuItem("'");
        opMenuItem.setToolTipText("Event / Walltime");
        opMenuItem.addActionListener(this);
        opMenu.add(opMenuItem);
    
        opMenuItem = new JMenuItem("\"");
        opMenuItem.setToolTipText("Event / CPUtime");
        opMenuItem.addActionListener(this);
        opMenu.add(opMenuItem);
    
        // we make the operation MenuBar
    
        opMenuBar = new JMenuBar();
        opMenuBar.setBorderPainted(true);
        opMenuBar.add(opMenu);
    
        label = new JLabel("Derived Event", SwingConstants.RIGHT);
        pane.add(label);
        pane.add(name);
    
        label = new JLabel("Event", SwingConstants.RIGHT);
        pane.add(label);
        pane.add(event);
    
        opLabel = new JLabel("'", SwingConstants.LEFT);
        pane.add(opMenuBar);
        pane.add(opLabel);
    
        pane.add(defineButton);
        pane.add(quitButton);
    
        setContentPane(pane);
        pack();
    
    } // DerivedFrame

    /**
     * {@inheritDoc}
     *
     * @see java.awt.event.ActionListener#actionPerformed(java.awt.event.ActionEvent)
     */
    public void actionPerformed(ActionEvent e) {

        Object actionSource = e.getSource();

        if (actionSource instanceof JMenuItem) {

            JMenuItem item = (JMenuItem) actionSource;

            opLabel.setText(item.getText());
            opLabel.setToolTipText(item.getToolTipText());

        }

        if (actionSource instanceof JButton) {

            JButton actionButton = (JButton) e.getSource();

            if (actionButton == quitButton) {
                setVisible(false);
            }

            if (actionButton == defineButton) {

                PerformanceEvent p;

                p = myPMConfiguration.findSupportedEvent((String) event.getSelectedItem());

                DerivedEvent derEvent = new DerivedEvent(name.getText(), p, opLabel.getText());

                boolean error = false;

                // try to add it to the table of events

                try {
                    
                    myPMConfiguration.addNewEvent(derEvent);
                    
                } catch (RuntimeException ex) {

                    // show error message

                    JOptionPane.showMessageDialog(null, ex.getMessage(), "Error for derived event", JOptionPane.ERROR_MESSAGE);

                    error = true; // in case of error we do not close

                } // catch

                if (!error) {
                    
                    setVisible(false);
                }

            } // B == Define

        } // actionPerformed on JButton

    } // actionPerformed

    // setEvent: might be called via the Event Panel

    /**
     * This routine is used to set an event that will be derived.
     * 
     * @param setEvent is the performance event that will be derived.
     */
    void setEvent(PerformanceEvent setEvent) {

        int index = myPMConfiguration.getEventIndex(setEvent);

        if (index < 0) {

            return;
        }

        event.setSelectedIndex(index);
    }

    /**
     * This method takes a derived event to fill the values of this frame.
     * 
     * @param derEvent is an existing derived event
     */
    void editEvent(DerivedEvent derEvent) {

        PerformanceEvent basicEvent = derEvent.getEvent();

        event.setSelectedIndex(myPMConfiguration.getEventIndex(basicEvent));

        opLabel.setText(derEvent.getKind());

        name.setText(derEvent.getName());

    } // editEvent
}
