/*
 * CTPanel.java
 * 
 * Panel for CalltreePM that contains a kill button for layout thread.
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

import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JPanel;

import org.apache.log4j.Logger;

/**
 * CTPanel is a panel for CalltreePM that contains a button that allows to
 * kill the current layout thread.
 *
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */
public class CTPanel extends JPanel implements ActionListener {


    /**
     * Logger for this class.
     */
    private static Logger logger = Logger.getLogger(CTPanel.class);

    /**
     * Button that is used to kill the layout thread.
     */
    private JButton myKillButton;

    /**
     * Reference to the layout thread that can be stopped.
     */
    private LayoutThread myThread;

    /**
     * Constructor for a new panel. The kill button is disabled by default.
     * 
     */
    public CTPanel() {
    
        setLayout(new GridLayout(1, 1, 10, 0));
    
        myKillButton = new JButton("Kill Layout Thread");
    
        myKillButton.addActionListener(this);
    
        add(myKillButton);
    
        myKillButton.setEnabled(false);
    
    } // constructor CTPanel

    /**
     * {@inheritDoc}
     *
     * @see java.awt.event.ActionListener#actionPerformed(java.awt.event.ActionEvent)
     */
    public void actionPerformed(ActionEvent e) {

        Object actionSource = (Object) e.getSource();

        logger.info("CTACtion: action performed");

        if (actionSource instanceof JButton) {

            JButton actionButton = (JButton) actionSource;

            if (actionButton == myKillButton) {

                // Attention: getId only available since 1.5
                
                if (myThread != null) {
                    
                    logger.info("cancel the layout thread " 
                                + myThread.getName()
                                + ", id =  x "); 

                    // myThread.getId()
                    
                    myThread.cancel();
                    
                    actionButton.setEnabled(false);
                }  
                
  
            }


        } // instanceof JButton 

    } // actionPerformed

    /**
     * This routine will enable the kill button that can be used to cancel
     * a thread.
     * 
     * @param theThread is the thread that will be canceled
     */
    public void enableKillButton(LayoutThread theThread) {

        myThread = theThread;
        myKillButton.setEnabled(true);

    }

    /**
     * The kill button will be disabled as the thread might have finished.
     * 
     */
    public void disableKillButton() {

        myKillButton.setEnabled(false);
    }

} // CTPanel

