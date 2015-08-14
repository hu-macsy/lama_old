/*
 * ValueButton.java
 *
 * File contains class for a new button class that contains editable value.
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
import java.io.File;

import javax.swing.JButton;
import javax.swing.JOptionPane;

import org.apache.log4j.Logger;

/**
 * The class ValueButton defines a button that contains a value of a certain type.
 *
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */

class ValueButton extends JButton implements ActionListener {
    
    /**
     * This is the logger variable for this class to use log4j.
     */
    private static Logger logger = Logger.getLogger(EditPM.class);
    
    /**
     * The variable "name" contains the name of the button for which
     * a value can be defined.
     */
    private String name = null;
    
    /**
     * This variable is the actual value for the button.
     */
    private Object value = null;
    
    /**
     * Construct a new value button with name, initial value and an 
     * additional ActionListener that is called after the value of the 
     * button has been changed.
     * 
     * @param buttonName is the name of the variable to set by button
     * @param initValue is the initial value
     * @param change is an addtional listener for updates
     */
    public ValueButton(String buttonName, Object initValue, ActionListener change) {
        
        super(buttonName);
        
        this.name = buttonName;
        value = initValue;
        
        setEnabled(false);
        
        // first action listener will be called after the second one
        
        addActionListener(change);
        
        // when button is selected its value can be changed
        
        addActionListener(this);
        
    }
    
    /**
     * {@inheritDoc}
     *
     * @see java.awt.event.ActionListener#actionPerformed(java.awt.event.ActionEvent)
     */
    public void actionPerformed(ActionEvent e) {
        
        ValueButton actionButton;
        
        actionButton = (ValueButton) e.getSource();
        
        String input;
        
        input = JOptionPane.showInputDialog(null, "Enter the " + name + ": ", value);
        
        // input can be cancelled, so do nothing
        
        if (input != null) {
            
            // the type of value decides how to make a new value
            
            if (value instanceof Integer) {
                
                try {
                    
                    value = new Integer(Integer.parseInt(input));
                    
                } catch (NumberFormatException excetpion) {
                    
                    JOptionPane.showMessageDialog(null, "Illegal integer value: " + input,
                                                  "Error for button " + name, JOptionPane.ERROR_MESSAGE);
                }
                
            } else if (value instanceof String) {
                
                value = input;
                
            } else if (value instanceof File) {
                
                File f = new File(input);
                
                if (!f.exists()) {
                    
                    JOptionPane.showMessageDialog(null, f + " does not exist", "RIF error", JOptionPane.ERROR_MESSAGE);
                }
                
                value = f;
                
            } else {
                
                logger.error("unsupported value for button " + name);
            }
            
            actionButton.setText(name + " = " + value);  
        }
    }
    
    /**
     * {@inheritDoc}
     *
     * @see java.awt.Component#setEnabled(boolean)
     */
    public void setEnabled(boolean flag) {
        
        logger.info("set enabled = " + flag + " for button " + name);
        
        String text = name;
        
        if (flag) {
            text += " = " + value;
        }
        
        setText(text);
        
        super.setEnabled(flag);
    }
    
    /**
     * Getter routine for access to the value of the button.
     * 
     * @return the value of the button
     */
    public Object getValue() {
        
        return value;
    }
    
    /**
     * Setter routine to set the value of the button and to enable it.
     * 
     * @param val becomes the actual value for the button
     */
    public void setValue(Object val) {
        
        logger.info("set value for button " + name + ", value = " + val);
        value = val;
        setEnabled(true);
    }
}

