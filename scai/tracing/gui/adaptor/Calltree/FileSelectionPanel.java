/*
 * FileSelectionPanel.java
 * 
 * FileSelectionPanel is an extension of JPanel that defines
 * a panel for the selection of different input files. Each 
 * file can be enabled or disabled.
 * 
 * Created: 2008-01-03 Thomas Brandes <thomas.brandes@scai.fraunhofer.de>
 * Changed:
 * 
 * $Id$
 * 
 * Copyright (C) 2008 Fraunhofer SCAI, Germany
 * 
 * All rights reserved
 *
 * http://www.scai.fhg.de/EP-CACHE/adaptor
 */

package adaptor.Calltree;

import java.awt.Color;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.BorderFactory;
import javax.swing.JCheckBox;
import javax.swing.JPanel;
import javax.swing.border.Border;

import org.apache.log4j.Logger;


/**
 * FileSelectionPanel is an extension of JPanel that defines
 * a panel for the selection of different input files. Each 
 * file can be enabled or disabled.
 * 
 *
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */
public class FileSelectionPanel extends JPanel implements ActionListener {

    /**
     * Logger for this class.
     */
    private static Logger logger = Logger.getLogger(FileSelectionPanel.class);
    
    /**
     * Array of check boxes for each available file.
     */
    private JCheckBox[] myFileBoxes;
    
    /**
     * Array of flags to take the current selection of enabled files.
     * It is valid: myFilesEnabled.length == myFileBoxes.length
     */
    private boolean[] myFilesEnabled;


    /**
     * Constructor for a file selection panel.
     * 
     * @param fileNames is an array of file names that can be enabled or disabled
     */
    public FileSelectionPanel(String[] fileNames) {
    
        super();
    
        int noFiles = fileNames.length;

        Border myBorder = BorderFactory.createLineBorder(Color.BLACK);
        
        myBorder = BorderFactory.createTitledBorder(myBorder, "File Selection");
        
        setBorder(myBorder);
        
        setLayout(new GridLayout(noFiles, 1, 10, 10));
    
        myFileBoxes = new JCheckBox[noFiles];
        myFilesEnabled = new boolean[noFiles];
    
        for (int i = 0; i < noFiles; i++) {
    
            myFilesEnabled[i] = true;
            myFileBoxes[i] = new JCheckBox(fileNames[i], null, myFilesEnabled[i]);
            myFileBoxes[i].addActionListener(this);
            add(myFileBoxes[i]);
    
        }
               
    }

    /**
     * This routine returns an boolean array to get for all
     * files the selected flag.
     * 
     * @return getSelectedFiles[i] is true if file[i] has been selected.
     */
    public boolean[] getSelectedFiles() {

        return myFilesEnabled;
    }
    
    /**
     * {@inheritDoc}
     * 
     * @see java.awt.event.ActionListener#actionPerformed(java.awt.event.ActionEvent)
     */
    public void actionPerformed(ActionEvent e) {

        Object actionSource = e.getSource();
        
        if (actionSource instanceof JCheckBox) {

            JCheckBox actionBox = (JCheckBox) actionSource;
            
            logger.info("file " + actionBox.getText() + " is selected = " + actionBox.isSelected());
            
            for (int i = 0; i < myFileBoxes.length; i++) {
                
                if (actionBox == myFileBoxes[i]) {
                    
                    myFilesEnabled[i] = myFileBoxes[i].isSelected();
                    
                }
            }


        } // Source was JCheckBox
    }

} // FileSelectionPanel
