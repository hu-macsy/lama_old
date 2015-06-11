/*
 * PMMenuBar.java
 *
 * Menu bar for PM editor.
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

import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.KeyStroke;

import adaptor.General.AboutAdaptorFrame;


/***************************************************************************************************
 * * PM Menu Bar * *
 **************************************************************************************************/

/**
 * The class PMMenuBar defines a menu for the PM configuration editor.
 *
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */

class PMMenuBar extends JMenuBar implements ActionListener {
    
    /**
     * Reference back to the editor, set by constructor.
     */
    private EditPM myEditPM;
    
    /**
     * Menu item to quit the EditPM application.
     */
    private JMenuItem quitItem = new JMenuItem("Quit");
    /**
     * Menu item to save the current file.
     */
    private JMenuItem saveItem = new JMenuItem("Save");
    /**
     * Menu item to save the file under a new name.
     */
    private JMenuItem saveAsItem = new JMenuItem("Save as");
    /**
     * Menu item to make a new configuration.
     */
    private JMenuItem newItem = new JMenuItem("New");
    /**
     * Menu item to load a configuration file.
     */
    private JMenuItem loadItem = new JMenuItem("Load");
    /**
     * Menu item to add a configuration file.
     */
    private JMenuItem addItem = new JMenuItem("Add");
      
    /**
     * Menu item to copy the a counter speicification.
     */
    private JMenuItem copyItem = new JMenuItem("Copy");
    
    /**
     * Menu item to delete the selected counter specification. 
     */
    private JMenuItem remItem = new JMenuItem("Delete");
    
    /**
     * Menu item to paste the saved counter to a new one. 
     */
    private JMenuItem pasteItem = new JMenuItem("Paste");
    
    /**
     * Menu item to delete all user counter specifications. 
     */
    private JMenuItem remAllItem = new JMenuItem("Delete All");
    
    /**
     * Menu item to show the frame with all events.
     */
    private JMenuItem showEvItem = new JMenuItem("Show Events");
    
    /**
     * Menu item to show the frame for defining composed events.
     */
    private JMenuItem newComItem = new JMenuItem("New Composed Event");
    
    /**
     * Menu item to show the frame for defining derived events. 
     */
    private JMenuItem newDerItem = new JMenuItem("New Derived Event");
    
    /**
     * Menu item to pop up information about this tool.
     */
    private JMenuItem aboutItem = new JMenuItem("About");
    
    /**
     * This frame will be generated/shown when "About" is selected.
     */
    private AboutAdaptorFrame myAbout = null;
    
    /**
     * Constructor for the menu bar in the PM editor.
     * 
     * @param editor is reference to the editor for which the menu bar is created
     */
    public PMMenuBar(EditPM editor) {
        
        myEditPM = editor;
        
        // JMenuBar HeaderMenu = new JMenuBar ();
        
        JMenu fileMenu = new JMenu("File");
        JMenu counterMenu = new JMenu("Table-Counter");
        JMenu viewMenu = new JMenu("Events");
        JMenu helpMeun = new JMenu("Help");
        
        fileMenu.setMnemonic('F');
        counterMenu.setMnemonic('T');
        viewMenu.setMnemonic('E');
        viewMenu.setToolTipText("Showing and editing events");
        
        fileMenu.add(newItem);
        fileMenu.add(loadItem);
        fileMenu.add(addItem);
        fileMenu.add(saveItem);
        fileMenu.add(saveAsItem);
        fileMenu.add(quitItem);
        
        counterMenu.add(remItem);
        counterMenu.add(copyItem);
        counterMenu.add(pasteItem);
        counterMenu.add(remAllItem);
        
        viewMenu.add(showEvItem);
        viewMenu.add(newComItem);
        viewMenu.add(newDerItem);
        
        helpMeun.add(aboutItem);
        
        add(fileMenu);
        add(counterMenu);
        add(viewMenu);
        add(helpMeun);
        
        quitItem.addActionListener(this);
        saveItem.addActionListener(this);
        saveAsItem.addActionListener(this);
        newItem.addActionListener(this);
        loadItem.addActionListener(this);
        addItem.addActionListener(this);
        remItem.addActionListener(this);
        copyItem.addActionListener(this);
        pasteItem.addActionListener(this);
        remAllItem.addActionListener(this);
        showEvItem.addActionListener(this);
        newComItem.addActionListener(this);
        newDerItem.addActionListener(this);
        aboutItem.addActionListener(this);
        
        showEvItem.setToolTipText("Showing all known performance events");
        
        KeyStroke stroke = KeyStroke.getKeyStroke("ctrl shift C");
        copyItem.setAccelerator(stroke);
        stroke = KeyStroke.getKeyStroke("ctrl shift V");
        pasteItem.setAccelerator(stroke);
        stroke = KeyStroke.getKeyStroke("ctrl shift X");
        remItem.setAccelerator(stroke);
        
    } // PMMenuBar

    /**
     * {@inheritDoc}
     *
     * @see java.awt.event.ActionListener#actionPerformed(java.awt.event.ActionEvent)
     */
    public void actionPerformed(ActionEvent e) {
        
        JMenuItem item;
        
        item = (JMenuItem) e.getSource();
        
        if (item == quitItem) {
            
            System.exit(0);
        }
        
        if (item == saveItem) {
            
            myEditPM.save();
            
        }
        
        if (item == saveAsItem) {
            
            myEditPM.saveAs();
            
        }
        
        if (item == loadItem) {
            
            myEditPM.load();
            
        }
        
        if (item == addItem) {
            
            myEditPM.add();
            
        }
        
        if (item == newItem) {
            
            myEditPM.clear();
        }
        
        if (item == remItem) {

            myEditPM.deleteCounter();
        }
      
        if (item == copyItem) {
            
            myEditPM.copyCounter();            
        }

        
        if (item == pasteItem) {
            
            myEditPM.pasteCounter();
        }

        
        if (item == remAllItem) {
            
            myEditPM.clearCounters();
        }

        
        if (item == showEvItem) {
            
            // show panel with all events

            myEditPM.showEventFrame();
        }

        // show frame for the composed events
        
        if (item == newComItem) {

            myEditPM.showComposedFrame();
        }

        
        // show frame for the derived events
        
        if (item == newDerItem) {

            myEditPM.showDerivedFrame();
        }
        
        // show About Window
        
        if (item == aboutItem) {
            
            if (myAbout == null) {
                
                String info = "Editor for Performance Monitoring";
                
                myAbout = new AboutAdaptorFrame("EditPM", info);
            }
            
            myAbout.setVisible(true);
        }
        
    } // actionPerformed
}
