/*
 * CTMenuBar.java
 * 
 * Menu bar for Calltree
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

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JPopupMenu;

import org.apache.log4j.Logger;

import adaptor.General.AboutAdaptorFrame;


/**
 * class CTMenuBar defines the menu for the calltree visualization tool.
 *
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */
public class CTMenuBar extends JMenuBar implements ActionListener {

    /**
     * Logger for this class.
     */
    private static Logger logger = Logger.getLogger(CTMenuBar.class);

    /**
     * myCalltree is reference to the Calltree where this menu is anchored.
     */
    private CTInterface myCalltree;

    /**
     * This frame will be generated/shown when "About" is selected.
     */
    private AboutAdaptorFrame myAbout = null;

    /**
     * Constructor to create a menu for the calltree. The 
     * calltree is needed for callbacks.
     * 
     * @param theCalltree is the calltree for which menu is generated
     */
    public CTMenuBar(CTInterface theCalltree) {
    
        logger.info("construct CTMenuBar");
    
        // set the Interface
    
        this.myCalltree = theCalltree;
    
        // We tell CTProperties to make all its menus
        
        CTProperties.createPropertyMenus();
           
        // ExportMenu
    
        String [] menuItems = theCalltree.getExportStrings();
        String    help      = "export the calltree to a file of a certain type";
        
        add(newMenu("Export", help, menuItems));
   
        String [] helpItems = { "About" };
        
        add(newMenu("Help", "more information", helpItems));
        
    } // constructor CTMenuBar

    /**
     * This routine returns for a menu or menu item the invoker.
     * 
     * @param item is an Object that might have an invoker
     * @return the invoker object if available
     */
    private Object getInvoker(Object item) {
        
        Object parent = null;
        Object invoker = null;
        
        if (item instanceof JMenuItem) {
            
            parent = ((JMenuItem) item).getParent();
            
        } else if (item instanceof JMenu) {
            
            parent = ((JMenu) item).getParent();
            
        } else {
            
            logger.error("no parent for this object " + item);
        }
        
        if (parent instanceof JPopupMenu) {
            
            invoker = ((JPopupMenu) parent).getInvoker();
            
        } else if (parent instanceof CTMenuBar) {
            
            invoker = null;
            
        } else {
            
            logger.error("unknown invoker for this parent " + parent);
            
        }
        
        return invoker;
        
    }
    
    /**
     * This routine returns the name of an invoker which is assumed
     * to be a JMenu.
     * 
     * @param invoker is the invoker Object
     * @return a string containing the name of the invoker
     */
    private String getInvokerName(Object invoker) {
        
        String name = "";
        
        if (invoker instanceof JMenu) {
            
            name = ((JMenu) invoker).getText();         
            
        } else {
            
            logger.error("no name for this invoker " + invoker);
        }
        
        return name;
    }
 
    /**
     * {@inheritDoc}
     *
     * @see java.awt.event.ActionListener#actionPerformed(java.awt.event.ActionEvent)
     */
    public void actionPerformed(ActionEvent e) {

        Object actionSource = (Object) e.getSource();

        if (!(actionSource instanceof JMenuItem)) {

            return;
        }

        String cmd = ((JMenuItem) actionSource).getText();

        // now get the name of the menu to which item belongs
        // (mandatory if there are items with the same name)

        Object parent = getInvoker(actionSource);
        
        String title  = getInvokerName(parent); 

        if (title.equals("Export")) {

            myCalltree.export(cmd);
            return;

        }

        if (title.equals("Help") && cmd.equals("About")) {
            
            if (myAbout == null) {
                
                String info = "Calltree for Performance Monitoring";
                
                myAbout = new AboutAdaptorFrame("CalltreePM", info);
            }
            
            myAbout.setVisible(true);
        }

    } // actionPerformed

    /**
     * 
     * Help routine to create a new menu with different items and tool tip
     * and adds this class as ActionListener for all entries.
     * 
     * @param title is the title of the new menu
     * @param tip is the tool tip for this menu
     * @param items are the items of this menu
     * @return a new menu 
     */
    private JMenu newMenu(String title, String tip, String[] items) {

        int i; // loop variable

        // CounterMenu

        JMenu menu = new JMenu(title);

        for (i = 0; i < items.length; i++) {

            JMenuItem numberItem = new JMenuItem(items[i]);
            numberItem.addActionListener(this);
            menu.add(numberItem);

        }

        menu.setToolTipText(tip);

        return menu;

    } // help routine NewMenu

} // CTMenuBar
