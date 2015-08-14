/*
 * ShowPM.java
 * 
 * GUI to show the results of Performance monitoring.
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

package adaptor.ShowPM;


import java.io.File;

import java.awt.Dimension;

import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JOptionPane;

import org.apache.log4j.Logger;

import adaptor.General.ArraySelection;
import adaptor.General.FileDescriptor;
import adaptor.General.OptionsAdaptor;
import adaptor.General.RegionDescriptor;
import adaptor.General.ShowFileFrame;
import adaptor.General.ShowInfo;

/**
 * The class ShowPM implements a GUI to show performance results of
 * measurements made by performance monitoring.
 *
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */
public class ShowPM implements SelectedInterface {

    /**
     * Preferred dimension for ShowPM, taken as default.
     */
    private static final Dimension SHOW_PM_DIMENSION = new Dimension(800, 400);

    /**
     * Constant for maximal number of data sets.
     */
    private static final int PM_MAX = 16;

    /**
     * Logger variable for this class.
     */
    private static Logger logger = Logger.getLogger(ShowPM.class);

    /**
     * Panel for the counter and region table.
     */
    private CounterDisplay myCounterDisplay;
    
    /**
     * This array contains all PM data sets that are used for comparison.
     * The value at position 0 is the main data set.
     */
    private PMData[] myPMDataSets = new PMData[PM_MAX]; 
    
    /**
     * This is the number of data sets currently in use.
     */
    private int numberDataSets;

    /**
     * myPMDataSets[selectedPMData] is the currently selected performance
     * data set. Currentlty, its value is always 0.
     */
    private int selectedPMData; // currently selected performance data
    
    /**
     * Currently selected processor for the TableDisplay.
     */
    private int selectedProcessor; 
    
    /**
     * Currently selected thread for the TableDisplay.
     */
    private int selectedThread;

    /**
     * Currently selected counter from the region counters. 
     * A value less 0 indicates that none is selected.
     */
    private int selectedRegionCounter;
    
    /**
     * Currently selected counter from the user counters. 
     * A value less 0 indicates that none is selected.
     */
    private int selectedUserCounter;
    
    /**
     * This frame is used to display a file that contains a
     * selected region.
     */
    private ShowFileFrame myShowFileFrame = null;
    
    /**
     * This frame is used to display counter values for all processors.
     */
    private ProcessDisplay myProcessDisplay = null;
    
    /**
     * This frame is used to display counter values from different data sets
     * to compare them.
     */
    private CompareDisplay myCompareDisplay = null;

    /**
     * This is the main frame of this GUI.
     */
    private JFrame myMainFrame;

    /**
     * A ShowPM editor is created by a filename for the file containg performance data.
     * 
     * @param pmFileName is the name of the file with performance data
     */
    public ShowPM(String pmFileName) {
    
        myPMDataSets[0] = new PMData(pmFileName);
        myPMDataSets[0].compressPM();
        numberDataSets = 1;
    
        selectedPMData = 0;
        selectedProcessor = 0;
        selectedThread = 0;
    
        //Schedule a job for the event-dispatching thread:
        //creating and showing this application's GUI.
    
        JFrame.setDefaultLookAndFeelDecorated(false);
    
        // create and set up the window
    
        myMainFrame = new JFrame("Performance results " + pmFileName);
    
        myMainFrame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    
        myCounterDisplay = new CounterDisplay(myPMDataSets[selectedPMData], this);
    
        setProcessor(0, 0);
    
        myShowFileFrame = new ShowFileFrame();
    
        myCounterDisplay.setOpaque(true);
        myMainFrame.setContentPane(myCounterDisplay);
        myMainFrame.pack();

        myMainFrame.setSize(SHOW_PM_DIMENSION);
    
        // frame.setSize (800, 400) would work here
    
        myMainFrame.setVisible(true);
        
        final double viewBoth = 0.5;
        
        setView(viewBoth);
    
    }

    /**
     * This routine shows a file in the show file frame.
     * 
     * @param fileDSP is the descriptor of the file to display.
     * @param start is the first line to highlight
     * @param stop is the last line to highlight
     */
    
    private void showFile(FileDescriptor fileDSP, int start, int stop) {
        
        myShowFileFrame.showFile(fileDSP, start, stop);
    }

    /**
     * This routine is used to show the display for all processes and
     * threads for a given counter.
     * 
     * @param isUser must be true if counter will be a user counter
     * @param counter is the integer index of the counter to be shown
     */
    public void showProcesses(boolean isUser, int counter) {

        if (myProcessDisplay == null) {

            // create a new Processor Frame for PM results

            myProcessDisplay = new ProcessDisplay(myPMDataSets[selectedPMData], this);  
            
        }

        // set the counter to be displayed in the Processor Frame

        myProcessDisplay.setCounter(isUser, counter);

        if (isUser) {

            selectedRegionCounter = ArraySelection.NO_SELECTION;
            selectedUserCounter = counter;

        } else {

            selectedRegionCounter = counter;
            selectedUserCounter = ArraySelection.NO_SELECTION;
        }

        myProcessDisplay.repaint();
        myProcessDisplay.setVisible(true);

    } // showProcesses

    /**
     * This routine is used to show the display for all data sets
     * for a given counter.
     * 
     * @param userFlag must be true if counter will be a user counter
     * @param counter is the integer index of the counter to be shown
     */
    public void showComparison(boolean userFlag, int counter) {
        
        // create the compare display if not available yet
        
        if (myCompareDisplay == null) {
            
            myCompareDisplay = new CompareDisplay(myPMDataSets, numberDataSets, this);            
        }

        myCompareDisplay.setCounter(userFlag, counter, selectedProcessor, selectedThread);

        if (userFlag) {

            selectedRegionCounter = ArraySelection.NO_SELECTION;
            selectedUserCounter = counter;

        } else {

            selectedRegionCounter = counter;
            selectedUserCounter = ArraySelection.NO_SELECTION;
        }

        myCompareDisplay.repaint();
        myCompareDisplay.setVisible(true);

    } // showProcesses

    /**
     * This routine is used to fix the selection of a certain processors 
     * (process, thread).
     * 
     * @param indexProc is the index of the process
     * @param indexThread is the index of the thread
     */
    public void setProcessor(int indexProc, int indexThread) {

        String procInfo = "";

        selectedProcessor = indexProc;
        selectedThread = indexThread;

        // PM_thread = PM[actual_PM].myPMCounterValues[selectedProcessor][selectedThread];

        int numberProcs = myPMDataSets[selectedPMData].getNumberProcesses();

        if (numberProcs > 1) {
            
            procInfo = "PID = " + indexProc + " of " + numberProcs;
        }

        int numberThreads = myPMDataSets[selectedPMData].getNumberThreads();

        if (numberThreads > 1) {
            
            if (procInfo.length() != 0) {
                
                procInfo = procInfo + ", ";               
            }
 
            procInfo = procInfo + "TID = " + indexThread + " of " + numberThreads;
        }

        myCounterDisplay.setTitleExtension(procInfo);
        myCounterDisplay.repaint();
    }

    /**
     * This routine selects the next processor (restarts with the
     * first one for the last processor).
     * 
     */
    public void nextProcess() {

        int indexProc = selectedProcessor + 1;
        
        if (indexProc >= myPMDataSets[selectedPMData].getNumberProcesses()) {
            
            indexProc = 0;            
        }

        setProcessor(indexProc, selectedThread);
    }

    /**
     * This routine selects the next thread (restarts with the first
     * one for the last thread).
     * 
     */
    public void nextThread() {

        int indexThread = selectedThread + 1;
        
        if (indexThread >= myPMDataSets[selectedPMData].getNumberThreads()) {
            
            indexThread = 0;      
        }

        setProcessor(selectedProcessor, indexThread);
    }

    /**
     * This routine allows to focus on one of the region or user table. A
     * small value focuses the user table, a larger value makes the region
     * table more visible.
     * 
     * @param divider must be a double value greater tan  0.0 and less than 1.0
     */
    public void setView(double divider) {

        myCounterDisplay.setView(divider);
        
    }

    /**
     * Routine asks for a new PM input file and puts it to the 
     * comparison data.
     * 
     */
    public void readCompareData() {

        // the user has to select a new PM file

        String workingDirectory = System.getProperty("user.dir");

        JFileChooser chooser = new JFileChooser(workingDirectory);

        chooser.setDialogTitle("Select PMS file");

        int result = chooser.showOpenDialog(myMainFrame);

        if (result == JFileChooser.CANCEL_OPTION) {
            
            return;     
        }

        File selectedFile = chooser.getSelectedFile();

        if (selectedFile == null) {
            
            return;       
        }

        if (numberDataSets >= PM_MAX) {

            String msg = "maximal number of PM files " + PM_MAX + " reached";
            
            showErrorMessage(msg);

            return;
        }

        String filename = selectedFile.getPath();

        logger.info("read comparison file " + filename);

        myPMDataSets[numberDataSets] = new PMData(filename);

        myPMDataSets[numberDataSets].compressPM();

        String compare = myPMDataSets[0].matchRegions(myPMDataSets[numberDataSets]);

        if (compare.equals("")) {

            // region matches

            numberDataSets += 1;

        } else {

            showErrorMessage(compare);

            // we force alignment with main PMData set

            myPMDataSets[numberDataSets].alignRegions(myPMDataSets[0]);

            numberDataSets += 1;
        }

    }

    /**
     * This routine is used to remove the latest data set that has been read.
     * 
     */
    public void removeCompareData() {

        logger.info("rm comparison file " + myPMDataSets[numberDataSets].getName());

        myPMDataSets[numberDataSets] = null;
        numberDataSets -= 1;

    }

    /**
     * {@inheritDoc}
     *
     * @see adaptor.ShowPM.SelectedInterface#getSelectedProcessor()
     */
    public int getSelectedProcess() {
        
        return selectedProcessor;
    }
    
    /**
     * {@inheritDoc}
     *
     * @see adaptor.ShowPM.SelectedInterface#getSelectedThread()
     */
    public int getSelectedThread() {
        
        return selectedThread;
    }
    
    /**
     * {@inheritDoc}
     *
     * @see adaptor.ShowPM.SelectedInterface#getSelectedUserCounter()
     */
    public int getSelectedUserCounter() {
        
        return selectedUserCounter;
    }
    
    /**
     * {@inheritDoc}
     *
     * @see adaptor.ShowPM.SelectedInterface#getSelectedRegion()
     */
    public int getSelectedRegionCounter() {
        
        return selectedRegionCounter;
    }

    /**
     * {@inheritDoc}
     *
     * @see adaptor.ShowPM.SelectedInterface#showFile(adaptor.ShowPM.PMData, int)
     */
    public void showFile(PMData data, int region) {
        
        RegionDescriptor regionDSP = data.getRegion(region);
        
        FileDescriptor fileDSP = regionDSP.getFile();

        logger.info("show file for region " + regionDSP.getName() 
                    + ": id = " + fileDSP.getFileId() + " name = " + fileDSP.getLongFileName());

        // might be better to show only one line

        showFile(fileDSP, regionDSP.getFirstLine(), regionDSP.getLastLine());
        
    }

    /**
     * {@inheritDoc}
     *
     * @see adaptor.ShowPM.SelectedInterface#infoFile(adaptor.ShowPM.PMData, int)
     */
    public void infoFile(PMData data, int regionIndex) {
        
        RegionDescriptor regionDSP = data.getRegion(regionIndex);
        
        StringBuffer info = new StringBuffer();
        
        regionDSP.appendInfo(info, "");
        
        ShowInfo.showInfoText(info.toString());
        
    }
    
    /**
     * This routine exports a region information file.
     * 
     */
    public void exportRIF() {
      
        if (numberDataSets == 1) {
            
            String currentDir = System.getProperty("user.dir");
            
            JFileChooser chooser = new JFileChooser(currentDir);
            
            chooser.setDialogTitle("Select RIF file to save:");
            
            int result = chooser.showSaveDialog(myMainFrame);
            
            if (result == JFileChooser.CANCEL_OPTION) {
                
                return;   
            }

            File selectedFile = chooser.getSelectedFile();
            
            // make sure that that really a file has been selected
            
            if (selectedFile == null) {
                
                return;         
            }

            myPMDataSets[0].writeRIF(selectedFile.getPath());
            
            
        } else {
            
            // show an error message that this is not allowed            
            
            JOptionPane.showMessageDialog(null, "Not allowed for multiple files", "Export RIF", JOptionPane.ERROR_MESSAGE);

        }
        
    }
    
    /**
     * {@inheritDoc}
     *
     * @see adaptor.ShowPM.SelectedInterface#counterSelection(boolean, int)
     */
    public void counterSelection(boolean isUser, int counter) {
        
        if (numberDataSets == 1) {
        
            // only one performance data set available

            showProcesses(isUser, counter);
  
        } else {
            
            showComparison(isUser, counter);
        }
     
    }

    /**
     * {@inheritDoc}
     *
     * @see adaptor.ShowPM.SelectedInterface#showErrorMessage(java.lang.String)
     */
    public void showErrorMessage(String msg) {
        
        JOptionPane.showMessageDialog(myMainFrame, msg, "ShowPM ERROR", JOptionPane.ERROR_MESSAGE);

    }
    
    /**
     * Main routine for ShowPM; allows to start this application.
     * 
     * @param arguments are the command line arguments
     */
    public static void main(String[] arguments) {
        
        final int errorExitCode = -1;
        
        String helpText =  "Results of Performance Monitoring";
        
        OptionsAdaptor theOptions = new OptionsAdaptor("ShowPM", helpText);
    
        // files can be specified without suffix
        
        // theOptions.setSuffix(".pms");
        
        // parse the arguments
        
        String[] inputFiles = null;
               
        try {
        
            inputFiles = theOptions.parseArguments(arguments);
                      
        } catch (IllegalArgumentException e) {
            
            // there might be illegal arguments; in this case we will stop
            
            logger.error(e.getMessage());
            
            theOptions.printUsage();
            
            System.exit(errorExitCode);
            
        }

        int nFiles = inputFiles.length;
        
        logger.info("main program starts, " + nFiles + " input files");
        
        for (int i = 0; i < nFiles; i++) {
            
            new ShowPM(inputFiles[i]);
            
        }
    
        logger.info("main program terminates, GUI thread remains");
    
    
    } // main 


} // class ShowPM

