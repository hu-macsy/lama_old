/*
 * FileData.java
 * 
 * An object of the class FileData contains all information about the available files.
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

package adaptor.General;

import org.apache.log4j.Logger;


/**
 * An object of the class FileTable contains all information about the available files.
 *
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */
public class FileTable {

    /**
     * This value will be the initial size of the file table.
     */
    private static final int INIT_NO_FILES = 10;
    
    /**
     * Logger for this class.
     */
    private static Logger logger = Logger.getLogger(FileTable.class);
    
    /**
     * This is the frame to show files. It will be reused for all files to have
     * not too many frames on the screen.
     */
    private static ShowFileFrame myShowFileFrame;

    /**
     * This is the array containing all regions. The region id is
     * used as an index for this array. The array itself can have holes
     * in it.
     */
    private FileDescriptor[] theFiles = null;
    
    /**
     * Constructor for FileTable with a given number of entries. This
     * constructor should be used if number of FileDescriptors is known.
     * 
     * @param noFiles is the number of entries.
     */
    public FileTable(int noFiles) {
        
        this.theFiles = new FileDescriptor[noFiles];
        
        for (int i = 0; i < noFiles; i++) {
            
            this.theFiles[i] = null;
        }
    }

    /**
     * Constructor for FileTable with a default number of entries.
     * 
     */
    public FileTable() {
    
        this(INIT_NO_FILES);
    }
    
    /**
     * This routine takes an already existing array of region descriptors.
     * 
     * @param inFiles is the array of known region descriptors.
     */
    public FileTable(FileDescriptor[] inFiles) {

        theFiles = inFiles;

    } // setRegions 

    /**
     * This routine reallocates the region data in such a way that an entry
     * for a certain region id can be made.
     * 
     * @param nr is the region id for which an entry in the array is needed
     */
    private void extendFiles(int nr) {

        int noFiles = theFiles.length;
        
        while (nr >= noFiles) {

            // logger.error("get_region: maximal number of regions reached");
            // logger.info("double region table to " + 2*MaxRegions  + " entries");

            // double number of Regions

            FileDescriptor[] newFiles = new FileDescriptor[2 * noFiles];

            for (int i = 0; i < noFiles; i++) {
                
                newFiles[i]             = theFiles[i];
                newFiles[i + noFiles] = null;
            }

            noFiles = 2 * noFiles;
            theFiles = newFiles;
        }

    } // extend_regions

    /**
     * This routine returns the region descriptor by the region id.
     * 
     * @param index is the index of the region
     * @return the descriptor of the region
     */
    public FileDescriptor getFile(int index) {

        extendFiles(index);

        return theFiles[index];
    }

    /**
     * This routine defines a file descriptor at a given position.
     * 
     * @param index is the index for the region descriptor in the table
     * @param file is the FileDescriptor to be defined
     */
    public void defineFile(int index, FileDescriptor file) {
        
        // make sure that the region array is big enough
        
        extendFiles(index);

        if (theFiles[index] == null) {
            
            theFiles[index] = file;
            
        }

        theFiles[index] = file;
                  
    }
    
    /**
     * This routine returns a file descriptor by the given 
     * file identification and the file name. If the file 
     * identification is already in use, the name of the files
     * must be identical.
     * 
     * @param nr is the fileId and position in file descriptor array
     * @param name is the name of the file
     * @return the file descriptor
     */
    public FileDescriptor getFile(int nr, String name) {

        extendFiles(nr);

        if (theFiles[nr] == null) {

            theFiles[nr] = new FileDescriptor(nr, name);
            return theFiles[nr];

        }

        // file is defined again but then it should have the same name

        String sname = FileDescriptor.getName(name);

        if (theFiles[nr].getShortFileName().equals(sname)) {
            
            return theFiles[nr];     
        }

        logger.error("The file " + nr + " already defined (was " + theFiles[nr].getShortFileName() + ", now: " + name);
        
        return null;

    }

    /**
     * Defines a new file descriptor at a free position.
     * 
     * @param fileName is the long name of the file
     * @return the exsiting or the new FileDescriptor
     */
    public FileDescriptor getFile(String fileName) {
        
        int pos = ArraySelection.NO_SELECTION;
        
        /* linear search for an existing file descriptor */
        
        for (int i = 0; i < theFiles.length; i++) {
            
            if (theFiles[i] != null) {
                
                if (theFiles[i].isSameFile(fileName)) {
                    
                    logger.info("File " + fileName + " has file id " + theFiles[i].getFileId());
                    
                    return theFiles[i];
                    
                }
                
            } else if (pos < 0) {
                
                // we have found a free position and we might use it
                
                pos = i;
            }
        }
        
        /* so we have a new file, make a new descriptor */
        
        if (pos < 0) {
            
            pos = theFiles.length;
            
        }

        logger.info("File " + fileName + " gets new file id " + pos);
        
        return getFile(pos, fileName);

    }
    
    /**
     * This routine return the highest index for which a FileDescriptor
     * is available. Please not that there might be null descriptors at
     * some positions before.
     * 
     * @return the highest index for which a FileDescriptor is available.
     */
    public int noFiles() {
        
        int n = 0;
        
        if (theFiles != null) {
            
            n = theFiles.length;
            
            // there might be unused regions at the end of the table
            
            for (int k = n - 1; k >= 0; k--) {
                
                if (theFiles[k] == null) {
                    
                    n = k;
                    
                } else {
                    
                    break;
                }
            }
        }
        
        return n;
    }
    
    /**
     * Getter routine for the array of files.
     * 
     * @return array with all file descriptors
     */
    public FileDescriptor[] getFiles() {

        return theFiles;
    }
    
    /**
     * Display a file in a FileFrame and highlight from line1 to line2.
     * 
     * @param dspFile is the file identification
     * @param line1 is the first line to highlight
     * @param line2 is the last line to highlight
     */
    public static void showFile(FileDescriptor dspFile, int line1, int line2) {

        if (myShowFileFrame == null) {
            
            myShowFileFrame = new ShowFileFrame();
        }

        myShowFileFrame.showFile(dspFile, line1, line2);

    } // showFile

}
