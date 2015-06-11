/*
 * FileDescriptor.java
 * 
 * FileDescriptor is a descriptor for a file
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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import org.apache.log4j.Logger;

/**
 * Class contains attributes and methods for description of a file.
 *
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */
public class FileDescriptor {

    /**
     * Constant value that determines the number of char used for the line
     * numbers. 
     */
    private static final int NO_CHARS_LINENUMBER = 6;

    /**
     * Logger for this class.
     */
    private static Logger logger = Logger.getLogger(FileDescriptor.class);

    /**
     * internal number of file in RIF in Calltree File.
     */
    private int myFileId;
    
    /**
     * This is the file denoted by this descriptor.
     */
    private File myFile = null;
    
    /**
     * content will be set if file has been opened for viewing.
     */
    private StringBuffer myContent = null; /* content of file if opened */
    
    /**
     * The array myLines is used to translate line numbers to character
     * positions in the file for which the descriptor stands.
     * myLines [i-1] is the position in the file where line i starts. 
     */
    private int[] myLines;
    
    /**
     * This constructor creates a file descriptor from an input line.
     * Its shape must be "file=id  name=name path=path.
     * 
     * @param theLine is the input line
     * @throws IOException in case of an illegal input line
     */
    public FileDescriptor(String theLine) throws IOException {
        
        try {
            
            myFileId = Utilities.readInt(theLine, "file");
            
            String name = Utilities.readString(theLine, "name");
            String path = Utilities.readString(theLine, "path");
            
            myFile = new File(path + File.separator + name);
            
        } catch (IOException e) {
            
            String msg = "Cannot create file descriptor from this line ("
                         + e.getMessage() + "): " + theLine;
            
            logger.error(msg);
            
            throw e;
            
        }
        
    } // FileDescriptor
    
    /**
     * construct a file descriptor with the full file name.
     * 
     * @param id is the file id for the file
     * @param name is the path name of the file
     */
    public FileDescriptor(int id, String name) {
        
        // what happens if the file does not exist ?
        
        myFileId   = id;
        myFile     = new File(name);

    } // FileDescriptor as used for ShowPM
    
    /**
     * Returns the name of the file denoted by this pathname. This is just the last name in the
     * pathname's name sequence.
     * 
     * @param pathname
     *            is the full pathname
     * @return the last name in the pathnames seqeuence
     */
    public static String getName(String pathname) {
        
        File testFile = new File(pathname);
        return testFile.getName();
        
    }
    
    /**
     * getter routine for the file identification.
     * 
     * @return the number of the file
     */
    public int getFileId() {
        
        return myFileId;
    }
    
    /**
     * This routine return the (restricted) path of the file descriptor.
     * Let /x1/x2/x3/x4 be the path, then getFilePath(1,1) returns x4,
     * getFilePath(1,2) returns x3/x4.
     * 
     * @param n1 is the last item in path counted from behind
     * @param n2 is the first item in path counted from behind
     * @return the (restricted) path
     */
    public String getFilePath(int n1, int n2) {
        
        String path = myFile.getParent();
        
        if ((n1 > 0) && (n1 <= n2)) {
            
            // this implies the selection of a sub path
            
            String[] pathItems = path.split(File.pathSeparator);
            
            int n = pathItems.length;
            
            int k1 = Math.min(n1, n);
            int k2 = Math.min(n2, n);
            
            path = pathItems[n - k2];
            
            for (int i = n - k2 + 1; i < n - k1; i++) {

                path += File.separator + pathItems[i];

            }             
        }        
        
        return path;
        
    }
    
    /**
     * This routine returns the short name of the file.
     * 
     * @return string containing the short file name
     */
    public String getShortFileName() {
        
        return myFile.getName();
    }
    
    /**
     * This routine returns the name of the file with its full path.
     * 
     * @return string containing the full file name
     */
    public String getLongFileName() {
        
        return myFile.getAbsolutePath();
    }
    
    /**
     * This routine tries to find a file with a given name in a directory.
     * 
     * @param dir is the existing directory
     * @param fname is the name of the file
     * @return the File if found and exists, null otherwise
     */
    private static File searchDirectory(File dir, String fname) {
        
        String fullName = dir.getPath() + File.separator + fname;
        
        logger.debug("try for " + fullName);
        
        File foundFile = new File(fullName);
        
        if (foundFile.exists()) {
            
            logger.info("file exists: " + foundFile.getPath());
            logger.info("absolute " + foundFile.getAbsolutePath());
            
        } else {
            
            // if not found in this directory we try other ones   
            
            foundFile = null;

            File[] dirEntries = dir.listFiles();
            
            for (int i = 0; i < dirEntries.length; i++) {
                
                File f = dirEntries[i];
                
                if (f.isDirectory()) {
                    
                    f = searchDirectory(f, fname);
                    
                    if (f != null) {
                        
                        foundFile = null;
                        break;             
                    }
                     
                } // if subdirectory
                
            } // for all entries in the directory
     
        }
        
        return foundFile;
        
    }
    
    /**
     * This routine tries to find a file with a given name in a directory.
     * 
     * @param dir is the directory
     * @param fname is the name of the file
     * @return the File if found and exists, null otherwise
     */
    public static File searchFile(File dir, String fname) {
        
        File foundFile = null;
    
        if (dir != null) {
            
            // okay, dir is relly a file
            
            if (dir.isDirectory()) {
                
                foundFile = searchDirectory(dir, fname);
                
            }
        }

        return foundFile;
        
    } // SearchFile
    
    /**
     * OpenFile returns File for this descriptor. It searches also the file
     * in all other subdirectories of the current subdirectory.
     * 
     * @return the File descriptor (can be null if not found)
     */
    File openFile() {
        
        if (myFile.exists()) {
            
            return myFile;            
        }
        
        logger.info("could not open file, try else");
        
        return searchDirectory(new File("."), myFile.getName());
        
    } // OpenFile
    
    /**
     * This routine provides a string with all descriptor data.
     * 
     * @return a line containing all file descriptor information
     */
    public String makeFileString() {
        
        return "file=" + myFileId + " name=" + myFile.getName() + " path=" + myFile.getParent();
    }

    /**
     * Missing file descriptors for output files can be generated by this routine.
     * 
     * @param fileId
     *            is the file identification for which no file is available
     * @return a string for a dummy file descriptor
     */
    public static String makeDummyString(int fileId) {

        return "file=" + fileId + " name=? path=?";

    }

    /**
     * This routine returns the content of the file in a long 
     * string.
     * 
     * @return String containing the full content of the file
     * @throws IOException in case of IO problems
     */
    public String getContent() throws IOException {
        
        final int maxLines = 65536;
        
        // return if content is already available
        
        if (myContent != null) {

            return myContent.toString();
            
        }
        
        String inputLine;
        
        int theLineCounter = 0; /* counts lines in file */
        
        int[] tmpLinePositions = new int[maxLines];
        
        File f = openFile();
        
        String fname = myFile.getCanonicalPath();
        
        if (f != null) {
            
            fname = f.getCanonicalPath();
            
        }
        
        FileReader file = new FileReader(fname);
        BufferedReader buff = new BufferedReader(file);
        
        logger.info("opened file " + fname + ", will get its content");
        
        myContent = new StringBuffer();
        
        StringBuffer buffLineNumber = new StringBuffer("");
        
        for (int i = 0; i < NO_CHARS_LINENUMBER; i++) {

            buffLineNumber.append(" ");
        }

        try {
            
            inputLine = buff.readLine();
            
        } catch (IOException e) {
            
            myContent = null;
            
            throw new IOException(fname + "(read error in line 1)");
        }
        
        while (inputLine != null) { 
            
            // put the line number in the string buffer
            
            String lineNr = "" + (theLineCounter + 1);
            
            int len = lineNr.length();
            int pos = NO_CHARS_LINENUMBER - 1;
            buffLineNumber.replace(pos - len, pos, lineNr);
            
            myContent.append(buffLineNumber.toString() + inputLine + "\n");
            
            tmpLinePositions[theLineCounter] = myContent.length();
            
            theLineCounter++;
            
            // sun.io.MalformedException can be possible in lines
            
            try {
                
                inputLine = buff.readLine();
                
            } catch (IOException e) {
                
                myContent = null;
                
                throw new IOException(fname + "(read error in line " + theLineCounter + ")");
            }
            
        }
        
        // logger.info ("file " + file_name + " has "
        // + no_lines + " Lines");
        
        myLines = new int[theLineCounter];
        
        for (int i = 0; i < theLineCounter; i++) {

            myLines[i] = tmpLinePositions[i];
            
        }
        
        return myContent.toString();
    }
    
    /**
     * This routine returns the file pos of a line in a ext file.
     * 
     * @param line is the number of line we search for the index
     * @return the line index, error value is -1
     */
    public int getLineIndex(int line) {
        
        boolean error = false;
        
        final int errorVal = -1;
        
        if (line < 1) {
            
            error = true;
            
        } else if (myLines == null) { 
            
            // logger.info ("getLineIndex: Lines == null");
            
            error = true;
            
        } else if (line > myLines.length) { 
            
            // logger.error ("Line " + line + " is out of range, file " +
            //                      file_name + " has only " + Lines.length + " lines");
            
            error = true;
        }
        
        if (error) {
            return errorVal;
        } else {
            return myLines[line - 1];
        }
        
    }
    
    /**
     * Verify if this file is identical to a file specified by a name.
     * 
     * @param name
     *            is the name of the other file
     * @return true if both files are the same
     */
    public boolean isSameFile(String name) {

        File f = new File(name);

        if (f == null) { 
            
            return false;
        }
        
        // logger.info("check if this File Descriptor (id=" + myFileId + ", name=" + myFileName + ",
        // path=" + myFilePath + ") is same as " + name);

        if (myFile.getAbsolutePath().equals(f.getAbsolutePath())) {

            return true;
        }

        return false;

    }


} // class FileDescriptor

