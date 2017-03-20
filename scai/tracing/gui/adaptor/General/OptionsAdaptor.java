/*
 * OptionsAdaptor.java
 *
 * Utilities to deal with command line arguments for the ADAPTOR GUI
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

import java.io.File;
import java.util.List;
import java.util.Vector;

import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.apache.log4j.PropertyConfigurator;

/***
 * OptionsAdaptor is a class that evaluates the command line arguments.
 *
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */
public class OptionsAdaptor
{

    /**
     * String with four blanks to fill infos.
     */
    private static final String FILL4 = "     ";

    /**
     *  Maximal number of characters for an option. This value is
     *  used to align all options when printing their usage.
     */
    private static final int MAX_ARG_LEN = 10;

    /**
     * Logger for this class.
     */
    private static Logger logger = Logger.getLogger( OptionsAdaptor.class );

    /**
     * This variable contains all flag options, e.g. -root, -rif, ...
     */
    private List<String> options1 = new Vector<String>();

    /**
     * This variable contains all options with an additional argument,
     * e.g. -output file, -log file, -group kind.
     */
    private List<String> options2 = new Vector<String>();

    /**
     * These are the descriptions for options1.
     */
    private List<String> help1    = new Vector<String>();

    /**
     * These are the descriptions for options2.
     */
    private List<String> help2    = new Vector<String>();

    /**
     * This variable contains for every possible option with an additional
     * argument a possibility to check it.
     */
    private List<String[]> vals2    = new Vector<String[]>();

    /**
     * After parsing the command line arguments this variable shows
     * which of the options has been set.
     */
    private boolean [] isSet1  = null;

    /**
     * This variable contains for all args in option2 the value
     * of this option after parsing the arguments.
     */
    private String  [] setVals = null;

    /**
     * This variable contains the program name used for printUsage.
     */
    private String programName = null;

    /**
     * This variable contains a description text for the program.
     */
    private String programDescription = null;

    /**
     * Integer variable to count the number of errors.
     */
    private int errors = 0;

    /**
     * String variable to which we append error messages.
     */
    private String errorString = null;

    /**
     * Create a new instance to allow parsing arguments.
     *
     * @param name is the name of the program
     * @param helpText is a description of the program
     */
    public OptionsAdaptor( String name, String helpText )
    {

        programName = name;
        programDescription = helpText;

        // These are the default options we have with each application using this class

        addValueOption( "-log4j", "configuration file for log4j" );

        addFlagOption( "-h", "prints this usage info" );
        addFlagOption( "-info", "set logger level to INFO" );
        addFlagOption( "-debug", "set logger level to DEBUG" );
        addFlagOption( "-nowarn", "set logger level to ERROR (no more WARNings)" );
    }

    /**
     * Help routine to print usage information. Currently this
     * is done via System.out, but we might change it here later
     * to a more convenient solution.
     *
     * @param msg is the line of usage to print
     */
    private static void myPrintUsage( String msg )
    {

        System.out.println( msg );
    }

    /**
     * This routine defines a new simple option for the program.
     *
     * @param option is the name of the option, e.g. -fast, -root
     * @param help is a description of the option
     */
    public void addFlagOption( String option, String help )
    {

        options1.add( option );
        help1.add( help );

    }

    /**
     * This routine defines a new option that follows a second value argument.
     *
     * @param option is the name of the option, e.g. -I xxx
     * @param help is a description of the option
     */
    public void addValueOption( String option, String help )
    {

        options2.add( option );
        help2.add( help );
        vals2.add( null );
    }

    /**
     * This routine defines a new option that follows a second argument.
     *
     * @param option is the name of the option, e.g. -I xxx
     * @param help is a description of the option
     * @param possibleVals is an array of strings for allowed values of xxx
     */
    public void addValueOption( String option, String help, String[] possibleVals )
    {

        options2.add( option );
        help2.add( help );
        vals2.add( possibleVals );
    }


    /**
     * Help routine to print a flag option.
     *
     * @param option is the name of the flag option
     * @param help is a description of the effect for the option
     */
    private static void printFlagOption( String option, String help )
    {

        String info = FILL4 + option + FILL4;

        for ( int i = 0; i < MAX_ARG_LEN - option.length(); i++ )
        {

            info += " ";
        }

        info += ": " + help;

        myPrintUsage( info );
    }


    /**
     * This routine prints info about a value option.
     *
     * @param option is the name of the option
     * @param help is the description
     * @param vals is an object describing restrictions
     */
    private static void printValueOption( String option, String help, Object vals )
    {

        String info = FILL4 + option + " arg ";

        String valHelp = "";

        if ( vals instanceof String [] )
        {

            String [] possibleValues = ( String [] ) vals;

            valHelp = "" + possibleValues[0];

            for ( int i = 1; i < possibleValues.length; i++ )
            {

                valHelp += "," + possibleValues[i];

            }

            valHelp = " [" + valHelp + "]";

        }

        for ( int i = 0; i < MAX_ARG_LEN - option.length(); i++ )
        {

            info += " ";

        }

        info += ": " + help + valHelp;

        myPrintUsage( info );
    }

    /**
     * This routine prints the usage of the program. It should be called after
     * all options have been defined.
     *
     */
    public void printUsage()
    {

        myPrintUsage( "Usage of " + programName + " (" + programDescription + ")" );
        myPrintUsage( "" );

        myPrintUsage( programName + " [options] files" );
        myPrintUsage( "" );
        myPrintUsage( "  These options are supported:" );
        myPrintUsage( "" );

        for ( int i = 0; i < options1.size(); i++ )
        {

            printFlagOption( ( String ) options1.get( i ), ( String ) help1.get( i ) );

        }

        myPrintUsage( "" );

        for ( int i = 0; i < options2.size(); i++ )
        {

            printValueOption( ( String ) options2.get( i ), ( String ) help2.get( i ), vals2.get( i ) );

        }

    }

    /**
     * This routine is used to get the value of an option with a second argument.
     *
     * @param option is the name of the option (used for error messages)
     * @param oldVal if not null this is an existing value
     * @param args is the argument list
     * @param index is the index where option appears
     * @return the value of the option
     */
    private String getOption2( String option, String oldVal, String[] args, int index )
    {

        String newVal = oldVal;

        if ( oldVal != null )
        {

            logger.error( "option " + option + " already used, val = " + oldVal );
        }

        // now make sure that the index is valid

        if ( index >= args.length )
        {

            logger.error( "missing second argument for " + option );

        }
        else
        {

            newVal = args[index];
        }

        return newVal;

    }

    /**
     * This routine checks whether an option with a value is a legal
     * option.
     *
     * @param option is the name of the option, e.g. "-file"
     * @param val is the value of the option (as string)
     * @param checkVal can be an array of strings or a file
     * @return true if val is a legal value for option
     */
    private boolean legalOption2( String option, String val, Object checkVal )
    {

        boolean isLegal = true;

        if ( checkVal instanceof String[] )
        {

            isLegal = false;

            String [] possibleValues = ( String [] ) checkVal;

            for ( int i = 0; i < possibleValues.length; i++ )
            {

                if ( val.equals( possibleValues[i] ) )
                {

                    isLegal = true;
                    break;
                }
            }
        }

        if ( checkVal instanceof File )
        {

            logger.info( "check if file exists" );

        }

        return isLegal;

    }

    /**
     * Routine resets all errors.
     *
     */
    private void resetErrors()
    {

        errors = 0;
        errorString = null;
    }

    /**
     * This routine increasses the error counter and adds an
     * error message to the error string.
     *
     * @param msg is the new error message.
     */
    private void addError( String msg )
    {

        errors++;

        if ( errorString == null )
        {

            errorString = msg;

        }
        else
        {

            errorString += ", " + msg;

        }
    }

    /**
     * After parsing the arguments this routine can be called to check
     * for such options that belong to this class itself.
     *
     */
    private void evalMyOptions()
    {

        if ( isOptionSet( "-h" ) )
        {

            printUsage();

            // we exit the application (not an error at all)

            System.exit( 0 );
        }

        if ( isOptionSet( "-debug" ) )
        {

            Logger.getRootLogger().setLevel( Level.DEBUG );

        }
        else if ( isOptionSet( "-info" ) )
        {

            Logger.getRootLogger().setLevel( Level.INFO );

        }
        else if ( isOptionSet( "-nowarn" ) )
        {

            Logger.getRootLogger().setLevel( Level.ERROR );
        }

        String logConfigFile = getOptionVal( "-log4j" );

        if ( logConfigFile != null )
        {

            logger.info( "configure now with file " + logConfigFile );
            PropertyConfigurator.configure( logConfigFile );
        }

    }

    /**
     * This routine parses all the arguments and checks for known
     * options. The known options are removed from the argument list
     * and can be asked for afterwards.
     *
     * @param args are the command line arguments
     * @return all arguments that are not options
     */
    public String [] parseArguments( String [] args )
    {

        List nonOptionArgs = new Vector();

        // we count the errors and then we might throw an Exception

        resetErrors();

        if ( isSet1 != null )
        {

            throw new IllegalArgumentException( "parseArguments should only be called once" );
        }

        isSet1 = new boolean [options1.size()];

        for ( int k = 0; k < isSet1.length; k++ )
        {
            isSet1[k] = false;
        }

        setVals = new String [options2.size()];

        for ( int k = 0; k < setVals.length; k++ )
        {

            setVals[k] = null;
        }

        boolean ignoreNext = false;

        for ( int i = 0; i < args.length; i++ )
        {

            if ( ignoreNext )
            {

                ignoreNext = false;

                continue;
            }

            // check for a simple flag option, e.g. -root, -rif

            boolean isOption = false;

            for ( int k = 0; k < options1.size(); k++ )
            {

                if ( args[i].startsWith( ( String ) options1.get( k ) ) )
                {

                    isOption = true;
                    isSet1[k] = true;
                    break;
                }
            }

            if ( isOption )
            {

                continue;

            }

            // check for an argument option, e.g. -e executable, -o outfile

            for ( int k = 0; k < options2.size(); k++ )
            {

                String option = ( String ) options2.get( k );

                if ( args[i].startsWith( option ) )
                {

                    isOption = true;

                    setVals[k] = getOption2( option, setVals[k], args, i + 1 );

                    if ( !legalOption2( option, setVals[k], vals2.get( k ) ) )
                    {
                        errors++;
                    }

                    break;
                }
            }

            if ( isOption )
            {

                ignoreNext = true;
                continue;

            }

            // argument was not an option, so we take it as input argument

            nonOptionArgs.add( args[i] );
        }

        // before we go on, we check for our own options

        evalMyOptions();

        String [] files = new String [nonOptionArgs.size()];

        for ( int i = 0; i < nonOptionArgs.size(); i++ )
        {

            files[i] = ( String ) nonOptionArgs.get( i );

            File testFile = new File( files[i] );

            if ( !testFile.exists() )
            {

                String msg = "argument " + files[i] + " is not an existing file";

                addError( msg );

            }
            else
            {

                logger.info( "argument " + files[i] + ": file exists" );
            }
        }

        if ( errorString != null )
        {

            throw new IllegalArgumentException( errorString );
        }

        return files;

    }

    /**
     * This routine checks if a single option has been set in command line.
     *
     * @param option is the name of the opiton
     * @return true if option has been used
     */
    public boolean isOptionSet( String option )
    {

        boolean isSet = false;
        boolean found = false;

        for ( int i = 0; i < options1.size(); i++ )
        {

            if ( options1.get( i ).equals( option ) )
            {

                found = true;
                isSet = isSet1[i];
                break;

            }
        }

        if ( !found )
        {

            logger.error( "illegal use of isOptionSet, " + option + " has never been defined." );
        }

        return isSet;
    }

    /**
     * Search for an value option the corresponding value.
     *
     * @param option is the option for which the value is searched
     * @return the value, is null if not found
     */
    public String getOptionVal( String option )
    {

        String val = null;

        boolean found = false;

        for ( int i = 0; i < options2.size(); i++ )
        {

            if ( options2.get( i ).equals( option ) )
            {

                val = setVals[i];
                found = true;
                break;

            }
        }

        if ( !found )
        {

            logger.error( "illegal use of getOptionVal, " + option + " has not been defined." );
        }

        return val;
    }

    static
    {

        // make a basic configuration of the logger

        BasicConfigurator.configure();

        // set the default level for logging to WARN (no INFO, no DEBUG)

        Logger.getRootLogger().setLevel( Level.WARN );
    }

} // class OptionsAdaptor

