/*
 * AboutAdaptorFrame.java
 *
 * Frame to display all relevant information about ADAPTOR GUI Project.
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

import java.io.IOException;
import java.net.URL;
import java.util.List;
import java.util.Vector;

import java.awt.image.BufferedImage;

import javax.imageio.ImageIO;

import org.apache.log4j.Logger;

import org.jfree.base.Library;
import org.jfree.ui.about.AboutFrame;
import org.jfree.ui.about.Contributor;
import org.jfree.ui.about.ProjectInfo;

/**
 * AboutAdaptorFrame is a general frame for all Java GUIs of the ADAPTOR project.
 * General information about all GUIs of the ADAPTOR project is kept here.
 *
 * @version $LastChangedRevision$
 * @author Thomas Brandes
 */
public class AboutAdaptorFrame
{

    /**
     * This is my logger variable to log infos via log4j.
     */
    private static Logger logger = Logger.getLogger( AboutAdaptorFrame.class );

    /**
     * Frame that contains all About information available.
     */
    private AboutFrame myFrame;

    /**
     * Constructor for the frame containg about information.
     *
     * @param name is the name of the project
     * @param info is the info about the project
     */
    public AboutAdaptorFrame( String name, String info )
    {

        final String gnuLesserLicence = "GNU Lesser General Public Licence 2.1";

        final String imgURL = "http://www.scai.fraunhofer.de/EP-CACHE/adaptor/www/pics/adaptor_logo.gif";

        // final String imgURL = "http://www.scai.fraunhofer.de/fileadmin/icons/scai_logo.gif";

        final String publicLicence = "Common Public Licence 1.0";

        ProjectInfo myInfo = new ProjectInfo();

        myInfo.setName( name );

        myInfo.setInfo( info );

        myInfo.setVersion( "1.0" );

        myInfo.setCopyright( "Copyright (C) 2006 Fraunhofer SCAI, Germany" );

        List<Contributor> authors = new Vector<Contributor>();

        Contributor aContributor = new Contributor( "Thomas Brandes", "Thomas.Brandes@scai.fhg.de" );

        authors.add( aContributor );

        myInfo.setContributors( authors );

        Library myL = new Library( "jcommon", "1.0.0", gnuLesserLicence, "General Utilities" );

        myInfo.addLibrary( myL );

        myL = new Library( "log4j", "1.2.13", "Apache Licence 2.0", "Logging for Java" );

        myInfo.addLibrary( myL );

        myL = new Library( "grappa", "1.2", publicLicence, "Graph Drawing Package" );

        myInfo.addLibrary( myL );

        myL = new Library( "jfreechart", "1.0.0", gnuLesserLicence, "Generation of Charts" );

        myInfo.addLibrary( myL );

        myInfo.setLicenceName( "Licence" );

        myInfo.setLicenceText( "Licence for ADAPTOR GUI has not been fixed yet. Please observe the copyright" );

        BufferedImage myLogo = null;

        try
        {

            URL myURL = new URL( imgURL );

            // other solution:
            // Toolkit t = Toolkit.getDefaultToolkit();
            // myLogo = t.createImage(myURL);

            myLogo = ImageIO.read( myURL );

        }
        catch ( IOException e )
        {

            logger.error( e );

        }

        if ( myLogo != null )
        {

            myInfo.setLogo( myLogo );

        }
        else
        {

            logger.error( "could not read logo (" + imgURL + ")" );

        }

        myFrame = new AboutFrame( name, myInfo );

        myFrame.setVisible( true );

    }

    /**
     * Shows or hides the AboutFrame.
     *
     * @param flag shows if true and hides if false
     */
    public void setVisible( boolean flag )
    {

        myFrame.setVisible( flag );
    }
}