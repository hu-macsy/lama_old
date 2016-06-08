#!/usr/bin/env python

#
#  @file PyRename.py
# 
#  @license
#  Copyright (c) 2009-2016
#  Fraunhofer Institute for Algorithms and Scientific Computing SCAI
#  for Fraunhofer-Gesellschaft
# 
#  This file is part of the Library of Accelerated Math Applications (LAMA).
# 
#  LAMA is free software: you can redistribute it and/or modify it under the
#  terms of the GNU Affero General Public License as published by the Free
#  Software Foundation, either version 3 of the License, or (at your option)
#  any later version.
# 
#  LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
#  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
#  FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for
#  more details.
# 
#  You should have received a copy of the GNU Affero General Public License
#  along with LAMA. If not, see <http://www.gnu.org/licenses/>.
#  @endlicense
# 
#  @brief Python script that updates all files with a new license file            
#  @author: Thomas Brandes
#  @date 29.03.2016
#

import os
import os.path
 
##  Configuration variables to set

##  File with current LAMA license to set in each file

LICENSE_FILE = "/home/brandes/LAMA_License.txt"

##  Files with one of these suffixes will never get a license header

SKIP_FILE_EXTENSIONS = [ "", ".rst" ]

##  Files with one of these suffixes will always get a license header

FORCE_FILE_EXTENSIONS = [ ".cpp", ".cu", ".hpp", ".sh" ]

## Note: all other files get a license header if @license is found as keyword

## Keywords that might appear in the header

LICENSE_KEYWORDS = [ "@brief", "@author", "@file", "@date", "@license", "Created on:", "Author:" ]

#   Read in the current license file as array of lines

def readLicense():

   license = []
   for line in open( LICENSE_FILE, 'r' ):
      license = license + [line] 
   return license

def endHeader( commChar ):

    if ( commChar == "*" ):
        return "*/" 
    else:
        return commChar + commChar + commChar 


def writeInd( outfile, commChar, txt ):

    if ( txt == "\n" ):
        outfile.write( " " + commChar + txt )
    else:
        outfile.write( " " + commChar + " " + txt )

def writeIndNL( outfile, commChar, txt ):
 
    if ( txt == "" ):
        outfile.write( " " + commChar + "\n" )
    else:
        outfile.write( " " + commChar + " " + txt + "\n" )

def writeIndBegin( outfile, commChar ):

    if ( commChar == "*" ):
        outfile.write( "/**\n" )
    else:
        outfile.write( commChar + commChar + commChar + "\n" )

def writeIndEnd( outfile, commChar ):

    if ( commChar == "*" ):
        outfile.write( " */\n" )
    else:
        outfile.write( commChar + commChar + commChar + "\n" )

def writeLicenseHeader( outfile, fileinfo ):

    commChar = fileinfo[ "commChar" ]

    writeIndBegin( outfile, commChar )
    writeIndNL( outfile, commChar, "@file " + fileinfo["@file"] )
    writeIndNL( outfile, commChar, "" )
    writeIndNL( outfile, commChar, "@license" )

    license = readLicense()

    for line in  license:
        writeInd( outfile, commChar, line )

    writeIndNL( outfile, commChar, "@endlicense" )
    writeIndNL( outfile, commChar, "" )

    brief = fileinfo[ "@brief" ].split( "|" )
    writeIndNL( outfile, commChar, "@brief " + brief[0] )
    for line in brief[1:]:
       writeIndNL( outfile, commChar, "       " + line )

    writeIndNL( outfile, commChar, "@author " + fileinfo["@author"] )
    writeIndNL( outfile, commChar, "@date " + fileinfo["@date"] )
    writeIndEnd( outfile, commChar )


monthIds = [ "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec" ]

def parseDate( datestring ):

    d = {}
    
    if datestring[2] == ".":  
       #  dd.mm.yyyy
       d["day"] = int( datestring[0:2] )
       d["month"] = int( datestring[3:5] )
       d["year"] = int( datestring[6:10] )
    elif datestring[4] == "-":
       #  yyyy-mm-dd
       d["day"] = int( datestring[8:10] )
       d["month"] = int( datestring[5:7] )
       d["year"] = int( datestring[0:4] )
    elif datestring[6] == ",":
       #  Mar 16, 2015
       d["day"] = int( datestring[4:6] )
       d["month"] = monthIds.index( datestring[0:3] ) + 1
       d["year"] = int( datestring[8:12] )
    else:
       print "Error: unknown date, ", datestring
       d["day"] = 16
       d["month"] = 3
       d["year"] = 2015

    # print "converted date ", datestring, " is " , d

    return d

def unparseDate( date ):
   
    if ( date["day"] < 10 ):
       day = "0" + str( date["day"] )
    else :
       day = str( date["day"] )

    if ( date["month"] < 10 ):
       month = "0" + str( date["month"] )
    else :
       month = str( date["month"] )

    year = str( date["year"] )

    # return day + "." + month + "." + year

    return monthIds[date["month"] - 1] + " " + day + ", " + year

#  Method that returns true if filenames are equal at the end

def nameMatch( name1, name2 ):

    if ( name1.endswith( name2 )):
       return True
    if ( name2.endswith( name1 )):
       return True

    return False

def checkFilename( fileinfo ):

    filename = fileinfo[ "filename" ]

    # remove leading ./

    if ( filename.startswith( "./" ) ):
       filename = filename[ 2: ]

    file_only, file_extension = os.path.splitext( filename )

    # check if file is a configurtion file  suffix starts with "in" 

    if ( file_extension.startswith( ".in" ) ):

       filename = file_only

    # check for correct filename

    if ( "@file" in fileinfo ):

       if ( not nameMatch( filename, fileinfo["@file"] ) ):

          print "ERROR: ", fileinfo[ "filename" ], ", illegal @file = ", fileinfo[ "@file" ], " no match to ", filename
          fileinfo[ "@file"] = filename
    else:

       print "WARN: ", fileinfo[ "filename" ], "set missing @file = ", filename
       fileinfo[ "@file"] = filename

def checkBrief( fileinfo ):

    if  ( "@brief" in fileinfo ):

       brief = fileinfo[ "@brief" ]

       if brief.endswith( "|" ):
          brief = brief[ : len(brief) - 1 ]
          print "WARN: ", fileinfo[ "filename" ], ", remove empty line in @brief"
          fileinfo[ "@brief" ] = brief

    else:

       print "WARN: ", fileinfo[ "filename" ], "missing @brief"
       fileinfo[ "@brief" ] = "ToDo: Missing description in " + fileinfo[ "filename" ]

def checkAuthor( fileinfo ):

    if "Author:" in fileinfo:

       if "@author" in fileinfo:

           print "ERROR: ", fileinfo[ "filename" ], "Author: and @author specified"

       else:
           fileinfo[ "@author" ] = fileinfo[ "Author:" ]

    if  ( not "@author" in fileinfo ):

        if "git_author" in fileinfo:
            fileinfo[ "@author" ] = fileinfo[ "git_author" ]
            print "INFO ", fileinfo[ "filename" ], "@author = ", fileinfo["@author"], "  from git"
        else:
            print "ERROR: ", fileinfo[ "filename" ], "no @author found"
            fileinfo[ "@author" ] = "LAMA Team, SCAI Fraunhofer"

def checkDate( fileinfo ):

    if "Created on:" in fileinfo:

       if "@date" in fileinfo:

           print "ERROR: ", fileinfo[ "filename" ], "Created on: and @date specified"

       else:
           fileinfo[ "@date" ] = fileinfo[ "Created on:" ]

    if  ( not "@date" in fileinfo ):

        if "git_date" in fileinfo:
            fileinfo[ "@date" ] = fileinfo[ "git_date" ]
            print "INFO ", fileinfo[ "filename" ], "@date = ", fileinfo["@date"], "  from git"
        else:
            print "ERROR: ", fileinfo[ "filename" ], "no @date found"
            fileinfo[ "@date" ] = "16.03.2016"

    fileinfo[ "@date" ] = unparseDate( parseDate( fileinfo[ "@date" ] ) ) 

def checkStartLine( fileinfo ):

    if  ( not "startLine" in fileinfo ):
        print "WARN: ", fileinfo[ "filename" ], "no start line found, assume no header lines"
        fileinfo[ "startLine" ] = 0
        return

    for key in LICENSE_KEYWORDS:
        if key in fileinfo: 
           return
 
    # no key found in the file then start line must be 0

    fileinfo[ "startLine" ] = 0

def checkCommentChar( fileinfo ):

    if  ( not "commChar" in fileinfo ):

        file_only, file_extension = os.path.splitext( fileinfo["@file"] )

        if file_extension in [ ".cpp", ".hpp", ".cu" ]:
            fileinfo[ "commChar" ] = "*"
        else:
            fileinfo[ "commChar" ] = "#"

        print "WARN: ", fileinfo[ "filename" ], "could not identify comment char, take ", fileinfo[ "commChar" ], ", start line = ", fileinfo[ "startLine" ]

def noHeaderFile( filename ):
 
    file_only, file_extension = os.path.splitext( filename )

    if file_extension in SKIP_FILE_EXTENSIONS:
        return True

    return False

def forceHeaderFile( filename ):
 
    file_only, file_extension = os.path.splitext( filename )

    if file_extension in FORCE_FILE_EXTENSIONS:
        return True

    return False

def gitInfo( filename, fileinfo ):

    tmp_filename = "blame.info"

    os.system( "git blame " + filename + " > " + tmp_filename )

    for  line in open( tmp_filename, 'r' ):

       if "git_author" in fileinfo:
          continue

       items = line.split()

       if len( items ) < 5: 
          continue

       if items[1].startswith( "(" ):
          fileinfo[ "git_author" ] = items[1][1:] + " " + items[2]
          fileinfo[ "git_date" ] = items[3]
       else:
          # item[1] stands for filename
          fileinfo[ "git_author" ] = items[2][1:] + " " + items[3]
          fileinfo[ "git_date" ] = items[4]

    os.system( "rm -f " + tmp_filename )

def parseInfo( filename ):

    # collect all relevant information in a Python dictionary

    fileinfo = {}

    # Set already filename and suffix

    fileinfo[ "filename" ] = filename

    lastKey = ""

    linenr = 0 

    for  line in open( filename, 'r' ):

        stripLine = line.strip();

        isHeaderLine = True

        if ( stripLine == "" ):

           # empty line indicates that header is not available
 
           isHeaderLine = False

        else:

           # no more comment -> first relevant line of file reached

           if ( "commChar" in fileinfo ):
               if ( stripLine[0] != fileinfo[ "commChar" ] ):
                   isHeaderLine = False

        # no more header line, set startLine if not done yet

        if ( not isHeaderLine and not "startLine" in fileinfo ) :
            print "INFO: ", filename, ", set start line = ", linenr
            fileinfo[ "startLine" ] = linenr

        newKey = ""

        for key in LICENSE_KEYWORDS:

            pos = line.find( key )

            if ( pos >= 0  and  not key in fileinfo ):

                first = pos + len( key )
                last  = len( line ) - 1
                fileinfo[ key ] = line[ first: last ].lstrip()
                newKey = key
                commChar = line.lstrip()[0]

                if ( not "commChar" in fileinfo ):
                   fileinfo[ "commChar" ] = commChar

        if ( newKey == "" and lastKey == "@brief" ):
           # continuation line for brief
           pos = line.find( commChar )
           if ( pos >= 0 ): 
               first = pos + 1
               last  = len( line ) - 1
               fileinfo[ "@brief" ] = fileinfo[ "@brief"] + "|" + line[ first:last ].lstrip()
        
        if ( line.find( "#include" ) == 0 ):
           fileinfo[ "isCPP" ] = True

        if ( line.find( "#pragma" ) == 0 ):
           fileinfo[ "isHPP" ] = True

        lastKey = newKey

        linenr = linenr + 1  # increase line counter
    
    gitInfo( filename, fileinfo )

    checkStartLine( fileinfo )
    checkFilename( fileinfo )
    checkBrief( fileinfo )
    checkAuthor( fileinfo )
    checkDate( fileinfo )
    checkCommentChar( fileinfo )

    return fileinfo
 
def newHeaderFile( filename, fileinfo ):

    outfile = open( filename + ".tmpxy", 'w')

    writeLicenseHeader( outfile, fileinfo )

    linenr = 0
    startLine = fileinfo[ "startLine" ]

    for line in  open( filename, 'r' ):
        if ( linenr >= startLine ):
            outfile.write( line )
        linenr = linenr + 1

    outfile.close()

    os.system( "mv " + filename + ".tmpxy " + filename )

#  Method called for each file (not for directory files)

def run_file( filename ):

    if ( not os.path.isfile( filename ) ):
        return

    fileinfo = parseInfo( filename )

    if ( noHeaderFile( filename ) ):
        print 'INFO: skip file ', filename
    elif ( forceHeaderFile( filename ) or "@license" in fileinfo ):
        print 'OKAY: set new license header', filename
        # newHeaderFile( filename, fileinfo)
    else:
        print 'INFO: skip file without @license', filename

#  Method to be passed for os.path.walk (called recursively for each (sub-)directory
#
#  @param x dummy argument (not needed here)
#  @param dir_name  name of the current directory
#  @param files list of all files in the current directory

def run_directory(x, dir_name, files):

    for filename in files:
       run_file( dir_name + '/' + filename )

os.path.walk( ".", run_directory, 0)
