###
 # @file CustomCommands/distclean.sh
 #
 # @license
 # Copyright (c) 2009-2016
 # Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 # for Fraunhofer-Gesellschaft
 #
 # This file is part of the SCAI framework LAMA.
 #
 # LAMA is free software: you can redistribute it and/or modify it under the
 # terms of the GNU Affero General Public License as published by the Free
 # Software Foundation, either version 3 of the License, or (at your option)
 # any later version.
 #
 # LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 # WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 # FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for
 # more details.
 #
 # You should have received a copy of the GNU Affero General Public License
 # along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 # @endlicense
 #
 # @brief ToDo: Missing description in ./CustomCommands/distclean.sh
 # @author Lauretta Schubert
 # @date 14.10.2015
###
#!/bin/sh

recursive_remove()
{
	#echo $PWD
	rm -f *.so *.cmake install_manifest.txt
	rm -rf tmp src CMakeFiles
	
    for dir in `ls`
    do
        if [ -d $dir ]
        then
	    	cd $dir
	    	recursive_remove
	    	cd ..
        fi
    done
}

recursive_remove $PWD