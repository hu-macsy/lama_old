###
 # @file scai_build.sh
 #
 # @license
 # Copyright (c) 2009-2017
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
 #
 # Other Usage
 # Alternatively, this file may be used in accordance with the terms and
 # conditions contained in a signed written agreement between you and
 # Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 # @endlicense
 #
 # @brief ToDo: Missing description in ./scai_build.sh
 # @author Thomas Brandes
 # @date 19.08.2015
###

###
 # @file scai_build.sh
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
 #
 # Other Usage
 # Alternatively, this file may be used in accordance with the terms and
 # conditions contained in a signed written agreement between you and
 # Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 # @endlicense
 #
 # @brief ToDo: Missing description in ./scai_build.sh
 # @author Thomas Brandes
 # @date 19.08.2015
###

# scai_build.sh clean
# scai_build.sh <install_dir>

# Skript to build all SCAI projects

function checkErrorValue( ) {
	$*
	if [ "$?" -ne 0 ];
	then
		echo "Build Failed. Aborting..."
		exit 1
	fi
}

if [ "$#" -ne 1 ]; then
   echo "Illegal call with $# arguments"
   echo "scai_build.sh clean"
   echo "scai_build.sh <install_prefix>"
   exit 1
fi

PROJECTS="common logging tracing tasking hmemo kregistry blaskernel utilskernel sparsekernel dmemo lama solver"

if [ "$1" == "clean" ]; then
    for project in $PROJECTS
    do
        echo "Clean SCAI project $project" 
        cd $project
        rm -rf build
        cd ..
    done
else
    for project in $PROJECTS
    do
        echo "Build SCAI project $project"
        cd $project
        mkdir -p build
        cd build
        # echo "cmake .. -DCMAKE_INSTALL_PREFIX=$1\""
        # optional:  -DCMAKE_CXX_COMPILER=icpc
        # optional:  -DCXX_SUPPORTS_C11=0
        # optional:  -DBoost_NO_BOOST_CMAKE=TRUE 
        checkErrorValue cmake .. -DCMAKE_INSTALL_PREFIX=$1 -DCMAKE_BUILD_TYPE=Release -DBUILD_TEST=0
        checkErrorValue make -j 8
        checkErrorValue make install
        
	cd ../..
    done
fi

