### set include directories ### 

include_directories ( . )


### set source files ###

set ( CXXSOURCES
        ${CXXSOURCES}
      
      ## criteria
        solver/criteria/Criterion.cpp
        solver/criteria/IterationCount.cpp
        solver/criteria/ResidualStagnation.cpp
        solver/criteria/ResidualThreshold.cpp
      
      ## logger
        solver/logger/CommonLogger.cpp
        solver/logger/FileLogger.cpp
        solver/logger/Logger.cpp
        solver/logger/OpenMPTimer.cpp
        solver/logger/Timer.cpp

      ## solvers
        solver/CG.cpp
        solver/DefaultJacobi.cpp
        solver/GMRES.cpp
        solver/InverseSolver.cpp
        #solver/LUSolver.cpp
        solver/SimpleAMG.cpp
        solver/SingleGridSetup.cpp
        solver/SOR.cpp
        #solver/SpecialLUSolver.cpp
        solver/SpecializedJacobi.cpp
        solver/TrivialPreconditioner.cpp
      
      ## solver abstract and utils
        solver/AMGSetup.cpp
        solver/AMGSetupFactory.cpp
        solver/AMGSetupManager.cpp
        solver/IterativeSolver.cpp
        solver/OmegaSolver.cpp
        solver/SolutionProxy.cpp
        solver/Solver.cpp
    )
    
if ( Boost_INCLUDE_DIR )
    message ( STATUS "Boost_INCLUDE_DIR = ${Boost_INCLUDE_DIR}" )
    
    get_filename_component ( Boost_PATH ${Boost_INCLUDE_DIR} PATH )
    message ( STATUS "Boost_PATH = ${Boost_PATH}" )
    
    find_file ( BOOST_SPIRIT_VERSION_HPP "boost/spirit/version.hpp" PATHS ${Boost_INCLUDE_DIR} NO_DEFAULT_PATH )
    find_file ( BOOST_SPIRIT_VERSION_HPP_2 "boost/spirit/include/version.hpp" PATHS ${Boost_INCLUDE_DIR} NO_DEFAULT_PATH )
    mark_as_advanced ( FORCE BOOST_SPIRIT_VERSION_HPP BOOST_SPIRIT_VERSION_HPP_2 )
    message ( STATUS "Boost spirit version file = ${BOOST_SPIRIT_VERSION_HPP}" )
    
### check boost spirit version ###
    if ( BOOST_SPIRIT_VERSION_HPP )
        file (READ ${BOOST_SPIRIT_VERSION_HPP} SOURCE)
        string (REGEX MATCH "^.*#define SPIRIT_VERSION 0x([0-9]+).*$" BOOST_SPIRIT_VERSION_MATCH ${SOURCE} )
        
        if ( BOOST_SPIRIT_VERSION_MATCH )
            set ( Boost_SPIRIT_VERSION ${CMAKE_MATCH_1} )
            message ( STATUS "Found Boost spirit version ${Boost_SPIRIT_VERSION}" )
            #  Boost_SPIRIT_VERSION % 100 is the sub-minor version             
            #  Boost_SPIRIT_VERSION / 100 % 1000 is the minor version         
            #  Boost_SPIRIT_VERSION / 100000 is the major version
        else ( BOOST_SPIRIT_VERSION_MATCH )
            if (BOOST_SPIRIT_VERSION_HPP_2)
                file ( READ ${BOOST_SPIRIT_VERSION_HPP_2} SOURCE )
                string ( REGEX MATCH "^.*#define SPIRIT_VERSION 0x([0-9]+).*$" BOOST_SPIRIT_VERSION_MATCH_2 ${SOURCE} )
                if ( BOOST_SPIRIT_VERSION_MATCH_2 )
                    set ( Boost_SPIRIT_VERSION ${CMAKE_MATCH_1} )
                    message ( STATUS "Found Boost spirit version ${Boost_SPIRIT_VERSION}" )
                endif ( BOOST_SPIRIT_VERSION_MATCH_2 )
            endif ( BOOST_SPIRIT_VERSION_HPP_2 )
        endif ( BOOST_SPIRIT_VERSION_MATCH )

    else ( BOOST_SPIRIT_VERSION_HPP )
        message ( SEND_ERROR "Found Boost but boost/spirit/version.hpp is missing." )
    endif ( BOOST_SPIRIT_VERSION_HPP )

else ( Boost_INCLUDE_DIR )
    message ( FATAL_ERROR "Boost (include directory) not found: give hint by environment variable BOOST_ROOT" ) 
endif ( Boost_INCLUDE_DIR )
    
if ( ${Boost_SPIRIT_VERSION} LESS "2010" )
    message ( WARNING "Not building MetaSolver because Boost Spirit is to old: ${Boost_Spirit_VERSION}." )
else ( ${Boost_SPIRIT_VERSION} LESS "2010" )

    set ( CXXSOURCES ${CXXSOURCES} 
            solver/creator/SolverFactory.cpp
            
            solver/creator/Creator.cpp
            solver/creator/SolverCreator.cpp

            solver/creator/InverseSolverCreator.cpp
            #solver/creator/LUSolverCreator.cpp

            solver/creator/CriteriaCreator.cpp
            solver/creator/LoggerCreator.cpp
            solver/creator/IterativeSolverCreator.cpp
            
            solver/creator/CGCreator.cpp
            solver/creator/SimpleAMGCreator.cpp
            solver/creator/GMRESCreator.cpp
            
            solver/creator/OmegaSolverCreator.cpp
            solver/creator/DefaultJacobiCreator.cpp
            solver/creator/SORCreator.cpp
            solver/creator/SpecializedJacobiCreator.cpp
            #solver/creator/SpecialLUSolverCreator.cpp
            
            solver/MetaSolver.cpp 
        )  
            
endif ( ${Boost_SPIRIT_VERSION} LESS "2010" )
    
# do this only once per file
set ( CXXSOURCES ${CXXSOURCES} PARENT_SCOPE )

### for external amg ###

file ( GLOB AMG_SUBDIRPRESENT "amg" )

if ( AMG_SUBDIRPRESENT AND MPI_FOUND )
    add_subdirectory ( ${AMG_SUBDIRPRESENT} )
endif ( AMG_SUBDIRPRESENT AND MPI_FOUND )
