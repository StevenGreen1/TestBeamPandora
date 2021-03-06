########################################################
# cmake file for building MarlinPandora
# @author Jan Engels, DESY
CMAKE_MINIMUM_REQUIRED(VERSION 2.6 FATAL_ERROR)
########################################################


# project name
PROJECT( CalicePandora )


# project version
SET( ${PROJECT_NAME}_VERSION_MAJOR 1 )
SET( ${PROJECT_NAME}_VERSION_MINOR 0 )
SET( ${PROJECT_NAME}_VERSION_PATCH 0 )



### DEPENDENCIES ############################################################

FIND_PACKAGE( ILCUTIL COMPONENTS ILCSOFT_CMAKE_MODULES REQUIRED )

# load default settings from ILCSOFT_CMAKE_MODULES
INCLUDE( ilcsoft_default_settings )


#FIND_PACKAGE( Marlin 1.0 REQUIRED )
#FIND_PACKAGE( MarlinUtil REQUIRED )
FIND_PACKAGE( Marlin REQUIRED )
FIND_PACKAGE( MarlinUtil REQUIRED )
FIND_PACKAGE( PandoraPFANew REQUIRED )
FIND_PACKAGE( LCIO REQUIRED )
FIND_PACKAGE( GEAR REQUIRED )
FIND_PACKAGE( CALICE_USERLIB REQUIRED )
FIND_PACKAGE( RAW2CALOHIT REQUIRED )
FIND_PACKAGE( ROOT REQUIRED )
FIND_PACKAGE( LCCD REQUIRED )


FOREACH( pkg Marlin MarlinUtil PandoraPFANew LCIO GEAR CALICE_USERLIB RAW2CALOHIT ROOT LCCD)
    IF( ${pkg}_FOUND )
        INCLUDE_DIRECTORIES( ${${pkg}_INCLUDE_DIRS} )
	MESSAGE(STATUS "bla=  ${${pkg}_LIBRARIES}")
        LINK_LIBRARIES( ${${pkg}_LIBRARIES} )
        ADD_DEFINITIONS ( ${${pkg}_DEFINITIONS} )
    ENDIF()
ENDFOREACH()



### DOCUMENTATION ############################################################
OPTION( INSTALL_DOC "Set to OFF to skip build/install Documentation" OFF )

IF( INSTALL_DOC AND EXISTS "${PROJECT_SOURCE_DIR}/doc/CMakeLists.txt" )
    ADD_SUBDIRECTORY( ./doc )
ELSE()
    MESSAGE( STATUS "INSTALL_DOC set to OFF" )
    SET( INSTALL_DOC OFF )
ENDIF()



### LIBRARY AND TOOLS #######################################################

#include directories
INCLUDE_DIRECTORIES( ./include )
INSTALL_DIRECTORY( ./include DESTINATION . FILES_MATCHING PATTERN "*.h" )

ADD_DEFINITIONS( "-Wall -ansi -pedantic" )
ADD_DEFINITIONS( "-Wno-long-long" )

AUX_SOURCE_DIRECTORY( ./src library_sources )
ADD_SHARED_LIBRARY( ${PROJECT_NAME} ${library_sources} )
INSTALL_SHARED_LIBRARY( ${PROJECT_NAME} DESTINATION lib )


# display some variables and write them to cache
DISPLAY_STD_VARIABLES()
