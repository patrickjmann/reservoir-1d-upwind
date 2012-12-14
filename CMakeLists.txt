#################################################################################
# $Source: /home/mann/Dropbox/research/1/ReservoirSimulation/1d/upwind1/RCS/CMakeLists.txt,v $
# $Date: 2012/09/19 21:47:05 $
# $Revision: 1.1 $
#################################################################################
SET( PROJECT_NAME reservoir_upwind1 )
SET( VERSION_MAJOR "1" )
SET( VERSION_MINOR "0" )
SET( PATCH_LEVEL "0" )
SET( VERSION "${VERSION_MAJOR}.${VERSION_MINOR}.${PATCH_LEVEL}" )
SET( AUTHOR "Patrick Mann" )
SET( SUMMARY "${PROJECT_NAME}: High-resolution Upwind Predictor-Corrector for 1d Porous Medium" )

PROJECT(${PROJECT_NAME})

MESSAGE( "\nConfiguring the build for ${PROJECT_NAME} v${VERSION}\n" )

# This module simply prints out cmake system information
##INCLUDE(CMakePrintSystemInformation)

#--------------------------------------
cmake_minimum_required(VERSION 2.6)
###cmake_policy(VERSION 2.6)  # Not necessary when minimum defined

#--------------------------------------
SET(SOURCE_DIR ${${PROJECT_NAME}_SOURCE_DIR})
SET(BINARY_DIR ${${PROJECT_NAME}_BINARY_DIR})

#-----------------------------------------
ADD_DEFINITIONS( -DDEBUG )
ADD_DEFINITIONS( -DVERSION="\\"${VERSION}\\"" -DPACKAGE="\\"${PROJECT_NAME}\\"" )

#-----------------------------------------
# Note: the "cmake -D.." puts the value in the cache.
#  -this is NOT overwritten when another cmake is run.
#  -so the old value is still there, and is therefore the default
#
#SET( CXX "g++" )

IF( DEFINED CXX )
ELSE( DEFINED CXX )  # NOT DEFINED
  MESSAGE( FATAL_ERROR "ERROR: Compiler must be specified: cmake -D CXX:STRING=\"g++\"" )
ENDIF( DEFINED CXX )

MESSAGE( STATUS "the compiler to be used is ${CXX}" )

#-----------------------------------------
IF( CXX STREQUAL "icpc" )
  MESSAGE( STATUS "in icpc section" )
  SET( CMAKE_CXX_COMPILER "icpc" )
  SET( BASE_DIR $ENV{HOME}/intel )
  ADD_DEFINITIONS( -O3 -static -opt_report -opt_report_file icc_opt.report)
  #### ADD_DEFINITIONS( -openmp -openmp_report2 )
  ADD_DEFINITIONS( -cxxlib-icc )
  ADD_DEFINITIONS( -wd424 -wd981 -wd111 -wd117 -wd810 -wd1419 -wd444 )
  ADD_DEFINITIONS( -Wall -Qoption,link,--noinhibit-exec )
ENDIF( CXX STREQUAL "icpc" )

IF( CXX STREQUAL "g++" )
  SET( CMAKE_CXX_COMPILER "g++" )
  SET( BASE_DIR $ENV{HOME}/gcc )

  ADD_DEFINITIONS( -O3 -Wall )
ENDIF( CXX STREQUAL "g++" )

MESSAGE( STATUS "BASE_DIR for includes, libraries and installs is \"${BASE_DIR}\"" )
MESSAGE( STATUS "CMAKE_CXX_COMPILER is \"${CMAKE_CXX_COMPILER}\"" )

ADD_DEFINITIONS( -g )

#-----------------------------------------
SET( INSTALL_DIR ${BASE_DIR}/bin )
INCLUDE_DIRECTORIES(
  ${SOURCE_DIR}/include
  ${BASE_DIR}/include
)

###SUBDIRS(evolve io main include)
add_subdirectory(evolve)
add_subdirectory(io)
add_subdirectory(include)
add_subdirectory(main)

INSTALL(PROGRAMS ${BINARY_DIR}/main/${PROJECT_NAME} DESTINATION ${INSTALL_DIR} )

#---------------------------------------
ADD_CUSTOM_TARGET( emacsclean
  COMMAND rm -f *~ */*~
  WORKING_DIRECTORY ${SOURCE_DIR}
)
#---------------------------------------
# RCS:  rcs targets should be in the subdirs

ADD_CUSTOM_TARGET( rcs
   COMMAND ci -l CMakeLists.txt LOG DESCRIPTION configure
   WORKING_DIRECTORY ${SOURCE_DIR}
)
add_dependencies(rcs evolve_RCS io_RCS include_RCS main_RCS )

# Handy when moving to another branch to be able to 
# clean out all the RCS files.  I should probably set this
# up automatically so all SUBDIRS are visited!

ADD_CUSTOM_TARGET( rcsclean
   COMMAND rm -rf RCS
   COMMAND rm -rf evolve/RCS
   COMMAND rm -rf io/RCS
   COMMAND rm -rf main/RCS
   COMMAND rm -rf include/RCS
   WORKING_DIRECTORY ${SOURCE_DIR}
)

#---------------------------------------
# NOTE: cpack copies source files into a sub-directory and
#  then does the usual tar, gzip.  The .tar.gz is then
#  copied back out.  THE TEMPORORARY AREAS ARE NOT REMOVED!
#  As they are part of the source tree they will be added
#  to the dist the next time cpack is run.  Very annoying.
#  So I have a target here which gets rid of the temps.

ADD_CUSTOM_TARGET( dist
   COMMAND cpack --config CPackSourceConfig.cmake
   COMMAND rm -rf _CPack_Packages
)

#---------------------------------------
# This generates CPack output to do the distribution

SET( CPACK_PACKAGE_VERSION_MAJOR ${VERSION_MAJOR} )
SET( CPACK_PACKAGE_VERSION_MINOR ${VERSION_MINOR} )
SET( CPACK_PACKAGE_VERSION_PATCH ${PATCH_LEVEL} )
SET( CPACK_PACKAGE_VENDOR ${AUTHOR} )
SET( CPACK_PACKAGE_DESCRIPTION_FILE "${SOURCE_DIR}/DESCRIPTION" )
SET( CPACK_PACKAGE_DESCRIPTION_SUMMARY "${SUMMARY}" )

INCLUDE(CPack)

# For some unknown reason the CPACK_GENERATOR has to be set
# after the INCLUDE(CPack).  Everything else comes before.
SET( CPACK_GENERATOR TGZ )

