#----------------------------------------------------------------------------
# Top level CMakeLists.txt file for SCDS
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.3)
PROJECT(SCDS)

#-----------------------------------------------------------------------------
# set cmake module path
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${SCDS_SOURCE_DIR}/CMake")

#-----------------------------------------------------------------------------
# set directory for compiled libraries and executables
SET (LIBRARY_OUTPUT_PATH
  ${PROJECT_BINARY_DIR}/bin
  CACHE PATH
  "Single Directory for all libraries.")
SET (EXECUTABLE_OUTPUT_PATH
  ${PROJECT_BINARY_DIR}/bin
  CACHE PATH
  "Single Directory for all executables.")
MARK_AS_ADVANCED (
  LIBRARY_OUTPUT_PATH
  EXECUTABLE_OUTPUT_PATH
)

#-----------------------------------------------------------------------------
# windows specific options

# OS dependent stuff
IF (WIN32)
  ADD_DEFINITIONS(-DPLAN_WINNT)
  #INCLUDE_DIRECTORIES(
  #${ComputationalOncology_SOURCE_DIR}/Libraries/xdr
  #)
ENDIF(WIN32)

#-----------------------------------------------------------------------------
# unix specific options
IF(UNIX)
  OPTION(UNIX_WALL_MODE "Display all warnings." ON)
  SET(COMPONC_UNIX_WALL_MODE ${UNIX_WALL_MODE})
ENDIF(UNIX)


#-----------------------------------------------------------------------------
# set unix specific variables
IF(UNIX)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated")
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-deprecated")
  IF(SCDS_UNIX_WALL_MODE)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall")
  ENDIF(SCDS_UNIX_WALL_MODE)

  IF(SCDS_USE_STATIC_LIBS)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -static")
  ENDIF(SCDS_USE_STATIC_LIBS)
ENDIF(UNIX)

#-----------------------------------------------------------------------------
# decide what to compile
OPTION(COMPILE_SCDS_LIBS "Compile SCDS libraries." ON)
SET(SCDS_COMPILE_SCDS_LIBS ${COMPILE_SCDS_LIBS})

OPTION(COMPILE_SCDS_APPLICATIONS "Compile SCDS applications." ON)
SET(SCDS_COMPILE_SCDS_APPLICATIONS ${COMPILE_SCDS_APPLICATIONS})

OPTION(COMPILE_SCDS_TESTING "Compile SCDS testing." OFF)
SET(SCDS_COMPILE_SCDS_TESTING ${COMPILE_SCDS_TESTING})

#-----------------------------------------------------------------------------
# load local macros
INCLUDE(${SCDS_SOURCE_DIR}/CMake/SCDSMacros.cmake)

INCLUDE_DIRECTORIES(
	${SCDS_SOURCE_DIR}/Libraries/planio
	${SCDS_SOURCE_DIR}/Libraries/DownsampleFilter
	${SCDS_SOURCE_DIR}/Libraries/ContourTiler
	${SCDS_SOURCE_DIR}/Libraries/Matrix
  ${SCDS_SOURCE_DIR}/Libraries/Algorithms
  ${SCDS_SOURCE_DIR}/Libraries/Base
  ${SCDS_SOURCE_DIR}/Libraries/DataTypes
  ${SCDS_SOURCE_DIR}/Libraries/Numerics
  ${SCDS_SOURCE_DIR}/Libraries/UtilitiesDataTypes
  ${SCDS_SOURCE_DIR}/Libraries/PowerCrust
  ${SCDS_SOURCE_DIR}/Libraries/tclap
)

SCDS_USE_LAPACK()

#-----------------------------------------------------------------------------
# recurse into subdirectories
IF(SCDS_COMPILE_SCDS_LIBS)
  SUBDIRS(${SCDS_SOURCE_DIR}/Libraries)
ENDIF(SCDS_COMPILE_SCDS_LIBS)


IF(SCDS_COMPILE_SCDS_TESTING)
  SUBDIRS(${SCDS_SOURCE_DIR}/Testing)
ENDIF(SCDS_COMPILE_SCDS_TESTING)

IF(SCDS_COMPILE_SCDS_APPLICATIONS)
  SUBDIRS(${SCDS_SOURCE_DIR}/Applications)
ENDIF(SCDS_COMPILE_SCDS_APPLICATIONS)






#-----------------------------------

IF(WIN32)
SET(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS})
ENDIF(WIN32)

# Setup ITK
# * include directories are added automatically
# * applications must add appropriate libraries to link libraries
MARK_AS_ADVANCED(USE_ITK)
OPTION(USE_ITK "Use ITK if found" ON)
FIND_PACKAGE(ITK)
IF(ITK_FOUND)
  INCLUDE (${ITK_USE_FILE})
ELSE(ITK_FOUND)
  MESSAGE("Could not find ITK directory.")
ENDIF(ITK_FOUND)

#-----------------------------------------------------------------------------
# Setup VTK
# * include directories are added automatically
# * applications must add appropriate libraries to link libraries
MARK_AS_ADVANCED(USE_VTK)
OPTION(USE_VTK "Use VTK if found" ON)
FIND_PACKAGE(VTK)
IF(VTK_FOUND)
  INCLUDE (${VTK_USE_FILE})
ELSE(VTK_FOUND)
  MESSAGE("Could not find VTK directory.")
ENDIF(VTK_FOUND)


MARK_AS_ADVANCED(USE_FFTWF)
OPTION(USE_FFTWF "Use single precision FFTW if found" ON)
FIND_PACKAGE(FFTW)


#pthread
if (WIN32)
FIND_PACKAGE(PTHREADS REQUIRED)
IF(PTHREADS_FOUND)
   INCLUDE_DIRECTORIES(${PTHREADS_INCLUDE_DIR})
ELSE(PTHREADS_FOUND)
    MESSAGE("Could not find pthread directory.")
ENDIF(PTHREADS_FOUND)
ENDIF(WIN32)
