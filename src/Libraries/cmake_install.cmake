# Install script for directory: /afs/radonc/home/sharonxx/public/ITK_based_Program/SCDS/code/Libraries

# Set the install prefix
IF(NOT DEFINED CMAKE_INSTALL_PREFIX)
  SET(CMAKE_INSTALL_PREFIX "/usr/local")
ENDIF(NOT DEFINED CMAKE_INSTALL_PREFIX)
STRING(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
IF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  IF(BUILD_TYPE)
    STRING(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  ELSE(BUILD_TYPE)
    SET(CMAKE_INSTALL_CONFIG_NAME "debug")
  ENDIF(BUILD_TYPE)
  MESSAGE(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
ENDIF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)

# Set the component getting installed.
IF(NOT CMAKE_INSTALL_COMPONENT)
  IF(COMPONENT)
    MESSAGE(STATUS "Install component: \"${COMPONENT}\"")
    SET(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  ELSE(COMPONENT)
    SET(CMAKE_INSTALL_COMPONENT)
  ENDIF(COMPONENT)
ENDIF(NOT CMAKE_INSTALL_COMPONENT)

# Install shared libraries without execute permission?
IF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  SET(CMAKE_INSTALL_SO_NO_EXE "0")
ENDIF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)

IF(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  INCLUDE("/afs/radonc/home/sharonxx/public/ITK_based_Program/SCDS/code/Libraries/Algorithms/cmake_install.cmake")
  INCLUDE("/afs/radonc/home/sharonxx/public/ITK_based_Program/SCDS/code/Libraries/Base/cmake_install.cmake")
  INCLUDE("/afs/radonc/home/sharonxx/public/ITK_based_Program/SCDS/code/Libraries/DataTypes/cmake_install.cmake")
  INCLUDE("/afs/radonc/home/sharonxx/public/ITK_based_Program/SCDS/code/Libraries/Numerics/cmake_install.cmake")
  INCLUDE("/afs/radonc/home/sharonxx/public/ITK_based_Program/SCDS/code/Libraries/UtilitiesDataTypes/cmake_install.cmake")
  INCLUDE("/afs/radonc/home/sharonxx/public/ITK_based_Program/SCDS/code/Libraries/ContourTiler/cmake_install.cmake")
  INCLUDE("/afs/radonc/home/sharonxx/public/ITK_based_Program/SCDS/code/Libraries/DownsampleFilter/cmake_install.cmake")
  INCLUDE("/afs/radonc/home/sharonxx/public/ITK_based_Program/SCDS/code/Libraries/Matrix/cmake_install.cmake")
  INCLUDE("/afs/radonc/home/sharonxx/public/ITK_based_Program/SCDS/code/Libraries/planio/cmake_install.cmake")
  INCLUDE("/afs/radonc/home/sharonxx/public/ITK_based_Program/SCDS/code/Libraries/PowerCrust/cmake_install.cmake")
  INCLUDE("/afs/radonc/home/sharonxx/public/ITK_based_Program/SCDS/code/Libraries/tclap/cmake_install.cmake")

ENDIF(NOT CMAKE_INSTALL_LOCAL_ONLY)

