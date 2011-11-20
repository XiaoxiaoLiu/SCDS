
IF(NOT VTK_DIR)
  FIND_PATH(VTK_DIR VTKConfig.cmake
    ${RADONC_LIBS_PATH}/VTK
    )
ENDIF(NOT VTK_DIR)

IF(VTK_DIR)
  # VTK_DIR must be set for this to work...
  INCLUDE (${CMAKE_ROOT}/Modules/FindVTK.cmake)
  IF(VTK_FOUND)
    # Load the compiler settings used for VTK.
    IF(VTK_BUILD_SETTINGS_FILE)
      INCLUDE(${CMAKE_ROOT}/Modules/CMakeImportBuildSettings.cmake)
      CMAKE_IMPORT_BUILD_SETTINGS(${VTK_BUILD_SETTINGS_FILE})
    ENDIF(VTK_BUILD_SETTINGS_FILE)

    # Add compiler flags needed to use VTK.
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${VTK_REQUIRED_C_FLAGS}")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${VTK_REQUIRED_CXX_FLAGS}")

    # Add include directories needed to use VTK.
    INCLUDE_DIRECTORIES(${VTK_INCLUDE_DIRS})

    # Add link directories needed to use VTK.
    LINK_DIRECTORIES(${VTK_LIBRARY_DIRS})

#    INCLUDE (${VTK_USE_FILE})
  ELSE(VTK_FOUND)
    MESSAGE("Could not find vtk directory.")
  ENDIF(VTK_FOUND)
ELSE(VTK_DIR)
  MESSAGE("Could not find vtk directory.")
ENDIF(VTK_DIR)

