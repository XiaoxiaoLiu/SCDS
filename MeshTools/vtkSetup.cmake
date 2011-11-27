
IF(NOT VTK_DIR)
  FIND_PATH(VTK_DIR VTKConfig.cmake
    ${NeuroLib_LIBS_PATH}/VTK-4.2.3-VC++
    ${NeuroLib_LIBS_PATH}/VTK-Linux 
    ${NeuroLib_LIBS_PATH}/VTK-Sun
    ${NeuroLib_LIBS_PATH}/VTK-VC++
    )
ENDIF(NOT VTK_DIR)

IF(VTK_DIR)

  # VTK_DIR must be set for this to work...
  FIND_PACKAGE(VTK REQUIRED)
  INCLUDE(${USE_VTK_FILE})
  IF(VTK_FOUND)
    LINK_LIBRARIES(${VTK_LIBRARIES})
    ADD_DEFINITIONS(-DUSE_VTK)
  ENDIF(VTK_FOUND)

ELSE(VTK_DIR)
   MESSAGE(FATAL_ERROR "VTK library not found!\n" "Please go to http://www.ia.unc.edu/dev/tutorials/InstallLib")
ENDIF(VTK_DIR)

