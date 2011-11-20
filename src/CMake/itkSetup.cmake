IF(NOT ITK_DIR)
  FIND_PATH(ITK_DIR ITKConfig.cmake
    ${RADONC_LIBS_PATH}/ITK
    )
ENDIF(NOT ITK_DIR)



IF(ITK_DIR)
  # ITK_DIR must be set for this to work...
  INCLUDE (${CMAKE_ROOT}/Modules/FindITK.cmake)
  IF(ITK_FOUND)
    INCLUDE (${ITK_USE_FILE})
  ELSE(ITK_FOUND)
    MESSAGE("Could not find itk directory.")
  ENDIF(ITK_FOUND)
ELSE(ITK_DIR)
  MESSAGE("Could not find itk directory.")
ENDIF(ITK_DIR)



