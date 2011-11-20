#-----------------------------------------------------------------------------
# macros used by SCDS

# this is not a macro but necessary for ADD_GUI_EXECUTABLE macro
IF(APPLE)
  # Apple applications need a special "resource" in order to function
  FIND_PROGRAM(APPLE_RESOURCE Rez /Developer/Tools)
  FIND_FILE(FLTK_RESOURCE mac.r /usr/local/include/FL)
  IF(NOT FLTK_RESOURCE)
    MESSAGE("Fltk resources not found, GUI application will not respond to mouse events")
  ENDIF(NOT FLTK_RESOURCE)
ENDIF(APPLE)

MACRO(ADD_GUI_EXECUTABLE name sources)
  ADD_EXECUTABLE(${name} ${sources})
  SET(EXEC_PATH ${EXECUTABLE_OUTPUT_PATH})
  IF(NOT EXEC_PATH)
    SET(EXEC_PATH ${CMAKE_CURRENT_BINARY_DIR})
  ENDIF(NOT EXEC_PATH)
  IF(APPLE_RESOURCE)
    ADD_CUSTOM_COMMAND(SOURCE ${name}
      COMMAND ${APPLE_RESOURCE}
      ARGS -t APPL ${FLTK_RESOURCE} -o
      ${EXEC_PATH}/${name}
      TARGET ${name})
  ENDIF(APPLE_RESOURCE)
ENDMACRO(ADD_GUI_EXECUTABLE)

MACRO(SCDS_USE_ITK)
INCLUDE(${SCDS_SOURCE_DIR}/CMake/itkSetup.cmake)  
ENDMACRO(SCDS_USE_ITK)

MACRO(SCDS_USE_VTK)
INCLUDE(${SCDS_SOURCE_DIR}/CMake/vtkSetup.cmake)  
ENDMACRO(SCDS_USE_VTK)

MACRO(SCDS_USE_FLTK)
INCLUDE(${SCDS_SOURCE_DIR}/CMake/fltkSetup.cmake)  
ENDMACRO(SCDS_USE_FLTK)

MACRO(SCDS_USE_FFTW)
INCLUDE(${SCDS_SOURCE_DIR}/CMake/fftwSetup.cmake)  
ENDMACRO(SCDS_USE_FFTW)

MACRO(SCDS_USE_FFTW3)
INCLUDE(${SCDS_SOURCE_DIR}/CMake/fftw3Setup.cmake)  
ENDMACRO(SCDS_USE_FFTW3)

MACRO(SCDS_USE_LAPACK)
INCLUDE(${SCDS_SOURCE_DIR}/CMake/lapackSetup.cmake)  
ENDMACRO(SCDS_USE_LAPACK)
