# this is FindFLTK.cmake with some modifications

#
# Find the native FLTK includes and library
#
# The following settings are defined
# FLTK_FLUID_EXECUTABLE, where to find the Fluid tool
# FLTK_WRAP_UI, This allows the FLTK_WRAP_UI command to work.
# FLTK_INCLUDE_DIR, where to find include files
# FLTK_LIBRARIES, list of fltk libraries
# FLTK_VERSION_1.0.11 Use this Version
# FLTK_VERSION_1.1 Use this Version
# FLTK_FOUND, Don't use FLTK if false.


# The following settings should not be used in general.
# FLTK_BASE_LIBRARY    = the full path to fltk.lib
# FLTK_GL_LIBRARY      = the full path to fltk_gl.lib
# FLTK_FORMS_LIBRARY   = the full path to fltk_forms.lib
# FLTK_IMAGES_LIBRARY  = the full path to fltk_images.lib

OPTION(FLTK_VERSION_1.1    "Use FLTK version 1.1"    1)
OPTION(FLTK_VERSION_1.0.11 "Use FLTK version 1.0.11" 0)

# Exclusion between the twFLTK_FLUID_EXECUTABLE o version

IF(FLTK_VERSION_1.1)
  SET(FLTK_VERSION_1.0.11 0)
ENDIF(FLTK_VERSION_1.1)

FIND_PATH(FLTK_BASE_DIR FL/Fl.H
  ${RADONC_LIBS_PATH}/fltk
)
SET(FLTK_INCLUDE_DIR ${FLTK_BASE_DIR})
SET (FLTK_LIBRARY_DIR ${FLTK_INCLUDE_DIR}/lib)

SET(FLTK_FLUID_EXECUTABLE ${FLTK_INCLUDE_DIR}/bin/fluid CACHE PATH "Path to fltk's ui builder")




# Platform dependent libraries required by FLTK

IF(WIN32)
  IF(NOT CYGWIN)
    IF(BORLAND)
      SET( FLTK_PLATFORM_DEPENDENT_LIBS import32 )
    ELSE(BORLAND)
      SET( FLTK_PLATFORM_DEPENDENT_LIBS wsock32 comctl32 )
    ENDIF(BORLAND)
  ENDIF(NOT CYGWIN)
ENDIF(WIN32)

IF(UNIX)
  INCLUDE(${CMAKE_ROOT}/Modules/FindX11.cmake)
  SET( FLTK_PLATFORM_DEPENDENT_LIBS ${X11_LIBRARIES} -lm)
ENDIF(UNIX)

IF(APPLE)
  # Should set the FindX11.cmake properly
  INCLUDE_DIRECTORIES(/usr/X11R6/include)
  SET( FLTK_PLATFORM_DEPENDENT_LIBS  "-framework Carbon -framework Cocoa -framework ApplicationServices -lz")
ENDIF(APPLE)


# Enable the Wrap UI command
IF (FLTK_FLUID_EXECUTABLE)
  SET ( FLTK_WRAP_UI 1 CACHE INTERNAL "Can we honor the FLTK_WRAP_UI command" )
ENDIF (FLTK_FLUID_EXECUTABLE)


MARK_AS_ADVANCED(
  FLTK_VERSION_1.0.11
  FLTK_VERSION_1.1
)

IF(APPLE)
  # FindFLTK does not add these which are commonly required
  SET (FLTK_PLATFORM_DEPENDENT_LIBS "${FLTK_PLATFORM_DEPENDENT_LIBS} -framework AGL -framework OpenGL")
ENDIF(APPLE)

IF (WIN32)

  SET (FLTK_LIBRARIES_RELEASE
    fltk
    fltkgl
    fltkforms
    fltkimages
    ${FLTK_PLATFORM_DEPENDENT_LIBS}
  )

  SET (FLTK_LIBRARIES_DEBUG
    fltkd
    fltkgld
    fltkformsd
    fltkimagesd
    ${FLTK_PLATFORM_DEPENDENT_LIBS}
  )

  OPTION(DEBUG_BUILD "Debug only build." 0)
  IF(DEBUG_BUILD)
    SET( FLTK_LIBRARIES ${FLTK_LIBRARIES_DEBUG} )
  ELSE(DEBUG_BUILD)
    SET( FLTK_LIBRARIES ${FLTK_LIBRARIES_RELEASE} )
  ENDIF(DEBUG_BUILD)

ENDIF (WIN32)

IF(UNIX)
    SET (FLTK_LIBRARIES
        fltk_gl
        fltk_forms
        fltk_images
        fltk
        ${FLTK_PLATFORM_DEPENDENT_LIBS}
    )
ENDIF(UNIX)

LINK_DIRECTORIES (${FLTK_LIBRARY_DIR})
INCLUDE_DIRECTORIES(${FLTK_INCLUDE_DIR})


