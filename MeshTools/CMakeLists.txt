CMAKE_MINIMUM_REQUIRED(VERSION 2.4)



## set up ITK
FIND_PACKAGE(ITK REQUIRED)
IF(ITK_FOUND)
  INCLUDE (${ITK_USE_FILE})
ELSE(ITK_FOUND)
  MESSAGE("Could not find ITK directory.")
ENDIF(ITK_FOUND)


## set up VTK
INCLUDE(vtkSetup.cmake)
IF(VTK_FOUND)
  INCLUDE (${VTK_USE_FILE})
ELSE(VTK_FOUND)
  MESSAGE("Could not find VTK directory.")
ENDIF(VTK_FOUND)

## set up fltk
include(fltkSetup.cmake)
INCLUDE_DIRECTORIES(${FLTK_INCLUDE_PATH})



# Shorthand to source and binary directories
SET(MY_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
SET(MY_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR})

# Include directories
INCLUDE_DIRECTORIES(
  ${MY_SOURCE_DIR}
  ${MY_BINARY_DIR}
 )


# Assign libraries to targets
Link_LIBRARIES( 
  ${FLTK_LIBRARIES}
  ITKIO
  ITKCommon
  ITKBasicFilters
  vtkRendering
  vtkGraphics
  vtkImaging
  vtkIO
  vtkCommon
  vtkHybrid
  vtkFiltering
  )

SET(SRC_VOL2SURF_CXX
vol2surf.cxx
vol2surfMain.cxx
)


SET(APPLICATION_FLUIDS
 vol2surfGUI.fl
)

# FLTK Wrapping
FLTK_WRAP_UI(Vol2surf ${APPLICATION_FLUIDS})
ADD_EXECUTABLE(Vol2surf ${SRC_VOL2SURF_CXX} ${Vol2surf_FLTK_UI_SRCS})


SET(SRC_MESHMATH_CXX MeshMath.cxx itkMeshTovtkPolyData.cxx vtkPolyDataToitkMesh.cxx)
ADD_EXECUTABLE(MeshMath ${SRC_MESHMATH_CXX})
TARGET_LINK_LIBRARIES(MeshMath ${ITK_LIBRARIES})


ADD_EXECUTABLE(NormalGrabber NormalGrabber.cxx)

ADD_EXECUTABLE(XYZgrabber XYZgrabber.cxx)
