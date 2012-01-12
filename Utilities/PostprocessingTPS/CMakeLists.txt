project(ParticleCorrespondencePostprocessingTPS)

cmake_minimum_required(VERSION 2.6)

# include a test for endianness (used to be done by default)
INCLUDE (${CMAKE_ROOT}/Modules/TestBigEndian.cmake)
TEST_BIG_ENDIAN(CMAKE_WORDS_BIGENDIAN)

find_package (VTK REQUIRED)
include (${VTK_USE_FILE})
find_package (ITK REQUIRED)
include (${ITK_USE_FILE})

CONFIGURE_FILE(
  ${CMAKE_CURRENT_SOURCE_DIR}/FreeSurferConfigure.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/FreeSurferConfigure.h
)

set(MyModule_SOURCE 
main.cxx 
MeshIO.cxx
itkMeshTovtkPolyData.cxx
vtkFSSurfaceReader.cxx 
vtkFSSurfaceScalarReader.cxx
vtkFSIO.cxx
KWMeshVisuAttributeIO.cxx
)

set (VTK_LIBRARIES
  vtkRendering
  vtkGraphics
  vtkImaging
  vtkIO
  vtkCommon
  vtkHybrid
  vtkFiltering 
) 

FIND_PACKAGE(GenerateCLP REQUIRED)
INCLUDE(${GenerateCLP_USE_FILE})

GENERATECLP(MyModule_SOURCE PostprocessingTPS.xml)
add_executable(ParticleCorrespondencePostprocessingTPS ${MyModule_SOURCE})
target_link_libraries( ParticleCorrespondencePostprocessingTPS ${ITK_LIBRARIES} ${VTK_LIBRARIES})
