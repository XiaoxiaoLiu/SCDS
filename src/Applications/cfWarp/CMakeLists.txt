CMAKE_MINIMUM_REQUIRED(VERSION 2.8.3)

SET(cfWarp_SRCS
  main.cxx
)

ADD_EXECUTABLE(cfWarp ${cfWarp_SRCS})

TARGET_LINK_LIBRARIES(cfWarp
  Algorithms
  UtilitiesDataTypes
  Numerics
  DataTypes
  Base
  ITKIO
  ITKCommon
  UtilDownsampleFilter
  ${FFTWF_LIB}
)

ENABLE_TESTING()
set(TESTING_OUTPUT_PATH "/home/xiaoxiao/work/tmp")
set(DATA_PATH "/home/xiaoxiao/work/src/SCDS/data")
ADD_TEST(cfWarpTest ${EXECUTABLE_OUTPUT_PATH}/cfWarp 
"--outputImageFilenamePrefix=${TESTING_OUTPUT_PATH}/deformed_" 
"--outputHFieldFilenamePrefix=${TESTING_OUTPUT_PATH}/hfield_" 
"--inputHFieldFile=${DATA_PATH}/hfield-cinePhase00.mhd" 
"--HWeights"  "1.0"  "1.0"  "1.0" 
"--scaleLevel=4 --numberOfIterations=100" 
"--scaleLevel=2 --numberOfIterations=50" 
"--scaleLevel=1 --numberOfIterations=25" 
"${DATA_PATH}/phase00.mhd" 
"${DATA_PATH}/phase50.mhd" 
)
#a pair of images

