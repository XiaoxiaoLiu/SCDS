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


