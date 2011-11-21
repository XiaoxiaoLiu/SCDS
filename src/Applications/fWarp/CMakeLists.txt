CMAKE_MINIMUM_REQUIRED(VERSION 2.8.3)

SET(fWarp_SRCS
  main.cxx
)

ADD_EXECUTABLE(fWarp ${fWarp_SRCS})

TARGET_LINK_LIBRARIES(fWarp
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


