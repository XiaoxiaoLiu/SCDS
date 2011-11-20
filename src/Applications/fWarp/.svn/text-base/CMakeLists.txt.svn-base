CMAKE_MINIMUM_REQUIRED(VERSION 2.4)
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
  # prefer ${LAPACK_LIBRARIES}
  #lapack
  #g2c
  ${FFTWF_LIB}
)


