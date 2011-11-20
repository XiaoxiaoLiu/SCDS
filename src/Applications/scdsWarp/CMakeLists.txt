SET(scdsWarp_SRCS
  main.cxx
 )
 
#SET(makeVectorImage_SRCS
#MakeVectorImage.cxx
#)

ADD_EXECUTABLE(scdsWarp "${scdsWarp_SRCS}")
#ADD_EXECUTABLE(makeVectorImage "${makeVectorImage_SRCS}")

TARGET_LINK_LIBRARIES(scdsWarp
  Algorithms
  UtilitiesDataTypes
  UtilDownsampleFilter
  DataTypes
  Base
  ITKIO
  ITKCommon
  ${FFTWF_LIB}
  ${PTHREADS_LIBRARY}
)


#TARGET_LINK_LIBRARIES(makeVectorImage
# Algorithms
# UtilitiesDataTypes
# DataTypes
# Base
# ITKIO
# ITKCommon
#)
