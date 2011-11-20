SET(AtlasWerks_SRCS
  main.cxx
)

ADD_EXECUTABLE(AtlasWerks "${AtlasWerks_SRCS}")

TARGET_LINK_LIBRARIES(AtlasWerks
  UtilDownsampleFilter
  Algorithms
  Base
  ITKIO
  ITKCommon
  ${FFTWF_LIB}
  ${PTHREADS_LIBRARY}
  )


