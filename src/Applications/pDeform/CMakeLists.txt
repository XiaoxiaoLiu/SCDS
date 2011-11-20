SET(pDeform_SRCS
  main.cxx
)

ADD_EXECUTABLE(pDeform "${pDeform_SRCS}")

TARGET_LINK_LIBRARIES(pDeform
  Algorithms
  DataTypes
  ITKIO
  ITKCommon
  ITKNumerics
  ITKAlgorithms
  ITKBasicFilters
  )


