SCDS_USE_LAPACK()

SET(txApply_SRCS
  main.cxx
)

ADD_EXECUTABLE(txApply "${txApply_SRCS}")



TARGET_LINK_LIBRARIES(txApply

UtilitiesDataTypes
  Numerics
  DataTypes
  Base
  UtilPowerCrust
  CTI
  ITKIO
  ITKCommon
  ${LAPACK_LIBRARIES}
  vtkGraphics
  vtkCommon 
 
  )

