ADD_EXECUTABLE(ShapeWorksGroom ShapeWorksGroom.cxx)
TARGET_LINK_LIBRARIES(ShapeWorksGroom ITKIO ITKNumerics ITKBasicFilters ITKCommon params)
INSTALL(TARGETS ShapeWorksGroom   RUNTIME DESTINATION .)
