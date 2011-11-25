FLTK_WRAP_UI(ShapeWorksView AnalyzeCorrespondenceGUI.fl)

ADD_EXECUTABLE(ShapeWorksView ShapeWorksView.cxx ShapeWorksViewApp.cxx ${ShapeWorksView_FLTK_UI_SRCS})
TARGET_LINK_LIBRARIES(ShapeWorksView VtkUtilities Utilities ITKParticleSystem ITKIO ITKNumerics ITKBasicFilters ITKCommon params ITKVtkFltk vtkIO vtkWidgets vtkRendering vtkImaging vtkFiltering vtkHybrid vtkGraphics vtkCommon ${FLTK_LIBRARIES})
INSTALL(TARGETS ShapeWorksView   RUNTIME DESTINATION .)
