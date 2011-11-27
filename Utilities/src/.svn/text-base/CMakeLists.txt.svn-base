cmake_minimum_required(VERSION 2.4)
PROJECT(Utilities)
	
INCLUDE (${CMAKE_ROOT}/Modules/FindITK.cmake) 

IF (USE_ITK_FILE) 
  INCLUDE(${USE_ITK_FILE}) 
ENDIF(USE_ITK_FILE) 


# Add the executable with our source file.

#ADD_EXECUTABLE(Demons_ITK_Filter_3d Demons_ITK_Filter_3d.cxx)
#TARGET_LINK_LIBRARIES(Demons_ITK_Filter_3d ITKAlgorithms ITKBasicFilters ITKCommon ITKIO)


ADD_EXECUTABLE(ResamplingImageBySpacing ResamplingImageBySpacing.cxx)
TARGET_LINK_LIBRARIES(ResamplingImageBySpacing ITKAlgorithms ITKBasicFilters ITKCommon ITKIO)

ADD_EXECUTABLE(ThresholdImage ThresholdImage.cxx)
TARGET_LINK_LIBRARIES(ThresholdImage ITKAlgorithms ITKBasicFilters ITKCommon ITKIO)

ADD_EXECUTABLE(BinaryErode BinaryErode.cxx)
TARGET_LINK_LIBRARIES(BinaryErode ITKAlgorithms ITKBasicFilters ITKCommon ITKIO)

ADD_EXECUTABLE(BinaryDilate BinaryDilate.cxx)
TARGET_LINK_LIBRARIES(BinaryDilate ITKAlgorithms ITKBasicFilters ITKCommon ITKIO)

ADD_EXECUTABLE(ConnectedThreshold ConnectedThreshold.cxx)
TARGET_LINK_LIBRARIES(ConnectedThreshold ITKAlgorithms ITKBasicFilters ITKCommon ITKIO)

ADD_EXECUTABLE(NeighborhoodConnectedThreshold NeighborhoodConnectedThreshold.cxx)
TARGET_LINK_LIBRARIES(NeighborhoodConnectedThreshold ITKAlgorithms ITKBasicFilters
ITKCommon ITKIO)

ADD_EXECUTABLE(Diffusion Diffusion.cxx)
TARGET_LINK_LIBRARIES(Diffusion ITKAlgorithms ITKBasicFilters ITKCommon ITKIO)

ADD_EXECUTABLE(HoleFilling HoleFilling.cxx)
TARGET_LINK_LIBRARIES(HoleFilling ITKAlgorithms ITKBasicFilters ITKCommon ITKIO)


ADD_EXECUTABLE(RollingBall RollingBall.cxx)
TARGET_LINK_LIBRARIES(RollingBall ITKAlgorithms ITKBasicFilters ITKCommon ITKIO)

ADD_EXECUTABLE(DistanceMap DistanceMap.cxx)
TARGET_LINK_LIBRARIES(DistanceMap ITKAlgorithms ITKBasicFilters ITKCommon ITKIO)

ADD_EXECUTABLE(SecondOrderGradient SecondOrderGradient.cxx)
TARGET_LINK_LIBRARIES(SecondOrderGradient  ITKBasicFilters ITKCommon ITKIO)

ADD_EXECUTABLE(TranslationRegistration TranslationRegistration.cxx)
TARGET_LINK_LIBRARIES(TranslationRegistration  ITKAlgorithms ITKBasicFilters ITKCommon ITKIO)

ADD_EXECUTABLE(CropAndPad CropAndPad.cxx)
TARGET_LINK_LIBRARIES(CropAndPad ITKAlgorithms ITKBasicFilters ITKCommon ITKIO)
