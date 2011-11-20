
#include <iostream>


#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkConnectedThresholdImageFilter.h"



int main( int argc, char *argv[] )
{
	if (argc < 7 )
	{
		std::cerr << "Missing Parameters " << std::endl;
		std::cerr << "Usage: " << argv[0];
		std::cerr << " inputImageFileName outputImageFileName lowThreshold_regionGrow[int]  highThreshold_regionGrow[int]  seed1_x seed1_y seed1_z [seed2_x seed2_y seed2_z]  " << std::endl;
		return 1;
	}

	//---------------------------------------------  DEFINE IMAGE TYPE -----------------------------------------------
	const unsigned int Dimension = 3;

	typedef itk::Image<float, Dimension > ImageType;
	typedef itk::Image<short,Dimension> OutImageType;

	//---------------------------------------------      READ  IMAGE ---------------------------------------------------

	typedef itk::ImageFileReader<ImageType> ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(argv[1]);

	try
	{
		reader->Update();
	//	reader->GetOutput()->Print(std::cout);
	}
	catch (itk::ExceptionObject &ex)
	{
		std::cout << ex;
		return EXIT_FAILURE;
	}

	//---------------------------------------------      PIPELINE --------------------------------------------------
	
	//region grow segmentation
	typedef itk::ConnectedThresholdImageFilter< ImageType,OutImageType > ConnectedFilterType;

	ConnectedFilterType::Pointer connectedThreshold = ConnectedFilterType::New();

	ImageType::PixelType lowerThreshold = atoi(argv[3]);	 //1100-1024 = 76    (76-(-350))/750*255 =  145
	ImageType::PixelType upperThreshold = atoi(argv[4]);   	 //1250-1024  =226   (226-(-350))/750*255 =  196 


	itk::Index<3> seedIndex = {{ atoi(argv[5]),  atoi(argv[6]),  atoi(argv[7])}};
	connectedThreshold->SetInput( reader->GetOutput() );
	connectedThreshold->SetLower( lowerThreshold );
	connectedThreshold->SetUpper( upperThreshold );
	connectedThreshold->SetReplaceValue(1);	
	connectedThreshold->AddSeed( seedIndex );
        std::cerr <<"seedPixel value ="<< reader->GetOutput() ->GetPixel(seedIndex) << std::endl;


	if (argc>8){//second seeds
        	itk::Index<3> seed2Index = {{ atoi(argv[8]),  atoi(argv[9]),  atoi(argv[10])}};
		connectedThreshold->AddSeed(seed2Index);
		 std::cerr <<"seed2Pixel value ="<< reader->GetOutput() ->GetPixel(seed2Index) << std::endl;

	}


	typedef itk::ImageFileWriter<OutImageType> regionGrowFileWriterType;
	regionGrowFileWriterType:: Pointer regionGrowWriter = regionGrowFileWriterType::New();
	regionGrowWriter->SetFileName(argv[2]);

	regionGrowWriter->SetInput(connectedThreshold->GetOutput());


	try 
	{
		regionGrowWriter->Update(); 
	}
	catch (itk::ExceptionObject &e)
	{
		std::cerr << e << std::endl;
	}



	return 0;

}






