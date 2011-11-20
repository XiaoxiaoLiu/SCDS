
#include <iostream>


#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkBinaryThresholdImageFilter.h"



int main( int argc, char *argv[] )
{
	if (argc < 5 )
	{
		std::cerr << "Missing Parameters " << std::endl;
		std::cerr << "Usage: " << argv[0];
		std::cerr << " imageDimension inputImageFileName outputImageFileName lowThreshold_regionGrow[int]  highThreshold_regionGrow[int] " << std::endl;
		return 1;
	}




	if (atoi(argv[1]) !=3 )
	{
		std::cerr<<"Implemented for 3D image right now!"<<std::endl;
		return 0;	
	}

	const unsigned int Dimension = 3;
	typedef unsigned short PixelType;

	typedef itk::Image<float, Dimension > InImageType;
	typedef itk::Image<PixelType,Dimension> OutImageType;



    typedef itk::ImageFileReader<InImageType> ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(argv[2]);

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





	typedef itk::BinaryThresholdImageFilter< InImageType,OutImageType > BinaryThresholdFilterType;

	BinaryThresholdFilterType::Pointer binaryThreshold = BinaryThresholdFilterType::New();

	InImageType::PixelType lowerThreshold = atoi(argv[4]);	
	InImageType::PixelType upperThreshold = atoi(argv[5]);   	 

	binaryThreshold->SetInput( reader->GetOutput() );	
	binaryThreshold->SetOutsideValue( 0 );
	binaryThreshold->SetInsideValue( 1 );


	binaryThreshold->SetLowerThreshold( lowerThreshold );
	binaryThreshold->SetUpperThreshold( upperThreshold );


	typedef itk::ImageFileWriter<OutImageType> imageFileWriterType;
	imageFileWriterType:: Pointer imageWriter = imageFileWriterType::New();
	imageWriter->SetFileName(argv[3]);

	imageWriter->SetInput(binaryThreshold->GetOutput());


	try 
	{
		imageWriter->Update(); 
	}
	catch (itk::ExceptionObject &e)
	{
		std::cerr << e << std::endl;
	}



	return 0;

}







