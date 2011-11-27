
#include <iostream>
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkVotingBinaryIterativeHoleFillingImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryDilateImageFilter.h"



int main( int argc, char *argv[] )
{

	if (argc < 9 )
	{
		std::cerr << "Missing Parameters " << std::endl;
		std::cerr << "Usage: " << argv[0];
		std::cerr << " inputImageFileName outputImageFileName erodeRadius dilateRadius"<<std::endl;
		std::cerr << " holeFillingRadius_x holeFillingRadius_y holeFillingRadius_z  majorityThreshold   maximumIterations"<<std::endl;
		std::cerr <<"[the number of foreground neighbors should be at  least (3x3 -1 )/2 + majorityThreshold.]"<<std::endl;

		return 1;
	}

	const unsigned int Dimension = 3;
	typedef unsigned short PixelType;
	typedef itk::Image< PixelType, Dimension > ImageType;

	typedef itk::ImageFileReader< ImageType > ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName( argv[1] );

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

	// dilate after erode: open
	typedef itk::BinaryBallStructuringElement<PixelType, Dimension > StructuringElementType;
	StructuringElementType structuringElement;

	typedef itk::BinaryErodeImageFilter< ImageType, ImageType, StructuringElementType > ErodeFilterType;
	ErodeFilterType::Pointer erode = ErodeFilterType::New();

	structuringElement.SetRadius( atoi(argv[3]) );
	structuringElement.CreateStructuringElement();
	erode->SetKernel( structuringElement );
	erode->SetInput(reader->GetOutput());
	erode->SetErodeValue(1);
	erode->Update();

	typedef itk::BinaryDilateImageFilter< ImageType, ImageType, StructuringElementType > DilateFilterType;
	DilateFilterType::Pointer dilate = DilateFilterType::New();
	StructuringElementType structuringElement2;

	structuringElement2.SetRadius( atoi(argv[4]) ); 
	structuringElement2.CreateStructuringElement();
	dilate->SetKernel( structuringElement2 );
	dilate->SetInput(erode->GetOutput());
	dilate->SetDilateValue(1);
	dilate->Update();



	typedef itk::ImageFileWriter<ImageType> segFileWriterType;
	segFileWriterType:: Pointer segWriter= segFileWriterType::New();
	segWriter->SetFileName(argv[2]);


	//holefilling
	//
	typedef itk::VotingBinaryIterativeHoleFillingImageFilter<ImageType> HoleFillingFilterType;
	HoleFillingFilterType::Pointer holeFillingFilter = HoleFillingFilterType::New();
	ImageType::SizeType holeRadius;
	holeRadius[0] = atoi(argv[5]);
	holeRadius[1] = atoi(argv[6]);
	holeRadius[2] = atoi(argv[7]);

	holeFillingFilter->SetRadius(holeRadius);
	holeFillingFilter->SetBackgroundValue(0);
	holeFillingFilter->SetForegroundValue(1);


	holeFillingFilter->SetMajorityThreshold(atoi(argv[8]));
	if (argc>9){
		holeFillingFilter->SetMaximumNumberOfIterations(atoi(argv[9]));
	}else{
		holeFillingFilter->SetMaximumNumberOfIterations(20);
	}


	holeFillingFilter->SetInput(dilate->GetOutput());
	holeFillingFilter->Update();

	if (atoi(argv[5])>0){
		segWriter->SetInput( holeFillingFilter->GetOutput());
	}else{
		segWriter->SetInput(dilate->GetOutput());
	}



	try 
	{
		segWriter->Update(); 
	}
	catch (itk::ExceptionObject &e)
	{
		std::cerr << e << std::endl;
	}


}

