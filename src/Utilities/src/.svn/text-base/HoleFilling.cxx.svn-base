
#include <iostream>
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkVotingBinaryIterativeHoleFillingImageFilter.h"


//marching cube contouring, get the binary volume

int main( int argc, char *argv[] )
{

	if (argc < 7 )
	{
		std::cerr << "Missing Parameters " << std::endl;
		std::cerr << "Usage: " << argv[0];
		std::cerr << " inputImageFileName outputImageFileName holeFillingRadius_x holeFillingRadius_y holeFillingRadius_z  majorityThreshold"<<std::endl;
		/*we set the majority value to 2, then we are requiring that the number of foreground neighbors should be at
		  leas/t (3x3 -1 )/2 + majority.*/

		return 1;
	}

	const unsigned int Dimension = 3;
	typedef float PixelType;
	typedef itk::Image< PixelType, Dimension > ImageType;

	typedef itk::ImageFileReader< ImageType > ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName( argv[1] );

	try
	{
		reader->Update();
		reader->GetOutput()->Print(std::cout);
	}
	catch (itk::ExceptionObject &ex)
	{	
		std::cout << ex;
		return EXIT_FAILURE;
	}

	typedef itk::VotingBinaryIterativeHoleFillingImageFilter<ImageType> HoleFillingFilterType;
	HoleFillingFilterType::Pointer holeFillingFilter = HoleFillingFilterType::New();
	ImageType::SizeType holeRadius;
	holeRadius[0] = atoi(argv[3]);
	holeRadius[1] = atoi(argv[4]);
	holeRadius[2] = atoi(argv[5]);

	holeFillingFilter->SetRadius(holeRadius);
	holeFillingFilter->SetBackgroundValue(0);
	holeFillingFilter->SetForegroundValue(1);

	holeFillingFilter->SetMajorityThreshold(atoi(argv[6]));
	holeFillingFilter->SetMaximumNumberOfIterations(20);


	holeFillingFilter->SetInput(reader->GetOutput());


	typedef itk::ImageFileWriter<ImageType> segFileWriterType;
	segFileWriterType:: Pointer segWriter= segFileWriterType::New();
	segWriter->SetFileName(argv[2]);
	segWriter->SetInput( holeFillingFilter->GetOutput());


	try 
	{
		segWriter->Update(); 
	}
	catch (itk::ExceptionObject &e)
	{
		std::cerr << e << std::endl;
	}


}

