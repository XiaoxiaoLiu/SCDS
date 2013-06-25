#include <iostream>

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkIntensityWindowingImageFilter.h"


int main( int argc, char *argv[] )
{
	if (argc < 6 )
	{
		std::cerr << "Missing Parameters " << std::endl;
		std::cerr << "Usage: " << argv[0];
		std::cerr << " inputImageFileName outputImageFileName input_lowerbound  input_upperbound output_lowerbound output_upperbound"<<std::endl; 
		return 1;
	}


	const unsigned int Dimension = 3;
	typedef float  PixelType;

	typedef itk::Image<float, Dimension > InImageType;
	typedef itk::Image<PixelType,Dimension> OutImageType;



    typedef itk::ImageFileReader<InImageType> ReaderType;
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


	typedef itk::IntensityWindowingImageFilter< InImageType,OutImageType > IntensityWindowingImageFilterType;

	IntensityWindowingImageFilterType::Pointer windowing = IntensityWindowingImageFilterType::New();

	InImageType::PixelType input_lowerBound = atoi(argv[3]);
	InImageType::PixelType input_upperBound = atoi(argv[4]);

	OutImageType::PixelType output_lowerBound = atoi(argv[5]);
	OutImageType::PixelType output_upperBound = atoi(argv[6]);

  windowing->SetInput(reader->GetOutput());
  windowing->SetWindowMinimum(input_lowerBound);
  windowing->SetWindowMaximum(input_upperBound);
  windowing->SetOutputMinimum(output_lowerBound);
  windowing->SetOutputMaximum(output_upperBound);
  windowing->Update();


	typedef itk::ImageFileWriter<OutImageType> imageFileWriterType;
	imageFileWriterType:: Pointer imageWriter = imageFileWriterType::New();
	imageWriter->SetFileName(argv[2]);

	imageWriter->SetInput(windowing->GetOutput());


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







