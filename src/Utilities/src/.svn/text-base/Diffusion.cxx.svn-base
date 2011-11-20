
/*=========================================================================


=========================================================================*/


#ifdef _MSC_VER
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif




/*NOTES:
  To RUN:
  
 */





#include <iostream>

#include "itkGradientAnisotropicDiffusionImageFilter.h"
#include "itkIntensityWindowingImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkImageRegionIterator.h"


int main( int argc, char *argv[] )
{

	if (argc < 5 )
	{
		std::cerr << "Missing Parameters " << std::endl;
		std::cerr << "Usage: " << argv[0];
		std::cerr << " inputImageFileName  outputImageFileName";
	    std::cerr << " conductanceTerm     diffusionIterations" << std::endl;

		return 1;
	}


	//---------------------------------------------  DEFINE IMAGE TYPE -----------------------------------------------
	const unsigned int Dimension = 3;

	typedef itk::Image<unsigned short, Dimension > ImageType;

	char  * fileName = new char[100];


	//read image

	typedef itk::ImageFileReader<ImageType> ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	const char * fn = argv[1];
	reader->SetFileName(fn);

	try
	{
		reader->Update();
		//reader->GetOutput()->Print(std::cout);
	}
	catch (itk::ExceptionObject &ex)
	{
		std::cout << ex;
		return EXIT_FAILURE;
	}



	/*
	   ----- anisotropic diffusion

	 */

	typedef itk::GradientAnisotropicDiffusionImageFilter<ImageType,ImageType> 
		DiffusionFilterType;
	DiffusionFilterType::Pointer diffusion = DiffusionFilterType::New();
	diffusion->SetConductanceParameter( atof(argv[4]) );
	diffusion->SetNumberOfIterations( atoi(argv[5]) );
	diffusion->SetTimeStep(0.125);

	typedef itk::ImageFileWriter<ImageType> DiffusedFileWriterType;
	DiffusedFileWriterType:: Pointer diffusedWriter = DiffusedFileWriterType::New();
	fileName = strcpy (fileName, argv[2]);
	diffusedWriter->SetFileName(fileName);



	diffusion->SetInput(reader->GetOutput());
	diffusedWriter->SetInput(diffusion->GetOutput());

	try 
	{
		diffusedWriter->Update(); 
		std::cerr <<"written:"<<fileName<<std::endl;

	}
	catch (itk::ExceptionObject &e)
	{
		std::cerr << e << std::endl;
	}


	return 0;

}

