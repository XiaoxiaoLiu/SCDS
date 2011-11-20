	
#include <iostream>

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryBallStructuringElement.h"



//binary erode operator, for 3D binary images
int main( int argc, char *argv[] )
{
	if (argc < 4 )
	{
		std::cerr << "Missing Parameters " << std::endl;
		std::cerr << "Usage: " << argv[0];
		std::cerr << " inputImageFileName outputImageFileName radius" << std::endl;
		return 1;
	}


	//---------------------------------------------  DEFINE IMAGE TYPE -----------------------------------------------
	const unsigned int Dimension = 3;

	typedef itk::Image<float, Dimension > ImageType;

	typedef itk::Image<unsigned short, Dimension> BinaryImageType;
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

	//--------------------------------------------- ERODE--------------------------------------------------------

	typedef itk::BinaryBallStructuringElement<ImageType::PixelType, Dimension > StructuringElementType;
	typedef itk::BinaryErodeImageFilter< ImageType, BinaryImageType, StructuringElementType > ErodeFilterType;
	ErodeFilterType::Pointer erode = ErodeFilterType::New();

	StructuringElementType structuringElement;
	structuringElement.SetRadius( atoi(argv[3]) ); // 3x3 structuring element

	structuringElement.CreateStructuringElement();
	erode->SetKernel( structuringElement );
	erode->SetInput(reader->GetOutput());
	erode->SetErodeValue(1);

	erode->Update();


	//--------------------------------------------    OUTPUT -------------------------------------------------
	typedef itk::ImageFileWriter<BinaryImageType> FileWriterType;
	FileWriterType:: Pointer fileWriter = FileWriterType::New();
	fileWriter->SetFileName(argv[2]);

	fileWriter->SetInput(erode->GetOutput());


	try 
	{
		fileWriter->Update(); 
	}
	catch (itk::ExceptionObject &e)
	{
		std::cerr << e << std::endl;
	}

return 0;
}
