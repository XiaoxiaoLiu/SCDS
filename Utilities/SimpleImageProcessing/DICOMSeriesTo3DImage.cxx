
#include <iostream>           
#include <fstream>       
#include <stdio.h>                    
#include "itkImage.h"                   
#include "itkImageFileWriter.h"

#include "itkDICOMImageIO2Factory.h"
#include "itkDICOMImageIO2.h"
#include "itkImageSeriesReader.h"
#include "itkDICOMSeriesFileNames.h"

#include "itkImageFileReader.h"     




int main(int argc, char *argv[])        
{

	if ( argc < 3 )     
	{ 
		std::cout << "Useage ex:  DICOMSeriesTo3DImage inputFilePath outputFileName(with file extension '.mhd/.hdr...')  " << std::endl;
		return 1;
	}           

	//---------------------------------------------  DEFINE IMAGE TYPE -----------------------------------------------
	const unsigned int Dimension = 3;

	typedef itk::Image<signed short, Dimension> InputImageType;


	//read 3D DICOM 
	typedef itk::ImageSeriesReader<InputImageType> ReaderType;
	itk::DICOMImageIO2::Pointer io = itk::DICOMImageIO2::New();

	//Get the DICOM filenames from the directory
	itk::DICOMSeriesFileNames::Pointer names = itk::DICOMSeriesFileNames::New();
	names->SetDirectory(argv[1]);

	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileNames(names->GetFileNames());
	reader->SetImageIO(io);

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


	typedef itk::ImageFileWriter<InputImageType> OutputFileWriterType;
	OutputFileWriterType:: Pointer OutputWriter = OutputFileWriterType::New();
	OutputWriter->SetFileName(argv[2]);

	OutputWriter->SetInput(reader->GetOutput());	
	try 
	{
		OutputWriter->Update(); 
		std::cerr <<"written:"<<argv[2]<<std::endl;
	}
	catch (itk::ExceptionObject &e)
	{
		std::cerr << e << std::endl;
	}


	return 0;
} 





