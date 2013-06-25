#include "itkImage.h"

#include "itkTranslationTransform.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkResampleImageFilter.h"


int main( int argc, char *argv[] )
{
	if( argc < 7 )
	{
		std::cerr << "Missing Parameters " << std::endl;
		std::cerr << "Usage: " << argv[0];
		std::cerr << " inputImageFile  ";
		std::cerr << " outputImagefile ";
		std::cerr << " tx ty tz  ResamplingDefaultPixelValue "<< std::endl;
		return EXIT_FAILURE;
	}

	const    unsigned int    Dimension = 3;
	typedef  float    PixelType;

	typedef itk::Image< PixelType, Dimension >  InputImageType;

	typedef itk::ImageFileReader< InputImageType  > InputImageReaderType;
  InputImageReaderType::Pointer  inputImageReader  = InputImageReaderType::New();
	inputImageReader->SetFileName(  argv[1] );
	inputImageReader->Update();
  InputImageType ::Pointer inputImage = inputImageReader->GetOutput();

	typedef itk::TranslationTransform< double,Dimension > TransformType;
	TransformType::Pointer  transform= TransformType::New();
  TransformType::OutputVectorType translation;
	translation[0]  = atof(argv[3]);
  translation[1]  = atof(argv[4]);
  translation[2]  = atof(argv[5]);
  transform->Translate(translation);


	typedef itk::ResampleImageFilter<
		InputImageType,
		InputImageType >    ResampleFilterType;
	ResampleFilterType::Pointer resample = ResampleFilterType::New();

	resample->SetTransform( transform);
	resample->SetInput( inputImage);

	resample->SetSize( inputImage->GetLargestPossibleRegion().GetSize() );
	resample->SetOutputOrigin(  inputImage->GetOrigin() );
	resample->SetOutputSpacing( inputImage->GetSpacing() );
	resample->SetOutputDirection( inputImage->GetDirection() );
	resample->SetDefaultPixelValue(atof(argv[6]));

	typedef itk::ImageFileWriter< InputImageType >  WriterType;

	WriterType::Pointer  writer =  WriterType::New();

	writer->SetFileName( argv[2] );

	writer->SetInput( resample->GetOutput()   );

	try
	{
		writer->Update();
	}
	catch( itk::ExceptionObject & excp )
	{
		std::cerr << "ExceptionObject while writing the resampled image !" << std::endl;
		std::cerr << excp << std::endl;
		return EXIT_FAILURE;
	}

  return EXIT_SUCCESS;
}


