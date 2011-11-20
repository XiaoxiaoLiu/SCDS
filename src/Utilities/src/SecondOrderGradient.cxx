
#include <iostream>
#include "itkGradientRecursiveGaussianImageFilter.h"
#include "itkGradientImageFilter.h"
#include "itkGradientToMagnitudeImageFilter.h"
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"

int main(int argc, char *argv[])
{


	if (argc < 3)
	{
		std::cerr << "Missing Parameters " << std::endl;
		std::cerr << "Usage: " << argv[0];
		std::cerr << " inputImageFileName outputImageFileName sigma  [firstGradientMagFileName]" << std::endl;
		return 1;
	}

	char * inputFileName = argv[1];
	char * outputFileName = argv[2];
	float m_SigmaGaussian = atof(argv[3]);


	typedef float	PixelType;
	typedef itk::Image<PixelType, 3>	ImageType;

	typedef itk::CovariantVector<float, 3> GradientType;
	typedef itk::Image<GradientType, 3>   GradientImageType;



	typedef itk::ImageFileReader <ImageType> ImageReaderType;
	ImageReaderType ::Pointer reader = ImageReaderType::New();
	reader->SetFileName(inputFileName);
	try {
		reader->Update();
	}
	catch(itk::ExceptionObject &err){
		std::cerr<<err<<std::endl;
	}

	ImageType::Pointer image = reader->GetOutput();



	std::cout<<"generate second derivative image..."<<std::endl;

	typedef itk::GradientRecursiveGaussianImageFilter <ImageType,GradientImageType>  grgFilterType;
	typedef itk::GradientImageFilter  <ImageType, float, float>            gFilterType;
	typedef itk::GradientToMagnitudeImageFilter <GradientImageType, ImageType>       g2mFilterType;

	grgFilterType ::Pointer grgfilter = grgFilterType::New();
	gFilterType :: Pointer gfilter = gFilterType::New();
	g2mFilterType:: Pointer g2mfilter = g2mFilterType::New();

	grgfilter->SetInput(image);
	grgfilter->SetSigma(m_SigmaGaussian);

	g2mfilter->SetInput(grgfilter->GetOutput());

	gfilter->SetInput(g2mfilter->GetOutput());
	GradientImageType::Pointer m_gradientImage = gfilter->GetOutput();
	gfilter->Update();
	std::cout<<"derivative image done"<<std::endl;


	if (argc>4){
	   typedef itk::ImageFileWriter <ImageType> WriterType;
	   WriterType ::Pointer writer = WriterType::New();
	   writer->SetFileName(argv[4]);
	   writer->SetInput(g2mfilter->GetOutput());
	   try {
	   writer->Update();
	   }
	   catch(itk::ExceptionObject &err){
	   std::cerr<<err<<std::endl;
	   }
	}

	typedef itk::ImageFileWriter <GradientImageType> WriterType2;
	WriterType2 ::Pointer writer2 = WriterType2::New();
	writer2->SetFileName(outputFileName);
	writer2->SetInput(m_gradientImage);
	try {
		writer2->Update();
	}
	catch(itk::ExceptionObject &err){
		std::cerr<<err<<std::endl;
	}

return 0;

}
