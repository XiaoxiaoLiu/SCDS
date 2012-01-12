
#include <iostream>
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkResampleImageFilter.h"
#include "itkAffineTransform.h"
#include "itkNearestNeighborInterpolateImageFunction.h"


int main( int argc, char *argv[] )
{
	if (argc < 7 )
	{
		std::cerr << "Missing Parameters " << std::endl;
		std::cerr << "Usage: " << argv[0];
		std::cerr << " imageDimension inputImageFileName outputImageFileName Spacing_x Spacing_y Spacing_z " << std::endl;
		return 1;
	}

	
	if (atoi(argv[1]) !=3 )
	{
		std::cerr<<"Implemented for 3D image only!"<<std::endl;
		return 0;	
	}

	const unsigned int Dimension = 3;
	typedef float PixelType;

	typedef itk::Image<float, Dimension > ImageType;


	typedef itk::ImageFileReader<ImageType> ReaderType;
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


	ImageType::Pointer image = reader->GetOutput();

	typedef ImageType::IndexType                   ImageIndexType;
	typedef ImageType::Pointer                     ImagePointerType;
	typedef ImageType::RegionType                  ImageRegionType;
	typedef ImageType::SizeType                    ImageSizeType;
	typedef ImageType::SpacingType                 ImageSpacingType;


	ImageSizeType size = image->GetLargestPossibleRegion().GetSize();
	ImageType::PointType origin = image->GetOrigin();
	ImageSpacingType spacing = image->GetSpacing();


	ImageSpacingType outSpacing;
	outSpacing[0]=atof(argv[4]);
	outSpacing[1]=atof(argv[5]);
	outSpacing[2]=atof(argv[6]);

	// Create an identity affine transformation

	typedef itk::AffineTransform<double,Dimension>   AffineTransformType;
	AffineTransformType::Pointer aff = AffineTransformType::New();
	typedef  itk::Vector<double,Dimension> VectorType;
	VectorType scale;
     scale[0] = spacing[0]/outSpacing[0]; 
	 scale[1] = spacing[1]/outSpacing[1]; 
	 scale[2]= spacing[2]/outSpacing[2]; 

	aff->Scale(scale);

	ImageSizeType outSize;
	outSize[0] = size[0]*spacing[0]/outSpacing[0];
	outSize[1] = size[1]*spacing[1]/outSpacing[1];
	outSize[2] = size[2]*spacing[2]/outSpacing[2];


	// Create a linear interpolation image function
	typedef itk::NearestNeighborInterpolateImageFunction<
		ImageType, double > InterpolatorType;
	InterpolatorType::Pointer interp = InterpolatorType::New();
	interp->SetInputImage(image);


	// Create and configure a resampling filter
	itk::ResampleImageFilter< ImageType, ImageType >::Pointer resample;
	resample = itk::ResampleImageFilter< ImageType, ImageType >::New();
	resample->SetInput(image);
	resample->SetDefaultPixelValue(0);

	resample->SetSize(outSize);
	resample->SetTransform(aff);
	resample->SetInterpolator(interp);

	ImageIndexType index;
	index.Fill( 0 );
	resample->SetOutputStartIndex( index );
	resample->SetOutputOrigin( origin );
	resample->SetOutputSpacing( outSpacing );

	
	try
	{
	// Run the resampling filter
	resample->Update();
	}
	catch (itk::ExceptionObject &e)
	{
		std::cerr << e << std::endl;
	}


	//output
	typedef itk::ImageFileWriter<ImageType> imageFileWriterType;
	imageFileWriterType:: Pointer imageWriter = imageFileWriterType::New();
	imageWriter->SetFileName(argv[3]);

	imageWriter->SetInput(resample->GetOutput());

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







