#include "itkVectorImage.h"
#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "vectorITKConversions.h"
#include <iostream>
#include <sstream>
#include <string>
#include <math.h>
#include <stdlib.h>
#include "HField3DUtils.h"
#include "HField3DIO.h"
#include "Array3D.h"
#include "Array3DIO.h"
#include "Array3DUtils.h"

using namespace itk;

//TODO: Intensity Windwoing
int main( int argc, char *argv[] ){
	if (argc < 3 )
	{       
		std::cerr << "Missing Parameters " << std::endl;
		std::cerr << "Usage: " << argv[0];
		std::cerr << " inputImageFileName hFieldFileName  outputVectorImageFileName[default=inputImageFIleName.withHField.mhd]" << std::endl;
		std::cerr << " This program is used to combine  3D image and its displacement vector field togather into a vector 3D image."<<std::endl;
		return 0;
	}

	//---------------------------------------------  DEFINE IMAGE TYPE -----------------------------------------------
	const unsigned int Dimension = 3;
	typedef itk::Image<float, Dimension > ImageType;
	typedef vnl_vector<float> VoxelType;

	char *  imageFileName = argv[1];
	std::string hFieldFileName = argv[2];
	std::string outputVectorImageFileName;
	if (argc >3){
		outputVectorImageFileName = argv[3];
	}else{
		outputVectorImageFileName = imageFileName;
		outputVectorImageFileName.insert(outputVectorImageFileName.size()-4,".withHField");

	}

	//load Image
	//TODO: Intensity Windowing
	typedef itk::ImageFileReader<ImageType> ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(imageFileName);

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
	ImageType::Pointer iImage = reader->GetOutput();


	// load hField
	Array3D<Vector3D<float> > hField;
	Vector3D<unsigned int> hSize;
	Vector3D<double> hOrigin, hSpacing;

	std::cerr << "Loading hField: " << hFieldFileName.c_str() << "...";
	HField3DIO::readMETAHeader(hSize,hOrigin,hSpacing,
		hFieldFileName.c_str());
	HField3DIO::readMETA(hField, hFieldFileName.c_str());
	std::cerr << "DONE" << std::endl;
	std::cerr << "  Dimensions: " << hSize << std::endl;
	std::cerr << "  Origin    : " << hOrigin << std::endl;
	std::cerr << "  Spacing   : " << hSpacing << std::endl;

	//combine the image with it's diffeo field
	typedef itk::VectorImage<float> ITKImageType;
	ITKImageType::Pointer vImage = ITKImageType::New();
	ITKImageType::SizeType size;
	size= iImage->GetBufferedRegion().GetSize();


	ITKImageType::RegionType region;
	region.SetSize(size);
	unsigned int vectorDim =4;
	vImage->SetVectorLength(vectorDim);
	vImage->SetRegions(region);
	vImage->Allocate();

	for (unsigned int k = 0; k < size[2]; k++){
		for (unsigned int j = 0; j < size[1]; j++){
			for (unsigned int i = 0; i < size[0]; i++)
			{
				ITKImageType::IndexType ind;
				ind[0] = i;
				ind[1] = j;
				ind[2] = k;
				ITKImageType::PixelType w(vectorDim);
				w[0]= iImage->GetPixel(ind);
				w[1] =hField(i, j, k).x;
				w[2] =hField(i, j, k).y;
				w[3] =hField(i, j, k).z;

				vImage->SetPixel(ind, w);
			}
		}
	}


	ITKImageType::SpacingType spacing;
	spacing = iImage->GetSpacing();

	vImage->SetSpacing(spacing);

	ITKImageType::PointType origin;
	origin = iImage->GetOrigin();

	vImage->SetOrigin(origin);


	//write the vector image
	typedef itk::VectorImage<float> ITKImageType;
	typedef itk::ImageFileWriter<ITKImageType> WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName(outputVectorImageFileName);
	writer->SetInput(vImage);
	writer->UseCompressionOn();
	writer->Update();

	return 0;
}
