#include "itkVectorImage.h"
#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkImageWriter.h"
#include "itkImageReader.h"
#include "vectorITKConversion.h"
#include <iostream>
#include <sstream>

#include <math.h>
#include <stdlib.h>
#include "HField3DUtils.h"
#include "HField3DIO.h"
#include "Array3D.h"
#include "Array3DIO.h"
#include "Array3DUtils.h"

using namespace itk;


int main( int argc, char *argv[] ){


	if (argc < 4 )
	{       
		std::cerr << "Missing Parameters " << std::endl;
		std::cerr << "Usage: " << argv[0];
		std::cerr << " inputImageFileName hFieldFileName  outputVectorImageFileName " << std::endl;
		std::cerr << " This program is used to combine  3D image and its displacement vector field togather into a vector 3D image.<<std::Endl;
		//	std::cerr <<"   The first input is the 3D image, the second input is its displacement vector field.	"
		return 1;
	}

	//---------------------------------------------  DEFINE IMAGE TYPE -----------------------------------------------
	const unsigned int Dimension = 3;
	typedef itk::Image<float, Dimension > ImageType;
	typedef vnl_vector<float> VoxelType;




	char *  imageFileName = argv[1];
	char * hFieldFileName = argv[2];


	//loadImage
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
	if (!doTranslation)
	{
		std::cerr << "Loading hField: " << hFieldFileName << "...";
		HField3DIO::readMETAHeader(hSize,hOrigin,hSpacing,
				hFieldFileName.c_str());
		HField3DIO::readMETA(hField, hFieldFileName.c_str());
		std::cerr << "DONE" << std::endl;
		std::cerr << "  Dimensions: " << hSize << std::endl;
		std::cerr << "  Origin    : " << hOrigin << std::endl;
		std::cerr << "  Spacing   : " << hSpacing << std::endl;


	}




	//combine
	//
	//vImage 
	typedef itk::VectorImage<float> ITKImageType;
	ITKImageType::Pointer vImage = ITKImageType::New();
	ITKImageType::SizeType size;
	size[0] = iImage->getSizeX();
	size[1] = iImage->getSizeY();
	size[2] = iImage->getSizeZ();
	unsigned int vectorDim =4;
	vImage->SetVectorLength(vectorDim);
	vImage->SetRegions(region);
	vImage->Allocate();

	for (unsigned int k = 0; k < size[2]; k++)
		for (unsigned int j = 0; j < size[1]; j++)
			for (unsigned int i = 0; i < size[0]; i++)
			{
				ITKImageType::IndexType ind;
				ind[0] = i;
				ind[1] = j;
				ind[2] = k;
				ITKImageType::PixelType w(vectorDim);
				w[0]= iImage->get(i, j, k);
				for (unsigned int d = 0; d < 3; d++)
					w[d] =hFiled->get(i, j, k,d-1);



				itkimg->SetPixel(ind, w);
			}

	ITKImageType::SpacingType spacing;
	spacing[0] = img->getSpacingX();
	spacing[1] = img->getSpacingY();
	spacing[2] = img->getSpacingZ();
	itkimg->SetSpacing(spacing);

	ITKImageType::PointType origin;
	origin[0] = img->getOriginX();
	origin[1] = img->getOriginY();
	origin[2] = img->getOriginZ();
	itkimg->SetOrigin(origin);










	//write the vector image
	typedef itk::VectorImage<float> ITKImageType;

	typedef itk::ImageFileWriter<ITKImageType> WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName(fn);
	writer->SetInput(vImage);
	writer->UseCompressionOn();
	writer->Update();




	return 0;
}
