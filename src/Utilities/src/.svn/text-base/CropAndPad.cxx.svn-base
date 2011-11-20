
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMetaImageIO.h"
#include "itkImageIOBase.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"


int main(int argc, char ** argv) {


	if(argc!=12) {
		std::cerr<<"Usage: 3D image cropping and padding"<<std::endl;
		std::cerr<<argv[0]<<" inputFile outputFile startX startY startZ  sizeX sizeY sizeZ paddingSizeX paddingSizeY paddingSizeZ"<<std::endl;
		return -1;
	}


	typedef itk::Image<float, 3>        ImageType;

	typedef itk::ImageFileReader< ImageType> ReaderType;


	typedef itk::ImageRegionConstIterator< ImageType > ConstIteratorType;
	typedef itk::ImageRegionIterator< ImageType> IteratorType;

	//read image
	std::string input  = argv[1];
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName( input.c_str() );
	reader->ReleaseDataFlagOn();
	try         {
		reader->Update();
	}     catch (itk::ExceptionObject &err)         {
		std::cout<< "Exception caught:" << std::endl;
		std::cout<< err <<std::endl;
		return -1;
	}



	//   find out the boundary: go throught the volumes, and find out the upper left corner and lower right corner

	ImageType::IndexType	inputStart;
	ImageType::SizeType   	inputSize;
	ImageType::IndexType     padWidth;

	inputStart[0]=atoi(argv[3]);//xStart;
	inputStart[1]=atoi(argv[4]);//yStart;
	inputStart[2]=atoi(argv[5]);//zStart;
	inputSize[0]=atoi(argv[6]);// sizeX;
	inputSize[1]=atoi(argv[7]);//sizeY;	
	inputSize[2]=atoi(argv[8]);// sizeZ;

	padWidth[0] =atoi(argv[9]); //padde width in x direciton
	padWidth[1] =atoi(argv[10]);
	padWidth[2] =atoi(argv[11]);





	//allocate final output image 
	ImageType::Pointer outputImage = ImageType::New();

	const ImageType::SpacingType& spacing = reader->GetOutput()->GetSpacing();
	const ImageType::PointType& inputOrigin = reader->GetOutput()->GetOrigin();

	ImageType::SizeType   	outputSize;
	ImageType::IndexType    outputStart;

	double outputOrigin[3 ];
	for(unsigned int i=0; i< 3; i++)
	{
		//pad image 
	//	outputStart[i]= inputStart[i]- padWidth[i];
		outputStart[i]=0;
		outputSize[i]= inputSize[i] + 2*padWidth[i];

		outputOrigin[i] = inputOrigin[i] + spacing[i] *(inputStart[i]- padWidth[i]);
	}
	outputImage->SetSpacing( spacing );
	outputImage->SetOrigin( outputOrigin );

	ImageType::RegionType outputRegion;
	outputRegion.SetSize(outputSize);
	outputRegion.SetIndex(outputStart);
       
	outputImage->SetRegions( outputRegion );
	outputImage->Allocate();

	//fill in with zeros
	IteratorType outputImIt( outputImage, outputRegion );
	for (  outputImIt.GoToBegin(); !outputImIt.IsAtEnd(); ++outputImIt)
	{
		outputImIt.Set(0.0 );
	}


	std::cerr<<"allocated output image: "<<outputSize[0]<<" "<<outputSize[1]<<" "<<outputSize[2]<<std::endl;
	std::cerr<<"new origin:"<<outputOrigin[0]<<" "<<outputOrigin[1]<<" "<<outputOrigin[2]<<std::endl;

	//copy
	ImageType::RegionType	CroppedRegion; 
	CroppedRegion.SetSize(inputSize);
	CroppedRegion.SetIndex(inputStart);
	ConstIteratorType inputIt( reader->GetOutput(), CroppedRegion );


	ImageType::IndexType subRegionIndex;
	for(unsigned int i=0; i< 3; i++)
	{
		subRegionIndex[i]= padWidth[i];

	}
	ImageType::RegionType outputSubRegion;
	outputSubRegion.SetSize(inputSize);
	outputSubRegion.SetIndex(subRegionIndex);

	IteratorType outputIt( outputImage, outputSubRegion );
	for ( inputIt.GoToBegin(), outputIt.GoToBegin(); !inputIt.IsAtEnd();
			++inputIt, ++outputIt)
	{
		outputIt.Set( inputIt.Get() );
	}



	std::string output = argv[2];
	typedef itk::ImageFileWriter< ImageType> WriterType;
	WriterType::Pointer writer = WriterType::New();

	writer->SetFileName( output.c_str() );
	writer->SetInput(outputImage);

	try {
		writer->Update();
	}     catch (itk::ExceptionObject &err)         {
		std::cout<< "Exception caught:" << std::endl;
		std::cout<< err <<std::endl;
		return -1;
	}

	return 0;
}


