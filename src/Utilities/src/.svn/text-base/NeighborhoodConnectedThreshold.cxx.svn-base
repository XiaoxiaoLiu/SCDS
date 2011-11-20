/*=========================================================================

Program:   Insight Segmentation & Registration Toolkit
Module:    $RCSfile:RegionGrowingSegmentation.cxx,v $
Language:  C++      
Date:      $Date: 2008/11/15 15:39:12 $
Version:   $Revision: 1.18 $

Copyright (c) 2002 Insight Consortium. All rights reserved.
See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

This software is distributed WITHOUT ANY WARRANTY; without even 
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#include <iostream>


#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkNeighborhoodConnectedImageFilter.h"


int main( int argc, char *argv[] )
{
	if (argc < 7 )
	{
		std::cerr << "Missing Parameters " << std::endl;
		std::cerr << "Usage: " << argv[0];
		std::cerr << " inputImageFileName outputImageFileName lowThreshold_regionGrow highThreshold_regionGrow seed_x seed_y seed_z   " << std::endl;
		return 1;
	}

	//---------------------------------------------  DEFINE IMAGE TYPE -----------------------------------------------
	const unsigned int Dimension = 3;

	typedef itk::Image<float, Dimension > ImageType;
	typedef itk::Image<short,Dimension> OutImageType;
	

	//---------------------------------------------      READ  IMAGE ---------------------------------------------------

	typedef itk::ImageFileReader<ImageType> ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(argv[1]);

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

	//---------------------------------------------      PIPELINE --------------------------------------------------
	
	//region grow segmentation
	typedef itk::NeighborhoodConnectedImageFilter< ImageType,OutImageType > ConnectedFilterType;
	ConnectedFilterType::Pointer connectedThreshold = ConnectedFilterType::New();

	ImageType::PixelType lowerThreshold = atoi(argv[3]);
	ImageType::PixelType upperThreshold = atoi(argv[4]); 


	itk::Index<3> seedIndex = {{ atoi(argv[5]),  atoi(argv[6]),  atoi(argv[7])}};
	connectedThreshold->SetInput( reader->GetOutput() );
	connectedThreshold->SetLower( lowerThreshold );
	connectedThreshold->SetUpper( upperThreshold );
	connectedThreshold->SetReplaceValue( 1 );	
	connectedThreshold->SetSeed( seedIndex );

    ImageType::SizeType radius;
    radius[0]= 1; radius[1] = 1; radius[2]=1;
 	connectedThreshold->SetRadius(radius);

	std::cerr <<"seedPixel value ="<< reader->GetOutput() ->GetPixel(seedIndex) << std::endl;


	typedef itk::ImageFileWriter<OutImageType> regionGrowFileWriterType;
	regionGrowFileWriterType:: Pointer regionGrowWriter = regionGrowFileWriterType::New();
	regionGrowWriter->SetFileName(argv[2]);

	regionGrowWriter->SetInput(connectedThreshold->GetOutput());


	try 
	{
		regionGrowWriter->Update(); 
	}
	catch (itk::ExceptionObject &e)
	{
		std::cerr << e << std::endl;
	}



	return 0;

}






