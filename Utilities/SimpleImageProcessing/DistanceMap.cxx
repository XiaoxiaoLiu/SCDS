#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkSignedDanielssonDistanceMapImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"

#include <fstream>
#include <iostream>
#include <vector>
using namespace std ;

typedef itk::Image < float, 3 > ImageType ;
typedef ImageType::Pointer ImagePointer ;
typedef itk::ImageFileWriter < ImageType > ImageWriterType ;
typedef itk::ImageFileReader < ImageType > ImageReaderType ;
typedef itk::SignedDanielssonDistanceMapImageFilter< ImageType, ImageType > DDMapType ;
typedef itk::DiscreteGaussianImageFilter < ImageType, ImageType > gaussFilterType ;

int main(int argc, const char **argv)
{
  std::cout << "Usage: " << argv[0] << " binaryInput output [smoothing]" << std::endl ;
  double smoothing = 0 ;
  if ( argc == 4 ) 
    smoothing = atof ( argv[3] ) ;
  std::cout << "smoothing: " << smoothing << std::endl ;

  // read the input volume
  ImageReaderType::Pointer reader = ImageReaderType::New () ;
  reader->SetFileName ( argv[1] ) ;
  reader->Update () ;
  std::cout << "image read" << std::endl ;

  ImagePointer fImage = reader->GetOutput () ;

  // now set up the Danielsson filter
  DDMapType::Pointer ddmap = DDMapType::New () ;
  ddmap->SetInput ( fImage ) ;
  ddmap->Update () ;
  std::cout << "distance map computed" << std::endl ;
  
  // gaussian smoothing
  gaussFilterType::Pointer smoothFilter ;
  if ( smoothing > 0 )
    {
      smoothFilter = gaussFilterType::New () ;
      smoothFilter->SetInput ( ddmap->GetDistanceMap () ) ;
      smoothFilter->SetVariance ( smoothing ) ;
      smoothFilter->Update () ;
      std::cout << "smoothing completed" << std::endl ;
    }

  // write out the result
  ImageWriterType::Pointer itkWriter = ImageWriterType::New() ;
  itkWriter->SetFileName ( argv[2] ) ; 
  if ( smoothing > 0 ) 
    itkWriter->SetInput ( smoothFilter->GetOutput() ) ;
  else
    itkWriter->SetInput ( ddmap->GetDistanceMap () ) ;
  itkWriter->Write() ;          

  return 0 ;
}
