#ifndef ApplicationUtils_txx
#define ApplicationUtils_txx

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include "Image.h"
#include <exception>
#include <itkImportImageFilter.h>

template <class VoxelType>
void 
ApplicationUtils::
LoadImageITK(const char *fileName, Image<VoxelType>& image)
{
  //
  // load itk image
  //
  typedef itk::Image<VoxelType, 3>        ImageType;
  typedef typename ImageType::Pointer     ImagePointer;
  typedef itk::ImageFileReader<ImageType> VolumeReaderType;
  typename VolumeReaderType::Pointer reader = VolumeReaderType::New();
  reader->SetFileName(fileName);
  ImagePointer imagePtr = reader->GetOutput();

  try
  {
    reader->Update();
  }
  catch(itk::ImageFileReaderException exc)
  {
    std::stringstream ss;
    ss << "ITK FileReaderException: " << exc;
    throw std::runtime_error(ss.str().c_str());
  }
  catch(itk::ExceptionObject exc)
  {
    std::stringstream ss;
    ss << "ITK Exception: " << exc;
    throw std::runtime_error(ss.str().c_str());
  }
  catch(...)
  {
    std::stringstream ss;
    ss << "Unknown Exception";
    throw std::runtime_error(ss.str().c_str());
  }

  //
  // copy data to Image
  //
  // instead we should just steal the pointer
  int sizeX = imagePtr->GetLargestPossibleRegion().GetSize()[0];
  int sizeY = imagePtr->GetLargestPossibleRegion().GetSize()[1];
  int sizeZ = imagePtr->GetLargestPossibleRegion().GetSize()[2];
  image.resize(sizeX, sizeY, sizeZ);
  memcpy(image.getDataPointer(),
         imagePtr->GetBufferPointer(),
         sizeX * sizeY * sizeZ * sizeof(VoxelType));

  //
  // copy origin and spacing
  //
  image.setOrigin(imagePtr->GetOrigin()[0],
                  imagePtr->GetOrigin()[1],
                  imagePtr->GetOrigin()[2]);

  image.setSpacing(imagePtr->GetSpacing()[0],
                   imagePtr->GetSpacing()[1],
                   imagePtr->GetSpacing()[2]);
}

template <class VoxelType>
void 
ApplicationUtils::
SaveImageITK(const char *fileName, const Image<VoxelType>& image)
{
  //
  // create an itk version of the image 
  //
  typedef itk::Image<VoxelType, 3>               ImageType;
  typedef itk::ImportImageFilter<VoxelType, 3>   ImportFilterType;
  typename ImportFilterType::Pointer importFilter = ImportFilterType::New();

  // construction image region (basically just the dimensions)
  typename ImportFilterType::IndexType start;
  start.Fill(0);
  typename ImportFilterType::SizeType size;
  size[0] = image.getSize()[0];
  size[1] = image.getSize()[1];
  size[2] = image.getSize()[2];

  typename ImportFilterType::RegionType region;
  region.SetSize(size);
  region.SetIndex(start);
  importFilter->SetRegion(region);

  // set origin and spacing 
  typename ImportFilterType::OriginType origin;
  origin[0] = image.getOrigin()[0];
  origin[1] = image.getOrigin()[1];
  origin[2] = image.getOrigin()[2];
  importFilter->SetOrigin(origin);

  typename ImportFilterType::SpacingType spacing;
  spacing[0] = image.getSpacing()[0];
  spacing[1] = image.getSpacing()[1];
  spacing[2] = image.getSpacing()[2];
  importFilter->SetSpacing(spacing);

  // copy the data into the itk image
  const bool importImageFilterWillOwnTheBuffer = false;
  importFilter->SetImportPointer((VoxelType*)image.getDataPointer(),
                                 image.getNumElements(),
                                 importImageFilterWillOwnTheBuffer);

  //
  // create writer
  //
  typedef itk::ImageFileWriter<ImageType> VolumeWriterType;
  typename VolumeWriterType::Pointer writer = VolumeWriterType::New();
  writer->SetFileName(fileName);
  writer->SetInput(importFilter->GetOutput());

  //
  // write image
  //
  try
  {
    writer->Update();
  }
  catch(itk::ImageFileWriterException exc)
  {
    std::stringstream ss;
    ss << "ITK FileWriterException: " << exc;
    throw std::runtime_error(ss.str().c_str());
  }
  catch(itk::ExceptionObject exc)
  {
    std::stringstream ss;
    ss << "ITK Exception: " << exc;
    throw std::runtime_error(ss.str().c_str());
  }
  catch(...)
  {
    std::stringstream ss;
    ss << "Unknown Exception";
    throw std::runtime_error(ss.str().c_str());
  }
}

#endif // ApplicationUtils_txx
