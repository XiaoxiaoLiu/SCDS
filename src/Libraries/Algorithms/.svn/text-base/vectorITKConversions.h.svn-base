
#ifndef _vectorITKConversions_h
#define _vectorITKConversions_h

#include "Image.h"

#include "itkVectorImage.h"

#include "vnl/vnl_vector.h"

// Convert Image to itk::VectorImage
itk::VectorImage<float>::Pointer
convertImageVolume(Image< vnl_vector<float> >* img);

// Convert itk::VectorImage to Image
Image< vnl_vector<float> >*
convertITKVolume(itk::VectorImage<float>* itkimg);

#endif
