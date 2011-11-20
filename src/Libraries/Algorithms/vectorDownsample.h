#ifndef _vectorDownsample_h
#define _vectorDownsample_h

#include "Image.h"

#include "vnl/vnl_vector.h"

Image< vnl_vector<float> >*
vectorDownsample(
  Image< vnl_vector<float> >* image,
  double factor,
  double sigma,
  double kernelSize
);

#endif
