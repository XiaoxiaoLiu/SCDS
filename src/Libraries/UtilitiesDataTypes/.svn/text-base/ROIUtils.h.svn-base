//////////////////////////////////////////////////////////////////////
//
// File: ROIUtils.h
//
// Supports basic operations against the ROI class
//
//////////////////////////////////////////////////////////////////////
#ifndef _ROI_UTILS_H_
#define _ROI_UTILS_H_

#include <Array3D.h>
#include <ROI.h>

template <class T>
class ROIUtils
{
public:
  // Creates a new 3D sub-volume based on ROI 
  static void extractROIFromArray3D(const Array3D<T>& inArray,
    Array3D<T>& outArray,
    const ROI<int, unsigned int>& roi) {
    outArray.resize(roi.getSize());
    
    typename Array3D<T>::IndexType start = roi.getStart();
    typename Array3D<T>::SizeType size = roi.getSize();
    
    for(unsigned int zIndex = 0; zIndex < size.z; zIndex++) {
      for(unsigned int yIndex = 0; yIndex < size.y; yIndex++) {
        for(unsigned int xIndex = 0; xIndex < size.x; xIndex++) {
          outArray.set(xIndex, yIndex, zIndex, inArray(start.x + xIndex,
            start.y + yIndex,
            start.z + zIndex));
        }
      }
    }
  }
}; // class ROIUtils

#endif // #ifndef _ROI_UTILS_H_
