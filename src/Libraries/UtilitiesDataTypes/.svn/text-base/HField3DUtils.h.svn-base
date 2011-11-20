#ifndef HFIELD3D_UTILS_H
#define HFIELD3D_UTILS_H

#include <cassert>
#include <limits>
#include <cmath>
#include <queue>

#include "Vector3D.h"
#include "Array3D.h"
#include "Array3DUtils.h"
#include "Image.h"
#include "ROI.h"
#include "Surface.h"
#include "Matrix3D.h"
#include "Array3DIO.h"

#ifdef max
#undef max
#endif

class HField3DUtils
{
public:
  enum BackgroundStrategy { BACKBACKGROUND_STRATEGY_PARTIAL_ID,
                            BACKGROUND_STRATEGY_ID,
                            BACKGROUND_STRATEGY_ZERO };

  //
  // set hfield to identity
  // i.e. h(x) = x
  //
  // davisb 2003
  //
  template <class T>
  static
  void 
  setToIdentity(Array3D<Vector3D<T> >& hField)
  {
    Vector3D<unsigned int> size = hField.getSize();
    for (unsigned int z = 0; z < size.z; ++z) {
      for (unsigned int y = 0; y < size.y; ++y) {
	for (unsigned int x = 0; x < size.x; ++x) {
	  hField(x, y, z).set(x, y, z);
	}
      }
    }
  }

  //
  // convert a velocity field to an h field
  // h(x) = x + v(x) * delta
  //
  // bcd 2004
  template <class T>
  static
  void 
  velocityToH(Array3D<Vector3D<T> >& hField, const T& delta)
  {
    Vector3D<unsigned int> size = hField.getSize();
    for (unsigned int z = 0; z < size.z; ++z) {
      for (unsigned int y = 0; y < size.y; ++y) {
	for (unsigned int x = 0; x < size.x; ++x) {
	  hField(x, y, z).x = x + hField(x,y,z).x * delta;
	  hField(x, y, z).y = y + hField(x,y,z).y * delta;
	  hField(x, y, z).z = z + hField(x,y,z).z * delta;
	}
      }
    }
  }

  //
  // compose two h fields using trilinear interpolation
  // h(x) = f(g(x))
  //
  // davisb 2003
  //
  template <class T>
  static 
  void
  compose(const Array3D<Vector3D<T> >& f,
          const Array3D<Vector3D<T> >& g,
          Array3D<Vector3D<T> >& h)
  {
    Vector3D<unsigned int> size = h.getSize();
    for (unsigned int z = 0; z < size.z; ++z) {
      for (unsigned int y = 0; y < size.y; ++y) {
	for (unsigned int x = 0; x < size.x; ++x) {
 	  trilerp(f, 
 		  g(x,y,z).x, g(x,y,z).y, g(x,y,z).z, 
 		  h(x,y,z).x, h(x,y,z).y, h(x,y,z).z);
	}
      }
    }    
  }

  //
  // compose a velocity and h field to get an hfield
  // h(x) = x + v(g(x))
  //
  // davisb 2007
  //
  template <class T>
  static 
  void
  composeVH(const Array3D<Vector3D<T> >& v,
            const Array3D<Vector3D<T> >& g,
            Array3D<Vector3D<T> >& h)
  {
    Vector3D<unsigned int> size = h.getSize();
    for (unsigned int z = 0; z < size.z; ++z) {
      for (unsigned int y = 0; y < size.y; ++y) {
	for (unsigned int x = 0; x < size.x; ++x) {
 	  trilerp(v, 
 		  g(x,y,z).x, g(x,y,z).y, g(x,y,z).z, 
 		  h(x,y,z).x, h(x,y,z).y, h(x,y,z).z);
          h(x,y,z).x += x;
          h(x,y,z).y += y;
          h(x,y,z).z += z;          
	}
      }
    }    
  }

  //
  // compose a h field and a velocify field to get an hfield
  // h(x) = g(x+v(x))
  //
  // davisb 2007
  //
  template <class T>
  static 
  void
  composeHV(const Array3D<Vector3D<T> >& g,
            const Array3D<Vector3D<T> >& v,
            Array3D<Vector3D<T> >& h)
  {
    Vector3D<unsigned int> size = h.getSize();
    for (unsigned int z = 0; z < size.z; ++z) {
      for (unsigned int y = 0; y < size.y; ++y) {
	for (unsigned int x = 0; x < size.x; ++x) {
 	  trilerp(g, 
 		  x+v(x,y,z).x, y+v(x,y,z).y, z+v(x,y,z).z, 
 		  h(x,y,z).x, h(x,y,z).y, h(x,y,z).z);
	}
      }
    }    
  }

  //
  // compose h field with translation
  // creating h(x) = f(x) + t
  //
  // davisb 2003
  //
  template <class T>
  static
  void
  composeTranslation(const Array3D<Vector3D<T> >& f,
		     const Vector3D<T>& t,
		     Array3D<Vector3D<T> >& h)
  {
    unsigned int numElements = f.getNumElements();
    for (unsigned int i = 0; i < numElements; ++i)
      {
        h(i) = f(i) + t;
      }
  }

  //
  // compose h field with translation
  // creating h(x) = f(x) + t
  //
  // davisb 2003
  //
  template <class T>
  static
  void
  composeTranslation(const Array3D<Vector3D<T> >& f,
		     const Vector3D<T>& t,
		     Array3D<Vector3D<T> >& h,
                     const ROI<int,unsigned int>& roi)
  {
    for (int z = roi.getStartZ(); z <= roi.getStopZ(); ++z) {
      for (int y = roi.getStartY(); y <= roi.getStopY(); ++y) {
        for (int x = roi.getStartX(); x <= roi.getStopX(); ++x) {
          h(x,y,z) = f(x,y,z) + t;
        }
      }
    }
  }

  //
  // precompose h field with translation
  // creating h(x) = f(x + t)
  //
  // davisb 2003
  //
  template <class T>
  static
  void
  preComposeTranslation(const Array3D<Vector3D<T> >& f,
			const Vector3D<T>& t,
			Array3D<Vector3D<T> >& h)
  {
    Vector3D<unsigned int> size = f.getSize();
    for (unsigned int z = 0; z < size.z; ++z) {
      for (unsigned int y = 0; y < size.y; ++y) {
	for (unsigned int x = 0; x < size.x; ++x) {
	  trilerp(f, 
		  x + t.x, y + t.y, z + t.z, 
		  h(x, y, z).x, h(x, y, z).y, h(x, y, z).z);
	}
      }
    }
  }

  //
  // approximate the inverse of an incremental h field using according
  // to the following derivation
  // 
  // hInv(x0) = x0 + d
  // x0 = h(x0 + d)  
  // x0 = h(x0) + d // order zero expansion
  // d  = x0 - h(x0)
  //
  // hInv(x0) = x0 + x0 - h(x0)
  //
  //
  // bcd 2004
  //
  template <class T>
  static
  void
  computeInverseZerothOrder(const Array3D<Vector3D<T> >& h,
			    Array3D<Vector3D<T> >& hinv)
  {
    assert(h.getSize() == hinv.getSize());

    Vector3D<unsigned int> size = h.getSize();
    for (unsigned int z = 0; z < size.z; ++z) {    
      for (unsigned int y = 0; y < size.y; ++y) {    
	for (unsigned int x = 0; x < size.x; ++x) {    
	  hinv(x,y,z).x = x + x - h(x,y,z).x;
	  hinv(x,y,z).y = y + y - h(x,y,z).y;
	  hinv(x,y,z).z = z + z - h(x,y,z).z;
	}
      }
    }
  }

  template <class T>
  static
  void
  computeInverseConsistencyError(const Array3D<Vector3D<T> >& h,
				 const Array3D<Vector3D<T> >& hinv,
				 double& minError, double& maxError,
				 double& meanError, double& stdDevError)
  {
    Vector3D<unsigned int> size = hinv.getSize();
    assert(h.getSize() == hinv.getSize());
    assert(h.getNumElements() > 0);

    // compute hinvh(x) = hinv(h(x))
    Array3D<Vector3D<T> > hinvh(size);
    compose(hinv, h, hinvh);

    // compute statistics on hinvh
    minError  = std::numeric_limits<double>::max();
    maxError  = 0;
    meanError = 0;

    //
    // compute min, max, and mean
    //
    //unsigned int numElements = hinvh.getNumElements();
    unsigned int z;
    for (z = 30; z < size.z-30; ++z) {    
      for (unsigned int y = 30; y < size.y-30; ++y) {    
	for (unsigned int x = 30; x < size.x-30; ++x) {    
	  double error = 
	    sqrt((hinvh(x,y,z).x - x) * (hinvh(x,y,z).x - x) +
		 (hinvh(x,y,z).y - y) * (hinvh(x,y,z).y - y) +
		 (hinvh(x,y,z).z - z) * (hinvh(x,y,z).z - z));
	  if (error > maxError) maxError = error;
	  if (error < minError) minError = error;
	  meanError += error;
	}
      }
    }
    meanError /= (size.x-60) * (size.y-60) * (size.z-60);

    //
    // now compute standard deviation
    // 
    stdDevError = 0;
    for (z = 30; z < size.z-30; ++z) {    
      for (unsigned int y = 30; y < size.y-30; ++y) {    
	for (unsigned int x = 30; x < size.x-30; ++x) {    
	  double error = 
	    sqrt((hinvh(x,y,z).x - x) * (hinvh(x,y,z).x - x) +
		 (hinvh(x,y,z).y - y) * (hinvh(x,y,z).y - y) +
		 (hinvh(x,y,z).z - z) * (hinvh(x,y,z).z - z));
	  stdDevError += (error - meanError) * (error - meanError);
	}
      }
    }
    stdDevError /= ((size.x-60) * (size.y-60) * (size.z-60) - 1);
    stdDevError = sqrt(stdDevError);
  }

  template <class T>
  static
  void
  computeInverseConsistencyError(const Array3D<Vector3D<T> >& h,
				 const Array3D<Vector3D<T> >& hinv,
				 double& hhinvMinError, 
				 double& hhinvMaxError,
				 double& hhinvMeanError, 
				 double& hhinvStdDevError,
				 double& hinvhMinError, 
				 double& hinvhMaxError,
				 double& hinvhMeanError, 
				 double& hinvhStdDevError)
  {
    computeInverseConsistencyError(hinv, h, hhinvMinError, hhinvMaxError,
				   hhinvMeanError, hhinvStdDevError);
    computeInverseConsistencyError(h, hinv, hinvhMinError, hinvhMaxError,
				   hinvhMeanError, hinvhStdDevError);
  }

  template <class T>
  static
  void
  reportInverseConsistencyError(const Array3D<Vector3D<T> >& h,
				const Array3D<Vector3D<T> >& hinv)
  {
    double 
      hhinvMinError    = 0, 
      hhinvMaxError    = 0, 
      hhinvMeanError   = 0, 
      hhinvStdDevError = 0;
    double 
      hinvhMinError    = 0, 
      hinvhMaxError    = 0, 
      hinvhMeanError   = 0, 
      hinvhStdDevError = 0;

    computeInverseConsistencyError(h,
				   hinv,
				   hhinvMinError, 
				   hhinvMaxError, 
				   hhinvMeanError, 
				   hhinvStdDevError,
				   hinvhMinError, 
				   hinvhMaxError,
				   hinvhMeanError,
				   hinvhStdDevError);

    std::cerr << "\t|h(hinv(x)) - x|L2" << std::endl
	      << "\tMin      : " << hhinvMinError << std::endl
	      << "\tMax      : " << hhinvMaxError << std::endl	
	      << "\tMean     : " << hhinvMeanError << std::endl	
	      << "\tSt. Dev. : " << hhinvStdDevError << std::endl	
	      << "\t|hinv(h(x)) - x|L2" << std::endl
	      << "\tMin      : " << hinvhMinError << std::endl
	      << "\tMax      : " << hinvhMaxError << std::endl	
	      << "\tMean     : " << hinvhMeanError << std::endl	
	      << "\tSt. Dev. : " << hinvhStdDevError << std::endl;
  }

  //
  // apply hField to an image 
  // defImage(x) = image(h(x))
  //
  // trilerp by default but will use nearest neighbor if flag is set
  // to true
  //
  // NOTE: this does not round for integer types
  //
  // davisb 2003
  //
  template <class T, class U>
  static
  void
  apply(const Array3D<T>& image,
	const Array3D<Vector3D<U> >& hField,
	Array3D<T>& defImage,
	const T& background = 0,
	bool useNearestNeighbor = false)
  {
    Vector3D<unsigned int> size = defImage.getSize();
    
    if (useNearestNeighbor)
      {
	for (unsigned int z = 0; z < size.z; ++z) {
	  for (unsigned int y = 0; y < size.y; ++y) {
	    for (unsigned int x = 0; x < size.x; ++x) {
	      defImage(x, y, z) = 
		Array3DUtils::nearestNeighbor(image,
					      hField(x, y, z),
					      background);
	    }
	  }
	}
      }
    else
      {
	for (unsigned int z = 0; z < size.z; ++z) {
	  for (unsigned int y = 0; y < size.y; ++y) {
	    for (unsigned int x = 0; x < size.x; ++x) {
	      defImage(x, y, z) = 
		static_cast<T>(Array3DUtils::trilerp(image,
						     hField(x, y, z),
						     background));
	    }
	  }
	}
      }
  }

  //
  // apply hField to the roi of an image 
  // defImage(x) = image(h(x))
  //
  // trilerp by default but will use nearest neighbor if flag is set
  // to true
  //
  // NOTE: this does not round for integer types
  //
  // davisb 2004
  //
  template <class T, class U>
  static
  void
  apply(const Array3D<T>& image,
	const Array3D<Vector3D<U> >& hField,
	Array3D<T>& defImage,
	int hFieldStartX,
	int hFieldStartY,
	int hFieldStartZ,
	const T& background = 0,
	bool useNearestNeighbor = false)
  {
    Vector3D<unsigned int> imageSize = image.getSize();
    Vector3D<unsigned int> hFieldSize = hField.getSize();

    // copy entire image
    defImage = image;

    // ammend in region of interest
    if (useNearestNeighbor)
      {
	for (unsigned int z = 0; 
	     z < hFieldSize.z && z + hFieldStartZ < imageSize.z; ++z) {
	  for (unsigned int y = 0; 
	       y < hFieldSize.y && y + hFieldStartY < imageSize.y; ++y) {
	    for (unsigned int x = 0; 
		 x < hFieldSize.x && x + hFieldStartX < imageSize.x; ++x) {
	      defImage(x+hFieldStartX, y+hFieldStartY, z+hFieldStartZ) = 
		Array3DUtils::
                nearestNeighbor(image,
                                hField(x, y, z).x + hFieldStartX,
                                hField(x, y, z).y + hFieldStartY,
                                hField(x, y, z).z + hFieldStartZ,
                                background);
	    }
	  }
	}
      }
    else
      {
	for (unsigned int z = 0; 
	     z < hFieldSize.z && z + hFieldStartZ < imageSize.z; ++z) {
	  for (unsigned int y = 0; 
	       y < hFieldSize.y && y + hFieldStartY < imageSize.y; ++y) {
	    for (unsigned int x = 0; 
		 x < hFieldSize.x && x + hFieldStartX < imageSize.x; ++x) {
	      defImage(x+hFieldStartX, y+hFieldStartY, z+hFieldStartZ) = 
		static_cast<T>(Array3DUtils::
                               trilerp(image,
                                       hField(x, y, z).x + hFieldStartX,
                                       hField(x, y, z).y + hFieldStartY,
                                       hField(x, y, z).z + hFieldStartZ,
                                       background));
	    }
	  }
	}
      }
  }

  //
  // Apply hField to an roi of an image, taking origin and spacing
  // into account.  This separately considers the origin and spacing
  // of the target image (specified in defImage), the moving image
  // ('image'), and the hField.  The roi is determined by hFieldOrigin
  // and hFieldSpacing, rather than an ROI object.  Use this if the
  // origin and spacing of the image being deformed are different from
  // those of the images used to generate the hfield.  The entries in
  // hField are assumed to be given in voxel coordinates, so they must
  // be scaled and shifted by hFieldSpacing and hFieldOrigin
  // respectively to get them in physical coordinates.
  //
  // defImage(x) = image(h(x))
  //
  // trilerp only
  //
  // NOTE: this does not round for T and U if they're integer types
  //
  // foskey 2004
  //
  template <class T, class U>
  static
  void 
  applyOldVersion(const Image<T>& image,
	const Array3D< Vector3D<U> >& hField,
	Image<T>& defImage,
        Vector3D<double> hFieldOrigin,
        Vector3D<double> hFieldSpacing,
	const T& background = 0)
  {
    Vector3D<unsigned int> defImageSize = defImage.getSize();
    Vector3D<unsigned int> hFieldSize = hField.getSize();

    Vector3D<double> imageScaleHtoI = hFieldSpacing / image.getSpacing();
    Vector3D<double> defImageScaleHtoI = hFieldSpacing / defImage.getSpacing();

    Vector3D<double> imageRoiStart = (hFieldOrigin - image.getOrigin()) /
                                     image.getSpacing();

    Vector3D<double> defImageRoiStart = (hFieldOrigin - defImage.getOrigin()) /
                                        defImage.getSpacing();

    Vector3D<double> defImageRoiSize = Vector3D<double>(hFieldSize) *
                                       defImageScaleHtoI;

    for (unsigned int i = 0; i < 3; ++i) {
      if (defImageRoiStart[i] < 0) defImageRoiStart[i] = 0;
      if (defImageRoiSize[i] + defImageRoiStart[i] + 0.5 > defImageSize[i]) {
        defImageRoiStart[i] = defImageSize[i] - defImageRoiStart[0] - 0.6;
      }
    }

    unsigned int hFieldStartX = (unsigned int)(defImageRoiStart[0] + 0.5);
    unsigned int hFieldStartY = (unsigned int)(defImageRoiStart[1] + 0.5);
    unsigned int hFieldStartZ = (unsigned int)(defImageRoiStart[2] + 0.5);

    unsigned int hFieldSizeX = (unsigned int)(defImageRoiSize[0] + 0.5);
    unsigned int hFieldSizeY = (unsigned int)(defImageRoiSize[1] + 0.5);
    unsigned int hFieldSizeZ = (unsigned int)(defImageRoiSize[2] + 0.5);

    for (unsigned int z = 0; z < hFieldSizeZ; ++z) {
      for (unsigned int y = 0; y < hFieldSizeY; ++y) {
        for (unsigned int x = 0; x < hFieldSizeX; ++x) {

          Vector3D<float> imagePoint = hField(x, y, z) * imageScaleHtoI + 
                                       imageRoiStart;

          defImage(x+hFieldStartX, y+hFieldStartY, z+hFieldStartZ) = 
            static_cast<T>(Array3DUtils::trilerp(image, imagePoint,
                                                 background));

        }
      }
    }
  }

  //
  // Apply hField to an roi of an image, taking origin and spacing
  // into account.  The roi is determined by hFieldOrigin and
  // hFieldSpacing, rather than an ROI object.  The size of the
  // resulting image is the size of the hField, which equals the size
  // of the ROI.  Use this if the origin and spacing of the image
  // being deformed are different from those of the images used to
  // generate the hfield.  The entries in hField are assumed to be
  // given in voxel coordinates, so they must be scaled and shifted by
  // hFieldSpacing and hFieldOrigin respectively to get them in
  // physical coordinates.
  //
  // defImage(x) = image(h(x))
  //
  // trilerp only
  //
  // NOTE: this does not round for T and U if they're integer types
  //
  // foskey 2004
  //
  template <class T, class U>
  static
  void 
  apply(const Image<T>& image,
	const Array3D< Vector3D<U> >& hField,
	Image<T>& defImage,
        Vector3D<double> hFieldOrigin,
        Vector3D<double> hFieldSpacing,
	const T& background = 0)
  {
    Vector3D<unsigned int> hFieldSize = hField.getSize();

    Vector3D<double> scaleHtoI = hFieldSpacing / image.getSpacing();

    Vector3D<double> roiStart = (hFieldOrigin - image.getOrigin()) /
                                     image.getSpacing();

    defImage.resize(hFieldSize);
    defImage.setOrigin(hFieldOrigin);
    defImage.setSpacing(hFieldSpacing);

    for (unsigned int z = 0; z < hFieldSize.z; ++z) {
      for (unsigned int y = 0; y < hFieldSize.y; ++y) {
        for (unsigned int x = 0; x < hFieldSize.x; ++x) {

          Vector3D<float> imagePoint = hField(x, y, z) * scaleHtoI + 
                                       roiStart;

          defImage(x, y, z) = static_cast<T>(
            Array3DUtils::trilerp(image, imagePoint, background));

        }
      }
    }
  }

  //
  // apply hField to an image by using a mask 
  // defImage(x) = image(h(x))
  //
  // The mask allows you to avoid applying the algorithm on certain
  // parts of the image, like for example the bone.  When the value of
  // the mask is false, the displacement is null.
  //
  // trilerp by default but will use nearest neighbor if flag is set
  // to true
  //
  // NOTE: this does not round for integer types
  //
  // dprigent 2004
  //
  template <class T, class U>
    static
    void
    applyWithMask(const Array3D<T>& image,
    Array3D<Vector3D<U> > hField,
    Array3D<T>& defImage,
    Array3D<bool>& mask,
    const T& background = 0,
    bool useNearestNeighbor = false)
  {
    Vector3D<unsigned int> size = image.getSize();
    
    if (useNearestNeighbor)
    {
      for (unsigned int z = 0; z < size.z; ++z) {
        for (unsigned int y = 0; y < size.y; ++y) {
          for (unsigned int x = 0; x < size.x; ++x) {
            if(mask(x,y,z))
            {
              defImage(x, y, z) = 
                Array3DUtils::nearestNeighbor(image,
                hField(x, y, z),
                background);
            }
            else
            {
              hField(x, y, z).x=x;
              hField(x, y, z).y=y;
              hField(x, y, z).z=z;
              
            }
          }
        }
      }
    }
  else
  {
    for (unsigned int z = 0; z < size.z; ++z) {
      for (unsigned int y = 0; y < size.y; ++y) {
        for (unsigned int x = 0; x < size.x; ++x) {
          if(mask(x,y,z))
          {
            defImage(x, y, z) = 
              static_cast<T>(Array3DUtils::trilerp(image,
              hField(x, y, z),
              background));
          }
          else
          {
            hField(x, y, z).x=x;
            hField(x, y, z).y=y;
            hField(x, y, z).z=z;
          }
        }
      }
    }
  }
  }

  //
  // note: you should not call this function with &image == &defImage
  //
  template <class T, class U>
  static
  void
  forwardApply(const Array3D<T>& image,
               const Array3D<Vector3D<U> >& hField,
               Array3D<T>& defImage,
               int hFieldStartX,
               int hFieldStartY,
               int hFieldStartZ,
               const T& background = 0)
  {
    Vector3D<unsigned int> imageSize = image.getSize();
    Vector3D<unsigned int> hFieldSize = hField.getSize();
    Array3D<float> count(defImage.getSize());
    count.fill(0);

    // copy entire image
    defImage = image;

    // fill roi with background
    for (unsigned int z = 0; 
         z < hFieldSize.z && z + hFieldStartZ < imageSize.z; ++z) {
      for (unsigned int y = 0; 
           y < hFieldSize.y && y + hFieldStartY < imageSize.y; ++y) {
        for (unsigned int x = 0; 
             x < hFieldSize.x && x + hFieldStartX < imageSize.x; ++x) {
          defImage(x+hFieldStartX,y+hFieldStartY,z+hFieldStartZ) = background;
        }
      }
    }

    // trilinear weights
    float w000, w001, w010, w011, w100, w101, w110, w111;

    T ix;
    
    // floor index and residuals
    U hx, hy, hz;
    int fx, fy, fz;
    float rx, ry, rz;

    for (unsigned int z = 0; 
         z < hFieldSize.z && z + hFieldStartZ < imageSize.z; ++z) {
      for (unsigned int y = 0; 
           y < hFieldSize.y && y + hFieldStartY < imageSize.y; ++y) {
        for (unsigned int x = 0; 
             x < hFieldSize.x && x + hFieldStartX < imageSize.x; ++x) {

	  // get intensity value from source image
	  ix = image(x + hFieldStartX,
                     y + hFieldStartY,
                     z + hFieldStartZ);
	  
          // get h field value---where this intensity should go
          hx = hField(x,y,z).x + hFieldStartX;
          hy = hField(x,y,z).y + hFieldStartY;
          hz = hField(x,y,z).z + hFieldStartZ;

          // a fast version of the floor function
	  fx = static_cast<int>(hx);
	  fy = static_cast<int>(hy);
	  fz = static_cast<int>(hz);
          if (hx < 0 && hx != static_cast<int>(hx)) --fx;
          if (hy < 0 && hy != static_cast<int>(hy)) --fy;
          if (hz < 0 && hz != static_cast<int>(hz)) --fz;

	  if (hx > -1 && hx < (int) imageSize.x && // inside vol w/ 1 px border
	      hy > -1 && hy < (int) imageSize.y &&
	      hz > -1 && hz < (int) imageSize.z)
	    {
              // compute trilinear weights
              rx = hx - fx;
              ry = hy - fy;
              rz = hz - fz;
              w000 = (1.0 - rx) * (1.0 - ry) * (1.0 - rz);
              w001 = (1.0 - rx) * (1.0 - ry) * (rz);
              w010 = (1.0 - rx) * (ry)       * (1.0 - rz);
              w011 = (1.0 - rx) * (ry)       * (rz);
              w100 = (rx)       * (1.0 - ry) * (1.0 - rz);
              w101 = (rx)       * (1.0 - ry) * (rz);
              w110 = (rx)       * (ry)       * (1.0 - rz);
              w111 = (rx)       * (ry)       * (rz);

              // see which corners of cube are valid
              bool
                floorXIn = (fx >= 0), ceilXIn = (fx < (int) imageSize.x - 1),
                floorYIn = (fy >= 0), ceilYIn = (fy < (int) imageSize.y - 1),
                floorZIn = (fz >= 0), ceilZIn = (fz < (int) imageSize.z - 1);

              if (floorXIn && floorYIn && floorZIn)
                {
                  defImage(fx, fy, fz)       += w000 * ix;
                  count(fx, fy, fz)          += w000;
                }
              if (floorXIn && floorYIn && ceilZIn)
                {
                  defImage(fx, fy, fz+1)     += w001 * ix;
                  count(fx, fy, fz+1)        += w001;
                }
              if (floorXIn && ceilYIn && floorZIn)
                {
                  defImage(fx, fy+1, fz)     += w010 * ix;
                  count(fx, fy+1, fz)        += w010;
                }
              if (floorXIn && ceilYIn && ceilZIn)
                {
                  defImage(fx, fy+1, fz+1)   += w011 * ix;
                  count(fx, fy+1, fz+1)      += w011;
                }
              if (ceilXIn && floorYIn && floorZIn)
                {
                  defImage(fx+1, fy, fz)     += w100 * ix;
                  count(fx+1, fy, fz)        += w100;
                }
              if (ceilXIn && floorYIn && ceilZIn)
                {
                  defImage(fx+1, fy, fz+1)   += w101 * ix;
                  count(fx+1, fy, fz+1)      += w101;
                }
              if (ceilXIn && ceilYIn && floorZIn)
                {
                  defImage(fx+1, fy+1, fz)   += w110 * ix;
                  count(fx+1, fy+1, fz)      += w110;
                }
              if (ceilXIn && ceilYIn && ceilZIn)
                {
                  defImage(fx+1, fy+1, fz+1) += w111 * ix;
                  count(fx+1, fy+1, fz+1)    += w111;
                }
	    }          
        }
      }
    }    
    
    // find holes
//     float zeta = 1.0/20.0;
//     std::queue<Vector3D<unsigned int> > holes;
//     for (unsigned int z = 0; 
//          z < hFieldSize.z && z + hFieldStartZ < imageSize.z; ++z) {
//       for (unsigned int y = 0; 
//            y < hFieldSize.y && y + hFieldStartY < imageSize.y; ++y) {
//         for (unsigned int x = 0; 
//              x < hFieldSize.x && x + hFieldStartX < imageSize.x; ++x) {
//           if (count(x+hFieldStartX, y+hFieldStartY, z+hFieldStartZ) < zeta) {
//             holes.push(Vector3D<unsigned int>(x+hFieldStartX,
//                                               y+hFieldStartY,
//                                               z+hFieldStartZ));
//           }
//         }
//       }
//     }

    // fill in holes with average of 6-neighbors
//     int maxIters = 10;
//     int iter = 0;
//     float tau = 1.0/3.0;
//     while (holes.size() && iter++ < maxIters) {
//       std::cerr << "filling " << holes.size() << " holes: " 
//                 << iter << std::endl;
//       for (int i = 0; i < (int) holes.size(); ++i) {
//         // pop coords from front of queue
//         unsigned int x = holes.front().x;
//         unsigned int y = holes.front().y;
//         unsigned int z = holes.front().z;
//         holes.pop();

//         if (x+1 < imageSize.x && count(x+1,y,z) > tau) {
//           defImage(x,y,z) += defImage(x+1,y,z);
//           count(x,y,z)    += count(x+1,y,z);
//         }

//         if (x-1 > 0 && count(x-1,y,z) > tau) { 
//           defImage(x,y,z) += defImage(x-1,y,z);
//           count(x,y,z)    += count(x-1,y,z);
//         }

//         if (y+1 < imageSize.y && count(x,y+1,z) > tau) {
//           defImage(x,y,z) += defImage(x,y+1,z);
//           count(x,y,z)    += count(x,y+1,z);
//         }

//         if (y-1 > 0 && count(x,y-1,z) > tau) { 
//           defImage(x,y,z) += defImage(x,y-1,z);
//           count(x,y,z)    += count(x,y-1,z);
//         }

//         if (z+1 < imageSize.z && count(x,y,z+1) > tau) {
//           defImage(x,y,z) += defImage(x,y,z+1);
//           count(x,y,z)    += count(x,y,z+1);
//         }

//         if (z-1 > 0 && count(x,y,z-1) > tau) { 
//           defImage(x,y,z) += defImage(x,y,z-1);
//           count(x,y,z)    += count(x,y,z-1);
//         }
        
//         // add back to queue if did'nt find any neighbors
//         if (count(x,y,z) < zeta) {
//           holes.push(Vector3D<unsigned int>(x,y,z));
//         }
//       }
//     }


      


      // normalize 
    for (unsigned int z = 0; 
           z < hFieldSize.z && z + hFieldStartZ < imageSize.z; ++z) {
        for (unsigned int y = 0; 
             y < hFieldSize.y && y + hFieldStartY < imageSize.y; ++y) {
          for (unsigned int x = 0; 
               x < hFieldSize.x && x + hFieldStartX < imageSize.x; ++x) {
            defImage(x+hFieldStartX, y+hFieldStartY, z+hFieldStartZ) =
              count(x+hFieldStartX, y+hFieldStartY, z+hFieldStartZ) > 0
              ? defImage(x+hFieldStartX, y+hFieldStartY, z+hFieldStartZ)
              / count(x+hFieldStartX, y+hFieldStartY, z+hFieldStartZ)
              : background;
          }
        }
      }
  }

  template <class T, class U>
  static
  void
  forwardApply(const Array3D<T>& image,
               const Array3D<Vector3D<U> >& hField,
               Array3D<T>& defImage,
               const T& background = 0)
  {
    Vector3D<unsigned int> size = image.getSize();
    defImage.fill(0);
    Array3D<float> count(defImage);

    // trilinear weights
    float w000, w001, w010, w011, w100, w101, w110, w111;

    T ix;
    
    // floor index and residuals
    U hx, hy, hz;
    int fx, fy, fz;
    float rx, ry, rz;
    
    for (unsigned int z = 0; z < size.z; ++z) {
      for (unsigned int y = 0; y < size.y; ++y) {
	for (unsigned int x = 0; x < size.x; ++x) {
	  // get intensity value
	  ix = image(x,y,z);
	  
          // get h field value---where this intensity should go
          hx = hField(x,y,z).x;
          hy = hField(x,y,z).y;
          hz = hField(x,y,z).z;

          // this is a fast version of the floor function
	  fx = static_cast<int>(hx);
	  fy = static_cast<int>(hy);
	  fz = static_cast<int>(hz);
          if (hx < 0 && hx != static_cast<int>(hx)) --fx;
          if (hy < 0 && hy != static_cast<int>(hy)) --fy;
          if (hz < 0 && hz != static_cast<int>(hz)) --fz;

	  if (hx > -1 && hx < (int) size.x &&  // inside vol with 1 pix border
	      hy > -1 && hy < (int) size.y &&
	      hz > -1 && hz < (int) size.z)
	    {
              // compute trilinear weights
              rx = hx - fx;
              ry = hy - fy;
              rz = hz - fz;
              w000 = (1.0 - rx) * (1.0 - ry) * (1.0 - rz);
              w001 = (1.0 - rx) * (1.0 - ry) * (rz);
              w010 = (1.0 - rx) * (ry)       * (1.0 - rz);
              w011 = (1.0 - rx) * (ry)       * (rz);
              w100 = (rx)       * (1.0 - ry) * (1.0 - rz);
              w101 = (rx)       * (1.0 - ry) * (rz);
              w110 = (rx)       * (ry)       * (1.0 - rz);
              w111 = (rx)       * (ry)       * (rz);

              // see which corners of cube are valid
              bool
                floorXIn = (fx >= 0), ceilXIn = (fx < (int) size.x - 1),
                floorYIn = (fy >= 0), ceilYIn = (fy < (int) size.y - 1),
                floorZIn = (fz >= 0), ceilZIn = (fz < (int) size.z - 1);

              if (floorXIn && floorYIn && floorZIn)
                {
                  defImage(fx, fy, fz)       += w000 * ix;
                  count(fx, fy, fz)          += w000;
                }
              if (floorXIn && floorYIn && ceilZIn)
                {
                  defImage(fx, fy, fz+1)     += w001 * ix;
                  count(fx, fy, fz+1)        += w001;
                }
              if (floorXIn && ceilYIn && floorZIn)
                {
                  defImage(fx, fy+1, fz)     += w010 * ix;
                  count(fx, fy+1, fz)        += w010;
                }
              if (floorXIn && ceilYIn && ceilZIn)
                {
                  defImage(fx, fy+1, fz+1)   += w011 * ix;
                  count(fx, fy+1, fz+1)      += w011;
                }
              if (ceilXIn && floorYIn && floorZIn)
                {
                  defImage(fx+1, fy, fz)     += w100 * ix;
                  count(fx+1, fy, fz)        += w100;
                }
              if (ceilXIn && floorYIn && ceilZIn)
                {
                  defImage(fx+1, fy, fz+1)   += w101 * ix;
                  count(fx+1, fy, fz+1)      += w101;
                }
              if (ceilXIn && ceilYIn && floorZIn)
                {
                  defImage(fx+1, fy+1, fz)   += w110 * ix;
                  count(fx+1, fy+1, fz)      += w110;
                }
              if (ceilXIn && ceilYIn && ceilZIn)
                {
                  defImage(fx+1, fy+1, fz+1) += w111 * ix;
                  count(fx+1, fy+1, fz+1)    += w111;
                }
	    }
	}
      }
    }    
    
    // normalize counts (NOTE: no rounding for integer types)
    unsigned int numElements = defImage.getNumElements();
    for (unsigned int e = 0; e < numElements; ++e)
      {
	defImage(e) = count(e) > 0 ? defImage(e) / count(e) : background;
      }
  }

  //
  // apply to a region of interest in an image
  // defImage(x+riox) = image(h(x))
  //
  // trilerp by default but will use nearest neighbor
  // if flag is set
  //
  //
  // I DONT THINK THIS IS WHAT WE WANT!!! bcd
  //
  //
  template <class T, class U>
  static
  void
  applyWithROI(const Array3D<T>& image,
  	       const ROI<int, unsigned int>& roi,
  	       const Array3D<Vector3D<U> >& hField,
  	       Array3D<T>& defImage,
  	       const T& background = 0,
  	       bool useNearestNeighbor = false)
  {
    if (useNearestNeighbor)
      {
  	for (unsigned int z = roi.getStartZ(), i=0; 
             z <= roi.getStopZ(); ++z, ++i)
  	  {
  	    for (unsigned int y = roi.getStartY(), j=0; 
                 y <= roi.getStopY(); ++y, ++j)
  	      {
  		for (unsigned int x = roi.getStartX(), k=0; 
                     x <= roi.getStopX(); ++x, ++k)
  		  {
  		    defImage(x, y, z) = 
  		      Array3DUtils::nearestNeighbor(image,
						    hField(k, j, i),
						    background);
  		  }
  	      }
  	  }
      }
    else
      {
  	for (unsigned int z = roi.getStartZ(), i=0; 
             z <= roi.getStopZ(); ++z, ++i)
  	  {
  	    for (unsigned int y = roi.getStartY(), j=0; 
                 y <= roi.getStopY(); ++y, ++j)
  	      {
  		for (unsigned int x = roi.getStartX(), k=0; 
                     x <= roi.getStopX(); ++x, ++k)
  		  {
  		    defImage(x, y, z) = 
  		      static_cast<T>(Array3DUtils::trilerp(image,
  							   hField(k, j, i),
  							   background));
  		  }
  	      }
  	  }
      }
  }

  //
  // apply to surface
  // vertexi = h(vertexi)
  //
  // assumes that surface is in image index coordinates
  //
  // davisb 2003
  //
  template <class T>
  static
  void
  apply(Surface& surface,
	const Array3D<Vector3D<T> >& h)
  {
    unsigned int numVertices = surface.numVertices();
    for (unsigned int i = 0; i < numVertices; ++i)
      {
        float hx, hy, hz;
	trilerp(h,
                (T) surface.vertices[i].x,
		(T) surface.vertices[i].y,
		(T) surface.vertices[i].z,
                hx, hy, hz);
        surface.vertices[i].x = hx;
        surface.vertices[i].y = hy;
        surface.vertices[i].z = hz;
      }
  }

  //
  // applyWithROI to surface
  // vertexi = h(vertexi)
  //
  // assumes that surface is in image index
  // coordinates
  //
  template <class T>
  static
  void
  applyWithROI(Surface& surface,
               const ROI<int, unsigned int>& roi,
	       const Array3D<Vector3D<T> >& h)
  {
    unsigned int numVertices = surface.numVertices();
    for (unsigned int i = 0; i < numVertices; ++i)
      {
	float x = surface.vertices[i].x - roi.getStartX();
	float y = surface.vertices[i].y - roi.getStartY();
	float z = surface.vertices[i].z - roi.getStartZ();
	float Nx,Ny,Nz;
	HField3DUtils::trilerp(h,x,y,z,Nx,Ny,Nz);
	surface.vertices[i].x = Nx + roi.getStartX();
	surface.vertices[i].y = Ny + roi.getStartY();
	surface.vertices[i].z = Nz + roi.getStartZ();
      }
  }

  //
  // inverse apply to surface
  // vertexi = hinv(vertexi)
  //
  // assumes that surface is in image index coordinates
  //
  // davisb 2004
  //
  template <class T>
  static
  void
  inverseApply(Surface& surface,
               const Array3D<Vector3D<T> >& h)
  {
    unsigned int numVertices = surface.numVertices();
    for (unsigned int i = 0; i < numVertices; ++i)
      {
        if (i % 10 == 0) std::cerr << i << std::endl;
        float hinvx, hinvy, hinvz;
        inverseOfPoint(h, 
                       (float) surface.vertices[i].x,
                       (float) surface.vertices[i].y,
                       (float) surface.vertices[i].z,
                       hinvx, hinvy, hinvz);
        surface.vertices[i].x = hinvx;
        surface.vertices[i].y = hinvy;
        surface.vertices[i].z = hinvz;
      }
  }

  template <class T>
  static
  void
  inverseApply(Surface& surface,
               const Array3D<Vector3D<T> >& h,
               const ROI<int, unsigned int>& roi)               
  {
    inverseApply(surface, h, 
                 roi.getStartX(), roi.getStartY(), roi.getStartZ());
  }
               
  template <class T>
  static
  void
  inverseApply(Surface& surface,
               const Array3D<Vector3D<T> >& h,
               int hFieldStartX,
               int hFieldStartY,
               int hFieldStartZ)
  {
    unsigned int numVertices = surface.numVertices();
    for (unsigned int i = 0; i < numVertices; ++i)
      {
	float x = surface.vertices[i].x - hFieldStartX;
	float y = surface.vertices[i].y - hFieldStartY;
	float z = surface.vertices[i].z - hFieldStartZ;
	float Nx,Ny,Nz;
        inverseOfPoint(h, x, y, z, Nx, Ny, Nz);
	surface.vertices[i].x = Nx + hFieldStartX;
	surface.vertices[i].y = Ny + hFieldStartY;
	surface.vertices[i].z = Nz + hFieldStartZ;
      }    
  }

  //
  // create an h field of another size
  // outHField becomes inHField scaled to 
  // outputSizeX x outputSizeY x outputSizeZ
  //
  // davisb 2003
  //
  template<class U>
  static 
  void 
  resample(const Array3D<Vector3D<U> >& inHField, 
	   Array3D<Vector3D<U> >& outHField,
	   const unsigned int outputSizeX,
	   const unsigned int outputSizeY,
	   const unsigned int outputSizeZ,
           BackgroundStrategy backgroundStrategy = BACKGROUND_STRATEGY_ID)
  {
    Vector3D<unsigned int> inSize = inHField.getSize();
    
    if((outputSizeX == inSize.x) &&
       (outputSizeY == inSize.y) &&
       (outputSizeZ == inSize.z))
      {
	// no need to resample, already the same size
	outHField = inHField;
      }
    else
      {
	outHField.resize(outputSizeX, outputSizeY, outputSizeZ);
	
	// scale factors to convert outputSize to inputSize
	double rX = static_cast<double>(inSize.x) 
	  / static_cast<double>(outputSizeX);
	double rY = static_cast<double>(inSize.y) 
	  / static_cast<double>(outputSizeY);
	double rZ = static_cast<double>(inSize.z) 
	  / static_cast<double>(outputSizeZ);
	
	Vector3D<float> hInOfX;
        Vector3D<double> inIndex;
	for (unsigned int z = 0; z < outputSizeZ; ++z){
	  inIndex.z = z * rZ;
	  for (unsigned int y = 0; y < outputSizeY; ++y){
	    inIndex.y = y * rY;
	    for (unsigned int x = 0; x < outputSizeX; ++x){
	      inIndex.x = x * rX;

	      // get vector for corresponding position in inHField
	      trilerp(inHField,
		      (float) inIndex.x, (float) inIndex.y, (float) inIndex.z,
		      hInOfX.x, hInOfX.y, hInOfX.z, 
                      backgroundStrategy);

	      // rescale vector before setting to outHField
	      outHField(x,y,z).set(hInOfX.x / rX,
				   hInOfX.y / rY,
				   hInOfX.z / rZ);
            }
          }
        }
      }
  }

  template<class U>
  static 
  void 
  resample(const Array3D<Vector3D<U> >& inHField, 
	   Array3D<Vector3D<U> >& outHField, 
	   const Vector3D<unsigned int>& outputSize,
           BackgroundStrategy backgroundStrategy = BACKGROUND_STRATEGY_ID)
  {
    resample(inHField, outHField, 
	     outputSize.x, outputSize.y, outputSize.z, backgroundStrategy);
  }

  //
  // trilerp into h field
  //
  // hx, hy, hz gets h(x, y, z)
  //
  // this code has been optimized
  // speed over beauty...
  //
  // davisb 2003
  //
  template <class T>
  static
  void
  trilerp(const Array3D<Vector3D<T> >& h,
	  const T& x, const T& y, const T& z,
	  T& hx, T& hy, T& hz,
//PP
	  BackgroundStrategy backgroundStrategy = BACKGROUND_STRATEGY_ID)
	  //BackgroundStrategy backgroundStrategy = BACKBACKGROUND_STRATEGY_PARTIAL_ID)
  {
    // a (much) faster version of the floor function
    int floorX = static_cast<int>(x);
    int floorY = static_cast<int>(y);
    int floorZ = static_cast<int>(z);
    if (x < 0 && x != static_cast<int>(x)) --floorX;
    if (y < 0 && y != static_cast<int>(y)) --floorY;
    if (z < 0 && z != static_cast<int>(z)) --floorZ;

    // this is not truly ceiling, but floor + 1, which is usually ceiling    
    int ceilX = floorX + 1;
    int ceilY = floorY + 1;
    int ceilZ = floorZ + 1;

    //
    // ^
    // |  v3   v2       -z->        v4   v5
    // y           --next slice-->      
    // |  v0   v1                   v7   v6
    //
    //      -x->
    //   

    T v0X, v0Y, v0Z;
    T v1X, v1Y, v1Z;
    T v2X, v2Y, v2Z;
    T v3X, v3Y, v3Z;
    T v4X, v4Y, v4Z;
    T v5X, v5Y, v5Z;
    T v6X, v6Y, v6Z;
    T v7X, v7Y, v7Z;

    int sizeX = h.getSizeX();
    int sizeY = h.getSizeY();
    int sizeZ = h.getSizeZ();
    if (floorX >= 0 && ceilX < sizeX &&
	floorY >= 0 && ceilY < sizeY &&
	floorZ >= 0 && ceilZ < sizeZ)
      {
	//
	// coordinate is inside volume, fill in 
	// eight corners of cube
	//
	v0X = h(floorX, floorY, floorZ).x;
	v0Y = h(floorX, floorY, floorZ).y;
	v0Z = h(floorX, floorY, floorZ).z;

	v1X = h(ceilX, floorY, floorZ).x;
	v1Y = h(ceilX, floorY, floorZ).y;
	v1Z = h(ceilX, floorY, floorZ).z;

	v2X = h(ceilX, ceilY, floorZ).x;
	v2Y = h(ceilX, ceilY, floorZ).y;
	v2Z = h(ceilX, ceilY, floorZ).z;

	v3X = h(floorX, ceilY, floorZ).x;
	v3Y = h(floorX, ceilY, floorZ).y;
	v3Z = h(floorX, ceilY, floorZ).z;

	v4X = h(floorX, ceilY, ceilZ).x;
	v4Y = h(floorX, ceilY, ceilZ).y;
	v4Z = h(floorX, ceilY, ceilZ).z;

	v5X = h(ceilX, ceilY, ceilZ).x;
	v5Y = h(ceilX, ceilY, ceilZ).y;
	v5Z = h(ceilX, ceilY, ceilZ).z;

	v6X = h(ceilX, floorY, ceilZ).x;
	v6Y = h(ceilX, floorY, ceilZ).y;
	v6Z = h(ceilX, floorY, ceilZ).z;

	v7X = h(floorX, floorY, ceilZ).x;
	v7Y = h(floorX, floorY, ceilZ).y;
	v7Z = h(floorX, floorY, ceilZ).z;
      }
    else if (backgroundStrategy == BACKGROUND_STRATEGY_ID)
      {
	//
	// coordinate is not inside volume, return identity
	//
	hx = x; hy = y; hz = z;
	return;
      }
    else if (backgroundStrategy == BACKGROUND_STRATEGY_ZERO)
      {
        hx = 0; hy = 0; hz = 0;
        return;
      }
    else
      {
	//
	// coordinate is not inside volume; initialize cube
	// corners to identity then set any corners of cube that
	// fall on volume boundary
	//
	v0X = x; v0Y = y; v0Z = z;
	v1X = x; v1Y = y; v1Z = z;
	v2X = x; v2Y = y; v2Z = z;
	v3X = x; v3Y = y; v3Z = z;
	v4X = x; v4Y = y; v4Z = z;
	v5X = x; v5Y = y; v5Z = z;
	v6X = x; v6Y = y; v6Z = z;
	v7X = x; v7Y = y; v7Z = z;
	
	bool floorXIn = floorX >= 0 && floorX < sizeX;
	bool floorYIn = floorY >= 0 && floorY < sizeY;
	bool floorZIn = floorZ >= 0 && floorZ < sizeZ;
	
	bool ceilXIn = ceilX >= 0 && ceilX < sizeX;
	bool ceilYIn = ceilY >= 0 && ceilY < sizeY;
	bool ceilZIn = ceilZ >= 0 && ceilZ < sizeZ;
	
	if (floorXIn && floorYIn && floorZIn)
	  {
	    v0X = h(floorX, floorY, floorZ).x;
	    v0Y = h(floorX, floorY, floorZ).y;
	    v0Z = h(floorX, floorY, floorZ).z;	  
	  }
	if (ceilXIn && floorYIn && floorZIn)
	  {
	    v1X = h(ceilX, floorY, floorZ).x;
	    v1Y = h(ceilX, floorY, floorZ).y;
	    v1Z = h(ceilX, floorY, floorZ).z;
	  }
	if (ceilXIn && ceilYIn && floorZIn)
	  {
	    v2X = h(ceilX, ceilY, floorZ).x;
	    v2Y = h(ceilX, ceilY, floorZ).y;
	    v2Z = h(ceilX, ceilY, floorZ).z;
	  }
	if (floorXIn && ceilYIn && floorZIn)
	  {
	    v3X = h(floorX, ceilY, floorZ).x;
	    v3Y = h(floorX, ceilY, floorZ).y;
	    v3Z = h(floorX, ceilY, floorZ).z;
	  }
	if (floorXIn && ceilYIn && ceilZIn)
	  {
	    v4X = h(floorX, ceilY, ceilZ).x;
	    v4Y = h(floorX, ceilY, ceilZ).y;
	    v4Z = h(floorX, ceilY, ceilZ).z;	  
	  }
	if (ceilXIn && ceilYIn && ceilZIn)
	  {
	    v5X = h(ceilX, ceilY, ceilZ).x;
	    v5Y = h(ceilX, ceilY, ceilZ).y;
	    v5Z = h(ceilX, ceilY, ceilZ).z;	  
	  }
	if (ceilXIn && floorYIn && ceilZIn)
	  {
	    v6X = h(ceilX, floorY, ceilZ).x;
	    v6Y = h(ceilX, floorY, ceilZ).y;
	    v6Z = h(ceilX, floorY, ceilZ).z;	  
	  }
	if (floorXIn && floorYIn && ceilZIn)
	  {
	    v7X = h(floorX, floorY, ceilZ).x;
	    v7Y = h(floorX, floorY, ceilZ).y;
	    v7Z = h(floorX, floorY, ceilZ).z;	  
	  }
      }

    //
    // do trilinear interpolation
    //
    const double t = x - floorX;
    const double u = y - floorY;
    const double v = z - floorZ;
    const double oneMinusT = 1.0 - t;
    const double oneMinusU = 1.0 - u;
    const double oneMinusV = 1.0 - v;
    
    //
    // this is the basic trilerp function...
    //
    //     h = 
    //       v0 * (1 - t) * (1 - u) * (1 - v) +
    //       v1 * t       * (1 - u) * (1 - v) +
    //       v2 * t       * u       * (1 - v) +
    //       v3 * (1 - t) * u       * (1 - v) +
    //       v4 * (1 - t) * u       * v       +
    //       v5 * t       * u       * v       +
    //       v6 * t       * (1 - u) * v       +
    //       v7 * (1 - t) * (1 - u) * v;
    //
    // the following nested version saves 30 multiplies.
    //
    
    hx = 
      oneMinusT * (oneMinusU * (v0X * oneMinusV + v7X * v)  +
		   u         * (v3X * oneMinusV + v4X * v)) +
      t         * (oneMinusU * (v1X * oneMinusV + v6X * v)  +
		   u         * (v2X * oneMinusV + v5X * v));
    
    hy = 
      oneMinusT * (oneMinusU * (v0Y * oneMinusV + v7Y * v)  +
		   u         * (v3Y * oneMinusV + v4Y * v)) +
      t         * (oneMinusU * (v1Y * oneMinusV + v6Y * v)  +
		   u         * (v2Y * oneMinusV + v5Y * v));
    
    hz = 
      oneMinusT * (oneMinusU * (v0Z * oneMinusV + v7Z * v)  +
		   u         * (v3Z * oneMinusV + v4Z * v)) +
      t         * (oneMinusU * (v1Z * oneMinusV + v6Z * v)  +
		   u         * (v2Z * oneMinusV + v5Z * v));
  }
  
  //
  // HField3DUtils::divergence
  //
  // Computes the divergence, "Nabla-Dot", of vector field
  // (Sometimes known as the trace of the Jacobi matrix)
  //
  //          del h1   del h2   del h3
  // div h =  ------ + ------ + ------
  //          del x1   del x2   del x3
  //
  // Where h(x) = [h1 h2 h3] with x = [x1 x2 x3] and 'del' 
  // indicates a partial derivative. 
  //
  // This is intended to be a diagnostic rather than real-time
  // method.
  //
  // Input
  //   hField          - 3D transformation field
  //
  // Input/Output
  //   divergenceImage - 3D scalar divergence image
  //
  // Output
  //   <none>
  //
  // P Lorenzen (2005)
  //
  template <class T, class U>
  static
  void
  divergence(const Array3D<Vector3D<U> >& hField, Array3D<T>& divergence)
  {
    Vector3D<unsigned int> size = hField.getSize();
    unsigned int numElements = hField.getNumElements();

    // resize divergence array if necessary
    if (divergence.getSize() != size)
      {
	divergence.resize(size);
      }

    // build scalar images h1, h2, and h3
    // from the transformation field
    Array3D<U> h1(size);
    Array3D<U> h2(size);
    Array3D<U> h3(size);
    unsigned int i; // stupid vc++
    for (i = 0; i < numElements; ++i)
      {
	h1(i) = hField(i).x;
	h2(i) = hField(i).y;
	h3(i) = hField(i).z;
      }
    
    // compute the gradients of h1, h2, and h3
    Array3D<Vector3D<U> > grad_h1(size);
    Array3D<Vector3D<U> > grad_h2(size);
    Array3D<Vector3D<U> > grad_h3(size);   
    Array3DUtils::computeGradient(h1, grad_h1);
    Array3DUtils::computeGradient(h2, grad_h2);
    Array3DUtils::computeGradient(h3, grad_h3);

    // compute the divergence
    T t1, t2, t3;
    for (i = 0; i < numElements; ++i)
      {
	t1 = static_cast<T>(grad_h1(i).x);
	t2 = static_cast<T>(grad_h2(i).y);
	t3 = static_cast<T>(grad_h3(i).z);

	divergence(i) = t1 + t2 + t3;
      }
  }
  

  //
  // HField3DUtils::jacobian
  //
  // The Jacobian (sometimes known as the determinant of the 
  // Jacobi matrix)
  //
  //        | del h1   del h1   del h1 |
  //        | ------   ------   ------ |
  //        | del x1   del x2   del x3 |
  //        |                          |
  //        | del h2   del h2   del h2 |
  // J(h) = | ------   ------   ------ |
  //        | del x1   del x2   del x3 |
  //        |                          |
  //        | del h3   del h3   del h3 |
  //        | ------   ------   ------ |
  //        | del x1   del x2   del x3 |
  //
  // Where h(x) = [h1 h2 h3] with x = [x1 x2 x3] and 'del' 
  // indicates a partial derivative. 
  //
  // This is intended to be a diagnostic rather than real-time
  // method.
  //
  // Input
  //   hField        - 3D transformation field
  //
  // Input/Output
  //   jacobianImage - 3D scalar jacobian image
  //
  // Output
  //   <none>
  //
  // P Lorenzen (2003)
  // davisb 2003
  //
  template <class T, class U>
  static
  void
  jacobian(const Array3D<Vector3D<U> >& hField,
	   Array3D<T>& jacobian)
  {
    Vector3D<unsigned int> size = hField.getSize();
    unsigned int numElements = hField.getNumElements();

    // resize jacobian array if necessary
    if (jacobian.getSize() != size)
      {
	jacobian.resize(size);
      }

    // build scalar images h1, h2, and h3
    // from the transformation field
    Array3D<U> h1(size);
    Array3D<U> h2(size);
    Array3D<U> h3(size);
    unsigned int i; // stupid vc++
    for (i = 0; i < numElements; ++i)
      {
	h1(i) = hField(i).x;
	h2(i) = hField(i).y;
	h3(i) = hField(i).z;
      }
    
    // compute the gradients of h1, h2, and h3
    Array3D<Vector3D<U> > grad_h1(size);
    Array3D<Vector3D<U> > grad_h2(size);
    Array3D<Vector3D<U> > grad_h3(size);   
    Array3DUtils::computeGradient(h1, grad_h1);
    Array3DUtils::computeGradient(h2, grad_h2);
    Array3DUtils::computeGradient(h3, grad_h3);

    // compute the jacobian
    T t1, t2, t3;
    for (i = 0; i < numElements; ++i)
      {
	t1 = static_cast<T>(grad_h1(i).x * (grad_h2(i).y * grad_h3(i).z -
					    grad_h2(i).z * grad_h3(i).y));
	t2 = static_cast<T>(grad_h1(i).y * (grad_h2(i).x * grad_h3(i).z -
					    grad_h2(i).z * grad_h3(i).x));
	t3 = static_cast<T>(grad_h1(i).z * (grad_h2(i).x * grad_h3(i).y -
					    grad_h2(i).y * grad_h3(i).x));
	jacobian(i) = t1 - t2 + t3;
      }
  }

  //
  // compute the minimum and maximum deformation distance
  //
  // bcd 2004
  //
  template <class T>
  static
  void
  minMaxDeformationL2Norm(const Array3D<Vector3D<T> >& hField,
			  double& min, double& max)
  {
    Vector3D<unsigned int> size = hField.getSize();

    //
    // compute min, max
    //
    min = std::numeric_limits<double>::max();
    max = 0;
    for (unsigned int z = 0; z < size.z; ++z) {    
      for (unsigned int y = 0; y < size.y; ++y) {    
	for (unsigned int x = 0; x < size.x; ++x) {    
	  double normL2 = 
	    sqrt((hField(x,y,z).x - x) * (hField(x,y,z).x - x) +
		 (hField(x,y,z).y - y) * (hField(x,y,z).y - y) +
		 (hField(x,y,z).z - z) * (hField(x,y,z).z - z));
	  if (normL2 < min) min = normL2;
	  if (normL2 > max) max = normL2;
	}
      }
    }
  }

  //
  // compute the minimum and maximum velocity
  //
  // bcd 2004
  //
  template <class T>
  static
  void
  minMaxVelocityL2Norm(const Array3D<Vector3D<T> >& hField,
			  double& min, double& max)
  {
    Vector3D<unsigned int> size = hField.getSize();

    //
    // compute min, max
    //
    min = std::numeric_limits<double>::max();
    max = 0;
    for (unsigned int z = 0; z < size.z; ++z) {    
      for (unsigned int y = 0; y < size.y; ++y) {    
	for (unsigned int x = 0; x < size.x; ++x) {    
	  double normL2 = 
	    sqrt((hField(x,y,z).x) * (hField(x,y,z).x) +
		 (hField(x,y,z).y) * (hField(x,y,z).y) +
		 (hField(x,y,z).z) * (hField(x,y,z).z));
	  if (normL2 < min) min = normL2;
	  if (normL2 > max) max = normL2;
	}
      }
    }
  }

  //
  // experimental, compute jacobian at non-grid point
  //
  // bcd 2004
  //
  template <class T>
  static
  void
  jacobian(const Array3D<Vector3D<T> >& h,
	   const T& x, const T& y, const T& z,
	   double* const J)
  {
    // a (much) faster version of the floor function
    int floorX = static_cast<int>(x);
    int floorY = static_cast<int>(y);
    int floorZ = static_cast<int>(z);
    if (x < 0 && x != static_cast<int>(x)) --floorX;
    if (y < 0 && y != static_cast<int>(y)) --floorY;
    if (z < 0 && z != static_cast<int>(z)) --floorZ;
    int ceilX = floorX + 1;
    int ceilY = floorY + 1;
    int ceilZ = floorZ + 1;

    //
    // ^
    // |  J3   J2       -z->        J4   J5
    // y           --next slice-->      
    // |  J0   J1                   J7   J6
    //
    //      -x->
    //   

    double J0[9]; double J1[9]; double J2[9]; double J3[9]; 
    double J4[9]; double J5[9]; double J6[9]; double J7[9];

    jacobianAtGridpoint(h, floorX, floorY, floorZ, J0);
    jacobianAtGridpoint(h, ceilX,  floorY, floorZ, J1);
    jacobianAtGridpoint(h, ceilX,  ceilY,  floorZ, J2);
    jacobianAtGridpoint(h, floorX, ceilY,  floorZ, J3);
    jacobianAtGridpoint(h, floorX, ceilY,  ceilZ,  J4);
    jacobianAtGridpoint(h, ceilX,  ceilY,  ceilZ,  J5);
    jacobianAtGridpoint(h, ceilX,  floorY, ceilZ,  J6);
    jacobianAtGridpoint(h, floorX, floorY, ceilZ,  J7);

    //
    // do trilinear interpolation
    //
    const double t = x - floorX;
    const double u = y - floorY;
    const double v = z - floorZ;
    const double oneMinusT = 1.0 - t;
    const double oneMinusU = 1.0 - u;
    const double oneMinusV = 1.0 - v;

    //
    // this is the basic trilerp function...
    //
    //     h = 
    //       v0 * (1 - t) * (1 - u) * (1 - v) +
    //       v1 * t       * (1 - u) * (1 - v) +
    //       v2 * t       * u       * (1 - v) +
    //       v3 * (1 - t) * u       * (1 - v) +
    //       v4 * (1 - t) * u       * v       +
    //       v5 * t       * u       * v       +
    //       v6 * t       * (1 - u) * v       +
    //       v7 * (1 - t) * (1 - u) * v;
    //
    // the following nested version saves 30 multiplies.
    //

    for (int i = 0; i < 9; ++i)
      {
	J[i] =
	  oneMinusT * (oneMinusU * (J0[i] * oneMinusV + J7[i] * v)  +
		       u         * (J3[i] * oneMinusV + J4[i] * v)) +
	  t         * (oneMinusU * (J1[i] * oneMinusV + J6[i] * v)  +
		       u         * (J2[i] * oneMinusV + J5[i] * v));
      }
  }
  
  //
  // compute the Jacobian of the transformation at a grid point
  //
  // bcd 2004
  //
  template <class T>
  static
  void
  jacobianAtGridpoint(const Array3D<Vector3D<T> >& h,
		      int x, int y, int z,
		      double* const J)
  {
    Vector3D<int> size = h.getSize();

    //
    // arbitrarily set Jacobian to identity if outside region
    //
    if (x < 0 || x >= size.x ||
	y < 0 || y >= size.y ||
	z < 0 || z >= size.z)
      {
	J[0] = J[4] = J[8] = 1;
	J[1] = J[2] = J[3] = J[5] = J[6] = J[7] = 0;
	return;
      }

    //
    // do symmetric difference unless on the edge, then do forward
    // difference
    //

    if (x > 0 && x < (size.x - 1))
      {
	// symmetric difference
	J[0] = (h(x+1,y,z).x - h(x-1,y,z).x) / 2.0;
	J[3] = (h(x+1,y,z).y - h(x-1,y,z).y) / 2.0;
	J[6] = (h(x+1,y,z).z - h(x-1,y,z).z) / 2.0;
      }
    else if (x == 0)
      {	  
	// forward difference
	J[0] = (h(x+1,y,z).x - h(x,y,z).x);
	J[3] = (h(x+1,y,z).y - h(x,y,z).y);
	J[6] = (h(x+1,y,z).z - h(x,y,z).z);
      }
    else
      {
	// backwards difference
	J[0] = (h(x,y,z).x - h(x-1,y,z).x);
	J[3] = (h(x,y,z).y - h(x-1,y,z).y);
	J[6] = (h(x,y,z).z - h(x-1,y,z).z);
      }

    if (y > 0 && y < (size.y - 1))
      {
	// symmetric difference
	J[1] = (h(x,y+1,z).x - h(x,y-1,z).x) / 2.0;
	J[4] = (h(x,y+1,z).y - h(x,y-1,z).y) / 2.0;
	J[7] = (h(x,y+1,z).z - h(x,y-1,z).z) / 2.0;
      }
    else if (y == 0)
      {
	// forward difference
	J[1] = (h(x,y+1,z).x - h(x,y,z).x);
	J[4] = (h(x,y+1,z).y - h(x,y,z).y);
	J[7] = (h(x,y+1,z).z - h(x,y,z).z);
      }
    else
      {
	// backwards difference
	J[1] = (h(x,y,z).x - h(x,y-1,z).x);
	J[4] = (h(x,y,z).y - h(x,y-1,z).y);
	J[7] = (h(x,y,z).z - h(x,y-1,z).z);
      }
    
    if (z > 0 && z < (size.z - 1))
      {
	// symmetric difference
	J[2] = (h(x,y,z+1).x - h(x,y,z-1).x) / 2.0;
	J[5] = (h(x,y,z+1).y - h(x,y,z-1).y) / 2.0;
	J[8] = (h(x,y,z+1).z - h(x,y,z-1).z) / 2.0;
      }
    else if (z == 0)
      {
	// forward difference
	J[2] = (h(x,y,z+1).x - h(x,y,z).x);
	J[5] = (h(x,y,z+1).y - h(x,y,z).y);
	J[8] = (h(x,y,z+1).z - h(x,y,z).z);
      }
    else
      {
	// backwards difference
	J[2] = (h(x,y,z).x - h(x,y,z-1).x);
	J[5] = (h(x,y,z).y - h(x,y,z-1).y);
	J[8] = (h(x,y,z).z - h(x,y,z-1).z);
      }
  }

  //
  // experimental
  // bcd 2004
  // dont worry about dxdy = dydx right now
  //
  template <class T>
  static
  void
  hessianAtGridpoint(const Array3D<Vector3D<T> >& h,
		     int x, int y, int z,
		     double* const H)
  {
    Vector3D<int> size = h.getSize();

    //
    // arbitrarily set hessian to zero if outside region
    //
    if (x < 1 || x >= size.x - 1 ||
	y < 1 || y >= size.y - 1 ||
	z < 1 || z >= size.z - 1)
      {
	for (int i = 0; i < 27; ++i)
	  {
	    H[i] = 0;
	  }
	return;
      }

    //
    // non mixed partials
    //
    // dx^2
    H[0]  = h(x+1,y,z).x - 2.0 * h(x,y,z).x + h(x-1,y,z).x;
    H[9]  = h(x+1,y,z).y - 2.0 * h(x,y,z).y + h(x-1,y,z).y;
    H[18] = h(x+1,y,z).z - 2.0 * h(x,y,z).z + h(x-1,y,z).z;

    // dy^2
    H[4]  = h(x,y+1,z).x - 2.0 * h(x,y,z).x + h(x,y-1,z).x;
    H[13] = h(x,y+1,z).y - 2.0 * h(x,y,z).y + h(x,y-1,z).y;
    H[22] = h(x,y+1,z).z - 2.0 * h(x,y,z).z + h(x,y-1,z).z;

    // dz^2
    H[8]  = h(x,y,z+1).x - 2.0 * h(x,y,z).x + h(x,y,z-1).x;
    H[17] = h(x,y,z+1).y - 2.0 * h(x,y,z).y + h(x,y,z-1).y;
    H[26] = h(x,y,z+1).z - 2.0 * h(x,y,z).z + h(x,y,z-1).z;
    
    //
    // mixed partials
    //
    // dxdy
    H[1]  = H[3]  =
      + 0                  + h(x,y-1,z).x/2.0 - h(x+1,y-1,z).x/2.0
      + h(x-1,y,z).x/2.0   - h(x,y,z).x       + h(x+1,y,z).x/2.0
      - h(x-1,y+1,z).x/2.0 + h(x,y+1,z).x/2.0 + 0;
    H[10] = H[12] = 
      + 0                  + h(x,y-1,z).y/2.0 - h(x+1,y-1,z).y/2.0
      + h(x-1,y,z).y/2.0   - h(x,y,z).y       + h(x+1,y,z).y/2.0
      - h(x-1,y+1,z).y/2.0 + h(x,y+1,z).y/2.0 + 0;
    H[19] = H[21] = 
      + 0                  + h(x,y-1,z).z/2.0 - h(x+1,y-1,z).z/2.0
      + h(x-1,y,z).z/2.0   - h(x,y,z).z       + h(x+1,y,z).z/2.0
      - h(x-1,y+1,z).z/2.0 + h(x,y+1,z).z/2.0 + 0;

    // dxdz
    H[2]  = H[6]  = 
      + 0                  + h(x,y,z-1).x/2.0 - h(x+1,y,z-1).x/2.0
      + h(x-1,y,z).x/2.0   - h(x,y,z).x       + h(x+1,y,z).x/2.0
      - h(x-1,y,z+1).x/2.0 + h(x,y,z+1).x/2.0 + 0;
    H[11] = H[15] = 
      + 0                  + h(x,y,z-1).y/2.0 - h(x+1,y,z-1).y/2.0
      + h(x-1,y,z).y/2.0   - h(x,y,z).y       + h(x+1,y,z).y/2.0
      - h(x-1,y,z+1).y/2.0 + h(x,y,z+1).y/2.0 + 0;
    H[20] = H[24] = 
      + 0                  + h(x,y,z-1).z/2.0 - h(x+1,y,z-1).z/2.0
      + h(x-1,y,z).z/2.0   - h(x,y,z).z       + h(x+1,y,z).z/2.0
      - h(x-1,y,z+1).z/2.0 + h(x,y,z+1).z/2.0 + 0;

    // dydz
    H[5]  = H[7]  = 
      + 0                  + h(x,y,z-1).x/2.0 - h(x,y+1,z-1).x/2.0
      + h(x,y-1,z).x/2.0   - h(x,y,z).x       + h(x,y+1,z).x/2.0
      - h(x,y-1,z+1).x/2.0 + h(x,y,z+1).x/2.0 + 0;
    H[14] = H[16] = 
      + 0                  + h(x,y,z-1).y/2.0 - h(x,y+1,z-1).y/2.0
      + h(x,y-1,z).y/2.0   - h(x,y,z).y       + h(x,y+1,z).y/2.0
      - h(x,y-1,z+1).y/2.0 + h(x,y,z+1).y/2.0 + 0;
    H[23] = H[25] = 
      + 0                  + h(x,y,z-1).z/2.0 - h(x,y+1,z-1).z/2.0
      + h(x,y-1,z).z/2.0   - h(x,y,z).z       + h(x,y+1,z).z/2.0
      - h(x,y-1,z+1).z/2.0 + h(x,y,z+1).z/2.0 + 0;
  }

  template <class T>
  static
  bool
  inverseOfPoint(const Array3D<Vector3D<T> >& h,
                 const T& x, const T& y, const T& z,
                 T& hinvx, T& hinvy, T& hinvz,
                 float thresholdDistance = 2.0f)
  {
    // get a coarse estimate
    int estimateX, estimateY, estimateZ;
    float dist = inverseClosestPoint(h, 
                        x, y, z,
                        estimateX, estimateY, estimateZ);

    if (dist > thresholdDistance) {
      hinvx = x;
      hinvy = y;
      hinvz = z;
      return false;
    }

    // refine estimate
    return inverseOfPointRefine(h, 
                                x, y, z,
                                estimateX, estimateY, estimateZ,
                                hinvx, hinvy, hinvz);
  }

  template <class T>
  static
  float 
  inverseClosestPoint(const Array3D<Vector3D<T> >& h,
                      const T& x, const T& y, const T& z,
                      int& hinvx, int& hinvy, int& hinvz)
  {
    if (h.getNumElements() == 0) return HUGE_VAL;

    Vector3D<T> xVec(x,y,z);
    Vector3D<int> guess(0,0,0);
    double minDistSq = h(0,0,0).distanceSquared(xVec);
    
    // search on a coarse grid
    int gridsize = 3;
    for (unsigned int z = 0; z < h.getSizeZ(); z += gridsize) {
      for (unsigned int y = 0; y < h.getSizeY(); y += gridsize) {
        for (unsigned int x = 0; x < h.getSizeX(); x += gridsize) {
          double currDistSq = h(x,y,z).distanceSquared(xVec);
          if (currDistSq < minDistSq) 
            {
              minDistSq = currDistSq;
              guess.x = x;
              guess.y = y;
              guess.z = z;
            }
        }
      }
    }

    int slop = 10;
    unsigned int xRegionMin = guess.x - slop > 0 ? guess.x - slop : 0;
    unsigned int xRegionMax = 
      guess.x + slop < (int) h.getSizeX() 
      ? guess.x + slop
      : h.getSizeX() - 1;
    unsigned int yRegionMin = guess.y - slop > 0 ? guess.y - slop : 0;
    unsigned int yRegionMax =
      guess.y + slop < (int) h.getSizeY() 
      ? guess.y + slop
      : h.getSizeY() - 1;
    unsigned int zRegionMin = guess.z - slop > 0 ? guess.z - slop : 0;
    unsigned int zRegionMax =
      guess.z + slop < (int) h.getSizeZ() 
      ? guess.z + slop
      : h.getSizeZ() - 1;
    
    // search around best coarse point
    for (unsigned int z = zRegionMin; z <= zRegionMax; ++z) {
      for (unsigned int y = yRegionMin; y <= yRegionMax; ++y) {
        for (unsigned int x = xRegionMin; x <= xRegionMax; ++x) {
          double currDistSq = h(x,y,z).distanceSquared(xVec);
          if (currDistSq < minDistSq) 
            {
              minDistSq = currDistSq;
              guess.x = x;
              guess.y = y;
              guess.z = z;
            }
        }
      }
    }

    hinvx = guess.x;
    hinvy = guess.y;
    hinvz = guess.z;
    return minDistSq;
  }

  //
  // find hinvx such that h(hinvx) = x.  near guess x0.
  //
  template <class T>
  static
  bool
  inverseOfPointRefine(const Array3D<Vector3D<T> >& h,
                       const T& x, const T& y, const T& z, 
                       const int& x0, const int& y0, const int& z0,
                       T& hinvx, T& hinvy, T& hinvz)
  {
    if (x0 < 0 || x0 >= (int) h.getSizeX() ||
        y0 < 0 || y0 >= (int) h.getSizeY() ||
        z0 < 0 || z0 >= (int) h.getSizeZ())
      {
        return false;
      }

    //
    // approximate the inverse using a first order Taylor expansion,
    // according to the following derivation
    // 
    // hInv(x) = x0 + d
    // h(x0 + d) = x
    // h(x0 + d) = h(x0) + Jh|x0 * d + h.o.t.
    // x = h(x0) + Jh|x0 * d + h.o.t.
    // d = JInvh|x0 * (x - h(x0))
    // 
    // hinv(x) = x0 + JInvh|x0 * (x - h(x0))
    //    
    double J[9];
    double JInv[9];
    jacobianAtGridpoint(h, x0, y0, z0, J);
    bool success = Matrix3D<double>::computeInverse(J, JInv);
    if (!success) 
      {
        // J(h(x)) was not invertable, arbitrarily go with identity
        std::cerr << "At (" << x << "," << y << "," << z 
                  << "): J(h(x)) is non-invertable." << std::endl;
        hinvx = x;
        hinvy = y;
        hinvz = z;
      }
    else
      {
        hinvx = x0 
          + JInv[0] * (x - h(x0,y0,z0).x) 
          + JInv[1] * (y - h(x0,y0,z0).y) 
          + JInv[2] * (z - h(x0,y0,z0).z); 
        hinvy = y0 
          + JInv[3] * (x - h(x0,y0,z0).x) 
          + JInv[4] * (y - h(x0,y0,z0).y) 
          + JInv[5] * (z - h(x0,y0,z0).z); 
        hinvz = z0 
          + JInv[6] * (x - h(x0,y0,z0).x) 
          + JInv[7] * (y - h(x0,y0,z0).y) 
          + JInv[8] * (z - h(x0,y0,z0).z); 
      }    
    return true;
  }

  // The function 'func' should take three coordinates representing a
  // position in world coordinates, and return a vector representing
  // an offset (essentially a velocity) in voxel coordinates.  The
  // function 'fillByFunction' produces an hfield 'h' that corresponds
  // to the given deformation formula and has the specified origin and
  // spacing.
  static void
  fillByFunction(Array3D< Vector3D<float> >& h,
                 Vector3D<double> origin,
                 Vector3D<double> spacing,
                 Vector3D<float> (*func)(double, double, double))
  {
    for (unsigned int z = 0; z < h.getSizeZ(); ++z) {
    for (unsigned int y = 0; y < h.getSizeY(); ++y) {
    for (unsigned int x = 0; x < h.getSizeX(); ++x) {
      Vector3D<double> point(x, y, z);
      Vector3D<double> pointWC = point * spacing + origin;
      Vector3D<float> offset = func(pointWC.x, pointWC.y, pointWC.z);
      Vector3D<float> hFieldVal = point + offset / spacing;
      h(x, y, z) = hFieldVal;
    }
    }
    }
  }  
};
#endif
