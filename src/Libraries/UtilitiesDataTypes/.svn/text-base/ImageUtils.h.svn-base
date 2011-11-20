#ifndef ImageUtils_h
#define ImageUtils_h

#include <float.h>

#include "Array3DUtils.h"
#include "Array3DIO.h"
#include "AffineTransform3D.h"
#include "Image.h"

class ImageUtils
{
public:

  //
  // return the intensity value for this image at a given position in
  // world coordinates
  //
  // bcd 2004
  //
  template <class T>
  static
  T
  trilerp(const Image<T>& image,
          const Vector3D<double>& worldCoordinates,
          const T& background)
  {
    Vector3D<double> voxelCoordinates;
    image.worldToImageIndexCoordinates(worldCoordinates,
                                       voxelCoordinates);
    return Array3DUtils::trilerp(image,
                                 voxelCoordinates,
                                 background);
  }

  //
  // translate this image.  note: this just changes the origin
  //
  // bcd 2004
  //
  template <class T>
  static
  void
  translate(Image<T>& image,
	    const double& tx, 
	    const double& ty, 
	    const double& tz)
  {
    typename Image<T>::ContinuousPointType origin = image.getOrigin();
    image.setOrigin(origin.x + tx, origin.y + ty, origin.z + tz);
  }

  //
  // translate this image.  note: this just changes the origin
  //
  // bcd 2004
  //
  template <class T>
  static
  void
  translate(Image<T>& image,
	    const Vector3D<double>& t)
  {
    typename Image<T>::ContinuousPointType origin = image.getOrigin();
    image.setOrigin(origin.x + t.x, origin.y + t.y, origin.z + t.z);
  }

  //
  // make this image have the given origin, spacing, and dimensions.
  // intensities should stay in the same place in world coordinates.
  //
  // bcd 2004
  //
  template <class T>
  static
  void
  resample(Image<T>& image,
           const Vector3D<double>& newOrigin,
           const Vector3D<double>& newSpacing,
           const Vector3D<unsigned int>& newDimensions)
  {
    Image<T> tmp(image);
    image.resize(newDimensions);
    image.setOrigin(newOrigin);
    image.setSpacing(newSpacing);
    resample(tmp, image);
  }

  //
  // fill in the destination image from the source image, taking
  // spacing, origin, and dimensions into account.  where overlapping,
  // source and dest will be the same in world coordinates.
  //
  // bcd 2004
  //
  template <class T>
  static
  void
  resample(const Image<T>& sourceImage,
	   Image<T>& destImage)
  {
    Vector3D<unsigned int> size = destImage.getSize();
    Vector3D<double> world;
    Vector3D<double> voxel;
    for (unsigned int z = 0; z < size.z; ++z) {
      for (unsigned int y = 0; y < size.y; ++y) {
	for (unsigned int x = 0; x < size.x; ++x) {
	  destImage.imageIndexToWorldCoordinates(x, y, z, 
                                                 world.x, world.y, world.z);

          sourceImage.worldToImageIndexCoordinates(world, voxel);

          destImage(x,y,z) = Array3DUtils::trilerp(sourceImage, voxel, 0.0F);
	}
      }
    }
  }  

  //
  // Fill in the destination image from the source image, taking
  // spacing, origin, and dimensions into account.  Where overlapping,
  // source and dest will be the same in world coordinates.  Where not
  // overlapping, destImage is left unchanged.  This requires that
  // sourceImage have only nonnegative intensities.
  //
  // foskey 2005
  //
  template <class T>
  static
  void
  resampleWithTransparency(const Image<T>& sourceImage,
                           Image<T>& destImage)
  {
    Vector3D<unsigned int> size = destImage.getSize();
    Vector3D<double> world;
    Vector3D<double> voxel;
    for (unsigned int z = 0; z < size.z; ++z) {
      for (unsigned int y = 0; y < size.y; ++y) {
	for (unsigned int x = 0; x < size.x; ++x) {
	  destImage.imageIndexToWorldCoordinates(x, y, z, 
                                                 world.x, world.y, world.z);

          sourceImage.worldToImageIndexCoordinates(world, voxel);

          float intensity = Array3DUtils::trilerp(sourceImage, voxel,
                                                  -FLT_MAX);
          if (intensity >= 0.0f) destImage(x,y,z) = intensity;
          
	}
      }
    }
  }  

  // Take a transform in world coordinates and return the corresponding
  // transform expressed in units of voxels.  
  //
  // For the input transform, the Spacing parameters give the
  // distances in world units between successive voxel centers, and
  // the Origin parameters give the location of the (0,0,0) voxel in
  // world coordinates.
  //
  // For the returned transform, the coordinate systems of the two
  // images have units of voxels, and the origins are the centers of
  // the (0,0,0) voxels.
  static
  AffineTransform3D<double>
  transformInIndexCoordinates(
    const Vector3D<double>& fixedImageOrigin,
    const Vector3D<double>& fixedImageSpacing,
    const Vector3D<double>& movingImageOrigin,
    const Vector3D<double>& movingImageSpacing,
    const AffineTransform3D<double>& transformInWorldCoordinates )
  {
    Matrix3D<double> S1inv;
    Matrix3D<double> S2;
    for (unsigned int i = 0; i < 3; ++i) {
      S1inv(i,i) = 1.0 / movingImageSpacing[i];
      S2(i,i) = fixedImageSpacing[i];
    }

    AffineTransform3D<double> result;
    const Matrix3D<double>& wM = transformInWorldCoordinates.matrix;
    const Vector3D<double>& wV = transformInWorldCoordinates.vector;
    result.matrix = S1inv * wM * S2;
    result.vector = S1inv * (wV - movingImageOrigin + wM * fixedImageOrigin);
    return result;
  }

  template <class T>
  static
  void
  applyAffine(Image<T>& image,
              const Vector3D<double>& newOrigin,
              const Vector3D<double>& newSpacing,
              const Vector3D<unsigned int>& newDimensions,
              const AffineTransform3D<double>& transformInWorldCoordinates, 
              const float& backgroundValue = 0)
  {
    Image<T> tmp(image);
    image.resize(newDimensions);
    image.setOrigin(newOrigin);
    image.setSpacing(newSpacing);
    applyAffine(tmp, image, transformInWorldCoordinates, backgroundValue);
  }

  // Transform sourceImage into the coordinate system of destImage.
  // After the call, destImage has data from sourceImage, transformed
  // by the inverse of 'transform'.  To find the intensity of a point
  // p in destImage, one looks to the point transform(p) in sourceImage.
  // If transform(p) is outside of sourceImage, backgroundValue is
  // used.  The origins of the respective coordinate systems are those
  // specified by getOrigin() for the images.
  template <class T>
  static
  void
  applyAffine(const Image<T>& sourceImage,
              Image<T>& destImage,
              const AffineTransform3D<double>& transformInWorldCoordinates,
              const float& backgroundValue = 0)
  {

    // Compute transform in index coordinates
    AffineTransform3D<double> newTransform = 
      transformInIndexCoordinates( destImage.getOrigin(),
                                   destImage.getSpacing(),
                                   sourceImage.getOrigin(),
                                   sourceImage.getSpacing(),
                                   transformInWorldCoordinates );

    // Get data pointer to write to quickly
    T* newImageDataPtr = destImage.getDataPointer();

    // Determine iteration bounds
    size_t zSize = destImage.getSizeZ();
    size_t ySize = destImage.getSizeY();
    size_t xSize = destImage.getSizeX();

    // Iterate through and interpolate the data
    double x, y, z;
    for( unsigned int zPos = 0; zPos < zSize; zPos++ ) {
      for( unsigned int yPos = 0; yPos < ySize; yPos++ ) {
        for( unsigned int xPos = 0; xPos < xSize; xPos++ ) {

          newTransform.transformCoordinates( xPos, yPos, zPos, x, y, z );

          // Store interpolated value
          *newImageDataPtr++ = Array3DUtils::trilerp(
            sourceImage, x, y, z, backgroundValue);
        }
      }
    }
  }

  //
  // make all voxel spacing the same.  the smallest current voxel
  // spacing is chosen as the spacing to use.
  //
  // bcd 2004
  //
  template <class T>
  static
  void
  makeVoxelsIsotropic(Image<T>& image)
  {
    Vector3D<double> spacingOld = image.getSpacing();    
    Vector3D<unsigned int> sizeOld = image.getSize();
    
    if (spacingOld.x == spacingOld.y && spacingOld.x == spacingOld.z)
      {
        return;
      }

    // 
    // set new spacing in all dimensions to min spacing
    //
    double minSpacing = std::min(spacingOld.x, 
                                 std::min(spacingOld.y, spacingOld.z));
    Vector3D<double> 
      spacingNew(minSpacing * spacingOld.x < 0 ? -minSpacing : minSpacing,
                 minSpacing * spacingOld.y < 0 ? -minSpacing : minSpacing,
                 minSpacing * spacingOld.z < 0 ? -minSpacing : minSpacing);
    
    //
    // origin stays same
    //
    Vector3D<double> originNew(image.getOrigin());

    //
    // determine new dimensions
    //
    Vector3D<unsigned int> 
      sizeNew((unsigned int)floor(double(sizeOld.x) * 
                                  spacingOld.x / spacingNew.x),
              (unsigned int)floor(double(sizeOld.y) * 
                                  spacingOld.y / spacingNew.y),
              (unsigned int)floor(double(sizeOld.z) * 
                                  spacingOld.z / spacingNew.z));
              
    //
    // resample
    //
    resample(image, originNew, spacingNew, sizeNew); 
  }

  //
  // bcd 2004
  //
  template <class T>
  static
  void
  resliceZMakeIsotropic(Image<T>& image)
  {
    Vector3D<double> spacingOld = image.getSpacing();    
    Vector3D<unsigned int> sizeOld = image.getSize();
    
    if (spacingOld.x == spacingOld.z)
      {
        return;
      }

    // 
    // set new z spacing to x spacing (maintain sign)
    //
    Vector3D<double> 
      spacingNew(spacingOld.x, spacingOld.y,
                 spacingOld.x * spacingOld.z < 0 
                 ? -spacingOld.x 
                 : spacingOld.x);
    
    //
    // origin stays same
    //
    Vector3D<double> originNew(image.getOrigin());

    //
    // determine new dimensions
    //
    Vector3D<unsigned int> 
      sizeNew(sizeOld.x, sizeOld.y,
              (unsigned int)floor(double(sizeOld.z) * 
                                  spacingOld.z / spacingNew.z) - 1);
              
    //
    // resample
    //
    resample(image, originNew, spacingNew, sizeNew); 
  }

  //
  // gaussian downsample an image.  spacing is updated appropriatly.
  //
  // bcd 2004
  //
  template <class T>
  static
  void
  gaussianDownsample(const Image<T>& imageIn,
                     Image<T>& imageOut,
                     const Vector3D<double>& factors,
                     const Vector3D<double>& sigma,
                     const Vector3D<double>& kernelSize)
  {
    Array3DUtils::gaussianDownsample(imageIn, imageOut,
                                     factors, 
                                     sigma,
                                     kernelSize);
    Vector3D<double> spacingOld = imageIn.getSpacing();
    imageOut.setSpacing(spacingOld.x * factors.x,
                        spacingOld.y * factors.y,
                        spacingOld.z * factors.z);
  }

  //
  // return the sum of voxelwise squared intensity difference between
  // two images.  a difference is calculated for each voxel of image1.
  //
  // bcd 2004
  //
  template <class T>
  static
  double
  squaredError(const Image<T>& image1,
               const Image<T>& image2)
  {
    Vector3D<unsigned int> size = image1.getSize();
    double squaredDifference = 0;

    T val1, val2;
    double d;
    Vector3D<double> scale(image1.getSpacing() / image2.getSpacing());
    Vector3D<double> offset((image1.getOrigin() - image2.getOrigin())
                            / image2.getSpacing());

    for (unsigned int z = 0; z < size.z; ++z) {
      for (unsigned int y = 0; y < size.y; ++y) {
	for (unsigned int x = 0; x < size.x; ++x) {
          val1 = image1(x,y,z);
          val2 = Array3DUtils::trilerp(image2,
                                       x * scale.x + offset.x,
                                       y * scale.y + offset.y,
                                       z * scale.z + offset.z,
                                       (T)0.0F);
	  d = ((double) val1) - ((double) val2); 
	  squaredDifference += d * d;
	}
      }
    }  
    return squaredDifference;    
  }

  // Centroid of the image, in world coordinates.
  template <class T>
  static
  Vector3D<double> computeCentroid(const Image<T>& image,
                                   const ROI<int, unsigned int> voxelROI)
  {
#if 0
    Vector3D<double> centroid = 
      Array3DUtils::computeCentroid( image, voxelROI.getStart(),
                                     voxelROI.getSpacing() );
    centroid *= image.getSpacing();
    centroid += image.getOrigin();
    return centroid;
#endif
    std::cerr << "[ImageUtils.h] INTERNAL ERROR: no definition of spacing "
	      << " in ROI class." << std::endl;
    Vector3D<double> centroid(0,0,0);
    return centroid;
  }

  template <class T>

  static
  Vector3D<double> computeCentroid(const Image<T>& image)
  {
    Vector3D<double> centroid = 
      Array3DUtils::computeCentroid( image, typename Array3D<T>::IndexType(0,0,0),
                                     image.getSize() );
    centroid *= image.getSpacing();
    centroid += image.getOrigin();
    return centroid;
  }

  //
  // write the META header file and data file for this image .mhd and
  // .raw are automatically appended to the filenamePrefix
  //
  // bcd 2004
  //
  template <class T>
  static
  void writeMETA(const Image<T>& image,
                 const char* filenamePrefix)
  {
    Array3DIO::writeMETAVolume(image, 
                               image.getOrigin(),
                               image.getSpacing(),
                               filenamePrefix);
  }  

  //
  // extract an roi (specified in voxel coordinates) from an image
  // 
  // bcd 2004
  //
  template <class T>
  static
  void extractROIVoxelCoordinates(const Image<T>& image,
                                  Image<T>& roiImage,
                                  const ROI<int, unsigned int>& voxelROI)
  {
    Array3DUtils::extractROI(image, roiImage, voxelROI);
    roiImage.setSpacing(image.getSpacing());
    roiImage.setOrigin(image.getOrigin().x
                       + voxelROI.getStart().x * image.getSpacing().x,
                       image.getOrigin().y
                       + voxelROI.getStart().y * image.getSpacing().y,
                       image.getOrigin().z
                       + voxelROI.getStart().z * image.getSpacing().z);
  }

  //
  // extract an roi (specified in world coordinates) from an image 
  // 
  // bcd 2004
  //
  template <class T>
  static
  void extractROIWorldCoordinates(const Image<T>& image,
                                  Image<T>& roiImage,
                                  const ROI<double, double>& worldROI)
  {
    Vector3D<double> newStart = 
      (worldROI.getStart() - image.getOrigin()) / image.getSpacing();
    Vector3D<double> newSize = worldROI.getSize() / image.getSpacing();
    // !! this should round instead of cast ??
    ROI<int, unsigned int> voxelROI((Vector3D<int>)newStart, 
                                    (Vector3D<unsigned int>)newSize );
    extractROIVoxelCoords(image, roiImage, voxelROI);
  }


  // This is a specialized function to handle CT images. These images
  // are often stored with a minimum pixel value of -1024 (consistent
  // with the definition of Hounsfield units).  Other times they are
  // shifted so that the minimum intensity is 0, which is what is
  // expected by ImMap and BeamLock.  This function performs that
  // shift.  (Note that PLUNC uses the other convention.)
  template <class T>
  static void
  makeImageUnsigned(Image<T>& image)
  {
    T min, max;
    Array3DUtils::getMinMax(image, min, max);
    if (min < 0) 
    {
      for (unsigned int i = 0; i < image.getNumElements(); ++i)
      {
        image(i) += 1024;
      }
    }
  }

};

#endif
