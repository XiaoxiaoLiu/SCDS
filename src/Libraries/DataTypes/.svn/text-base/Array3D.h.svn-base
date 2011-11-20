#ifndef ARRAY3D_H
#define ARRAY3D_H

#include "Vector3D.h"
#include <iosfwd>
#include <stdexcept>

// NOTE: The data layout is such that the X coordinate is the one that
// changes most rapidly, and Z the least.  So, when iterating via 3
// nested loops, the innermost loop should be controlled by the
// changing X index.

template <class T>
class Array3D
{
public:
  typedef Vector3D<unsigned int>   IndexType;
  typedef Vector3D<unsigned int>   SizeType;

  Array3D();
  Array3D(unsigned int xSize,
	  unsigned int ySize, 
	  unsigned int zSize);

  Array3D(const SizeType& size);
  Array3D(const Array3D<T>& rhs);
  virtual ~Array3D();
  
  Array3D<T>& operator=(const Array3D<T>& rhs);

  const SizeType& getSize() const;
  unsigned int getSizeX() const;
  unsigned int getSizeY() const;
  unsigned int getSizeZ() const;

  bool isEmpty() const;

  // does not maintain data
  void resize(const SizeType& size);
  void resize(unsigned int xSize, 
	      unsigned int ySize,
	      unsigned int zSize);
  
 
  unsigned int getNumElements() const;
  unsigned int getSizeBytes() const;

  void set(unsigned int xIndex,
	   unsigned int yIndex,
	   unsigned int zIndex, 
	   const T& item);
  T& get(unsigned int xIndex,
	 unsigned int yIndex,
	 unsigned int zIndex);

  T& operator()(const IndexType& index);
  T& operator()(unsigned int xIndex, 
		unsigned int yIndex,
		unsigned int zIndex);
  T& operator()(unsigned int elementIndex);

  bool isValidIndex(const IndexType& index);

  const T& operator()(const IndexType& index) const;
  const T& operator()(unsigned int xIndex, 
		      unsigned int yIndex,
		      unsigned int zIndex) const;
  const T& operator()(unsigned int elementIndex) const;

  void fill(const T& fillValue);

  // Multiply all values by d
  void scale(const double& d);
  void addScalar(const double& scalar);

  void pointwiseMultiplyBy(const Array3D<T>& rhs);

  bool operator==(const Array3D<T>& rhs) const;
  bool operator!=(const Array3D<T>& rhs) const;

  T* getDataPointer(const IndexType& index);
  T* getDataPointer(unsigned int xIndex = 0, 
		    unsigned int yIndex = 0, 
		    unsigned int zIndex = 0);

  const T* getDataPointer(const IndexType& index) const;
  const T* getDataPointer(unsigned int xIndex = 0, 
			  unsigned int yIndex = 0, 
			  unsigned int zIndex = 0) const;
  
  void setData(const Array3D<T>& rhs);
  void copyData(const void* const dataPtr);

protected:
  SizeType     _size;
  unsigned int _xySize;
  unsigned int _yzSize;
  unsigned int _zxSize;
  unsigned int _xyzSize;
  T ***_slicePtrs;
  T  **_rowPtrs;
  T   *_dataPtr;

  void _allocateData();
  void _deallocateData();
};

#include "Array3D.txx"

#endif
