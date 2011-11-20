#ifndef ARRAY2D_H
#define ARRAY2D_H

#include "Vector2D.h"
#include <iosfwd>
#include <stdexcept>

template <class T>
class Array2D
{
public:

  typedef Vector2D<unsigned int>   IndexType;
  typedef Vector2D<unsigned int>   SizeType;


  Array2D();
  Array2D(unsigned int xSize,
		  unsigned int ySize);

  Array2D(const SizeType& size);
  Array2D(const Array2D<T>& rhs);
  ~Array2D();
  
  Array2D<T>& operator=(const Array2D<T>& rhs);

  SizeType getSize() const;
  unsigned int getSizeX() const;
  unsigned int getSizeY() const;

  bool isEmpty() const;

  // does not maintain data
  void resize(const SizeType& size);
  void resize(unsigned int xSize, 
			  unsigned int ySize);

  unsigned int getNumElements() const;
  unsigned int getSizeBytes() const;

  void set(unsigned int xIndex,
	       unsigned int yIndex, 
		   T item);
  T& get(unsigned int xIndex,
	     unsigned int yIndex);

  T& operator()(const IndexType& index);
  T& operator()(unsigned int xIndex, 
				unsigned int yIndex);

  T& operator()(unsigned int elementIndex);


  const T& operator()(const IndexType& index) const;
  const T& operator()(unsigned int xIndex, 
					  unsigned int yIndex) const;

  const T& operator()(unsigned int elementIndex) const;


  void fill(const T& fillValue);
  void scale(const double& d);
  void outputPNG(const std::string outputFileName);
  void outputPNGcolor(const std::string outputFileName);


  bool operator==(const Array2D<T>& rhs) const;
  bool operator!=(const Array2D<T>& rhs) const;

  T* getDataPointer(const IndexType& index);
  T* getDataPointer(unsigned int xIndex = 0, 
					unsigned int yIndex = 0);

  const T* getDataPointer(const IndexType& index) const;
  const T* getDataPointer(unsigned int xIndex = 0, 
						  unsigned int yIndex = 0) const;

protected:
  SizeType _size;
  unsigned int _xySize;
  T  **_rowPtrs;
  T   *_dataPtr;

  void _allocateData();
  void _deallocateData();
};

template <class T>
std::ostream& 
operator<<(std::ostream& output, const Array2D<T>& array);

#include "Array2D.txx"

#endif
