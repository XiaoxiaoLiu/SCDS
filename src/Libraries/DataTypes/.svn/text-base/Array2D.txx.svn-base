#ifndef ARRAY2D_TXX
#define ARRAY2D_TXX

#include <iostream>
#include <fstream>
#include "Array2D.h"
//////////////////
// constructor  //
//////////////////

template <class T>
Array2D<T>
::Array2D()
  : _size(0, 0),
    _xySize(0),
    _rowPtrs(0),
    _dataPtr(0)
{}

////////////////////////////////////
// constructor : two parameters //
////////////////////////////////////

template <class T>
Array2D<T>
::Array2D(unsigned int xSize,
	      unsigned int ySize)
  : _size(xSize, ySize),
    _xySize(xSize * ySize),
    _rowPtrs(0),
    _dataPtr(0)
{
  _allocateData();
}

////////////////////////////
// constructor : vector2D //
////////////////////////////

template <class T>
Array2D<T>
::Array2D(const SizeType& size)
{
  Array2D(size.x, size.y);
}

//////////////////////
// copy constructor //
//////////////////////


template <class T>
Array2D<T>
::Array2D(const Array2D<T>& rhs)
  : _size(rhs._size),
    _xySize(rhs._xySize),
    _rowPtrs(0),
    _dataPtr(0)
{
  _allocateData();
  memcpy(_dataPtr, rhs._dataPtr, _xySize * sizeof(T));
}

////////////////
// destructor //
////////////////

template <class T>
Array2D<T>
::~Array2D()
{
  _deallocateData();
}

///////////////
// operator= //
///////////////

template <class T>
Array2D<T>&
Array2D<T>
::operator=(const Array2D<T>& rhs)
{
  if (this == &rhs) return *this;

  if (_size != rhs._size)
    {
      _size    = rhs._size;
      _xySize  = rhs._xySize;
      _deallocateData();
      _allocateData();
    }
  memcpy(_dataPtr, rhs._dataPtr, _xySize * sizeof(T));
  return *this;
}

/////////////
// getSize //
/////////////

template <class T>
inline
typename Array2D<T>::SizeType
Array2D<T>
::getSize() const 
{ 
  return _size; 
}


//////////////
// getSizeX //
//////////////

template <class T>
inline
unsigned int
Array2D<T>
::getSizeX() const 
{ 
  return _size.x; 
}

//////////////
// getSizeY //
//////////////

template <class T>
inline
unsigned int
Array2D<T>
::getSizeY() const 
{ 
  return _size.y; 
}


/////////////
// isEmpty //
/////////////

template <class T>
inline
bool
Array2D<T>
::isEmpty() const 
{
return ( _rowPtrs== NULL);
}


////////////
// resize //
////////////

template <class T>
inline
void 
Array2D<T>
::resize(const SizeType& size) 
{ 
  resize(size.x, size.y); 
}


////////////
// resize //
////////////

template <class T>
void
Array2D<T>
::resize(unsigned int xSize, unsigned int ySize)
{
  if (_size.x == xSize && _size.y == ySize) return;

  _size.set(xSize, ySize);
  _xySize  = xSize * ySize;
  _deallocateData();
  _allocateData();
}


////////////////////
// getNumElements //
////////////////////

template <class T>
inline
unsigned int 
Array2D<T>
::getNumElements() const 
{ 
  return _xySize; 
}


//////////////////
// getSizeBytes //
//////////////////

template <class T>
inline
unsigned int 
Array2D<T>
::getSizeBytes() const 
{ 
  return _xySize * sizeof(T); 
}

/////////
// set //
/////////

template <class T>
inline
void
Array2D<T>
::set(unsigned int xIndex,unsigned int yIndex, T item)
{
  if(_rowPtrs == NULL){
	  throw std::runtime_error("Data is NULL");

    }
    if(xIndex >= 0 && xIndex < _size.x &&
       yIndex >= 0 && yIndex < _size.y ){
        _rowPtrs[yIndex][xIndex] = item;
    }else{
        throw std::out_of_range("Index is out of range");
    }
 }

/////////
// get //
/////////

template <class T>
inline
T&
Array2D<T>
::get(unsigned int xIndex,unsigned int yIndex)
{
  if(_rowPtrs == NULL){
        throw std::runtime_error("Data is NULL");

    }
    if(xIndex >= 0 && xIndex < _size.x &&
       yIndex >= 0 && yIndex < _size.y){
        return _rowPtrs[yIndex][xIndex];
    }else{
        throw std::out_of_range("Index is out of range");
    }
 }




////////////////
// operator() //
////////////////

template <class T>
inline
T&
Array2D<T>
::operator()(unsigned int xIndex, 
	     unsigned int yIndex)
{
  return _rowPtrs[yIndex][xIndex];
}


//////////////////////
// const operator() //
//////////////////////

template <class T>
inline
const T&
Array2D<T>
::operator()(unsigned int xIndex, 
	     unsigned int yIndex) const
{
  return _rowPtrs[yIndex][xIndex];
}


///////////////////////////
// operator() : Vector2D //
///////////////////////////

template <class T>
inline
T&
Array2D<T>
::operator()(const IndexType& index)
{
  return _rowPtrs[index.y][index.x];
}


////////////////////////////////
// const operator() : vector2D//
////////////////////////////////

template <class T>
inline
const T&
Array2D<T>
::operator()(const IndexType& index) const
{
  return _rowPtrs[index.y][index.x];
}


//////////
// fill //
//////////

template <class T>
inline
T&
Array2D<T>
::operator()(unsigned int elementIndex)
{
  return _dataPtr[elementIndex];
}

template <class T>
inline
const T&
Array2D<T>
::operator()(unsigned int elementIndex) const
{
  return _dataPtr[elementIndex];
}

template <class T>
void
Array2D<T>
::fill(const T& fillValue)
{
  for (unsigned int i = 0; i < _xySize; ++i)
    {
      _dataPtr[i] = fillValue;
    }
}


////////////////
// operator== //
////////////////

template <class T>
void
Array2D<T>
::scale(const double& d)
{
  for (unsigned int i = 0; i < _xySize; ++i)
    {
      _dataPtr[i] *= d;
    }
}

template <class T>
bool
Array2D<T>
::operator==(const Array2D<T>& rhs) const
{
  return dimensionsMatch(rhs) &&
    memcmp(_dataPtr, rhs._dataPtr, _xySize * sizeof(T));
}

////////////////
// operator!= //
////////////////

template <class T>
bool
Array2D<T>
::operator!=(const Array2D<T>& rhs) const
{
  return !(*this == rhs);
}

////////////////////
// getDataPointer //
////////////////////

template <class T>
inline
T* 
Array2D<T>
::getDataPointer(unsigned int xIndex, 
		 unsigned int yIndex)
{
  return &_rowPtrs[yIndex][xIndex];
}

///////////////////////////////
// getDataPointer : vector2D //
///////////////////////////////

template <class T>
inline
T* 
Array2D<T>
::getDataPointer(const IndexType& index)
{
  return &_rowPtrs[index.y][index.x];
}

//////////////////////////
// const getDataPointer //
//////////////////////////

template <class T>
inline
const T* 
Array2D<T>
::getDataPointer(unsigned int xIndex, 
		 unsigned int yIndex) const
{
  return &_rowPtrs[yIndex][xIndex];
}

/////////////////////////////////////
// const getDataPointer : vector2D //
/////////////////////////////////////

template <class T>
inline
const T* 
Array2D<T>
::getDataPointer(const IndexType& index) const
{
  return &_rowPtrs[index.y][index.x];
}

////////////////////////////////
// operator<<() [text output] //
////////////////////////////////

template <class T>
std::ostream& 
operator<<(std::ostream& output, const Array2D<T>& array)
{
  for (unsigned int y = 0; y < array.getSizeY(); ++y) {
    for (unsigned int x = 0; x < array.getSizeX(); ++x) {
      output << array(x, y) << " ";
    }
    output << '\n';
  }
  return output;
}

///////////////////
// _allocateData //
///////////////////

template <class T>
void
Array2D<T>
::_allocateData()
{
  _rowPtrs   = new T  *[_size.y];
  _dataPtr   = new T   [_xySize];

  for (unsigned int rowIndex = 0; rowIndex < _size.y; ++rowIndex) 
    {
      _rowPtrs[rowIndex] = &_dataPtr[rowIndex * _size.x];
    }
}

/////////////////////
// _deallocateData //
/////////////////////

template <class T>
void
Array2D<T>
::_deallocateData()
{
  if (!_rowPtrs)   delete [] _rowPtrs;
  if (!_dataPtr)   delete [] _dataPtr;

  _rowPtrs   = 0;
  _dataPtr   = 0;
}

#endif
