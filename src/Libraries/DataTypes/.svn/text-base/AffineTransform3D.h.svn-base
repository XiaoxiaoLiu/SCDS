// Mark Foskey 5/04

#ifndef AFFINETRANSFORM_H
#define AFFINETRANSFORM_H

#include <stdexcept>

#include "Matrix3D.h"
#include "Vector3D.h"

template< class T >
class AffineTransform3D
{
public:

  typedef T CoordinateType;
  typedef AffineTransform3D< CoordinateType > Self;
  typedef Matrix3D< CoordinateType > MatrixType;
  typedef Vector3D< CoordinateType > VectorType; 

  MatrixType matrix;
  VectorType vector;

  AffineTransform3D() : matrix(), vector() {} // Identity.
  AffineTransform3D(const MatrixType& m, const VectorType& v ) 
    : matrix(m), vector(v) {}
  explicit AffineTransform3D( const MatrixType& m ) : matrix(m), vector() {}
  explicit AffineTransform3D( const VectorType& v ) : matrix(), vector(v) {}

  template< class U >
  AffineTransform3D( const AffineTransform3D<U>& t )
    : matrix( t.matrix ), vector( t.vector ) {}

  template< class U >
  Self& operator=( const AffineTransform3D<U>& rhs )
  { 
    this->matrix = rhs.matrix; 
    vector = rhs.vector; 
    return *this; 
  }

  void eye() { *this = AffineTransform3D<T>(); } // Sets to identity
  bool invert();
  Self& operator*=( const Self& rhs );
    
  Self operator*( const Self& other ) const;
  VectorType operator*( const VectorType& v ) const;
  void transformVector( const VectorType& in, VectorType& out ) const;
  void transformVector( VectorType& v ) const;

  template <class U, class V>
  void transformCoordinates( const U& xIn, const U& yIn, const U& zIn, 
                             V& xOut, V& yOut, V& zOut ) const
  {
    xOut =(V)(xIn*matrix.a[0] + yIn*matrix.a[1] + zIn*matrix.a[2] + vector[0]);
    yOut =(V)(xIn*matrix.a[3] + yIn*matrix.a[4] + zIn*matrix.a[5] + vector[1]);
    zOut =(V)(xIn*matrix.a[6] + yIn*matrix.a[7] + zIn*matrix.a[8] + vector[2]);
  }

  bool operator==( const Self t ) const;

  // No other operators; e.g., '+' isn't really a sensible op.

  void writePLUNCStyle(const char* filename) const;
  void writePLUNCStyle(const std::string& filename) const
    { writePLUNCStyle(filename.c_str()); }
  void readPLUNCStyle(const char* filename);
  void readPLUNCStyle(const std::string& filename)
    { readPLUNCStyle(filename.c_str()); }

private:

  AffineTransform3D( const Self& first, const Self& second );

};

template <typename T>
std::ostream& 
operator<<(std::ostream& output, const AffineTransform3D<T>& t)
{
  return output << '\n' << t.matrix << "\n\n" << t.vector[0] << " "
                << t.vector[1] << " " << t.vector[2];
}


#include "AffineTransform3D.txx"

#endif  // ndef AFFINETRANSFORM_H
