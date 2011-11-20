#ifndef AFFINETRANSFORM_TXX
#define AFFINETRANSFORM_TXX

// Included at end of AffineTransform3D.h

#include <exception>
#include <fstream>
#include <iostream>
#include <stdexcept>

//------------------------------------------------------------------------
// Public methods

template<class T>
AffineTransform3D<T>& AffineTransform3D<T>::
operator*=( const AffineTransform3D<T>& rhs )
{
    matrix *= rhs.matrix;
    vector = *this * rhs.vector;
    return *this;
}

// If w = Mv + T, then M^-1 w - M^-1 T = v.
template<class T> 
bool AffineTransform3D<T>::
invert()
{
    if (!matrix.invert()) return false;
    vector *= static_cast<T>(-1);
    vector = matrix * vector;
    return true;
}

// Composes two trasformations.  (A*B)*v = A*(B*v).
template<class T> 
AffineTransform3D<T> AffineTransform3D<T>::
operator*( const AffineTransform3D<T>& other ) const
{ return AffineTransform3D<T>(*this, other); }

// Note that v*A is undefined.
template<class T> 
Vector3D<T>
AffineTransform3D<T>::operator*( const Vector3D<T>& v ) const
{
    Vector3D<T> result;
    transformVector( v, result );
    return result;
}

template<class T> 
inline void
AffineTransform3D<T>::transformVector( Vector3D<T>& v ) const
{
    Vector3D<T> in = v;
    transformVector( in, v );
}

// Apply the transformation (*this) to 'in', yielding 'out'.  Avoids
// copying the data structure.
template<class T> 
inline void
AffineTransform3D<T>::transformVector( const Vector3D<T>& in,
                                       Vector3D<T>& out ) const
{
    transformCoordinates( in[0], in[1], in[2], out[0], out[1], out[2] );
}

// template<class T> 
// inline void AffineTransform3D<T>::
// transformCoordinates( const T& xIn, const T& yIn, const T& zIn, 
//                       T& xOut, T& yOut, T& zOut ) const 
// {
//     xOut = xIn*matrix.a[0] + yIn*matrix.a[1] + zIn*matrix.a[2] + vector[0];
//     yOut = xIn*matrix.a[3] + yIn*matrix.a[4] + zIn*matrix.a[5] + vector[1];
//     zOut = xIn*matrix.a[6] + yIn*matrix.a[7] + zIn*matrix.a[8] + vector[2];
// }

template<class T> bool AffineTransform3D<T>::
operator==( const AffineTransform3D<T> t ) const
{
    return ( matrix == t.matrix && vector == t.vector );
}

template <class T>
void
AffineTransform3D<T>
::writePLUNCStyle(const char* filename) const
{
    std::ofstream output(filename);
    if (output.fail() || output.bad())
    {
        throw std::runtime_error("Could not open file.");
    }
    output << "matrix" << std::endl;
    output << matrix(0) << " " << matrix(3) << " " << matrix(6) << " " << 0 
           << std::endl;
    output << matrix(1) << " " << matrix(4) << " " << matrix(7) << " " << 0 
           << std::endl;
    output << matrix(2) << " " << matrix(5) << " " << matrix(8) << " " << 0 
           << std::endl;
    output << vector.x  << " " << vector.y  << " " << vector.z  << " " << 1
           << std::endl;
  
    if (output.bad())
    {
        throw std::runtime_error("Could not write file.");
    }

    output.close();
}


template <class T>
void
AffineTransform3D<T>
::readPLUNCStyle(const char* filename)
{
    std::ifstream input(filename);
    if (input.fail() || input.bad())
    {
        throw std::runtime_error("Could not open file.");
    }

    std::string key;
    std::getline(input, key, '\n');
    key.erase(key.begin()+6, key.end());
    if (key != ("matrix"))
    {
        throw std::runtime_error("Invalid matrix file format.");
    }

    T dummy1, dummy2, dummy3, dummy4;

    input >> matrix(0) >> matrix(3) >> matrix(6) >> dummy1;
    input >> matrix(1) >> matrix(4) >> matrix(7) >> dummy2;
    input >> matrix(2) >> matrix(5) >> matrix(8) >> dummy3;
    input >> vector.x  >> vector.y  >> vector.z  >> dummy4;
  
    if (input.fail())
    {
        throw std::runtime_error("Could not read file.");
    }

    input.close();
}



//------------------------------------------------------------------------
// Private methods

// Constructs by composing two others.  Makes it easier for the *
// operator to optimize away the copy when the value is returned (I
// think).  Not really appropriate for user code; just use the *
// operator.
template<class T> AffineTransform3D<T>::
AffineTransform3D( const AffineTransform3D<T>& first,
		   const AffineTransform3D<T>& second )
  : matrix( first.matrix * second.matrix ), vector( first * second.vector ) {}

#endif  // ndef AFFINETRANSFORM_TXX
