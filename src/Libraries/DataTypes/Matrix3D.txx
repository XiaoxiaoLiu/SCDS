#ifndef MATRIX3D_TXX
#define MATRIX3D_TXX

#include <cstring>
#include <iostream>
#include <algorithm>

template <typename T>
Matrix3D<T>
::Matrix3D() 
{
  eye();
}

template <typename T>
inline T&
Matrix3D<T>
::operator()(unsigned int elementIndex)
{
  return a[elementIndex];
}

template <typename T>
inline const T&
Matrix3D<T>
::operator()(unsigned int elementIndex) const
{
  return a[elementIndex];
}

template <typename T>
inline T&
Matrix3D<T>
::operator()(unsigned int rowIndex,
	     unsigned int columnIndex)
{
  return a[rowIndex * 3 + columnIndex];
}

template <typename T>
inline const T&
Matrix3D<T>
::operator()(unsigned int rowIndex,
	     unsigned int columnIndex) const
{
  return a[rowIndex * 3 + columnIndex];
}

template <typename T>
void 
Matrix3D<T>
::eye()
{
  a[0] = a[4] = a[8] = 1;
  a[1] = a[2] = a[3] = a[5] = a[6] = a[7] = 0;
}

template <typename T>
std::ostream& 
Matrix3D<T>
::writeASCII(std::ostream& output) const
{
  output 
    << a[0] << " " << a[1] << " " << a[2] << '\n'  // Don't want to use endl
    << a[3] << " " << a[4] << " " << a[5] << '\n'  // because that flushes the
    << a[6] << " " << a[7] << " " << a[8];         // buffer.
  return output;
}

template <typename T>
std::ostream& 
operator<<(std::ostream& output, const Matrix3D<T>& matrix)
{
  return matrix.writeASCII(output);
}

//
// compute the inverse of the matrix a, store the result in ainv.
// true is returned if the inverse computation was successful, false
// if there was a problem (e.g. a is singular).  a and ainv may point
// to the same memory, both must be initialized.  If false is
// returned, ainv is unchanged.
//
// uses the clapack functions dgetrf_ and dgetri_ to compute the LU
// decomposition and then the inverse.
//
// bcd 2004
//
template <typename T>
bool
Matrix3D<T>
::computeInverse(const T* const a, T* const ainv)
{
  // setup dgetrf parameters
  const integer M = 3;     // in: num rows
  const integer N = 3;     // in: num columns
  double A[9];   // in: the matrix, out: L and U w/o unit diagonal of L  
  integer LDA = 3;   // in: the leading dimension of the array
  integer IPIV[M*N]; // out: pivot indices
  integer INFO;      // out: 0 success, <0 arg problem, >0 singular matrix

  // copy a' into A
  _transposeMatrixData(a, A);

  // call dgetrf for LU decomposition
  dgetrf_(&M, &N, A, &LDA, IPIV, &INFO);
  if (INFO != 0) return false;

  // setup dgetri parameters
  const integer LWORK = 192;    // experimentally the optimal size for 3x3
  double WORK[LWORK];           // scratch room  

  // call dgetri to get the inverse
  dgetri_(&M, A, &LDA, IPIV, WORK, &LWORK, &INFO);
  if (INFO != 0) return false;

  // now take the transpose and we're done
  _transposeMatrixData(A, ainv);

  return true;
}

//  Compute the singular value decomposition (SVD) of the matrix
//  determined by the coordinates in a.  In the 3-by-3 case, the SVD
//  is written
//
//       A = U * SIGMA * transpose(V)
//
//  where SIGMA is a 3-by-3 matrix which is zero except for its
//  diagonal elements, and U and V are 3-by-3 orthogonal matrices.
//  The diagonal elements of SIGMA are the singular values of A; they
//  are real and non-negative, and are returned in descending order.
//  The columns of U and V are the left and right singular vectors of
//  A.  The effect is to represent A as a rotation, followed by a
//  possibly non-uniform scaling, followed by another rotation.
//
//  Note that the routine returns the transpose of V, not V.  The
//  results are stored in row-major order (the order used by
//  Matrix3D).
template <typename T>
bool
Matrix3D<T>
::computeSVD(const T* const A, T* const U, T* const Vt, T* const sv)
{

  // Input to svd method (get all rows/cols of u,vt in their own mem
  // locations).  
  char jobu = 'A';                // in: return all of U in own space
  char jobvt = 'A';               // in: return all of Vt in own space
  const integer m = 3;            // in: num rows
  const integer n = 3;            // in: num columns
  double a[9];                    // in: input matrix
  const integer lda = n;          // in: Leading dimension of a
  double s[3];                    // out: singular vals, largest to smallest
  double u[9];                    // out: U array
  const integer ldu = m;          // in: leading dimension of U
  double vt[9];                   // out: Vt array
  const integer ldvt = n;         // in: leading dimension of Vt
  const integer lwork = 6*(m+n);  // in: size of workspace
  double work[lwork];             // workspace array
  integer info;                   // out: 0 success, <0 bad args, >0 other

  // copy a' into A
  _transposeMatrixData(A, a);

  // Call svd method
  dgesvd_(&jobu,&jobvt,&m,&n,a,&lda,s,u,&ldu,vt,&ldvt,work,&lwork,&info);

  // Stop here if it didn't work
  if (info != 0) return false;

  // transpose matrix results into output
  _transposeMatrixData(u, U);
  _transposeMatrixData(vt, Vt);

  // Copy singular values
  sv[0] = static_cast<T>(s[0]);
  sv[1] = static_cast<T>(s[1]);
  sv[2] = static_cast<T>(s[2]);

  return true;

}

// If matrix is singular, returns false and leaves matrix unchanged.
template <typename T>
bool
Matrix3D<T>
::invert()
{
  return Matrix3D<T>::computeInverse(this->a, this->a);
}

template <typename T, typename U>
Vector3D<T> 
operator*(const Matrix3D<T>& m, const Vector3D<U>& v)
{
  return Vector3D<T>(v.x * m.a[0] +
		     v.y * m.a[1] +
		     v.z * m.a[2],
		     v.x * m.a[3] +
		     v.y * m.a[4] +
		     v.z * m.a[5],
		     v.x * m.a[6] +
		     v.y * m.a[7] +
		     v.z * m.a[8]);
}

template <typename T, typename U>
Matrix3D<T>& 
operator*=(Matrix3D<T>& lhs, const Matrix3D<U>& rhs)
{
  double tmp1, tmp2;
  tmp1     = lhs.a[0] * rhs.a[0] + lhs.a[1] * rhs.a[3] + lhs.a[2] * rhs.a[6];
  tmp2     = lhs.a[0] * rhs.a[1] + lhs.a[1] * rhs.a[4] + lhs.a[2] * rhs.a[7];
  lhs.a[2] = lhs.a[0] * rhs.a[2] + lhs.a[1] * rhs.a[5] + lhs.a[2] * rhs.a[8];
  lhs.a[0] = tmp1;
  lhs.a[1] = tmp2;

  tmp1     = lhs.a[3] * rhs.a[0] + lhs.a[4] * rhs.a[3] + lhs.a[5] * rhs.a[6];
  tmp2     = lhs.a[3] * rhs.a[1] + lhs.a[4] * rhs.a[4] + lhs.a[5] * rhs.a[7];
  lhs.a[5] = lhs.a[3] * rhs.a[2] + lhs.a[4] * rhs.a[5] + lhs.a[5] * rhs.a[8];
  lhs.a[3] = tmp1;
  lhs.a[4] = tmp2;

  tmp1     = lhs.a[6] * rhs.a[0] + lhs.a[7] * rhs.a[3] + lhs.a[8] * rhs.a[6];
  tmp2     = lhs.a[6] * rhs.a[1] + lhs.a[7] * rhs.a[4] + lhs.a[8] * rhs.a[7];
  lhs.a[8] = lhs.a[6] * rhs.a[2] + lhs.a[7] * rhs.a[5] + lhs.a[8] * rhs.a[8];
  lhs.a[6] = tmp1;
  lhs.a[7] = tmp2;

  return lhs;
}

template <typename T, typename U>
Matrix3D<T>
operator*(const Matrix3D<T>& lhs, const Matrix3D<U>& rhs)
{
  Matrix3D<T> tmp(lhs);
  tmp *= rhs;
  return tmp;
}

#endif
