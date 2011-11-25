/*=========================================================================
  Program:   ShapeWorks: Particle-based Shape Correspondence & Visualization
  Module:    $RCSfile: pca.txx,v $
  Date:      $Date: 2009/05/06 21:49:15 $
  Version:   $Revision: 1.1.1.1 $
  Author:    $Author: cates $

  Copyright (c) 2009 Scientific Computing and Imaging Institute.
  See ShapeWorksLicense.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.
=========================================================================*/
#ifndef __st_pca_txx
#define __st_pca_txx

#include "pca.h"
#include "itkImageRegionIterator.h"

namespace shapetools
{

template <class T, unsigned int D>
pca<T,D>::pca(param::parameterFile &pf)
{
  // Set some parameters.
  bool ok = true;
  PARAMSET(pf, m_background, "background", 0, ok, static_cast<T>(1));

  if (ok == false)
    {  throw param::Exception("pca:: missing parameters"); }
}

template <class T, unsigned int D> 
void pca<T,D>::operator()(typename image_type::Pointer img)
{
  // Compute the sample mean.
  itk::ImageRegionIteratorWithIndex<image_type>
    it(img, img->GetRequestedRegion());
  vnl_vector<double> mean(D);
  mean.fill(0.0);
  double count = 0.0;
  itk::Point<double, D> point;

  for (it.GoToBegin(); ! it.IsAtEnd(); ++it)
    {
    if (it.Get() != m_background)
      {
      // Get the physical index from the image index.
      img->TransformIndexToPhysicalPoint(it.GetIndex(), point);
      for (unsigned int i = 0; i < D; i++) { mean[i] += point[i]; }
      count += 1.0;
      }
    }

  mean /= count;

  // Cache for output.
  m_mean = mean;
  
  // Compute the covariance matrix.  
  vnl_matrix<double> pts_minus_mean(D, static_cast<int>(count));
  int j = 0;
  for (it.GoToBegin(); ! it.IsAtEnd(); ++it)
    {
    if (it.Get() != m_background)
      {
      // Get the physical index from the image index.
      img->TransformIndexToPhysicalPoint(it.GetIndex(), point);
      for (unsigned int i = 0; i < D; i++)
        { pts_minus_mean(i, j) = point[i] - mean[i]; }
      j++;
      }
    }
  vnl_matrix<double> covar = (pts_minus_mean * pts_minus_mean.transpose())
    * (1.0/(count - 1.0));;

  // Compute the eigenvalues & eigenvectors of the covariance matrix.
  vnl_symmetric_eigensystem<double> symEigen(covar);   
  m_eigenvectors = symEigen.V;
  m_eigenvalues  = symEigen.D;
}

} // end namespace

#endif
