/*=========================================================================
  Program:   ShapeWorks: Particle-based Shape Correspondence & Visualization
  Module:    $RCSfile: itkParticleSystem.cxx,v $
  Date:      $Date: 2009/05/06 21:49:15 $
  Version:   $Revision: 1.1.1.1 $
  Author:    $Author: cates $

  Copyright (c) 2009 Scientific Computing and Imaging Institute.
  See ShapeWorksLicense.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.
=========================================================================*/
#include "itkParticleSystem.h"

namespace itk {

/** For efficiency, we specialize for 3D and 2D */
template<>
ParticleSystem<3>::PointType
ParticleSystem<3>
::TransformPoint(const PointType &p, const TransformType &T) const
{
  PointType ans;

  ans[0] = T[0][0] * p[0] + T[0][1] * p[1] + T[0][2] * p[2] + T[0][3];
  ans[1] = T[1][0] * p[0] + T[1][1] * p[1] + T[1][2] * p[2] + T[1][3];
  ans[2] = T[2][0] * p[0] + T[2][1] * p[1] + T[2][2] * p[2] + T[2][3];

  return ans;
}

/** For efficiency, we specialize for 3D and 2D */
template<>
ParticleSystem<2>::PointType
ParticleSystem<2>
::TransformPoint(const PointType &p, const TransformType &T) const
{
  PointType ans;

  ans[0] = T[0][0] * p[0] + T[0][1] * p[1] + T[0][2];
  ans[1] = T[1][0] * p[0] + T[1][1] * p[1] + T[1][2];

  return ans;
}


/** For efficiency, we specialize for 3D and 2D */
template<>
ParticleSystem<3>::VectorType
ParticleSystem<3>
::TransformVector(const VectorType &p, const TransformType &T) const
{
  VectorType ans;

  ans[0] = T[0][0] * p[0] + T[0][1] * p[1] + T[0][2] * p[2];
  ans[1] = T[1][0] * p[0] + T[1][1] * p[1] + T[1][2] * p[2];
  ans[2] = T[2][0] * p[0] + T[2][1] * p[1] + T[2][2] * p[2];

  return ans;
}

/** For efficiency, we specialize for 3D and 2D */
template<>
ParticleSystem<2>::VectorType
ParticleSystem<2>
::TransformVector(const VectorType &p, const TransformType &T) const
{
  VectorType ans;

  ans[0] = T[0][0] * p[0] + T[0][1] * p[1];
  ans[1] = T[1][0] * p[0] + T[1][1] * p[1];

  return ans;
}

} // end namespace
