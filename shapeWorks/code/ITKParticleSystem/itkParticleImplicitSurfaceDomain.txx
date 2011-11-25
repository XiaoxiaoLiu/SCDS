/*=========================================================================
  Program:   ShapeWorks: Particle-based Shape Correspondence & Visualization
  Module:    $RCSfile: itkParticleImplicitSurfaceDomain.txx,v $
  Date:      $Date: 2009/05/06 21:49:15 $
  Version:   $Revision: 1.1.1.1 $
  Author:    $Author: cates $

  Copyright (c) 2009 Scientific Computing and Imaging Institute.
  See ShapeWorksLicense.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.
=========================================================================*/
#ifndef __itkParticleImplicitSurfaceDomain_txx
#define __itkParticleImplicitSurfaceDomain_txx

#include "vnl/vnl_math.h"
#include "vnl/vnl_cross.h"
#define PARTICLE_DEBUG_FLAG 0

namespace itk
{

template<class T, unsigned int VDimension>
void
ParticleImplicitSurfaceDomain<T, VDimension>::
SetCuttingPlane(const vnl_vector<double> &a, const vnl_vector<double> &b,
                const vnl_vector<double> &c)
{
  // See http://mathworld.wolfram.com/Plane.html, for example
  vnl_vector<double> q;
  if (VDimension == 3)  q = vnl_cross_3d((b-a),(c-a));
  else if (VDimension == 2)  q = vnl_cross_2d((b-a),(c-a));
  
  m_CuttingPlaneNormal = q / q.magnitude();
  m_CuttingPlanePoint = a;
  m_UseCuttingPlane = true;
}

template<class T, unsigned int VDimension>
bool
ParticleImplicitSurfaceDomain<T, VDimension>::
ApplyVectorConstraints(vnl_vector_fixed<double, VDimension> &gradE,
                       const PointType &pos,
                       double maxtimestep) const
{

  // NOTE --- DISABLED
  return Superclass::ApplyVectorConstraints(gradE,pos,maxtimestep);
  // END --- DISABLED
  
  if (this->m_UseCuttingPlane == true)
    {    
    // See http://mathworld.wolfram.com/Point-PlaneDistance.html, for example
    vnl_vector_fixed<double, 3> x;
    vnl_vector_fixed<T, VDimension> grad = this->SampleGradientVnl(pos);
    for (unsigned int i = 0; i < VDimension; i++)
      { x[i] = pos[i]; }
    const double D = dot_product(m_CuttingPlaneNormal, x- m_CuttingPlanePoint);
    
    //    x = m_CuttingPlaneNormal * fabs(1.0 / (D + 1.0e-3));

    // x = m_CuttingPlaneNormal * lambda * exp(-lambda * fabs(D));
    
    // Gradient of simple force 1/D^2 to push away from cutting plane
    vnl_vector_fixed<double, VDimension> df;
    const double k = (-2.0 / (D * D * D));
    df[0] = k * grad[0] * m_CuttingPlaneNormal[0];
    df[1] = k * grad[1] * m_CuttingPlaneNormal[1];
    df[2] = k * grad[2] * m_CuttingPlaneNormal[2];

    gradE = gradE + df;

    // Make sure force is not huge relative to other forces.
    if (gradE.magnitude() > maxtimestep)
      {
      gradE.normalize();
      gradE = gradE * maxtimestep;
      }    
    }
  
  return Superclass::ApplyVectorConstraints(gradE,pos,maxtimestep);
}



template<class T, unsigned int VDimension>
bool
ParticleImplicitSurfaceDomain<T, VDimension>::ApplyConstraints(PointType &p) const
{
  // First apply and constraints imposed by superclasses.  This will
  // guarantee the point starts in the correct image domain.
  bool flag = Superclass::ApplyConstraints(p);

  if (this->m_ConstraintsEnabled == true)
    {
  
    unsigned int k = 0;
    double mult = 1.0;
    
    const T epsilon = m_Tolerance * 0.001;
    T f = this->Sample(p);
    
    T gradmag = 1.0;
    while ( fabs(f) > (m_Tolerance * mult) || gradmag < epsilon)
      //  while ( fabs(f) > m_Tolerance || gradmag < epsilon)
      {
      vnl_vector_fixed<T, VDimension> grad = this->SampleGradientVnl(p);
      
      gradmag = grad.magnitude();
      vnl_vector_fixed<T, VDimension> vec   =  grad  * ( f / (gradmag + epsilon) );
      for (unsigned int i = 0; i < VDimension; i++)
        {
        p[i] -= vec[i];
        }
      
#ifdef  PARTICLE_DEBUG_FLAG
      if ( ! this->IsInsideBuffer(p) )
        {
      itkExceptionMacro("A Point, " << p << ", was projected outside the given image domain." );
        }
#endif  
      
      f = this->Sample(p);
      
#ifdef  PARTICLE_DEBUG_FLAG
      if ( gradmag < epsilon && fabs(f) > m_Tolerance)
        {
        itkExceptionMacro("Newton-Raphson iteration failed to find the zero level-set.  Gradient is zero, but f = "  <<  f );
        }
#endif
      
      // Raise the tolerance if we have done too many iterations.
      k++;
      if (k > 10000)
        {
        mult *= 2.0;
        k = 0;
        }
      } // end while
    } // end if m_ConstraintsEnabled == true

  return flag; 
}

template <class T, unsigned int VDimension>
double
ParticleImplicitSurfaceDomain<T, VDimension>::Distance(const PointType &a, const PointType &b) const
{
  // Check to see if the normals are >= 90 degrees apart.
  if ( dot_product(this->SampleGradientVnl(a), this->SampleGradientVnl(b) ) <= 0.0  )
    {
    return -1.0;
    }
  else  // Return Euclidean distance
    {
    double sum = 0.0;
    for (unsigned int i = 0; i < VDimension; i++)
      {      sum += (b[i]-a[i]) * (b[i]-a[i]);      }
    return sqrt(sum);
    }
}
} // end namespace

#endif
