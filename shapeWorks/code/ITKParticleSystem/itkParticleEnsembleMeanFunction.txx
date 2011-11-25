/*=========================================================================
  Program:   ShapeWorks: Particle-based Shape Correspondence & Visualization
  Module:    $RCSfile: itkParticleEnsembleMeanFunction.txx,v $
  Date:      $Date: 2009/05/06 21:49:15 $
  Version:   $Revision: 1.1.1.1 $
  Author:    $Author: cates $

  Copyright (c) 2009 Scientific Computing and Imaging Institute.
  See ShapeWorksLicense.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.
=========================================================================*/
#ifndef __itkParticleEnsembleMeanFunction_txx
#define __itkParticleEnsembleMeanFunction_txx

#include "vnl/vnl_vector_fixed.h"
#include "itkParticleImageDomainWithGradients.h"

namespace itk
{

template <unsigned int VDimension>
typename ParticleEnsembleMeanFunction<VDimension>::VectorType
ParticleEnsembleMeanFunction<VDimension>
::Evaluate(unsigned int idx, unsigned int d, const ParticleSystemType * system,
           double &maxmove, double &energy) const
{
  const ParticleImageDomainWithGradients<float, VDimension> * domain
    = static_cast<const ParticleImageDomainWithGradients<float, VDimension> *>(system->GetDomain(d));
  
  const double epsilon = 1.0e-8;
  const double N = (double)(system->GetNumberOfDomains() / m_DomainsPerShape);

  // Get the position for which we are computing the gradient.
  PointType pos = system->GetTransformedPosition(idx, d);

  // Find the mean position for the ENSEMBLE neighborhood.
  VectorType gradE;
  for (unsigned int n = 0; n < VDimension; n++)
    {
    gradE[n] = pos[n];
    }
  for (unsigned int i = d % m_DomainsPerShape; i < system->GetNumberOfDomains();
       i += m_DomainsPerShape)
    {
    if (i != d)
      {
      PointType neighpos = system->GetTransformedPosition(idx, i);

      for (unsigned int n = 0; n < VDimension; n++)
        {
        gradE[n] += neighpos[n];
        }
      }
    }
  for (unsigned int n = 0; n < VDimension; n++)
    {
    gradE[n] = pos[n] - (gradE[n] / N);
    }
  

  energy = gradE.magnitude();
  maxmove = domain->GetImage()->GetSpacing()[0];
  //  maxmove = energy * 0.5;
  
  //  Transform the gradient according to the transform of the given domain and
  //  return.
  return system->TransformVector(gradE, system->GetInversePrefixTransform(d)
                                 * system->GetInverseTransform(d));
}



} // end namespace
#endif
