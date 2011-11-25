/*=========================================================================
  Program:   ShapeWorks: Particle-based Shape Correspondence & Visualization
  Module:    $RCSfile: itkParticlePointIndexPair.h,v $
  Date:      $Date: 2009/05/06 21:49:15 $
  Version:   $Revision: 1.1.1.1 $
  Author:    $Author: cates $

  Copyright (c) 2009 Scientific Computing and Imaging Institute.
  See ShapeWorksLicense.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.
=========================================================================*/
#ifndef __itkParticlePointIndexPair_h
#define __itkParticlePointIndexPair_h

#include "itkPoint.h"
namespace itk {
/** 
    Struct containing a Point and an index value associated with a point.  This
    object is used mainly by itkParticleNeighborhood*/
template <unsigned int VDimension>
struct ParticlePointIndexPair
{
  typedef Point<double, VDimension> PointType;
  
  ParticlePointIndexPair() {}
  ParticlePointIndexPair(const PointType &p, unsigned int i=0)
  {
    Point = p;
    Index = i;
  }
  ParticlePointIndexPair(const ParticlePointIndexPair &o)
  { this->operator=(o); }
  const ParticlePointIndexPair &operator=(const ParticlePointIndexPair &o)
  {
    Point = o.Point;
    Index = o.Index;
	return *this;
  }
  PointType Point;
  unsigned int Index;
};

}
#endif
