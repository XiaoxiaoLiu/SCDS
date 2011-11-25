/*=========================================================================
  Program:   ShapeWorks: Particle-based Shape Correspondence & Visualization
  Module:    $RCSfile: fastmarching.h,v $
  Date:      $Date: 2009/05/06 21:49:15 $
  Version:   $Revision: 1.1.1.1 $
  Author:    $Author: cates $

  Copyright (c) 2009 Scientific Computing and Imaging Institute.
  See ShapeWorksLicense.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.
=========================================================================*/
#ifndef __st_fastmarching_h
#define __st_fastmarching_h

#include "itkImage.h"
#include "param.h"
#include "tool.h"

namespace shapetools
{
/**
 * \class fastmarching
 *
 *
 */
template <class T, unsigned int D> 
class fastmarching : public tool<T, D>
{
public:
  typedef T pixel_type;
  typedef itk::Image<T, D> image_type;
  
  fastmarching(param::parameterFile &);
  fastmarching() { m_levelset_value = 0.0; }
  virtual ~fastmarching() {}
  
  virtual void operator()(typename image_type::Pointer);

  
private: 
  double m_levelset_value;
  
};

 
} // end namespace 
#endif

#ifndef ST_MANUAL_INSTANTIATION
#include "fastmarching.txx"
#endif
