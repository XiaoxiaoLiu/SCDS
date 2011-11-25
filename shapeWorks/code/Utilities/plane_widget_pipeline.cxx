/*=========================================================================
  Program:   ShapeWorks: Particle-based Shape Correspondence & Visualization
  Module:    $RCSfile: plane_widget_pipeline.cxx,v $
  Date:      $Date: 2009/05/06 21:49:16 $
  Version:   $Revision: 1.1.1.1 $
  Author:    $Author: cates $

  Copyright (c) 2009 Scientific Computing and Imaging Institute.
  See ShapeWorksLicense.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.
=========================================================================*/
#include "plane_widget_pipeline.h"

void plane_widget_pipeline::SetTransformCallback(itk::Object *o, const itk::EventObject &e)
{
  // NOTE: Ignoring scale
  const itk::ParticleTransformSetEvent &event = dynamic_cast<const itk::ParticleTransformSetEvent &>(e);
  const itk::ParticleSystem<3> *ps= dynamic_cast<const itk::ParticleSystem<3> *>(o);

  unsigned int d = event.GetDomainIndex();
  if (d != m_MyDomain)
    {
    return;
    }

  vtkMatrix4x4 *m1;
  vtkTransform *t;
  //  vtkMatrix4x4 *m;
  //  vtkMatrix4x4 *m0;
  
  //  m  = vtkMatrix4x4::New();
  //  m0 = vtkMatrix4x4::New();
  m1 = vtkMatrix4x4::New();
  t  = vtkTransform::New();

  itk::ParticleSystem<3>::TransformType T = ps->GetTransform(d) * ps->GetPrefixTransform(d);

  m1->SetElement(0, 0, T(0,0)); // * transforms[i].scale);
  m1->SetElement(1, 0, T(1,0)); // * transforms[i].scale);
  m1->SetElement(2, 0, T(2,0)); // * transforms[i].scale);
  m1->SetElement(3, 0, 0.0);

  m1->SetElement(0, 1, T(0,1)); // * transforms[i].scale);
  m1->SetElement(1, 1, T(1,1)); // * transforms[i].scale);
  m1->SetElement(2, 1, T(2,1)); // * transforms[i].scale);
  m1->SetElement(3, 1, 0.0);
    
  m1->SetElement(0, 2, T(0,2)); // * transforms[i].scale);
  m1->SetElement(1, 2, T(1,2)); // * transforms[i].scale);
  m1->SetElement(2, 2, T(2,2)); // * transforms[i].scale);
  m1->SetElement(3, 2, 0.0);

  m1->SetElement(0, 3, T(0,3)); //transforms.translation(0));
  m1->SetElement(1, 3, T(1,3)); //transforms.translation(1));
  m1->SetElement(2, 3, T(2,3)); //transforms.translation(2));
  m1->SetElement(3, 3, 1.0);

  //  m0 =reinterpret_cast<vtkTransform *>(m_transformer->GetTransform())->GetMatrix();
      
  //  vtkMatrix4x4::Multiply4x4(m1, m0, m);
  t->SetMatrix(m1);
    
  //print final matrix
//    std::cout << "Domain " << d << std::endl;
//    for (unsigned int r = 0; r < 4; r++)
//      {
//      for (unsigned int c = 0; c < 4; c++)
//        {
//        std::cout << m1->GetElement(r,c) << " ";
//        }
//      std::cout << std::endl;
//      }
  
  m_transformer->SetTransform(t);
  //  this->output()->SetUserTransform(t);
}

void plane_widget_pipeline::SetPrefixTransformCallback(itk::Object *o, const itk::EventObject &e)
{
  // NOTE: Ignoring scale
  const itk::ParticlePrefixTransformSetEvent &event = dynamic_cast<const itk::ParticlePrefixTransformSetEvent &>(e);
  const itk::ParticleSystem<3> *ps= dynamic_cast<const itk::ParticleSystem<3> *>(o);

  unsigned int d = event.GetDomainIndex();
  if (d != m_MyDomain)
    {
    return;
    }

  vtkMatrix4x4 *m1;
  vtkTransform *t;
  //  vtkMatrix4x4 *m;
  //  vtkMatrix4x4 *m0;
  
  //  m  = vtkMatrix4x4::New();
  //  m0 = vtkMatrix4x4::New();
  m1 = vtkMatrix4x4::New();
  t  = vtkTransform::New();

  itk::ParticleSystem<3>::TransformType T = ps->GetTransform(d) * ps->GetPrefixTransform(d);

  m1->SetElement(0, 0, T(0,0)); // * transforms[i].scale);
  m1->SetElement(1, 0, T(1,0)); // * transforms[i].scale);
  m1->SetElement(2, 0, T(2,0)); // * transforms[i].scale);
  m1->SetElement(3, 0, 0.0);

  m1->SetElement(0, 1, T(0,1)); // * transforms[i].scale);
  m1->SetElement(1, 1, T(1,1)); // * transforms[i].scale);
  m1->SetElement(2, 1, T(2,1)); // * transforms[i].scale);
  m1->SetElement(3, 1, 0.0);
    
  m1->SetElement(0, 2, T(0,2)); // * transforms[i].scale);
  m1->SetElement(1, 2, T(1,2)); // * transforms[i].scale);
  m1->SetElement(2, 2, T(2,2)); // * transforms[i].scale);
  m1->SetElement(3, 2, 0.0);

  m1->SetElement(0, 3, T(0,3)); //transforms.translation(0));
  m1->SetElement(1, 3, T(1,3)); //transforms.translation(1));
  m1->SetElement(2, 3, T(2,3)); //transforms.translation(2));
  m1->SetElement(3, 3, 1.0);

  //  m0 =reinterpret_cast<vtkTransform *>(m_transformer->GetTransform())->GetMatrix();
      
  //  vtkMatrix4x4::Multiply4x4(m1, m0, m);
  t->SetMatrix(m1);
    
  //print final matrix
//    std::cout << "Domain " << d << std::endl;
//    for (unsigned int r = 0; r < 4; r++)
//      {
//      for (unsigned int c = 0; c < 4; c++)
//        {
//        std::cout << m1->GetElement(r,c) << " ";
//        }
//      std::cout << std::endl;
//      }
  
  m_transformer->SetTransform(t);
  //  this->output()->SetUserTransform(t);
}

plane_widget_pipeline::plane_widget_pipeline()
{
  //  m_plane_widget = vtkPlaneWidget::New();
  m_plane_source = vtkPlaneSource::New();
    
  m_transformer= vtkTransformPolyDataFilter::New();

  vtkTransform *tmp = vtkTransform::New();
  tmp->Identity();
  m_transformer->SetTransform(tmp);
  m_transformer->SetInput(m_plane_source->GetOutput());
  
  m_mapper = vtkPolyDataMapper::New();
  
  m_mapper->SetInput(m_transformer->GetOutput());
  // m_mapper->StaticOn();
  // m_mapper->ImmediateModeRenderingOff();
  m_mapper->ScalarVisibilityOff();
  
  m_actor = vtkActor::New();    
  //   m_actor->GetProperty()->SetSpecularColor(1.0, 1.0, 1.0);
  //   m_actor->GetProperty()->SetDiffuse(0.0);
  //   m_actor->GetProperty()->SetSpecular(0.8);
  //   m_actor->GetProperty()->SetSpecularPower(5.0);
  m_actor->SetMapper(m_mapper);

}

plane_widget_pipeline::~plane_widget_pipeline()
{
  m_mapper->Delete();
  m_actor->Delete();
  m_plane_source->Delete();
}


