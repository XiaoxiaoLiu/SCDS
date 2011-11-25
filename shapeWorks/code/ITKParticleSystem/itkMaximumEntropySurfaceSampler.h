/*=========================================================================
  Program:   ShapeWorks: Particle-based Shape Correspondence & Visualization
  Module:    $RCSfile: itkMaximumEntropySurfaceSampler.h,v $
  Date:      $Date: 2009/05/06 21:49:15 $
  Version:   $Revision: 1.1.1.1 $
  Author:    $Author: cates $

  Copyright (c) 2009 Scientific Computing and Imaging Institute.
  See ShapeWorksLicense.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.
=========================================================================*/
#ifndef __itkMaximumEntropySurfaceSampler_h
#define __itkMaximumEntropySurfaceSampler_h

#include "itkParticleSystem.h"
#include "itkParticleGradientDescentPositionOptimizer.h"
#include "itkParticleEntropyGradientFunction.h"
#include "itkParticleQualifierEntropyGradientFunction.h"
#include "itkParticleImplicitSurfaceDomain.h"
#include "itkInPlaceImageFilter.h"
#include "itkParticleContainerArrayAttribute.h"
#include "itkParticleCurvatureEntropyGradientFunction.h"
#include "itkParticleMeanCurvatureAttribute.h"
#include "itkParticleSurfaceNeighborhood.h"
#include "itkParticleOmegaGradientFunction.h"

namespace itk
{
  
/** \class MaximumEntropySurfaceSampler
 *
 * 
 *
 */
template <class TImage>
class ITK_EXPORT MaximumEntropySurfaceSampler
  : public InPlaceImageFilter<TImage,TImage> 
{
public:
  /** Standard class typedefs. */
  typedef MaximumEntropySurfaceSampler  Self;
  typedef InPlaceImageFilter<TImage,TImage>  Superclass;
  typedef SmartPointer<Self>   Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);
  
  /** Run-time type information (and related methods). */
  itkTypeMacro(MaximumEntropySurfaceSampler, InPlaceImageFilter);
  
  /**Expose the image dimension. */
  itkStaticConstMacro(Dimension, unsigned int, TImage::ImageDimension);
  
  /** Convenient typedef for storing cutting plane information */
  struct CuttingPlaneType
  {
    vnl_vector_fixed<double,Dimension> a;
    vnl_vector_fixed<double,Dimension> b;
    vnl_vector_fixed<double,Dimension> c;
  };

  /** Convenient typedef for storing sphere information */
  struct SphereType
  {
    vnl_vector_fixed<double,Dimension> center;
    double radius;
  };
  
  /** Type of the input/output image. */
  typedef TImage ImageType;
  typedef ParticleGradientDescentPositionOptimizer<typename ImageType::PixelType, Dimension> OptimizerType;

  /**
   * THIS IS A HACK UNTIL I CAN FIGURE OUT HOW TO ALLOCATE THE APPROPRIATE
   * NUMBER OF OUTPUTS AND GRAFT THEM TO THE INPUTS.
   */
  virtual void AllocateWorkingImages();
  
  void SetInput(const TImage *image)
  { this->SetInput(0, image);  }
  
  
  /**
   * Override parent classes to expand input list on new inputs.
   */
  void SetInput( unsigned int index, const TImage * image ) 
  {
    if (this->GetNumberOfInputs() < index+1)
      {
      this->SetNumberOfRequiredInputs(index+1);
      }
    
    this->ProcessObject::SetNthInput(index, const_cast< TImage *>( image ) );
  }
  
  /** Returns the particle system used in the surface sampling. */
  itkGetObjectMacro(ParticleSystem, ParticleSystem<Dimension>);
  itkGetConstObjectMacro(ParticleSystem, ParticleSystem<Dimension>);
  
  /** Returns a pointer to the gradient function used. */
  ParticleEntropyGradientFunction<typename ImageType::PixelType, Dimension>
  *GetGradientFunction()
  {
    return m_GradientFunction;
  }
     
  ParticleQualifierEntropyGradientFunction<typename ImageType::PixelType, Dimension>
  *GetQualifierGradientFunction()
  {
    return m_QualifierGradientFunction;
  }
  ParticleCurvatureEntropyGradientFunction<typename ImageType::PixelType, Dimension>
  *GetCurvatureGradientFunction()
  {
    return m_CurvatureGradientFunction;
  }

    ParticleOmegaGradientFunction<typename ImageType::PixelType, Dimension>
  *GetOmegaGradientFunction()
  {
    return m_OmegaGradientFunction;
  }  
  
  /** Return a pointer to the optimizer object. */
  itkGetObjectMacro(Optimizer, OptimizerType);
  itkGetConstObjectMacro(Optimizer, OptimizerType);
  
  /**Optionally provide a filename for an initial point set.*/
  void SetPointsFile(unsigned int i, const std::string &s)
  {
    if (m_PointsFiles.size() < i+1)
      {
      m_PointsFiles.resize(i+1);
      }
    m_PointsFiles[i] = s;
  }
  void SetPointsFile(const std::string &s)
  {
    this->SetPointsFile(0,s);
  }

  /** Optionally supply a cutting plane that will be set as a particle
      optimization constraint in the image domains. */
  void SetCuttingPlane(unsigned int i,
                       const vnl_vector_fixed<double,Dimension> &va,
                       const vnl_vector_fixed<double,Dimension> &vb,
                       const vnl_vector_fixed<double,Dimension> &vc)
  {
    if (m_CuttingPlanes.size() < i+1)
      {
      m_CuttingPlanes.resize(i+1);
      }
    
    m_CuttingPlanes[i].a = va;
    m_CuttingPlanes[i].b = vb;
    m_CuttingPlanes[i].c = vc;

    if (m_Initialized == true)
      {
      m_DomainList[i]->SetCuttingPlane(va,vb,vc);
      }
  }


  /** Optionally add spheres that may be used as constraints to the domain. */

  void AddSphere(unsigned int i, vnl_vector_fixed<double,Dimension> &c, double r)
  {
    if (m_Spheres.size() < i+1)
      {
      m_Spheres.resize(i+1);
      }

    m_Spheres[i].push_back( SphereType() );
    m_Spheres[i][m_Spheres[i].size()-1].center = c;
    m_Spheres[i][m_Spheres[i].size()-1].radius = r;
    
    if (m_Initialized == true)
      {
      m_DomainList[i]->AddSphere(c,r);
      }
  }
                 
  
  /** This method sets the optimization function for the sampling.
      mode 0 = no adaptivity
      mode 1 = isotropic adaptivity
      mode 2 = anisotropic adaptivity
  */
  virtual void SetAdaptivityMode(int mode)
  {
    if (mode == 0)
      {
      m_Optimizer->SetGradientFunction(m_CurvatureGradientFunction);
      }
    else if (mode ==1)
      {
      m_Optimizer->SetGradientFunction(m_GradientFunction);
      }
    else if (mode ==2)
      {
      m_Optimizer->SetGradientFunction(m_QualifierGradientFunction);
      }
    else if (mode ==3)
      {
      m_Optimizer->SetGradientFunction(m_OmegaGradientFunction);
      }
    
    m_AdaptivityMode = mode;
    this->Modified();
    
  }
  int GetAdaptivityMode() const
  { return m_AdaptivityMode; }

  void SetTransformFile(const std::string& s)
  { m_TransformFile = s; }
  void SetTransformFile(const char *s)
  { m_TransformFile = std::string(s); }

  void SetPrefixTransformFile(const std::string& s)
  { m_PrefixTransformFile = s; }
  void SetPrefixTransformFile(const char *s)
  { m_PrefixTransformFile = std::string(s); }
  
  void ReadTransforms();
  void ReadPointsFiles();
  virtual void AllocateDataCaches();
  virtual void AllocateDomainsAndNeighborhoods();
  virtual void InitializeOptimizationFunctions();
  
  /** */
  virtual void Initialize()
  {
    this->m_Initializing = true;
    this->Update();
    this->m_Initializing = false;
  }
  
protected:
  MaximumEntropySurfaceSampler();
  virtual ~MaximumEntropySurfaceSampler() {};

  void PrintSelf(std::ostream& os, Indent indent) const
  {
    Superclass::PrintSelf(os, indent);
  }

  void GenerateData();

  itkSetMacro(Initialized, bool);
  itkGetMacro(Initialized, bool);
  itkSetMacro(Initializing, bool);
  itkGetMacro(Initializing, bool);
  
  bool m_Initialized;
  int m_AdaptivityMode;
  bool m_Initializing;

  std::vector<typename TImage::Pointer> m_WorkingImages;
  
  typename OptimizerType::Pointer m_Optimizer;

  typename ParticleEntropyGradientFunction<typename ImageType::PixelType,  Dimension>
  ::Pointer m_GradientFunction;
  typename ParticleQualifierEntropyGradientFunction<typename ImageType::PixelType,  Dimension>
  ::Pointer m_QualifierGradientFunction;
  typename ParticleCurvatureEntropyGradientFunction<typename ImageType::PixelType, Dimension>
  ::Pointer m_CurvatureGradientFunction;

  typename ParticleOmegaGradientFunction<typename ImageType::PixelType, Dimension>
  ::Pointer m_OmegaGradientFunction;
  
  typename ParticleContainerArrayAttribute<double, Dimension>::Pointer m_Sigma1Cache;
  typename ParticleContainerArrayAttribute<double, Dimension>::Pointer m_Sigma2Cache;
  
  typename ParticleMeanCurvatureAttribute<typename ImageType::PixelType, Dimension>
  ::Pointer m_MeanCurvatureCache;
  
  typename ParticleSystem<Dimension>::Pointer m_ParticleSystem;
  
  std::vector<typename ParticleImplicitSurfaceDomain<typename
                                 ImageType::PixelType, Dimension>::Pointer> m_DomainList;
  //    std::vector<typename ParticleRegionNeighborhood<Dimension>::Pointer>
  //    m_NeighborhoodList;
    std::vector<typename ParticleSurfaceNeighborhood<ImageType>::Pointer> m_NeighborhoodList;
  
private:
  MaximumEntropySurfaceSampler(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  std::vector<std::string> m_PointsFiles;
  std::string m_TransformFile;
  std::string m_PrefixTransformFile;
  std::vector< CuttingPlaneType> m_CuttingPlanes;
  std::vector< std::vector<SphereType> > m_Spheres;
};

} // end namespace itk

#if ITK_TEMPLATE_EXPLICIT
#include "Templates/itkMaximumEntropySurfaceSampler+-.h"
#endif

#if ITK_TEMPLATE_TXX
#include "itkMaximumEntropySurfaceSampler.txx"
#endif

#endif


