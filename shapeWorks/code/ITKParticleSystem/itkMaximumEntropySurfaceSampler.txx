/*=========================================================================
  Program:   ShapeWorks: Particle-based Shape Correspondence & Visualization
  Module:    $RCSfile: itkMaximumEntropySurfaceSampler.txx,v $
  Date:      $Date: 2009/05/06 21:49:15 $
  Version:   $Revision: 1.1.1.1 $
  Author:    $Author: cates $

  Copyright (c) 2009 Scientific Computing and Imaging Institute.
  See ShapeWorksLicense.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.
=========================================================================*/
#ifndef __itkMaximumEntropySurfaceSampler_txx
#define __itkMaximumEntropySurfaceSampler_txx

#include "itkParticlePositionReader.h"
#include "itkImageRegionIterator.h"
#include "object_reader.h"

namespace itk
{

template <class TImage>
MaximumEntropySurfaceSampler<TImage>::MaximumEntropySurfaceSampler()
{
  m_AdaptivityMode = 0;
  m_Initializing = false;


  m_PrefixTransformFile = "";
  m_TransformFile = "";
  
  // Allocate the particle system members.
  m_ParticleSystem = ParticleSystem<Dimension>::New();

  m_GradientFunction
    = ParticleEntropyGradientFunction<typename ImageType::PixelType, Dimension>::New();
  m_QualifierGradientFunction
    = ParticleQualifierEntropyGradientFunction<typename ImageType::PixelType, Dimension>::New();
  m_CurvatureGradientFunction
    = ParticleCurvatureEntropyGradientFunction<typename ImageType::PixelType, Dimension>::New();
  m_OmegaGradientFunction
    = ParticleOmegaGradientFunction<typename ImageType::PixelType, Dimension>::New();
  
  // Allocate some optimization code.
  m_Optimizer = OptimizerType::New();

  m_Initialized = false;
  m_PointsFiles.push_back("");
}

template <class TImage>
void
MaximumEntropySurfaceSampler<TImage>::AllocateWorkingImages()
{
  m_WorkingImages.resize(this->GetNumberOfInputs());
  for (unsigned int i = 0; i < this->GetNumberOfInputs(); i++)
     {
     m_WorkingImages[i] = const_cast<TImage *>(this->GetInput(i));
     }
}

template <class TImage>
void
MaximumEntropySurfaceSampler<TImage>::AllocateDataCaches()
{
  // Set up the various data caches that the optimization functions will use.
  m_Sigma1Cache = ParticleContainerArrayAttribute<double, Dimension>::New();
  m_ParticleSystem->RegisterAttribute(m_Sigma1Cache);
  m_GradientFunction->SetSpatialSigmaCache(m_Sigma1Cache);
  m_QualifierGradientFunction->SetSpatialSigmaCache(m_Sigma1Cache);
  m_CurvatureGradientFunction->SetSpatialSigmaCache(m_Sigma1Cache);
  m_OmegaGradientFunction->SetSpatialSigmaCache(m_Sigma1Cache);

  m_Sigma2Cache = ParticleContainerArrayAttribute<double, Dimension>::New();
  m_ParticleSystem->RegisterAttribute(m_Sigma2Cache);
  
  m_MeanCurvatureCache = ParticleMeanCurvatureAttribute<typename ImageType::PixelType, Dimension>::New();
  m_CurvatureGradientFunction->SetMeanCurvatureCache(m_MeanCurvatureCache);
  m_OmegaGradientFunction->SetMeanCurvatureCache(m_MeanCurvatureCache);
  m_ParticleSystem->RegisterAttribute(m_MeanCurvatureCache);
}

template <class TImage>
void
MaximumEntropySurfaceSampler<TImage>::AllocateDomainsAndNeighborhoods()
{
  // Allocate all the necessary domains and neighborhoods. This must be done
  // *after* registering the attributes to the particle system since some of
  // them respond to AddDomain.
  for (int i = 0; i < this->GetNumberOfInputs(); i++)
    { 
    m_DomainList.push_back( ParticleImplicitSurfaceDomain<typename
                            ImageType::PixelType, Dimension>::New() );
    //    m_NeighborhoodList.push_back(ParticleRegionNeighborhood<Dimension>::New());
    m_NeighborhoodList.push_back( ParticleSurfaceNeighborhood<ImageType>::New() );

    m_DomainList[i]->SetSigma(m_WorkingImages[i]->GetSpacing()[0] * 2.0);
    
    m_DomainList[i]->SetImage(m_WorkingImages[i]);

    if (m_CuttingPlanes.size() > i)
      {        
      m_DomainList[i]->SetCuttingPlane(m_CuttingPlanes[i].a,
                                       m_CuttingPlanes[i].b,
                                       m_CuttingPlanes[i].c);
      }

    if (m_Spheres.size() > i)
      {
      for (unsigned int j = 0; j < m_Spheres[i].size();j++)
        {
        m_DomainList[i]->AddSphere(m_Spheres[i][j].center,m_Spheres[i][j].radius);
        }
      }
    
      // END TEST CUTTING PLANE
    
    m_ParticleSystem->AddDomain(m_DomainList[i]);
    m_ParticleSystem->SetNeighborhood(i, m_NeighborhoodList[i]);
    }
}

template <class TImage>
void
MaximumEntropySurfaceSampler<TImage>::ReadPointsFiles()
{
  // If points file names have been specified, then read the initial points.
  for (unsigned int i = 0; i < m_PointsFiles.size(); i++)
    {
    if (m_PointsFiles[i] != "")
      {
      itk::ParticlePositionReader<3>::Pointer reader
        = itk::ParticlePositionReader<3>::New();
      reader->SetFileName(m_PointsFiles[i].c_str());
      reader->Update();
      this->GetParticleSystem()->AddPositionList(reader->GetOutput(), i);
      }
    }

  // Push position information out to all observers (necessary to correctly
  // fill out the shape matrix).
  this->GetParticleSystem()->SynchronizePositions();
}

template <class TImage>
void
MaximumEntropySurfaceSampler<TImage>::InitializeOptimizationFunctions()
{
  // Set the minimum neighborhood radius and maximum sigma based on the
  // domain of the 1st input image.
  unsigned int maxdim = 0;
  double maxradius = 0.0;
  double spacing = this->GetInput()->GetSpacing()[0];
  for (unsigned int i = 0; i < TImage::ImageDimension; i++)
    {
    if (this->GetInput()->GetRequestedRegion().GetSize()[i] > maxdim)
      {
      maxdim = this->GetInput()->GetRequestedRegion().GetSize()[i];
      maxradius = maxdim * this->GetInput()->GetSpacing()[i];
      }
    }
  
  // Initialize member variables of the optimization functions.
  //  m_GradientFunction->SetMinimumNeighborhoodRadius(maxradius / 3.0);
  m_GradientFunction->SetMinimumNeighborhoodRadius(spacing * 5.0);
  m_GradientFunction->SetMaximumNeighborhoodRadius(maxradius);

  m_QualifierGradientFunction->SetMinimumNeighborhoodRadius(spacing * 5.0);
  m_QualifierGradientFunction->SetMaximumNeighborhoodRadius(maxradius);
  
  m_CurvatureGradientFunction->SetMinimumNeighborhoodRadius(spacing * 5.0);
  m_CurvatureGradientFunction->SetMaximumNeighborhoodRadius(maxradius);
  m_CurvatureGradientFunction->SetParticleSystem(this->GetParticleSystem());
  m_CurvatureGradientFunction->SetDomainNumber(0);

  m_OmegaGradientFunction->SetMinimumNeighborhoodRadius(spacing * 5.0);
  m_OmegaGradientFunction->SetMaximumNeighborhoodRadius(maxradius);
  m_OmegaGradientFunction->SetParticleSystem(this->GetParticleSystem());
  m_OmegaGradientFunction->SetDomainNumber(0);
}
  
template <class TImage>
void
MaximumEntropySurfaceSampler<TImage>::GenerateData()
{
  this->SetInPlace(false); // this is required so that we don't release our inputs
  if (m_Initialized == false)
    {
    // TEMPORARY HACK.  REALLY WHAT I WANT TO DO IS GRAFT THE INPUT IMAGES TO
    // AN ARRAY OF OUTPUT IMAGES AND USE THOSE.
    this->AllocateWorkingImages();
    this->AllocateDataCaches();
    this->SetAdaptivityMode(m_AdaptivityMode);
    this->AllocateDomainsAndNeighborhoods();

    // Point the optimizer to the particle system.
    m_Optimizer->SetParticleSystem(m_ParticleSystem);
    this->ReadTransforms();
    this->ReadPointsFiles();
    this->InitializeOptimizationFunctions();
    m_Initialized = true;
    }

  if (m_Initializing == true) return;
  m_Optimizer->StartOptimization();
}


template <class TImage>
void
MaximumEntropySurfaceSampler<TImage>::ReadTransforms()
{
  if (m_TransformFile != "")
    {
    object_reader< itk::ParticleSystem<3>::TransformType > reader;
    reader.SetFileName(m_TransformFile.c_str());
    reader.Update();
    
    for (unsigned int i = 0; i <this->GetParticleSystem()->GetNumberOfDomains();
         i++)
      {
      std::cout << "Transform " << i << std::endl << reader.GetOutput()[i] << std::endl;
      this->GetParticleSystem()->SetTransform(i, reader.GetOutput()[i]);
      }
    }
// }
// void CorrespondenceApp::ReadPrefixTransformFile(const std::string &fn)
// {

   if (m_PrefixTransformFile != "")
     {
     object_reader< itk::ParticleSystem<3>::TransformType > reader;
     reader.SetFileName(m_PrefixTransformFile.c_str());
     reader.Update();
     
     for (unsigned int i = 0; i < this->GetParticleSystem()->GetNumberOfDomains();
          i++)
       {
       std::cout << "Prefix Transform " << i << std::endl << reader.GetOutput()[i] << std::endl;
       this->GetParticleSystem()->SetPrefixTransform(i, reader.GetOutput()[i]);
       }
     }

   // }
}


} // end namespace

#endif
