/*=========================================================================
  Program:   ShapeWorks: Particle-based Shape Correspondence & Visualization
  Module:    $RCSfile: itkParticleGradientDescentPositionOptimizer.txx,v $
  Date:      $Date: 2009/05/06 21:49:15 $
  Version:   $Revision: 1.1.1.1 $
  Author:    $Author: cates $

  Copyright (c) 2009 Scientific Computing and Imaging Institute.
  See ShapeWorksLicense.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.
=========================================================================*/
#ifndef __itkParticleGradientDescentPositionOptimizer_txx
#define __itkParticleGradientDescentPositionOptimizer_txx

namespace itk
{
template <class TGradientNumericType, unsigned int VDimension>
ParticleGradientDescentPositionOptimizer<TGradientNumericType, VDimension>
::ParticleGradientDescentPositionOptimizer()
{
  m_StopOptimization = false;
  m_NumberOfIterations = 0;
  m_MaximumNumberOfIterations = 0;
  m_Tolerance = 0.0;
  m_TimeStep = 1.0;
  m_OptimizationMode = 0;
}

/*** ADAPTIVE GAUSS SEIDEL ***/
template <class TGradientNumericType, unsigned int VDimension>
void
ParticleGradientDescentPositionOptimizer<TGradientNumericType, VDimension>
::StartAdaptiveGaussSeidelOptimization()
{
  const double factor = 1.1;//1.1;
  //  const double epsilon = 1.0e-4;

  // NOTE: THIS METHOD WILL NOT WORK AS WRITTEN IF PARTICLES ARE
  // ADDED TO THE SYSTEM DURING OPTIMIZATION.
  m_StopOptimization = false;
  m_GradientFunction->SetParticleSystem(m_ParticleSystem);

  bool reset = false;
  // Make sure the time step vector is the right size
  while (m_TimeSteps.size() != m_ParticleSystem->GetNumberOfDomains() )
    {
    reset = true;
    std::vector<double> tmp;
    m_TimeSteps.push_back( tmp );
    }

  for (unsigned int i = 0; i < m_ParticleSystem->GetNumberOfDomains(); i++)
    {
    unsigned int np = m_ParticleSystem->GetPositions(i)->GetSize();
    if (m_TimeSteps[i].size() != np)
      {
      // resize and initialize everything to 1.0
      m_TimeSteps[i].resize(np);
      for (unsigned int j = 0; j < np; j++)
        {
        m_TimeSteps[i][j] = 1.0;
        }
      reset = true;
      }
    }

  //  if (reset == true) m_GradientFunction->ResetBuffers();
  
  VectorType gradient;
  VectorType original_gradient;
  PointType newpoint;

  unsigned int numdomains = m_ParticleSystem->GetNumberOfDomains();
  std::vector<double> meantime(numdomains);
  std::vector<double> maxtime(numdomains);
  std::vector<double> mintime(numdomains);

  for (unsigned int q = 0; q < numdomains; q++)
    {
    meantime[q] = 0.0;
    maxtime[q]  = 1.0e30;
    mintime[q]  = 1.0;
    }

  double maxchange = 0.0;
  while (m_StopOptimization == false)
    {
    m_GradientFunction->BeforeIteration();
    double maxdt;

    // Iterate over each domain
    for (unsigned int dom = 0; dom < numdomains; dom++)
      {
      meantime[dom] = 0.0;
      // skip any flagged domains
      if (m_ParticleSystem->GetDomainFlag(dom) == false)
        {
        // Tell function which domain we are working on.
        m_GradientFunction->SetDomainNumber(dom);
     
        // Iterate over each particle position
        unsigned int k = 0;
        maxchange = 0.0;
        typename ParticleSystemType::PointContainerType::ConstIterator endit =
          m_ParticleSystem->GetPositions(dom)->GetEnd();
        for (typename ParticleSystemType::PointContainerType::ConstIterator it
               = m_ParticleSystem->GetPositions(dom)->GetBegin(); it != endit; it++, k++)
          {
          bool done = false;

          // Compute gradient update.
          double energy = 0.0;
          m_GradientFunction->BeforeEvaluate(it.GetIndex(), dom, m_ParticleSystem);
          original_gradient = m_GradientFunction->Evaluate(it.GetIndex(), dom, m_ParticleSystem,
                                                           maxdt, energy);
          PointType pt = *it;

          double newenergy, gradmag;
          while ( !done )
            {
            gradient = original_gradient * m_TimeSteps[dom][k];

            dynamic_cast<DomainType *>(m_ParticleSystem->GetDomain(dom))
              ->ApplyVectorConstraints(gradient,
                                       m_ParticleSystem->GetPosition(it.GetIndex(), dom), maxdt);
            gradmag = gradient.magnitude();
            
            // Prevent a move which is too large
            if (gradmag * m_TimeSteps[dom][k] > maxdt)
              {
              m_TimeSteps[dom][k] /= factor;
              }
            else // Move is not too large
              {
              // Make a move and compute new energy
              for (unsigned int i = 0; i < VDimension; i++)
                {  newpoint[i] = pt[i] - gradient[i]; }
              m_ParticleSystem->SetPosition(newpoint, it.GetIndex(), dom);
              newenergy = m_GradientFunction->Energy(it.GetIndex(), dom, m_ParticleSystem);
              
              if (newenergy < energy) // good move, increase timestep for next time
                {
                meantime[dom] += m_TimeSteps[dom][k];
                m_TimeSteps[dom][k] *= factor;
                if (m_TimeSteps[dom][k] > maxtime[dom]) m_TimeSteps[dom][k] = maxtime[dom];
                if (gradmag > maxchange) maxchange = gradmag;
                done = true;
                }
              else 
                {// bad move, reset point position and back off on timestep            
                if (m_TimeSteps[dom][k] > mintime[dom])
                  {
                  m_ParticleSystem->SetPosition(pt, it.GetIndex(), dom);

                  m_TimeSteps[dom][k] /= factor;
                  }
                else // keep the move with timestep 1.0 anyway
                  {
                  meantime[dom] += m_TimeSteps[dom][k];
                  m_TimeSteps[dom][k] = mintime[dom];
                  if (gradmag > maxchange) maxchange = gradmag;
                  done = true;
                  }
                }
              }
            } // end while not done
          } // for each particle
        
        // Compute mean time step
        meantime[dom] /= static_cast<double>(k);

        if (meantime[dom] < 1.0) meantime[dom] = 1.0;
        //        std::cout << "meantime = " << meantime[dom] << std::endl;
        maxtime[dom] = meantime[dom] + meantime[dom] * 0.2;
        mintime[dom] = meantime[dom] - meantime[dom] * 0.1;
        } // if not flagged
      }// for each domain
    
    m_NumberOfIterations++;
    m_GradientFunction->AfterIteration();
    this->InvokeEvent(itk::IterationEvent());
    
    // Check for convergence.  Optimization is considered to have converged if
    // max number of iterations is reached or maximum distance moved by any
    // particle is less than the specified precision.
    //    std::cout << "maxchange = " << maxchange << std::endl;
    if ((m_NumberOfIterations == m_MaximumNumberOfIterations)
      || (m_Tolerance > 0.0 &&  maxchange <  m_Tolerance) )
      {
      m_StopOptimization = true;
      }
    
    } // end while stop optimization
  
}

/*** GAUSS SEIDEL ***/
template <class TGradientNumericType, unsigned int VDimension>
void
ParticleGradientDescentPositionOptimizer<TGradientNumericType, VDimension>
::StartGaussSeidelOptimization()
{
  // NOTE: THIS METHOD WILL NOT WORK AS WRITTEN IF PARTICLES ARE
  // ADDED TO THE SYSTEM DURING OPTIMIZATION.
  m_StopOptimization = false;
  m_GradientFunction->SetParticleSystem(m_ParticleSystem);
  
  VectorType gradient;
  PointType newpoint;
  
  while (m_StopOptimization == false)
    {
    m_GradientFunction->BeforeIteration();
    double maxdt;
    
    // Iterate over each domain
    for (unsigned int dom = 0; dom < m_ParticleSystem->GetNumberOfDomains(); dom++)
      {
      // skip any flagged domains
      if (m_ParticleSystem->GetDomainFlag(dom) == false)
        {
        // Tell function which domain we are working on.
        m_GradientFunction->SetDomainNumber(dom);
        
        // Iterate over each particle position
        unsigned int k = 0;
        typename ParticleSystemType::PointContainerType::ConstIterator endit =
          m_ParticleSystem->GetPositions(dom)->GetEnd();
        for (typename ParticleSystemType::PointContainerType::ConstIterator it
               = m_ParticleSystem->GetPositions(dom)->GetBegin(); it != endit; it++, k++)
          {
          // Compute gradient update.
          m_GradientFunction->BeforeEvaluate(it.GetIndex(), dom, m_ParticleSystem);
          gradient = m_GradientFunction->Evaluate(it.GetIndex(), dom, m_ParticleSystem,
                                                  maxdt);
          
          // May modify gradient
          dynamic_cast<DomainType *>(m_ParticleSystem->GetDomain(dom))
            ->ApplyVectorConstraints(gradient,
                                     m_ParticleSystem->GetPosition(it.GetIndex(), dom), maxdt);
          
          // Hack to avoid blowing up under certain conditions.
          if (gradient.magnitude() > maxdt)
            { gradient = (gradient / gradient.magnitude()) * maxdt; }
          
          // Compute particle move based on update.
           for (unsigned int i = 0; i < VDimension; i++)
             {
             newpoint[i] = (*it)[i] - gradient[i] * m_TimeStep;
             }
           
           // Apply update
           m_ParticleSystem->SetPosition(newpoint, it.GetIndex(), dom);
           } // for each particle
         } // if not flagged
       }// for each domain

     m_NumberOfIterations++;
     m_GradientFunction->AfterIteration();
     this->InvokeEvent(itk::IterationEvent());
  
     // Check for convergence.  Optimization is considered to have converged if
     // max number of iterations is reached or maximum distance moved by any
     // particle is less than the specified precision.
     if (m_NumberOfIterations == m_MaximumNumberOfIterations)
       {
       m_StopOptimization = true;
       }
  
     } // end while

 }

template <class TGradientNumericType, unsigned int VDimension>
void
ParticleGradientDescentPositionOptimizer<TGradientNumericType, VDimension>
::StartJacobiOptimization()
{
  // NOTE: THIS METHOD WILL NOT WORK AS WRITTEN IF PARTICLES ARE
  // ADDED TO THE SYSTEM DURING OPTIMIZATION.
  m_StopOptimization = false;

  m_GradientFunction->SetParticleSystem(m_ParticleSystem);
  
  // Vector for storing updates.
  std::vector<std::vector<PointType> > updates(m_ParticleSystem->GetNumberOfDomains());

  for (unsigned int d = 0; d < m_ParticleSystem->GetNumberOfDomains(); d++)
    {
    updates[d].resize(m_ParticleSystem->GetPositions(d)->GetSize());
    }
  VectorType gradient;
  PointType  newpoint;

  while (m_StopOptimization == false)
    {
    m_GradientFunction->BeforeIteration();
    double maxdt;
    
    // Iterate over each domain
    for (unsigned int dom = 0; dom < m_ParticleSystem->GetNumberOfDomains(); dom++)
      {
      // skip any flagged domains
      if (m_ParticleSystem->GetDomainFlag(dom) == false)
        {
        // Tell function which domain we are working on.
        m_GradientFunction->SetDomainNumber(dom);
        
        // Iterate over each particle position
        unsigned int k = 0;
        typename ParticleSystemType::PointContainerType::ConstIterator endit =
          m_ParticleSystem->GetPositions(dom)->GetEnd();
        for (typename ParticleSystemType::PointContainerType::ConstIterator it
               = m_ParticleSystem->GetPositions(dom)->GetBegin(); it != endit; it++, k++)
          {
          // Compute gradient update.
          m_GradientFunction->BeforeEvaluate(it.GetIndex(), dom, m_ParticleSystem);
          gradient = m_GradientFunction->Evaluate(it.GetIndex(), dom, m_ParticleSystem,
                                                  maxdt);
          // May modify gradient
          dynamic_cast<DomainType *>(m_ParticleSystem->GetDomain(dom))
            ->ApplyVectorConstraints(gradient, m_ParticleSystem->GetPosition(it.GetIndex(), dom), maxdt);
          
          // Hack to avoid blowing up under certain conditions.
          if (gradient.magnitude() > maxdt)
            {
            gradient = (gradient / gradient.magnitude()) * maxdt;
            }
          
          // Compute particle move based on update.
          for (unsigned int i = 0; i < VDimension; i++)
            {
            newpoint[i] = (*it)[i] - gradient[i] * m_TimeStep;
            }
          
          // Store update
          updates[dom][k] = newpoint;
          } // for each particle
        } // if not flagged
      }// for each domain
    
    for (unsigned int dom = 0; dom < m_ParticleSystem->GetNumberOfDomains(); dom++)
      {
      // skip any flagged domains
      if (m_ParticleSystem->GetDomainFlag(dom) == false)
        {
        unsigned int k = 0;
        typename ParticleSystemType::PointContainerType::ConstIterator endit =
          m_ParticleSystem->GetPositions(dom)->GetEnd();
        for (typename ParticleSystemType::PointContainerType::ConstIterator it
               = m_ParticleSystem->GetPositions(dom)->GetBegin();  it != endit; it++, k++)
          {
          m_ParticleSystem->SetPosition(updates[dom][k], it.GetIndex(), dom);
          } // for each particle
        } // if not flagged
      } // for each domain
    
    m_NumberOfIterations++;
    m_GradientFunction->AfterIteration();
    this->InvokeEvent(itk::IterationEvent());
    
    // Check for convergence.  Optimization is considered to have converged if
    // max number of iterations is reached or maximum distance moved by any
    // particle is less than the specified precision.
    if (m_NumberOfIterations == m_MaximumNumberOfIterations)
      {
      m_StopOptimization = true;
      }
    
    } // end while
  
}

} // end namespace

#endif
