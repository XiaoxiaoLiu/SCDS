//
// Implemenation of STAPLE algorithm
// Simultaneous Truth and Performance Level Estimation
// See Warfield et al. IEEE TMI Vol 23 No 7 July 2004
// See also Rholfing et al. IPMI 2003
//
// B. Davis 2004
//

#ifndef STAPLE_h
#define STAPLE_h

#include <vector>
#include "Array3D.h"
#include "Vector3D.h"
#include <iostream>

template <class LabelType, class ProbabilityType>
class STAPLE
{
public:
  typedef Array3D<LabelType>        LabelImageType;
  typedef Array3D<ProbabilityType>  ProbabilityImageType;

  static
  void
  estimateGroundTruth(const unsigned int numExperts,
		      const ProbabilityType& g,
		      LabelImageType** expertSegmentations,
		      ProbabilityType* p,
		      ProbabilityType* q,
		      ProbabilityImageType& probabilisticSegmentation,
		      unsigned int maxIterations,
		      const double& minWChange,
		      bool restrictToDisputedVoxels,
		      bool showOutput)
  {
    if (numExperts < 1)
      {
	std::cerr << "STAPLE: No Expert Segmentations."  << std::endl;
	return;
      }
    
    //
    // get input voxels (only disputed if desired)
    //
    std::vector<std::vector<LabelType> > d;
    std::vector<Vector3D<unsigned int> > disputedVoxelLocations;
    if (restrictToDisputedVoxels)
      {
	// side effect: this also sets probabilisticSegmentation for
	// non-disputed voxles
	_getDisputedVoxels(numExperts, expertSegmentations, 
			   d, disputedVoxelLocations,
			   probabilisticSegmentation);
      }
    else
      {
	_getAllVoxels(numExperts, expertSegmentations, d);
      }
    unsigned int numVoxels = d[0].size();
    std::vector<ProbabilityType> w(numVoxels);
    
    std::cerr << "STAPLE: numExperts = " << numExperts << std::endl;  
    std::cerr << "STAPLE: numVoxels  = " << numVoxels << std::endl;
    
    //
    // main EM loop
    //
    double sumForegroundProbabilities;
    double sumBackgroundProbabilities;
    double* sumExpertForegroundProbabilities = new double[numExperts];
    double* sumExpertBackgroundProbabilities = new double[numExperts];
    bool done = false;
    ProbabilityType maxCurrentWChange = 0;
    unsigned int expertIndex; // stupid vcc
    for (unsigned int iter = 0; !done && (iter < maxIterations); ++iter)
      {
	maxCurrentWChange = 0;
	
	//
	// reset global and expert specific foreground/background
	// probabilities
	//
	sumForegroundProbabilities = 0;
	sumBackgroundProbabilities = 0;
	for (expertIndex = 0; 
	     expertIndex < numExperts; ++expertIndex)
	  {
	    sumExpertForegroundProbabilities[expertIndex] = 0;
	    sumExpertBackgroundProbabilities[expertIndex] = 0;
	  }      
	
	//
	// compute probabilistic segmentation from prior and current
	// sensitivities and specifities
	//
	for (unsigned int voxelIndex = 0; voxelIndex < numVoxels; ++voxelIndex)
	  {
	    //
	    // compute alpha and beta for this voxel index
	    //
	    double alpha = 1;
	    double beta  = 1;
	    for (expertIndex = 0; 
		 expertIndex < numExperts; ++expertIndex)
	      {
		if (d[expertIndex][voxelIndex] > 0)
		  {
		    alpha *= p[expertIndex];
		    beta  *= 1 - q[expertIndex];
		  }
		else
		  {
		    alpha *= 1 - p[expertIndex];
		    beta  *= q[expertIndex];
		  }
	      }
	    
	    //
	    // update the probability of foreground for this pixel based
	    // on alpha, beta, and the prior
	    //
	    ProbabilityType lastProbability = w[voxelIndex];
	    w[voxelIndex] = 
	      static_cast<ProbabilityType>((g * alpha) / 
					   (g * alpha + (1-g) * beta)); 
	    ProbabilityType wChange = 
	      fabs(w[voxelIndex] - lastProbability); 
	    if (wChange > maxCurrentWChange)
	      {
		maxCurrentWChange = wChange;
	      }
	    
	    //
	    // increment sum of foreground/background probabilities
	    // globaly and for each expert
	    //
	    sumForegroundProbabilities += w[voxelIndex];
	    sumBackgroundProbabilities += 1 - w[voxelIndex];
	    
	    for (expertIndex = 0; expertIndex < numExperts; ++expertIndex)
	      {
		if (d[expertIndex][voxelIndex] > 0)
		  {
		    sumExpertForegroundProbabilities[expertIndex] 
		      += w[voxelIndex];
		  }
		else
		  {
		    sumExpertBackgroundProbabilities[expertIndex] 
		      += 1 - w[voxelIndex];
		  }
	      }
	  }
	
	//
	// guard against divide by zero
	//
	if (sumForegroundProbabilities == 0)
	  {
	    std::cout << "STAPLE: Stopping: sumForegroundProbabilities == 0"
		      << std::endl;
	    done = true;
	    break;
	  }
	if (sumBackgroundProbabilities == 0)
	  {
	    std::cout << "STAPLE: Stopping: sumBackgroundProbabilities == 0"
		      << std::endl;
	    done = true;
	    break;
	  }
	
	//
	// update expert sensitivities and specifities based on new
	// probabilistic segmentation
	//
	for (expertIndex = 0; 
	     expertIndex < numExperts; ++expertIndex)
	  {
	    p[expertIndex] = 
	      sumExpertForegroundProbabilities[expertIndex]
	      / sumForegroundProbabilities;
	    q[expertIndex] = 
	      sumExpertBackgroundProbabilities[expertIndex]
	      / sumBackgroundProbabilities; 
	    
	    if (showOutput) 
	      {
		if (iter == 0 && expertIndex == 0)
		  {
		    std::cout << "STAPLE: [iteration, expert index, "
			      << "p, q, max w change]" << std::endl;
		  }
		std::cout << "STAPLE: [" << iter << ", " << expertIndex << ", "
			  << p[expertIndex] << ", "
			  << q[expertIndex] << ", "
			  << maxCurrentWChange << "]" << std::endl;
	      }
	  }
	
	//
	// check for minWChange convergence
	//
	if (maxCurrentWChange < minWChange)
	  {
	    std::cout << "STAPLE: Stopping: maxCurrentWChange < minWChange"
		      << std::endl;
	    break;	  
	  }
      }
    
    //
    // fill probability image
    //
    if (restrictToDisputedVoxels)
      {
	_replaceDisputedVoxels(w,
			       disputedVoxelLocations,
			       probabilisticSegmentation);      
      }
    else
      {
	_replaceAllVoxels(w,probabilisticSegmentation);
      }

    //
    // clean up memory
    //
    delete [] sumExpertForegroundProbabilities;
    delete [] sumExpertBackgroundProbabilities;
  }
  
private:
  static
  void
  _getAllVoxels(const unsigned int numExperts,
		LabelImageType** expertSegmentations,
		std::vector<std::vector<LabelType> >& voxels)
  {
    voxels.resize(numExperts);
    unsigned int numVoxels = (*expertSegmentations[0]).getNumElements();
    for (unsigned int i = 0; i < numExperts; ++i)
      {
	voxels[i].resize(numVoxels);
	memcpy(&voxels[i][0],
	       (*expertSegmentations[i]).getDataPointer(),
	       (*expertSegmentations[i]).getSizeBytes());
      }
  }

  static
  void
  _getDisputedVoxels(const unsigned int numExperts,
		     LabelImageType** expertSegmentations,
		     std::vector<std::vector<LabelType> >& disputedVoxels,
		     std::vector<Vector3D<unsigned int> >&
		     disputedVoxelLocations,
		     ProbabilityImageType& probabilisticSegmentation)
  {
    disputedVoxelLocations.clear();
    disputedVoxels.clear();
    disputedVoxels.resize(numExperts);
    
    Vector3D<unsigned int> size = expertSegmentations[0]->getSize();
    for (unsigned int z = 0; z < size.z; ++z) {
      for (unsigned int y = 0; y < size.y; ++y) {
	for (unsigned int x = 0; x < size.x; ++x) {
	  LabelType zeroLabel = (*expertSegmentations[0])(x,y,z);
	  bool disputed = false;
	  for (unsigned int i = 1; i < numExperts; ++i) 
	    {
	      if ((*expertSegmentations[i])(x,y,z) != zeroLabel)
		{
		  disputed = true;
		  break;
		}
	    }
	  if (disputed)
	    {
	      disputedVoxelLocations.push_back(Vector3D<unsigned int>(x,y,z));
	      for (unsigned int i = 0; i < numExperts; ++i) {
		disputedVoxels[i].push_back((*expertSegmentations[i])(x,y,z));
	      }	    
	    }
	  else 
	    {
	      probabilisticSegmentation(x,y,z) = (zeroLabel > 0 ? 1 : 0);
	    }
	}
      }
    }
  }

  static
  void
  _replaceAllVoxels(std::vector<ProbabilityType> w,
		    ProbabilityImageType& image)
  {
    memcpy(image.getDataPointer(), &w[0], image.getSizeBytes());
  }

  static
  void
  _replaceDisputedVoxels(std::vector<ProbabilityType> w,
			 std::vector<Vector3D<unsigned int> >& 
			 disputedVoxelLocations,
			 ProbabilityImageType& image)
  {
    unsigned int numDisputedVoxels = disputedVoxelLocations.size();
    for (unsigned int i = 0; i < numDisputedVoxels; ++i)
      {
	image(disputedVoxelLocations[i]) = w[i];
      }
  }
};
#endif
