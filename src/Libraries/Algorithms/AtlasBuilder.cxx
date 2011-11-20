// todo: 
// save memory in update h
// trade off scratch: per thread or per image, what whatever is smaller

#include <AtlasBuilder.h>
#include <strstream>
#include <HField3DUtils.h>
#include <assert.h>
#include <Timer.h>
#include <iomanip>
#define _USE_MATH_DEFINES
#include <math.h>


struct ThreadInfo
{
  AtlasBuilder* atlasBuilder;
  unsigned int  threadIndex;
};

AtlasBuilder::
AtlasBuilder()
  : 
  _IterationDataOutputStream(std::cerr)
{
  this->_NumberOfImages               = 0;
  this->_AverageImagePointer          = NULL;
  this->_NumberOfThreads              = 1;

  this->_FFTWNumberOfThreads          = 1;
  this->_FFTWMeasure                  = true;

  this->_DeltaSelectionMethod         = DELTA_USE_INDIVIDUAL;
  this->_UpdateAfterEverySubIteration = true;
  this->_ComputeInverseDeformations   = false;
  this->_MeanComputationStrategy      = 
    new ArithmeticMeanComputationStrategy<AtlasBuilder::VoxelType>();

  this->_MeanSquaredError             = 0;
  this->_Iteration                    = 0;
  this->_NextImageToProcess           = 0;

  pthread_mutex_init(&this->_NextImageToProcessMutex, NULL);
  pthread_mutex_init(&this->_AverageImageMutex, NULL);
}

AtlasBuilder::
~AtlasBuilder()
{
  this->DeleteScratchMemory();
}

void
AtlasBuilder::
SetLogOutputStream(std::ostream& out)
{
  //this->_IterationDataOutputStream = out;
}

void
AtlasBuilder::
SetNumberOfThreads(unsigned int numThreads)
{
  this->_NumberOfThreads = numThreads;
}

unsigned int 
AtlasBuilder::
GetNumberOfThreads() const
{
  return this->_NumberOfThreads;
}

void
AtlasBuilder::
SetFFTWNumberOfThreads(unsigned int numThreads)
{
  this->_FFTWNumberOfThreads = numThreads;
}

unsigned int 
AtlasBuilder::
GetFFTWNumberOfThreads() const
{
  return this->_FFTWNumberOfThreads;
}

void 
AtlasBuilder::
SetFFTWMeasure(bool b)
{
  this->_FFTWMeasure = b;
}

void 
AtlasBuilder::
SetFFTWMeasureOn()
{
  this->SetFFTWMeasure(true);
}

void 
AtlasBuilder::
SetFFTWMeasureOff()
{
  this->SetFFTWMeasure(false);
}

bool 
AtlasBuilder::
GetFFTWMeasure() const
{
  return this->_FFTWMeasure;
}

void 
AtlasBuilder::
SetUpdateAverageEverySubIterationOn()
{
  this->SetUpdateAverageEverySubIteration(true);
}

void 
AtlasBuilder::
SetUpdateAverageEverySubIterationOff()
{
  this->SetUpdateAverageEverySubIteration(false);
}

void 
AtlasBuilder::
SetUpdateAverageEverySubIteration(bool b)
{
  this->_UpdateAfterEverySubIteration = b;
}

bool
AtlasBuilder::
GetUpdateAverageEverySubIteration() const
{
  return this->_UpdateAfterEverySubIteration;
}

void 
AtlasBuilder::
SetComputeInverseDeformationsOn()
{
  this->SetComputeInverseDeformations(true);
}

void 
AtlasBuilder::
SetComputeInverseDeformationsOff()
{
  this->SetComputeInverseDeformations(false);
}

void 
AtlasBuilder::
SetComputeInverseDeformations(bool b)
{
  this->_ComputeInverseDeformations = b;
}

bool
AtlasBuilder::
GetComputeInverseDeformations() const
{
  return this->_ComputeInverseDeformations;
}

void
AtlasBuilder::
SetFluidWarpParameters(const FluidWarpParameters& fluidParams)
{
  this->_FluidWarpParameters = fluidParams;
}

FluidWarpParameters&
AtlasBuilder::
GetFluidWarpParameters()
{
  return this->_FluidWarpParameters;
}

void
AtlasBuilder::
SetMeanComputationStrategy(AtlasBuilder::MeanComputationStrategyType* s)
{
  this->_MeanComputationStrategy = s;
}

AtlasBuilder::MeanComputationStrategyType*
AtlasBuilder::
GetMeanComputationStrategy() const
{
  return this->_MeanComputationStrategy;
}

void 
AtlasBuilder::
SetNumberOfInputImages(unsigned int n)
{
  this->_NumberOfImages = n;

  this->_ImagePointers.resize(n);
  this->_DeformationFieldPointers.resize(n);
  this->_DeformationFieldInversePointers.resize(n);
  this->_Delta.assign(n, 0.0);
}

unsigned int
AtlasBuilder::
GetNumberOfInputImages() const
{
  return this->_NumberOfImages;
}

void
AtlasBuilder::
SetNthInputImage(unsigned int n, ImageType* imagePointer)
{
  this->_ImagePointers[n] = imagePointer;
}

AtlasBuilder::ImageType*
AtlasBuilder::
GetNthInputImage(unsigned int n) const
{
  return this->_ImagePointers[n];
}

AtlasBuilder::ImageType*
AtlasBuilder::
GetNthDeformedImage(unsigned int n) const
{
  if (this->_DeformedImagePointers.size() <= n)
    {
      return NULL;
    }
  return this->_DeformedImagePointers[n];
}

void
AtlasBuilder::
SetNthDeformationField(unsigned int n, VectorFieldType* fieldPointer)
{
  this->_DeformationFieldPointers[n] = fieldPointer;
}

AtlasBuilder::VectorFieldType*
AtlasBuilder::
GetNthDeformationField(unsigned int n) const
{
  return this->_DeformationFieldPointers[n];
}

void
AtlasBuilder::
SetNthDeformationFieldInverse(unsigned int n, VectorFieldType* fieldPointer)
{
  this->_DeformationFieldInversePointers[n] = fieldPointer;
}

AtlasBuilder::VectorFieldType*
AtlasBuilder::
GetNthDeformationFieldInverse(unsigned int n) const
{
  return this->_DeformationFieldInversePointers[n];
}

void
AtlasBuilder::
SetAverageImage(ImageType* imagePointer)
{
  this->_AverageImagePointer = imagePointer;
}

AtlasBuilder::ImageType*
AtlasBuilder::
GetAverageImage() const
{
  return this->_AverageImagePointer;
}

AtlasBuilder::IterationData
AtlasBuilder::
GetIterationData(unsigned int iteration) const
{
  return this->_IterationDataLog[iteration];
}

void
AtlasBuilder::
CheckInputs()
{
  //
  // make sure that a valid mean strategy is set
  //
  if (this->_MeanComputationStrategy == NULL)
  {
    std::runtime_error("Mean computation strategy is null.");    
  }

  //
  // make sure that the mean image is there there
  //
  if (this->_AverageImagePointer == NULL)
  {
    std::runtime_error("Average Image is null.");
  }

  // use the size of the average as the standard
  Vector3D<unsigned int> imageSize = this->_AverageImagePointer->getSize();

  //
  // check the data for all the images
  //
  for (unsigned int i = 0; i < this->_NumberOfImages; ++i)
  {
    // check that the image is there
    if (this->_ImagePointers[i] == NULL)
    {
      std::strstream ss;
      ss << "Input image is null: " << i;
      throw std::runtime_error(ss.str());
    }

    // check that is is the right size
    if (imageSize != this->_ImagePointers[i]->getSize())
    {
      std::strstream ss;
      ss << "Incompatible image size for image " << i;
      throw std::runtime_error(ss.str());
    }

    // check that the deformation field is there
    if (this->_DeformationFieldPointers[i] == NULL)
    {
      std::strstream ss;
      ss << "Deformation field is null: " << i;
      throw std::runtime_error(ss.str());
    }

    // check that the deformation field is the correct size
    if (imageSize != this->_DeformationFieldPointers[i]->getSize())
    {
      std::strstream ss;
      ss << "Incompatible deformation field size for image " << i;
      throw std::runtime_error(ss.str());
    }

    // check that the inverse deformation field is there
    if (this->_ComputeInverseDeformations &&
        this->_DeformationFieldInversePointers[i] == NULL)
    {
      std::strstream ss;
      ss << "Inverse deformation field is null: " << i;
      throw std::runtime_error(ss.str());
    }

    // check inverse deformation field is the correct size
    if (this->_ComputeInverseDeformations &&
        imageSize != this->_DeformationFieldInversePointers[i]->getSize())
    {
      std::strstream ss;
      ss << "Incompatible inverse deformation field size for image " << i;
      throw std::runtime_error(ss.str());
    }
  }
}

void
AtlasBuilder::
InitializeScratchMemory()
{
  //
  // this function assumes that CheckInput has been called without error.
  //

  //
  // resize memory holders
  this->_DeformedImagePointers.resize(this->_NumberOfImages);
  this->_MaxL2Displacements.resize(this->_NumberOfImages, 0.0);
  this->_MaxVelocityL2Displacements.resize(this->_NumberOfImages, 0.0);

  this->_ScratchVectorFieldPointers.resize(this->_NumberOfThreads);
  this->_FFTWForwardPlans.resize(this->_NumberOfThreads);
  this->_FFTWBackwardPlans.resize(this->_NumberOfThreads);


  // use the size of the average as the standard
  Vector3D<unsigned int> imageSize = this->_AverageImagePointer->getSize();

  // this multipurpose array is used to hold gradient, body force, and
  // velocity arrays (in that order) this saves LOTS of memory
  //
  // VERY IMPORTANT: however, because it is used to hold the result of an fft
  // it needs to be padded just a little bit, see www.fftw.org
  // this wont affect its use as gradient and body force arrays
  // as long as logical size (image size) is used for this array
  // and access into this array is done **ONLY** via the (x, y, z) operator
  // and not via incremented pointers.  I REPEAT, dont access the vector
  // field via incremented pointers unless you *know* what you are doing.
  unsigned xSizeFFT = 2 * (imageSize.x / 2 + 1);

  for (unsigned int i = 0; i < this->_NumberOfThreads; ++i)
  {
    // allocate scratch vector field
    this->_ScratchVectorFieldPointers[i] = 
      new VectorFieldType(xSizeFFT, imageSize.y, imageSize.z);
  }

  for (unsigned int i = 0; i < this->_NumberOfImages; ++i)
  {
    // allocate deformed image
    this->_DeformedImagePointers[i] = new ImageType(imageSize);

    // generate the deformed image
    HField3DUtils::apply(*this->_ImagePointers[i], 
                         *this->_DeformationFieldPointers[i],
                         *this->_DeformedImagePointers[i]);
  }
}

void
AtlasBuilder::
DeleteScratchMemory()
{
  // order is important
  this->DeleteFFTWPlans();

  for (unsigned int i = 0; i < this->_NumberOfThreads; ++i)
  {
    if (this->_ScratchVectorFieldPointers[i])
    {
      delete this->_ScratchVectorFieldPointers[i];
      this->_ScratchVectorFieldPointers[i] = NULL;
    }
  }

  for (unsigned int i = 0; i < this->_NumberOfImages; ++i)
  {
    if (this->_DeformedImagePointers[i])
    {
      delete this->_DeformedImagePointers[i];
      this->_DeformedImagePointers[i] = NULL;
    }
  }
}

void
AtlasBuilder::
InitializeFFTWPlans()
{
 assert(this->_AverageImagePointer != NULL);

  // use the size of the average as the standard
  Vector3D<unsigned int> imageSize = this->_AverageImagePointer->getSize();  

  for (unsigned int i = 0; i < this->_NumberOfThreads; ++i)
  {
    int rank = 3;
    int logicalSizeParam[3];
    logicalSizeParam[0] = imageSize.z;
    logicalSizeParam[1] = imageSize.y;
    logicalSizeParam[2] = imageSize.x;

    int howMany = 3;
    int stride  = 3;
    int dist    = 1;
    float *dataPtr = 
      (CoordinateType*) this->_ScratchVectorFieldPointers[i]->getDataPointer();
    
    fftwf_plan_with_nthreads(this->_FFTWNumberOfThreads);
    
    this->_FFTWForwardPlans[i] = 
      fftwf_plan_many_dft_r2c(rank, logicalSizeParam, howMany, 
                              dataPtr, 
                              0, stride, dist, 
                              (fftwf_complex*) dataPtr, 
                              0, stride, dist,
                              this->_FFTWMeasure ? FFTW_MEASURE : 
                              FFTW_ESTIMATE);

    this->_FFTWBackwardPlans[i] = 
      fftwf_plan_many_dft_c2r(rank, logicalSizeParam, howMany, 
                              (fftwf_complex*) dataPtr,
                              0, stride, dist, 
                              dataPtr,
                              0, stride, dist,
                              this->_FFTWMeasure ? FFTW_MEASURE : 
                              FFTW_ESTIMATE);
    
    if (!_FFTWForwardPlans[i])
    {
      throw std::runtime_error("FFTW forward plan failed to initialize");
    }
    if (!_FFTWBackwardPlans[i])
    {
      throw std::runtime_error("FFTW backward plan failed to initialize");
    }  
  }
}
  
void
AtlasBuilder::
DeleteFFTWPlans()
{
  for (unsigned int i = 0; i < this->_NumberOfThreads; ++i)
  {
    fftwf_destroy_plan(this->_FFTWForwardPlans[i]);
    fftwf_destroy_plan(this->_FFTWBackwardPlans[i]);
  }
}

void 
AtlasBuilder::
InitializeOperatorLookupTable()
{
  Vector3D<unsigned int> imageSize = this->_AverageImagePointer->getSize();
  AtlasBuilder::LUT lut(imageSize.x, imageSize.y, imageSize.z);

  // hardcode these for now
  double deltaX = 1.0;
  double deltaY = 1.0;
  double deltaZ = 1.0;

  //
  // precompute some values
  //

  double sX = deltaX * 2.0 * M_PI / imageSize.x; 
  double sY = deltaY * 2.0 * M_PI / imageSize.y; 
  double sZ = deltaZ * 2.0 * M_PI / imageSize.z; 

  double deltaXSq = deltaX * deltaX;
  double deltaYSq = deltaY * deltaY;
  double deltaZSq = deltaZ * deltaZ;

  //
  // fill in luts
  //
  for (unsigned int x = 0; x < lut.cosWX.size(); ++x) 
    {
      lut.cosWX[x] = (2.0 * cos(sX * static_cast<float>(x)) - 2.0) / deltaXSq;
      lut.sinWX[x] = sin(sX * static_cast<float>(x)) / deltaX;
    }
  for (unsigned int y = 0; y < lut.cosWY.size(); ++y)
    {
      lut.cosWY[y] = (2.0 * cos(sY * static_cast<float>(y)) - 2.0) / deltaYSq;
      lut.sinWY[y] = sin(sY * static_cast<float>(y)) / deltaY;
    }
  for (unsigned int z = 0; z < lut.cosWZ.size(); ++z)
    {
      lut.cosWZ[z] = (2.0 * cos(sZ * static_cast<float>(z)) - 2.0) / deltaZSq;
      lut.sinWZ[z] = sin(sZ * static_cast<float>(z)) / deltaZ;
    }  

  for (unsigned int x = 0; x < lut.cosWX.size(); ++x)
    for (unsigned int y = 0; y < lut.cosWY.size(); ++y)
      for (unsigned int z = 0; z < lut.cosWZ.size(); ++z)
        lut.nsq(x,y,z) = lut.sinWX[x]*lut.sinWX[x]
          + lut.sinWY[y]*lut.sinWY[y]
          + lut.sinWZ[z]*lut.sinWZ[z]; // norm squared of projection
                                       // vector for incompressible flow

  //
  // copy values to the ivar
  //
  this->_OperatorLookupTable = lut;
}

void
AtlasBuilder::
GenerateAverage()
{
  //
  // start the total timer
  this->_TotalTimer.restart();

  //
  // check that input data is consistent: an exception will be thrown
  // if there is a problem
  std::cerr << "Checking inputs...";
  this->CheckInputs();
  std::cerr << "DONE" << std::endl;  

  //
  // set num threads to min(num threads, num images)
  if (this->_NumberOfThreads > this->_NumberOfImages)
  {
    std::cerr << "WARNING: More threads than images." << std::endl
              << "WARNING: Setting number of threads to " 
              << this->_NumberOfImages << std::endl;
    this->_NumberOfThreads = this->_NumberOfImages;
  }

  //
  // initialize memory: order is important
  std::cerr << "Initializing memory...";
  std::cerr << "Scratch memory...";
  this->InitializeScratchMemory();
  std::cerr << "FFTW Plans...";
  this->InitializeFFTWPlans();
  std::cerr << "LUT...";
  this->InitializeOperatorLookupTable();
  std::cerr << "DONE" << std::endl;

  //
  // start algorithm
  //
  this->RunAlgorithm();
}

unsigned int 
AtlasBuilder::
GetJobImageIndex()
{
  pthread_mutex_lock(&this->_NextImageToProcessMutex);
  unsigned int imageIndex = this->_NextImageToProcess++;
  pthread_mutex_unlock(&this->_NextImageToProcessMutex);
  return imageIndex;
}

void 
AtlasBuilder::
LockMeanImage()
{
  pthread_mutex_lock(&this->_AverageImageMutex);
}

void 
AtlasBuilder::
UnlockMeanImage()
{
  pthread_mutex_unlock(&this->_AverageImageMutex);
}

pthread_t
AtlasBuilder::
GetThreadID()
{
  return pthread_self();
}

void*
AtlasBuilder::
ThreadedUpdateImages(void* arg)
{
  ThreadInfo* threadInfoPtr = static_cast<ThreadInfo*>(arg);
  AtlasBuilder* ptr = threadInfoPtr->atlasBuilder;
  unsigned int threadIndex = threadInfoPtr->threadIndex;

  // loop until all jobs are taken
  while (true)
  {
    //
    // get an image index to process
    unsigned int imageIndex = ptr->GetJobImageIndex();
    //std::cerr << "{" << threadIndex << " / " << ptr->GetThreadID()
    //<< "} Got assignment: " 
    //<< imageIndex << std::endl;    
    if (imageIndex >= ptr->_NumberOfImages)
    {
      // all the jobs are taken
      break;
    }

    //
    // start the timer
    Timer iterationTimer;
    iterationTimer.start();

    //
    // compute gradient of deformed image
    //std::cerr << "{" << threadIndex << " / " << ptr->GetThreadID()
    //<< "} Gradient: " 
    //<< imageIndex << std::endl;    
    ptr->UpdateGradient(imageIndex, threadIndex);

    // 
    // compute body force 
    //
    // Note: the body force computation requires the average image,
    // make sure we don't access it while it is being updated
    if (ptr->_UpdateAfterEverySubIteration) ptr->LockMeanImage();
    //std::cerr << "{" << threadIndex << " / " << ptr->GetThreadID()
    //<< "} Body force: " 
    //<< imageIndex << std::endl;    
    
    ptr->UpdateBodyForce(imageIndex, threadIndex);

    if (ptr->_UpdateAfterEverySubIteration) ptr->UnlockMeanImage();


	
    //
    // update velocity field according to Euler-Lagrange equation
    //std::cerr << "{" << threadIndex << " / " << ptr->GetThreadID()
    //<< "} Velocity: " 
    //<< imageIndex << std::endl;    
    ptr->UpdateVelocityField(imageIndex, threadIndex);

    //
    // compose velocity to generate new deformation field
    //std::cerr << "{" << threadIndex << " / " << ptr->GetThreadID()
    //<< "} Deformation fields: " 
    //<< imageIndex << std::endl;    
    ptr->UpdateDeformationFields(imageIndex, threadIndex);

    //
    // measure maximum displacement
    //std::cerr << "{" << threadIndex << " / " << ptr->GetThreadID()
    //<< "} Min/max L2: " 
    //<< imageIndex << std::endl;    
    double minL2, maxL2;
    HField3DUtils::
      minMaxDeformationL2Norm(*ptr->_DeformationFieldPointers[imageIndex], 
                              minL2, maxL2);
    ptr->_MaxL2Displacements[imageIndex] = maxL2;
    //std::cerr << "{" << threadIndex << " / " << imageIndex
    //<< "} Min/max L2: " << maxL2 << std::endl;    
    
    //
    // deform the image according to the new deformation field
    //std::cerr << "{" << threadIndex << " / " << ptr->GetThreadID()
    //<< "} Deformed image: " 
    //<< imageIndex << std::endl;    
    ptr->UpdateDeformedImage(imageIndex, threadIndex);

    //
    // if the average should be updated after each sub-iteration, do
    // that now and report the results
    //
    if (ptr->_UpdateAfterEverySubIteration)
    {
      ptr->LockMeanImage();
      //std::cerr << "{" << threadIndex << " / " << ptr->GetThreadID()
      //<< "} Average: " 
      //<< imageIndex << std::endl;    
      ptr->UpdateAverageImage();
      //std::cerr << "{" << threadIndex << " / " << ptr->GetThreadID()
      //<< "} Error: " 
      //<< imageIndex << std::endl;    
      ptr->UpdateError();
      //std::cerr << "{" << threadIndex << " / " << ptr->GetThreadID()
      //<< "} Log: " 
      //<< imageIndex << std::endl;    
      ptr->LogIterationData(imageIndex, threadIndex, 
                             iterationTimer.getSeconds());
      ptr->UnlockMeanImage();
    }
  }  
 return 0;
}

void
AtlasBuilder::
UpdateDelta()
{
  if (this->_DeltaSelectionMethod == DELTA_USE_MEAN)
  {
    double l2Sum = std::accumulate(this->_MaxVelocityL2Displacements.begin(),
                                   this->_MaxVelocityL2Displacements.end(),
                                   0.0);
    this->_Delta.assign(this->_NumberOfImages, 
                        this->_FluidWarpParameters.maxPerturbation /
                        (this->_NumberOfImages * l2Sum));
  }
  else if (this->_DeltaSelectionMethod == DELTA_USE_INDIVIDUAL)
  {
    this->_Delta.resize(this->_NumberOfImages);
    for (unsigned int i = 0; i < this->_NumberOfImages; ++i)
    {
		if (this->_MaxVelocityL2Displacements[i]>0){
      this->_Delta[i] = this->_FluidWarpParameters.maxPerturbation / 
		  this->_MaxVelocityL2Displacements[i];}
		else {
			 this->_Delta[i] =0.0;
		}
    }
  }
  else
  {
    std::runtime_error("Unknown delta update method.");    
  }
}

void 
AtlasBuilder::
UpdateGradient(unsigned int imageIndex, unsigned int threadIndex)
{
  Array3DUtils::
    computeGradient(*this->_DeformedImagePointers[imageIndex],
                    *this->_ScratchVectorFieldPointers[threadIndex]);
}

void 
AtlasBuilder::
UpdateBodyForce(unsigned int imageIndex, unsigned int threadIndex)
{
  Vector3D<unsigned int> size = this->_AverageImagePointer->getSize();
  double di;
  for (unsigned int z = 0; z < size.z; ++z) {
    for (unsigned int y = 0; y < size.y; ++y) {
      for (unsigned int x = 0; x < size.x; ++x) {
        di = 
          (*this->_AverageImagePointer)(x, y, z) - 
          (*this->_DeformedImagePointers[imageIndex])(x, y, z);
        (*this->_ScratchVectorFieldPointers[threadIndex])(x,y,z) *= di;
      }
    }
  }    
}

void 
AtlasBuilder::
UpdateVelocityField(unsigned int imageIndex, unsigned int threadIndex)
{
  Vector3D<unsigned int> logicalSize = this->_AverageImagePointer->getSize();

  // forward fft (scale array, then compute fft)
  this->_ScratchVectorFieldPointers[threadIndex]->
    scale(1.0 / logicalSize.productOfElements());
  fftwf_execute(this->_FFTWForwardPlans[threadIndex]);

  // apply operator
  double lambda;
  float L00;
  float L10, L11;
  float L20, L21, L22;
  double alpha = this->_FluidWarpParameters.alpha;
  double beta  = this->_FluidWarpParameters.beta;
  double gamma = this->_FluidWarpParameters.gamma;

  unsigned int xFFTMax = logicalSize.x / 2 + 1;
  for (unsigned int z = 0; z < logicalSize.z; ++z)
    {
      for (unsigned int y = 0; y < logicalSize.y; ++y)
	{
	  for (unsigned int x = 0; x < xFFTMax; ++x)
	    {
	      //
	      // compute L (it is symmetric, only need lower triangular part)
	      //
	      
	      // maybe lambda should be stored in a lut
	      // it would reduce computation but may cause cache misses
	      lambda = - alpha 
		* (this->_OperatorLookupTable.cosWX[x] + 
                   this->_OperatorLookupTable.cosWY[y] + 
                   this->_OperatorLookupTable.cosWZ[z]) 
		+ gamma;	      
	      
	      L00 = lambda - beta * this->_OperatorLookupTable.cosWX[x];
	      L11 = lambda - beta * this->_OperatorLookupTable.cosWY[y];
	      L22 = lambda - beta * this->_OperatorLookupTable.cosWZ[z];
	      L10 = beta * this->_OperatorLookupTable.sinWX[x] * 
                this->_OperatorLookupTable.sinWY[y];
	      L20 = beta * this->_OperatorLookupTable.sinWX[x] * 
                this->_OperatorLookupTable.sinWZ[z];
	      L21 = beta * this->_OperatorLookupTable.sinWY[y] * 
                this->_OperatorLookupTable.sinWZ[z];

	      //
	      // compute V = Linv F (for real and imaginary parts)
	      //
              CoordinateType* complexPtr =
                &(*this->
                  _ScratchVectorFieldPointers[threadIndex])(x * 2, y, z).x; 
	      AtlasBuilder::
                InverseOperatorMultiply(complexPtr,
                                        L00,
                                        L10, L11,
                                        L20, L21, L22);

              if (this->_FluidWarpParameters.divergenceFree && x ^ y ^ z != 0)
                {
                // Project onto incompressible field
                AtlasBuilder::ProjectIncomp(complexPtr,x,y,z);
                }
	    }
	}
    }

  // backward fft
  fftwf_execute(this->_FFTWBackwardPlans[threadIndex]);

  //
  // compute max L2 norm for delta computation
  double maxL2 = 0;
  for (unsigned int z = 0; z < logicalSize.z; ++z)
  {
    for (unsigned int y = 0; y < logicalSize.y; ++y)
    {
      for (unsigned int x = 0; x < logicalSize.x; ++x)
      {
        double L2sq = 
          (*this->_ScratchVectorFieldPointers[threadIndex])(x,y,z).
          lengthSquared();
        if (maxL2 < L2sq)
        {
          maxL2 = L2sq;
        }
      }
    }
  }
  this->_MaxVelocityL2Displacements[imageIndex] = sqrt(maxL2);  
}

void
AtlasBuilder::
ProjectIncomp(CoordinateType* complexPtr, unsigned int x, unsigned int y, unsigned int z)
{
  // in Fourier space we project onto (-i*sin(u),-i*sin(v),-i*sin(w)) and remove that component
  // 2008 jdh

  float bRealX = complexPtr[0];
  float bRealY = complexPtr[2];
  float bRealZ = complexPtr[4];

  float bImagX = complexPtr[1];
  float bImagY = complexPtr[3];
  float bImagZ = complexPtr[5];

  float& vRealX = complexPtr[0];
  float& vRealY = complexPtr[2];
  float& vRealZ = complexPtr[4];

  float& vImagX = complexPtr[1];
  float& vImagY = complexPtr[3];
  float& vImagZ = complexPtr[5];

  AtlasBuilder::LUT *lut = &this->_OperatorLookupTable;

  // This is now in LUT
  //   nsq = lut->sinWX[x]*lut->sinWX[x]
  //       + lut->sinWY[y]*lut->sinWY[y]
  //       + lut->sinWZ[z]*lut->sinWZ[z]; // norm squared of projection vector

  // S=(sinwx,sinwy,sinwz)
  // Real part of S dot V in Fourier
  double ReSdotV = ( bRealX*lut->sinWX[x]
                     +bRealY*lut->sinWY[y]
                     +bRealZ*lut->sinWZ[z]);
  // Imag part of S dot V in Fourier
  double ImSdotV = ( bImagX*lut->sinWX[x]
                     +bImagY*lut->sinWY[y]
                     +bImagZ*lut->sinWZ[z]);
      
  // Subtract S dot V (normalizing S)
  vRealX = bRealX - ReSdotV*lut->sinWX[x]/lut->nsq(x,y,z);
  vRealY = bRealY - ReSdotV*lut->sinWY[y]/lut->nsq(x,y,z);
  vRealZ = bRealZ - ReSdotV*lut->sinWZ[z]/lut->nsq(x,y,z);

  vImagX = bImagX - ImSdotV*lut->sinWX[x]/lut->nsq(x,y,z);
  vImagY = bImagY - ImSdotV*lut->sinWY[y]/lut->nsq(x,y,z);
  vImagZ = bImagZ - ImSdotV*lut->sinWZ[z]/lut->nsq(x,y,z);
}

// TODO: we could probably rearrange this to save some memory...
void 
AtlasBuilder::
UpdateDeformationFields(unsigned int imageIndex, unsigned int threadIndex)
{
  Vector3D<unsigned int> size = this->_AverageImagePointer->getSize();

//add by xx, 2010
    //adaptively change the delta, to make sure each step is small enough, prevent divergence near the optimum
	if(this->_MaxVelocityL2Displacements[imageIndex]>0){
		  if (this->_MaxVelocityL2Displacements[imageIndex] 
		                      < this->_FluidWarpParameters.maxPerturbation)
		  {
			  //converged, maximum delta=1.0
			  this->_Delta[imageIndex]= 1.0;
		  }else{
		      this->_Delta[imageIndex]= this->_FluidWarpParameters.maxPerturbation
		                                        / this->_MaxVelocityL2Displacements[imageIndex];
		  }
	}else{
		this->_Delta[imageIndex]=0.0;
	}


  //
  // compute hIncremental(x) = x + velocity(x) * delta
  //

  VectorFieldType hIncremental(size);
  for (unsigned int z = 0; z < size.z; ++z) {
    for (unsigned int y = 0; y < size.y; ++y) {
      for (unsigned int x = 0; x < size.x; ++x) {
	hIncremental(x,y,z).x = x + 
          (*this->_ScratchVectorFieldPointers[threadIndex])(x,y,z).x * 
          this->_Delta[imageIndex];
	hIncremental(x,y,z).y = y + 
          (*this->_ScratchVectorFieldPointers[threadIndex])(x,y,z).y * 
          this->_Delta[imageIndex];
	hIncremental(x,y,z).z = z + 
          (*this->_ScratchVectorFieldPointers[threadIndex])(x,y,z).z * 
          this->_Delta[imageIndex];

      }
    }
  }

  //
  // compute h(x) = h(hIncremental(x))
  //
  VectorFieldType oldDeformation(*this->_DeformationFieldPointers[imageIndex]);
  HField3DUtils::compose(oldDeformation, hIncremental, 
                         *this->_DeformationFieldPointers[imageIndex]);

  // update inverse deformation field
  if (this->_ComputeInverseDeformations)
  {
    //
    // compute 
    // hIncrementalInv(x) = x + x - hIncremental(x)
    //
    VectorFieldType& hIncrementalInv = oldDeformation; // reuse memory here
    HField3DUtils::computeInverseZerothOrder(hIncremental, hIncrementalInv);
    
    //
    // compute hInv(x) = hIncrementalInv(hInv(x))
    //
    HField3DUtils::
      compose(hIncrementalInv, 
              *this->_DeformationFieldInversePointers[imageIndex], 
              *this->_DeformationFieldInversePointers[imageIndex]);
  }
}

void 
AtlasBuilder::
UpdateDeformedImage(unsigned int imageIndex, unsigned int threadIndex)
{
  HField3DUtils::apply(*this->_ImagePointers[imageIndex],
                       *this->_DeformationFieldPointers[imageIndex],
                       *this->_DeformedImagePointers[imageIndex]);
}

void
AtlasBuilder::
RunAlgorithm()
{
  //
  // compute average image
  this->ThreadedUpdateAverageImage();

  for (unsigned int iter = 0; iter < this->_FluidWarpParameters.numIterations;
       ++iter)
  {
    Timer iterationTimer;
    iterationTimer.start();

    this->_Iteration = iter;

	 
	 //
    // start threads which will update each image
    std::vector<pthread_t> threads(this->_NumberOfThreads);

    // reset job counter
    this->_NextImageToProcess = 0;

    // info passed to threads
    std::vector<ThreadInfo> threadInfo(this->_NumberOfThreads);

    // create threads that will update images
    for (unsigned int threadIndex = 0; threadIndex < this->_NumberOfThreads; 
         ++threadIndex)
    {
      threadInfo[threadIndex].atlasBuilder = this;
      threadInfo[threadIndex].threadIndex  = threadIndex;
  
      int rv = pthread_create(&threads[threadIndex], NULL,
                              &AtlasBuilder::ThreadedUpdateImages, 
                              &threadInfo[threadIndex]);
      if (rv != 0)
      {
        throw std::runtime_error("Error creating thread.");
      }
    }

    // join threads
    for (int threadIndex = 0; threadIndex < this->_NumberOfThreads; 
         ++threadIndex)
    {
      int rv = pthread_join(threads[threadIndex], NULL);
      if (rv != 0)
      {
        throw std::runtime_error("Error joining thread.");
      }
    }

    //
    // delta will will be average deltas for individual images.  A
    // delta for an individual image means that maxpert will be
    // realized for that image.
    if (this->_Iteration == 0)
    {
      this->UpdateDelta();
    }

    //
    // update the average if it is only done once per iteration.
    // otherwise, it will be updated by the threads.
    if (this->_UpdateAfterEverySubIteration == false)
    {
      this->ThreadedUpdateAverageImage();
      this->UpdateError();
      this->LogIterationData(this->_NumberOfImages-1, 0,
                             iterationTimer.getSeconds());
    }
  }
}

void
AtlasBuilder::
UpdateAverageImage()
{
  unsigned int numElements = this->_AverageImagePointer->getNumElements();
  std::vector<VoxelType> voxelData(this->_NumberOfImages);

  for (unsigned int i = 0; i < numElements; ++i)
  {
    // get data from each image for this voxel
    for (unsigned int j = 0; j < this->_NumberOfImages; ++j)
    {
      voxelData[j] = (*this->_DeformedImagePointers[j])(i);
    }

    // compute the mean from the list of voxels
    (*this->_AverageImagePointer)(i) = 
      this->_MeanComputationStrategy->ComputeMean(this->_NumberOfImages,
                                                  &voxelData[0]);
  }
}

void*
AtlasBuilder::
ThreadedUpdateAverageImage(void* arg)
{
  ThreadInfo* threadInfoPtr = static_cast<ThreadInfo*>(arg);
  AtlasBuilder* ptr = threadInfoPtr->atlasBuilder;
  unsigned int threadIndex = threadInfoPtr->threadIndex;
  unsigned int numberOfThreads  = ptr->_NumberOfThreads;
  unsigned int numberOfElements = ptr->_AverageImagePointer->getNumElements();
  unsigned int thisThreadBegin  = 
    threadIndex * (numberOfElements/numberOfThreads);
  unsigned int thisThreadEnd    = 
    (threadIndex + 1) * (numberOfElements/numberOfThreads);
  
  if (threadIndex == numberOfThreads-1)
  {
    thisThreadEnd = numberOfElements;
  }

  std::vector<VoxelType> voxelData(ptr->_NumberOfImages);
  for (unsigned int i = thisThreadBegin; i < thisThreadEnd; ++i)
  {
    // get data from each image for this voxel
    for (unsigned int j = 0; j < ptr->_NumberOfImages; ++j)
    {
      voxelData[j] = (*ptr->_DeformedImagePointers[j])(i);
    }

    // compute the mean from the list of voxels
    (*ptr->_AverageImagePointer)(i) = 
      ptr->_MeanComputationStrategy->ComputeMean(ptr->_NumberOfImages,
                                                 &voxelData[0]);    
  }
  return 0;
}

void
AtlasBuilder::
ThreadedUpdateAverageImage()
{
  if (this->_AverageImagePointer->getNumElements() < this->_NumberOfThreads)
  {
    // unlikely, but possible
    this->UpdateAverageImage();
    return;
  }

  //
  // start threads which will update the average
  std::vector<pthread_t> threads(this->_NumberOfThreads);

  // info passed to threads
  std::vector<ThreadInfo> threadInfo(this->_NumberOfThreads);
  
  // create threads that will update the average
  for (unsigned int threadIndex = 0; threadIndex < this->_NumberOfThreads; 
       ++threadIndex)
  {
    threadInfo[threadIndex].atlasBuilder = this;
    threadInfo[threadIndex].threadIndex  = threadIndex;
    
    int rv = pthread_create(&threads[threadIndex], NULL,
                            &AtlasBuilder::ThreadedUpdateAverageImage, 
                            &threadInfo[threadIndex]);
    if (rv != 0)
    {
      throw std::runtime_error("Error creating thread.");
    }
  }
  
  // join threads
  for (int threadIndex = 0; threadIndex < this->_NumberOfThreads; 
       ++threadIndex)
  {
    int rv = pthread_join(threads[threadIndex], NULL);
    if (rv != 0)
    {
      throw std::runtime_error("Error joining thread.");
    }
  }
}

void
AtlasBuilder::
UpdateError()
{
  unsigned int numElements = this->_AverageImagePointer->getNumElements();
  double MSE = 0;

  for (unsigned int i = 0; i < this->_NumberOfImages; ++i)
    {
      for (unsigned int j = 0; j < numElements; ++j)
      {
        double diff = (*this->_DeformedImagePointers[i])(j) - 
          (*this->_AverageImagePointer)(j);
        MSE += diff * diff;
      }
    }
  MSE /= (this->_NumberOfImages * numElements);
  this->_MeanSquaredError = MSE;
}

void
AtlasBuilder::
LogIterationData(unsigned int imageNumber, 
                 int threadID, 
                 double iterationTime)
{
  IterationData newData;
  newData.IterationNumber                = this->_Iteration;
  newData.ImageNumber                    = imageNumber;
  newData.IterationEllapsedTimeInSeconds = iterationTime;
  newData.TotalEllapsedTimeInSeconds     = this->_TotalTimer.getSeconds();
  newData.ProcessingThreadID             = threadID;
  newData.MeanSquaredError               = this->_MeanSquaredError;
  newData.RootMeanSquaredError           = sqrt(this->_MeanSquaredError);
  newData.Delta                          = this->_Delta[imageNumber];
  newData.MaxL2Displacement              = 
    *std::max_element(this->_MaxL2Displacements.begin(), 
                      this->_MaxL2Displacements.end());

  double MSEPercent  = 100.0;
  double RMSEPercent = 100.0;
  bool ValueDidImprove = true;
	static int continue_error_count = 0;
  if (newData.IterationNumber > 0)
  {
    double initialMSE = this->_IterationDataLog[0].MeanSquaredError;
    MSEPercent = 100.0 * newData.MeanSquaredError / initialMSE;

    double initialRMSE = this->_IterationDataLog[0].RootMeanSquaredError;
    RMSEPercent = 100.0 * newData.RootMeanSquaredError / initialRMSE;

    if (newData.MeanSquaredError > 
		    this->_IterationDataLog.back().MeanSquaredError)
    {
	    ValueDidImprove = false;
	    ++ continue_error_count;
    }else{
	    continue_error_count =0; // clean up 
    }
  }

  this->_IterationDataLog.push_back(newData);

  this->_IterationDataOutputStream 
	  << "["  << newData.IterationNumber
	  << "/"  << this->_FluidWarpParameters.numIterations
	  << " " << newData.ImageNumber
	  << "/" << this->_NumberOfImages << "] "
	  << "Time=" << newData.IterationEllapsedTimeInSeconds 
	  << "|" << newData.TotalEllapsedTimeInSeconds
	  << " " << "Thd=" << newData.ProcessingThreadID;

  if (this->_DeltaSelectionMethod == DELTA_USE_MEAN)
  {
	  this->_IterationDataOutputStream 
		  << " " 
		  << "D=M" << std::setprecision(4) << newData.Delta;
  }
  else
  {
	  this->_IterationDataOutputStream 
		  << " " << "D=IND";
  }
  this->_IterationDataOutputStream 
	  << " " << "MaxL2=" << newData.MaxL2Displacement
	  << " " << "RMSE=" << newData.RootMeanSquaredError 
	  << " (" << RMSEPercent << "%)";

  if (!ValueDidImprove)
  {

	  this->_IterationDataOutputStream 
		  << " <<--<<--<<-- MSE Increased! --<<--<<--<<";
	  if ( continue_error_count > 10 ){
		  //stopping condition  added by xx, 2010
		  this->_FluidWarpParameters.maxPerturbation /= 2.0; 
		  this->_IterationDataOutputStream 
			  << " Reduce maxPerturbation to  "<< this->_FluidWarpParameters.maxPerturbation;
		  if (  this->_FluidWarpParameters.maxPerturbation < 0.01){
			  this->_IterationDataOutputStream 
				  << " <<--<<--<<-- Stopped --<<--<<--<<";
			  this->_FluidWarpParameters.numIterations=0;
		  }
	  }
  }
  this->_IterationDataOutputStream << std::endl;
}


void
AtlasBuilder::
InverseOperatorMultiply(CoordinateType* complexPtr,
		float& L00,
		float& L10, float& L11,
		float& L20, float& L21, float& L22)
{
	float G00;
	float G10, G11;
	float G20, G21, G22;
	float y0, y1, y2;
	//
	// Given that A is pos-def symetric matrix, solve Ax=b by finding
	// cholesky decomposition GG'=A
	// and then performing 2 back-solves, Gy=b and then G'x=y to get x.
	// 

	// 1. find cholesky decomposition by finding G such that GG'=A.
	//    A must be positive definite symetric (we assume that here)
	//    G is then lower triangular, see algorithm 4.2.1 p142-3
	//    in Golub and VanLoan
	// Note: these are in matlab notation 1:3
	// [ G(1,1)   0      0    ]   [ G(1,1) G(2,1) G(3,1) ]   
	// [ G(2,1) G(2,2)   0    ] * [   0    G(2,2) G(3,2) ] = Amatrix
	// [ G(3,1) G(3,2) G(3,3) ]   [   0      0    G(3,3) ]

	float bRealX = complexPtr[0];
	float bRealY = complexPtr[2];
	float bRealZ = complexPtr[4];

	float bImagX = complexPtr[1];
	float bImagY = complexPtr[3];
	float bImagZ = complexPtr[5];

	float& vRealX = complexPtr[0];
	float& vRealY = complexPtr[2];
	float& vRealZ = complexPtr[4];

	float& vImagX = complexPtr[1];
	float& vImagY = complexPtr[3];
	float& vImagZ = complexPtr[5];

	G00 = sqrt(L00);
	G10 = L10 / G00;
	G20 = L20 / G00;

	G11 = L11 - G10 * G10;
	G21 = L21 - G20 * G10;
	G11 = sqrt(G11);
	G21 = G21 / G11;

	G22 = L22 - (G20*G20 + G21*G21);
	G22 = sqrt(G22);

	// back-solve Gy=b to get a temporary vector y
	// back-solve G'x=y to get answer in x
	//
	// Note: these are in matlab notation 1:3
	// [ G(1,1)   0      0    ]   [ y(1) ] = b(1)
	// [ G(2,1) G(2,2)   0    ] * [ y(2) ] = b(2)
	// [ G(3,1) G(3,2) G(3,3) ]   [ y(3) ] = b(3)
	//
	// [ G(1,1) G(2,1) G(3,1) ]   [ x(1) ] = y(1)
	// [   0    G(2,2) G(3,2) ] * [ x(2) ] = y(2)
	// [   0      0    G(3,3) ]   [ x(3) ] = y(3)
	y0 = bRealX / G00;
	y1 = (bRealY - G10*y0) / G11;
	y2 = (bRealZ - G20*y0 - G21*y1) / G22;

	vRealZ = y2 / G22;
	vRealY = (y1 - G21*vRealZ) / G11;
	vRealX = (y0 - G10*vRealY - G20*vRealZ) / G00;

	y0 = bImagX / G00;
	y1 = (bImagY - G10*y0) / G11;
	y2 = (bImagZ - G20*y0 - G21*y1) / G22;

	vImagZ = y2 / G22;
	vImagY = (y1 - G21*vImagZ) / G11;
	vImagX = (y0 - G10*vImagY - G20*vImagZ) / G00;
}
