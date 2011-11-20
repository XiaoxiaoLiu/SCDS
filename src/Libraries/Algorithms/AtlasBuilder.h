#ifndef AtlasBuilder_h
#define AtlasBuilder_h

#include <vector>
#include <Array3D.h>
#include <FluidWarpParameters.h>
#include <fftw3.h>
#include <Timer.h>
#include <ostream>
#include <map>
#include <pthread.h>


template <class VoxelType>
class MeanComputationStrategyBase
{
 public:
  virtual VoxelType ComputeMean(unsigned int numValues, VoxelType* values)=0;
};

template <class VoxelType>
class ArithmeticMeanComputationStrategy:
public MeanComputationStrategyBase<VoxelType>
{
 public:
  void SetNumberOfElements(unsigned int n)
    {
      this->_weights.resize(n);
      this->SetWeightsToEqual();
    }

  unsigned int GetNumberOfElements() const
    {
      return this->_weights.size();
    }

  void SetNthWeight(unsigned int n, double weight)
    {
      this->_weights[n] = weight;
    }

  double GetNthWeight(unsigned int n) const
    {
      return this->_weights[n];
    }
  void SetWeightsToEqual()
    {
      std::fill(this->_weights.begin(), this->_weights.end(),
                1.0/this->_weights.size());
    }

  virtual VoxelType ComputeMean(unsigned int numValues, VoxelType* values) 
    {
      double mean = 0;
      for (unsigned int i = 0; i < numValues; ++i)
      {
        mean += values[i] * _weights[i];
      }
      return static_cast<VoxelType>(mean);
    } 

 private:
  std::vector<double> _weights;
};

class AtlasBuilder
{

public:
  typedef float                                  VoxelType;
  typedef float                                  CoordinateType;
  typedef Array3D<VoxelType>                     ImageType;
  typedef Array3D<Vector3D<CoordinateType> >     VectorFieldType;
  typedef MeanComputationStrategyBase<VoxelType> MeanComputationStrategyType;

  struct IterationData
  {
    unsigned int        IterationNumber;
    unsigned int        ImageNumber;
    double              IterationEllapsedTimeInSeconds;
    double              TotalEllapsedTimeInSeconds;
    unsigned int        ProcessingThreadID;
    double              MeanSquaredError;
    double              MaxL2Displacement;
    double              RootMeanSquaredError;
    double              Delta;
	double              DVFSquareError;
	double              RootMeanDVFSquareError;
	double              RootMeanCompDVFSquareError;
	double              DVFError;
	double              MaxDVFError;
	double              TotalPenalty;
  };

  //
  // constructors/destructors
  //
  AtlasBuilder();
  ~AtlasBuilder();

 
  //
  // fftw interface
  //
  void         SetFFTWNumberOfThreads(unsigned int numThreads);
  unsigned int GetFFTWNumberOfThreads() const;
  
  void         SetFFTWMeasureOn();
  void         SetFFTWMeasureOff();
  void         SetFFTWMeasure(bool b);
  bool         GetFFTWMeasure() const;

  //
  // ouput options
  //
  void         SetLogOutputStream(std::ostream& ostream);

  //
  // algorithm options
  //
  void         SetNumberOfThreads(unsigned int numThreads);
  unsigned int GetNumberOfThreads() const;

  void         SetUpdateAverageEverySubIterationOn();
  void         SetUpdateAverageEverySubIterationOff();
  void         SetUpdateAverageEverySubIteration(bool b);
  bool         GetUpdateAverageEverySubIteration() const;

  void         SetComputeInverseDeformationsOn();
  void         SetComputeInverseDeformationsOff();
  void         SetComputeInverseDeformations(bool b);
  bool         GetComputeInverseDeformations() const;

  void         SetFluidWarpParameters(const FluidWarpParameters& fluidParams);
  FluidWarpParameters& GetFluidWarpParameters();

  void         SetMeanComputationStrategy(MeanComputationStrategyType* s);
  MeanComputationStrategyType* GetMeanComputationStrategy() const;

  enum DeltaSelectionType { DELTA_USE_MEAN,
                            DELTA_USE_INDIVIDUAL
  };
  void               SetDeltaSelectionToIndividual()
  { this->_DeltaSelectionMethod = DELTA_USE_INDIVIDUAL; }
  void               SetDeltaSelectionToMean()
  { this->_DeltaSelectionMethod = DELTA_USE_MEAN; }
  DeltaSelectionType GetDeltaSelectionType() const
  { return this->_DeltaSelectionMethod; }

  //
  // set inputs
  //
  virtual void     SetNumberOfInputImages(unsigned int n);
  unsigned int     GetNumberOfInputImages() const;

  void             SetNthInputImage(unsigned int n, ImageType* imagePointer);
  ImageType*       GetNthInputImage(unsigned int n) const;

  ImageType*       GetNthDeformedImage(unsigned int n) const;

  void             SetNthDeformationField(unsigned int n, 
                                          VectorFieldType* fieldPointer);
  VectorFieldType* GetNthDeformationField(unsigned int n) const;

  void             SetNthDeformationFieldInverse(unsigned int n, 
                                                 VectorFieldType* 
                                                 fieldPointer);
  VectorFieldType* GetNthDeformationFieldInverse(unsigned int n) const;

  void             SetAverageImage(ImageType* imagePointer);
  ImageType*       GetAverageImage() const;

  //
  // run it
  //
  void             GenerateAverage();

  IterationData    GetIterationData(unsigned int iteration) const;

 protected:
  void             RunAlgorithm();
  static void*     ThreadedUpdateImages(void* arg);
  static void*     ThreadedUpdateAverageImage(void* arg);

  void             ProjectIncomp(CoordinateType* complexPtr, unsigned int x, unsigned int y, unsigned int z);

  void             InitializeScratchMemory();
  void             DeleteScratchMemory();
  void             CheckInputs();
  void             InitializeFFTWPlans();
  void             DeleteFFTWPlans();
  void             InitializeOperatorLookupTable();
  
  void             UpdateAverageImage();
  void             ThreadedUpdateAverageImage();
  virtual void     UpdateError();
  void             UpdateGradient         (unsigned int imageIndex, 
                                           unsigned int threadIndex);
  virtual void     UpdateBodyForce        (unsigned int imageIndex, 
                                           unsigned int threadIndex);
  void             UpdateVelocityField    (unsigned int imageIndex, 
                                           unsigned int threadIndex);
  virtual void     UpdateDeformationFields(unsigned int imageIndex, 
                                           unsigned int threadIndex);
  virtual  void    UpdateDeformedImage    (unsigned int imageIndex, 
                                           unsigned int threadIndex);
  virtual void     UpdateDelta();

  virtual void LogIterationData(unsigned int imageNumber,  int threadID, double iterationTime);

 //virtual  void             LogIterationData(unsigned int imageIndex,
   //                                 int threadID,
     //                               double iterationEllapsedSeconds);

  void             LockMeanImage();
  void             UnlockMeanImage();
  unsigned int     GetJobImageIndex();

  static void      InverseOperatorMultiply(CoordinateType* complexPtr,
                                           float& L00,
                                           float& L10, float& L11,
                                           float& L20, float& L21, float& L22);

  static pthread_t GetThreadID();
  long GetThreadIndex();
  
  //
  // look up table to for Linv computation
  //
  struct LUT
  {
    std::vector<CoordinateType> cosWX, cosWY, cosWZ;
    std::vector<CoordinateType> sinWX, sinWY, sinWZ;
    Array3D<float> nsq;

    LUT()
      : cosWX(0), cosWY(0), cosWZ(0),
	sinWX(0), sinWY(0), sinWZ(0)
    {}

    LUT(unsigned int xSize, 
	unsigned int ySize, 
	unsigned int zSize)
      : cosWX(xSize / 2 + 1), cosWY(ySize), cosWZ(zSize),
	sinWX(xSize / 2 + 1), sinWY(ySize), sinWZ(zSize), nsq(xSize,ySize,zSize)
    {}
  };

  AtlasBuilder(const AtlasBuilder& rhs);             // not implemented
  AtlasBuilder& operator =(const AtlasBuilder& rhs); // not implemented

  //
  // how many images will be processed
  unsigned int                     _NumberOfImages;

  //
  // input data managed by the user
  std::vector<ImageType*>          _ImagePointers;
  ImageType*                       _AverageImagePointer;

  std::vector<VectorFieldType*>    _DeformationFieldPointers;
  std::vector<VectorFieldType*>    _DeformationFieldInversePointers;

  // 
  // scratch data managed by this class
  std::vector<VectorFieldType*>    _ScratchVectorFieldPointers;
  std::vector<ImageType*>          _DeformedImagePointers;

  //
  // FFT parameters
  bool                             _FFTWMeasure;
  unsigned int                     _FFTWNumberOfThreads;
  std::vector<fftwf_plan>          _FFTWForwardPlans;
  std::vector<fftwf_plan>          _FFTWBackwardPlans;

  //
  // algorithm parameters
  bool                             _UpdateAfterEverySubIteration;
  bool                             _ComputeInverseDeformations;  
  MeanComputationStrategyType*     _MeanComputationStrategy;
  FluidWarpParameters              _FluidWarpParameters;
  LUT                              _OperatorLookupTable;
  
  //
  // processing current state
  std::vector<double>              _Delta;
  double                           _MeanSquaredError;
  unsigned int                     _Iteration;
  std::vector<IterationData>       _IterationDataLog;
  Timer                            _TotalTimer;
  std::vector<double>              _MaxL2Displacements;
  std::vector<double>              _MaxVelocityL2Displacements;

  //
  // processing options
  DeltaSelectionType               _DeltaSelectionMethod;

  //
  // used for multithreading
  unsigned int                     _NumberOfThreads;
  unsigned int                     _NextImageToProcess;
  pthread_mutex_t                  _NextImageToProcessMutex;
  pthread_mutex_t                  _AverageImageMutex;

  //
  // output settings
  std::ostream&                    _IterationDataOutputStream;
};

#endif
