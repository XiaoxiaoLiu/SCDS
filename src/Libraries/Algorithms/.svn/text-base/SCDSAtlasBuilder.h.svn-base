#ifndef SCDSAtlasBuilder_h
#define SCDSAtlasBuilder_h

#include <AtlasBuilder.h>


class SCDSAtlasBuilder : public AtlasBuilder {

  //
  // constructors/destructors
  //
public:
  SCDSAtlasBuilder();
  ~SCDSAtlasBuilder();
  
   void     SetSigmas( double * sigmas);
   double * GetSigmas(){return _sigmas;}
   void     SetShapeScores (double * shapeScores);
   double * GetShapeScores() {return _shapeScores;}
   void     SetNumberOfInputImages(unsigned int n);
   void     SetNthInputHField(unsigned int n,  VectorFieldType* HFieldPointer);
   void    SetNthInputCompHField(unsigned int n, VectorFieldType* CompHFieldPointer );
   void    SetCompareHFieldOn(){COMPARE_HFIELD_ON=true;}
   void     SetWriteXSlice(int sliceNum){_writeXSlice = sliceNum;}
   void     SetWriteYSlice(int sliceNum){_writeYSlice = sliceNum;}
   void     SetWriteZSlice(int sliceNum){_writeZSlice = sliceNum;}
   void     WriteVolumeOn(){_writeVolume = true ;}
   void     WriteVolumeOff(){_writeVolume = false ;}

   int     GetWriteXSlice(){return _writeXSlice;}
   int     GetWriteYSlice(){return _writeYSlice; }
   int     GetWriteZSlice(){return _writeZSlice; }
   bool    GetWriteVolumeOn(){return _writeVolume ;}
   void    SetOutputFileNamePrefix(std::string outputFileNamePrefix){
	   _outputFileNamePrefix = outputFileNamePrefix;}


protected:
	void     writeXSlice();
	void	 writeYSlice();
	void	 writeZSlice();
    void     writeVolume();
	void     LogIterationData(unsigned int imageNumber,  int threadID, double iterationTime);
    void     UpdateError();
	void     UpdateBodyForce(unsigned int imageIndex,  unsigned int threadIndex);
	void     UpdateDeformedImage(unsigned int imageIndex, unsigned int threadIndex);
	void     UpdateHFieldBodyForce (unsigned int imageIndex,unsigned int threadIndex);
	void     UpdateDeformationFields(unsigned int imageIndex, unsigned int threadIndex);

	void     UpdateDelta();
   
	double     _DVFSquareError;//DVF Sqaure Magnitude  Error  among  all voxels from all images
	double      _CompDVFSquareError;// For comparison purpose, compare to the grount truth diffeo
	double     _MaxDVFError; //MaxiMum DVF  Magnitude Error  among  all voxels from all images
	double     _DVFError;//DVF  Magnitude Error  among  all voxels from all images
    double *    _sigmas; // weights for each channel of the body force( intensity + displacement vector field)
    double *    _shapeScores;
    std::vector<VectorFieldType*>   _inputHFieldPointers;
	std::vector<VectorFieldType*>   _inputCompHFieldPointers;// For comparison purpose, compare to the grount truth diffeo

	bool _writeVolume ;
	int _writeXSlice;
	int _writeYSlice;
	int _writeZSlice;
	std::string _outputFileNamePrefix;
	bool COMPARE_HFIELD_ON;
};

#endif
