#include <SCDSAtlasBuilder.h>
#include <strstream>
#include <HField3DUtils.h>
#include <assert.h>
#include <Timer.h>
#include <iomanip>
#define _USE_MATH_DEFINES
#include <math.h>
#include <HField3DIO.h>
#include <ImageUtils.h>




SCDSAtlasBuilder::
SCDSAtlasBuilder(){
	//set up the equall weights
	_sigmas = new double[4];
	_sigmas[0]=0.25;
	_sigmas[1]=0.25;
	_sigmas[2]=0.25;
	_sigmas[3]=0.25;
	_shapeScores = NULL;
	_writeVolume = false;
	_writeXSlice = -1;
	_writeYSlice = -1;
	_writeZSlice = -1;
	COMPARE_HFIELD_ON=false;
	_outputFileNamePrefix =std::string("");
	_MaxDVFError = 0;
	_DVFSquareError=0;
	_DVFError= 0;
}


SCDSAtlasBuilder::
~SCDSAtlasBuilder(){
}



void SCDSAtlasBuilder::
SetSigmas(double * sigmas){
	_sigmas[0] = sigmas[0];
	_sigmas[1] = sigmas[1];
	_sigmas[2] = sigmas[2];
	_sigmas[3] = sigmas[3];
}

void SCDSAtlasBuilder::
SetShapeScores (double * shapeScores){
	_shapeScores = shapeScores;

}
void 
SCDSAtlasBuilder::
SetNumberOfInputImages(unsigned int n)
{
	this->_NumberOfImages = n;

	this->_ImagePointers.resize(n);
	this->_DeformationFieldPointers.resize(n);
	this->_DeformationFieldInversePointers.resize(n);
	this->_Delta.assign(n, 0.0);

	this->_inputHFieldPointers.resize(n);
	this->_inputCompHFieldPointers.resize(n);
	if(_shapeScores==NULL){
	_shapeScores = new double[n];
			for (unsigned int i = 0; i < n; ++i)
				_shapeScores[i]=1.0;
	}
	if (_shapeScores == NULL){
		_shapeScores = new double[_NumberOfImages];
		for (unsigned int i = 0; i < this->_NumberOfImages; ++i)
			_shapeScores[i]=1.0;
		}

}


void
SCDSAtlasBuilder::
SetNthInputHField(unsigned int n, VectorFieldType* HFieldPointer)
{
	this->_inputHFieldPointers[n] = HFieldPointer;
	
}
void
SCDSAtlasBuilder::
SetNthInputCompHField(unsigned int n, VectorFieldType* CompHFieldPointer )
{
  this->_inputCompHFieldPointers[n] = CompHFieldPointer;
}

void
SCDSAtlasBuilder:: UpdateError() {
	unsigned int numElements = this->_AverageImagePointer->getNumElements();
	double MSE = 0;
	double DVFSquareError=0;
	double CompDVFSquareError=0;
	double DVFError=0;
	double MaxDVFError=0;
	for (unsigned int i = 0; i < this->_NumberOfImages; ++i)
	{
		for (unsigned int j = 0; j < numElements; ++j)
		{
			double diff = (*this->_DeformedImagePointers[i])(j) - 
				(*this->_AverageImagePointer)(j);
			MSE += diff * diff;

			Vector3D<double> diff1 = (*this->_DeformationFieldPointers[i])(j) - 
				(*this->_inputHFieldPointers[i])(j);
			DVFError += diff1.length();

			if ( MaxDVFError < diff1.length())  
			{
				MaxDVFError = diff1.length();
			}

			DVFSquareError += diff1.length()* diff1.length();
          

			if (COMPARE_HFIELD_ON)
		   { Vector3D<double> diff2 = (*this->_DeformationFieldPointers[i])(j) - 
				(*this->_inputCompHFieldPointers[i])(j);
        	CompDVFSquareError+= diff2.length()* diff2.length();
		   }
		}
	}
	MSE /= (this->_NumberOfImages * numElements);
	DVFSquareError /= (this->_NumberOfImages * numElements);
	CompDVFSquareError /=(this->_NumberOfImages * numElements);
	DVFError /= (this->_NumberOfImages * numElements);

	this->_MeanSquaredError = MSE;
	this->_MaxDVFError = MaxDVFError;
	this->_DVFSquareError = DVFSquareError;
	this->_CompDVFSquareError = CompDVFSquareError;
	this->_DVFError = DVFError;
}


void   
SCDSAtlasBuilder::
LogIterationData(unsigned int imageNumber,  int threadID, double iterationTime)
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
	newData.MaxL2Displacement              = *std::max_element(this->_MaxL2Displacements.begin(), 
			this->_MaxL2Displacements.end());
	newData.DVFSquareError = this->_DVFSquareError;
	newData.RootMeanDVFSquareError =sqrt(this->_DVFSquareError);
    newData.RootMeanCompDVFSquareError =sqrt(this->_CompDVFSquareError); newData.DVFError = this->_DVFError; newData.MaxDVFError = this->_MaxDVFError;

	double IntensityErrorPercent  = 100.0;
	double DiffeoDisPercent = 100.0;
	double CompDiffeoDisPercent = 100.0;
	newData.TotalPenalty = 200.0;  //initial total percentage
	bool ValueDidImprove = true;
	static int continue_error_count = 0;


	if (newData.IterationNumber == 0){
	//print out notes
		std::cout<<"image: RMSE of Intensity;  diffeo: RMSE of DVF magnitude"<<std::endl;
	}

	if (newData.IterationNumber > 0)
	{

		double initialTotalPenalty = this->_IterationDataLog[0].TotalPenalty;
		IntensityErrorPercent =100.0 * newData.RootMeanSquaredError/ this->_IterationDataLog[0].RootMeanSquaredError;
		DiffeoDisPercent = 100.0 * newData.RootMeanDVFSquareError/ this->_IterationDataLog[0].RootMeanDVFSquareError;
		if (COMPARE_HFIELD_ON){ CompDiffeoDisPercent= 100.0 * newData.RootMeanCompDVFSquareError/ this->_IterationDataLog[0].RootMeanCompDVFSquareError;}
		//this condition should be removed later
		if (_sigmas[1] == 0)
		{//no deformation constrains are actually used, the inputHfields just used as evaluation purpose
			newData.TotalPenalty = IntensityErrorPercent;}
		else{
			newData.TotalPenalty = IntensityErrorPercent +  DiffeoDisPercent;
		}


		if (newData.TotalPenalty > this->_IterationDataLog.back().TotalPenalty)
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
		<< "/"  << this->_FluidWarpParameters.numIterations;
	// << " " << newData.ImageNumber
	//  << "/" << this->_NumberOfImages << "] "
	//  << "Time=" << newData.IterationEllapsedTimeInSeconds 
	//  << "|" << newData.TotalEllapsedTimeInSeconds
	// << " " << "Thd=" << newData.ProcessingThreadID;

	if (this->_DeltaSelectionMethod == DELTA_USE_MEAN)
	{
		this->_IterationDataOutputStream 
			<< " " 
			<< "D=M" << std::setprecision(3) << newData.Delta;
	}
	else
	{
		this->_IterationDataOutputStream 
			<< " " << "D=IND";
	}
	this->_IterationDataOutputStream 
		<< " " << "MaxL2=" <<  std::setprecision(3) << newData.MaxL2Displacement
		<<" " << "MaxDVFError=" <<  std::setprecision(3) << newData.MaxDVFError
		<< " " << "Image= " << std::setprecision(3)<< newData.RootMeanSquaredError 
		<< " (" << std::setprecision(3)<<IntensityErrorPercent << "%)"
		<< " " << "Diffeo= "<< std::setprecision(3)<<newData.RootMeanDVFSquareError
		<< " (" <<std::setprecision(3)<< DiffeoDisPercent << "%)";
	if (COMPARE_HFIELD_ON){
		this->_IterationDataOutputStream 	   
			<<" " << "CompareDiffeo= "<< std::setprecision(3)<<newData.RootMeanCompDVFSquareError;
	}
	this->_IterationDataOutputStream 
		<< " (" <<std::setprecision(3)<< CompDiffeoDisPercent << "%)"
		<< " " << "Total= "<< std::setprecision(3)<<newData.TotalPenalty <<std::endl;

	//stopping condition
	if (!ValueDidImprove)
	{

		this->_FluidWarpParameters.maxPerturbation /= 2.0; 
		this->_IterationDataOutputStream 
			<< " Reduce maxPerturbation to  "<< this->_FluidWarpParameters.maxPerturbation<<" --<<--<<--<<";

		
		//debug, save the deformation field
		//	 Vector3D<float> origin, spacing; 
		//	origin = Vector3D<float>(0,0,0);
		//	spacing = Vector3D<float>(1,1,1);

		//	HField3DIO::writeMETA(*this->_DeformationFieldPointers[0], origin, spacing, "test_hfield");
		//	std::cout << "DONE" << std::endl; 


		if ( continue_error_count > 10 ){
			this->_IterationDataOutputStream 
				<< " <<--<<--<<-- Total Penalty Increased!";

			this->_FluidWarpParameters.maxPerturbation /= 2.0; 
			this->_IterationDataOutputStream 
				<< " Reduce maxPerturbation to  "<< this->_FluidWarpParameters.maxPerturbation<<" --<<--<<--<<";

			if( this->_FluidWarpParameters.maxPerturbation < 0.001){
		//stopping the iterations, exit
				
				this->_IterationDataOutputStream 
					<< " <<--<<--<<-- Stopped --<<--<<--<<";
				this->_FluidWarpParameters.numIterations=0;
			}
		}

	}
	this->_IterationDataOutputStream << std::endl;				

	//TO to removed: save the middle results  for debug
	if (_writeXSlice >=0){writeXSlice();}

	if (_writeYSlice >=0){writeYSlice();}
	if (_writeZSlice >=0){writeZSlice();}

	if (_writeVolume){writeVolume();}

}


void SCDSAtlasBuilder::
writeXSlice()
{
	int  _xSliceNumber =_writeXSlice;
	std::ostringstream oss;
	oss << _outputFileNamePrefix <<"x_"<< this->_Iteration;  

	std::cout << "writing x slice: " << oss.str().c_str() << "...";
	Array3DIO::writeMETASliceX<VoxelType>(*this->_AverageImagePointer, _xSliceNumber, oss.str().c_str());
	std::cout << "DONE" << std::endl;
}


void SCDSAtlasBuilder::
writeYSlice()
{
	int  _ySliceNumber =_writeYSlice;
	std::ostringstream oss;
	oss << _outputFileNamePrefix <<"y_"<< this->_Iteration;  

	std::cout << "writing y slice: " << oss.str().c_str() << "...";
	Array3DIO::writeMETASliceY<VoxelType>(*this->_AverageImagePointer, _ySliceNumber, oss.str().c_str());
	std::cout << "DONE" << std::endl;
}


void SCDSAtlasBuilder::
writeZSlice()
{
	int  _zSliceNumber =_writeZSlice;
	std::ostringstream oss;
	oss << _outputFileNamePrefix <<"z_"<< this->_Iteration;  

	std::cout << "writing z slice: " << oss.str().c_str() << "...";
	Array3DIO::writeMETASliceZ<VoxelType>(*this->_AverageImagePointer, _zSliceNumber, oss.str().c_str());
	std::cout << "DONE" << std::endl;
}


void SCDSAtlasBuilder::
writeVolume()
{
	std::ostringstream oss;
	oss << _outputFileNamePrefix<< this->_Iteration;  

	std::cout << "writing volume: "<< oss.str().c_str()<<"...";
	ImageUtils::writeMETA<VoxelType>(*(this->_AverageImagePointer), oss.str().c_str());
	std::cout << "DONE" << std::endl;
}



void 
SCDSAtlasBuilder::
UpdateBodyForce(unsigned int imageIndex, unsigned int threadIndex)
{
	Vector3D<unsigned int> size = this->_AverageImagePointer->getSize();
	double di;
	for (unsigned int z = 0; z < size.z; ++z) {
		for (unsigned int y = 0; y < size.y; ++y) {
			for (unsigned int x = 0; x < size.x; ++x) {

				di =  (*(this->_AverageImagePointer))(x, y, z)-
					(*this->_DeformedImagePointers[imageIndex])(x, y, z);

				(*this->_ScratchVectorFieldPointers[threadIndex])(x,y,z) *=  _sigmas[0]*di ;

			}
		}
	}
	//combine the forces together
	UpdateHFieldBodyForce(imageIndex,threadIndex);


} 


void 
SCDSAtlasBuilder::
UpdateDeformedImage(unsigned int imageIndex, unsigned int threadIndex)
{
	HField3DUtils::apply(*this->_ImagePointers[imageIndex],
			*this->_DeformationFieldPointers[imageIndex],
			*this->_DeformedImagePointers[imageIndex]);
}



void 
SCDSAtlasBuilder::
UpdateHFieldBodyForce(unsigned int imageIndex, unsigned int threadIndex)
{
	// allocate scratch vector field
	Vector3D<unsigned int> imageSize = this->_AverageImagePointer->getSize();

	unsigned xSizeFFT =  2 * (imageSize.x / 2 + 1);
	VectorFieldType * *  tmpScratchVectorFieldPointer = new VectorFieldType*[3];


	// Zero out the vectors in scratch memory before accumulating forces
	Vector3D<VoxelType> zerov;
	zerov[0] = 0;
	zerov[1] = 0;
	zerov[2] = 0;

	unsigned int numChannels = 3;

	Vector3D<unsigned int> size = this->_AverageImagePointer->getSize();
	// put xyz channels of the hfield into three images
	ImageType ** hFieldImages = new ImageType*[3];
	hFieldImages[0] = new ImageType(size.x, size.y, size.z);
	hFieldImages[1] = new ImageType(size.x, size.y, size.z);
	hFieldImages[2] = new ImageType(size.x, size.y, size.z);
	for (unsigned int z = 0; z < size.z; ++z) {
		for (unsigned int y = 0; y < size.y; ++y) {
			for (unsigned int x = 0; x < size.x; ++x) {
				hFieldImages[0]->set(x,y,z,  this->_DeformationFieldPointers[imageIndex]->get(x, y, z).x);
				hFieldImages[1]->set(x,y,z,  this->_DeformationFieldPointers[imageIndex]->get(x, y, z).y);
				hFieldImages[2]->set(x,y,z,  this->_DeformationFieldPointers[imageIndex]->get(x, y, z).z);
			}
		}
	}

	// no smoothing, original method
	VectorFieldType * inputH =_inputHFieldPointers[imageIndex];

	for (int i =0; i<3; i++){

		tmpScratchVectorFieldPointer[i] = new VectorFieldType(xSizeFFT, imageSize.y, imageSize.z);
		tmpScratchVectorFieldPointer[i]->fill(zerov);  //store  the gradient



		Array3DUtils::computeGradient(*hFieldImages[i],*tmpScratchVectorFieldPointer[i]);


	}


	/* Turns out to be too much, over kill the problem, converge too earlier!!!
	//smooth version

	VectorFieldType * inputH =_inputHFieldPointers[imageIndex];

	//compute the gradient of hfield images
	Image<VoxelType> *	tmpHFieldChanhnelImage = new Image<VoxelType> (*hFieldImages[0]);
	for (int i =0; i<3; i++){

	tmpScratchVectorFieldPointer[i] = new VectorFieldType(xSizeFFT, imageSize.y, imageSize.z);
	tmpScratchVectorFieldPointer[i]->fill(zerov);  //store  the gradient

	//smooth the hField first, not downsizing	
	ImageUtils::gaussianDownsample<VoxelType>(*hFieldImages[i],*tmpHFieldChanhnelImage,
	Vector3D<int>(1,1,1),//factor,original size 
	Vector3D<double>(2, 2, 2),//sigmas
	Vector3D<int>(4, 4, 4));//kernal size

	Array3DUtils::computeGradient(*tmpHFieldChanhnelImage,*tmpScratchVectorFieldPointer[i]);

	 *hFieldImages[i]=*tmpHFieldChanhnelImage;

	 }
	 delete tmpHFieldChanhnelImage;

*/


	for (int i =0; i<3; i++){
		for (unsigned int z = 0; z < size.z; ++z) {
			for (unsigned int y = 0; y < size.y; ++y) {
				for (unsigned int x = 0; x < size.x; ++x) {
					(*tmpScratchVectorFieldPointer[i])(x,y,z)  *=  inputH->get(x,y,z)[i] - hFieldImages[i]->get(x,y,z)  ;

					(*this->_ScratchVectorFieldPointers[threadIndex])(x,y,z) += _sigmas[i+1]* (*tmpScratchVectorFieldPointer[i])(x,y,z);
				}
			}
		}

	}


	//clean up
	for (int i = 0; i < 3;++i)
	{
		delete tmpScratchVectorFieldPointer[i];
		delete hFieldImages[i];
	}

	delete  []  tmpScratchVectorFieldPointer;
	delete  [] hFieldImages;
}




void
SCDSAtlasBuilder::
UpdateDelta()
{
	if (this->_DeltaSelectionMethod == DELTA_USE_INDIVIDUAL)
	{
		this->_Delta.resize(this->_NumberOfImages);

		for (unsigned int i = 0; i < this->_NumberOfImages; ++i)
		{
			if(this->_MaxVelocityL2Displacements[i]>0)
			{this->_Delta[i] = _shapeScores[i]*  this->_FluidWarpParameters.maxPerturbation / this->_MaxVelocityL2Displacements[i];}
			else
			{this->_Delta[i]=0;}
		}
	}
	else
	{
		std::runtime_error("Unknown delta update method.");    
	}
}


void
SCDSAtlasBuilder::
UpdateDeformationFields(unsigned int imageIndex, unsigned int threadIndex)
{
	Vector3D<unsigned int> size = this->_AverageImagePointer->getSize();

	//
	// compute hIncremental(x) = x + velocity(x) * delta
	//


    //adaptively change the delta, to make sure each step is small enough, prevent divergence near the optimum
	if(this->_MaxVelocityL2Displacements[imageIndex]>0){
		  if (this->_MaxVelocityL2Displacements[imageIndex] 
		                      <(_shapeScores[imageIndex]*  this->_FluidWarpParameters.maxPerturbation))
		  {
			  //converged, maximum delta=1.0
			  this->_Delta[imageIndex]= 1.0;
		  }else{
		      this->_Delta[imageIndex]=_shapeScores[imageIndex]*  this->_FluidWarpParameters.maxPerturbation
		                                        / this->_MaxVelocityL2Displacements[imageIndex];
		  }
	}else{
		this->_Delta[imageIndex]=0.0;
	}

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


