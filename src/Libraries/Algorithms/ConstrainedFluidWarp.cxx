#include "ConstrainedFluidWarp.h"
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <algorithm>
#include <sstream>

#include "Timer.h"
#include "Array3DUtils.h"
#include "HField3DUtils.h"
#include <time.h>

#include <cfloat> // old g++ compiler cant find limits
#include "Array3DIO.h"
#include "HField3DIO.h"
#include <sstream>
#define _USE_MATH_DEFINES
#include <math.h>


ConstrainedFluidWarp::ConstrainedFluidWarp()
{

  DVFSigma[0] = 1.0;
  DVFSigma[1] = 1.0;
  DVFSigma[2] = 1.0;

}


ConstrainedFluidWarp::~ConstrainedFluidWarp()
{


}


void ConstrainedFluidWarp::
computeHFieldAsymmetric(const ImageType& fixed,
			  const ImageType& moving,
      VectorFieldType* h_constrain,
			Parameters& parameters,
			  VectorFieldType*h)
{
  if (parameters.numIterations == 0)
    {
      return;
    }

  Timer totalTimer;
  totalTimer.start();

  //
  // make sure images are the same size, and that h and hinv are the
  // same size as the images
  //
  Vector3D<unsigned int> imageSize = fixed.getSize();
  if (imageSize != moving.getSize())
    {
      _report("Incompatible array sizes.\n", FW_OUTPUT_TYPE_ERROR);
      throw std::runtime_error("incompatible image array sizes");
    }
  if (imageSize != h->getSize())
    {
      _report("Incompatible h and image array sizes.\n", FW_OUTPUT_TYPE_ERROR);
      throw std::runtime_error("incompatible h and image array sizes");
    }


  //
  // allocate memory for velocity field and deformed image
  //
  _time("allocating memory...");

  // this multipurpose array is used to hold
  // gradient, body force, and velocity arrays (in that order)
  // this saves LOTS of memory
  //
  // VERY IMPORTANT: however, because it is used to hold the result of an fft
  // it needs to be padded just a little bit, see www.fftw.org
  // this wont affect its use as gradient and body force arrays
  // as long as logical size (image size) is used for this array
  // and access into this array is done **ONLY** via the (x, y, z) operator
  // and not via incremented pointers.  I REPEAT, don't access the vector
  // field via incremented pointers unless you *know* what you are doing.
  unsigned xSizeFFT = 2 * (imageSize.x / 2 + 1);
  VectorFieldType vf(xSizeFFT, imageSize.y, imageSize.z);

  // allocate deformed image (I(hk(x)))
  ImageType def(imageSize);
  //VectorFieldType hhinv(imageSize);
  _stopTime();

  //
  // initialize fftw plans
  //
  std::stringstream ss;
  ss << "initializing fftw plans (" 
     << this->_FFTWNumberOfThreads 
     << (this->_FFTWNumberOfThreads > 1 ? " threads)..." : " thread)...");
  _time(ss.str().c_str());
  fftwf_plan fftwForwardPlan, fftwBackwardPlan;
  _createFFTWPlans(vf, imageSize, fftwForwardPlan, fftwBackwardPlan);
  _stopTime();

  //
  // precompute info into luts
  //
  _time("precomputing luts...");
  LUT lut(imageSize.x, imageSize.y, imageSize.z);
  _initializeLUT(lut, imageSize);
  _stopTime();

  float delta = 0;
  bool haveDelta = false;
  double lastSquaredError = FLT_MAX;
  _openSpreadsheetOfstream();
  for (unsigned int iter = 0; iter < parameters.numIterations; ++iter)
    {
      //
      // deform the moving image
      //
      _time("deforming image...");
      HField3DUtils::apply(moving, *h, def, (float)0.0, false);
      _stopTime();

      //
      // debug write images
      //
//       std::ostringstream oss;
//       oss << "Def" << iter;
//       Array3DIO::writeMETAVolume(def, oss.str().c_str());
      //

      if( parameters.jacobianScale ) 
        {
          std::cerr << "Scaling by Jac." << std::endl;
          ImageType jacobian;
          HField3DUtils::jacobian(*h, jacobian);
          float min, max;
          Array3DUtils::getMinMax( jacobian, min, max );
          std::cout << "min, max = " << min << ", " << max << std::endl;

          ImageType::SizeType size = def.getSize();
          for(unsigned int k = 0; k < size[2]; ++k) {
            for(unsigned int j = 0; j < size[1]; ++j) {
              for(unsigned int i = 0; i < size[0]; ++i) {
                def(i, j, k) *= jacobian(i, j, k);
              }
            }
          }
        }

      double intensitySquaredError;
      double dvfSquaredError;
      if( parameters.jacobianScale ) 
        {
          _time("generating body force...");
          // squared error gets set here
          _generateBodyForceJacobianScale(fixed, def, vf, intensitySquaredError);
          _stopTime();
	} 
      else 
        {
          //
          // compute the gradient of the deformed image
          //
          _time("computing gradient...");
          Array3DUtils::computeGradient(def, vf);
          _stopTime();

          //
          // debug: write gradient for display
          //
//           Array3D<CoordinateType> gy(imageSize);
//           for (unsigned int z = 0; z < imageSize.z; ++z) {
//             for (unsigned int y = 0; y < imageSize.y; ++y) {
//               for (unsigned int x = 0; x < imageSize.x; ++x) {
//                 gy(x,y,z) = vf(x,y,z).normL2();
//               }
//             }
//           }
//           //Array3DUtils::rescaleElements(gy, 0.0F, 1.0F);
//           std::ostringstream tmpnamegy;
//           tmpnamegy << "GradY" << iter;
//           Array3DIO::writeMETAVolume(gy, tmpnamegy.str().c_str());      
          //

          //
          // generate the body force
          //
          _time("generating body force...");
          // squared error gets set here
          _generateBodyForce(fixed, def,h_constrain, h,vf, intensitySquaredError,dvfSquaredError);
          _stopTime(); 

          //
          // debug: write body force for display
          //
//            Array3D<CoordinateType> by(imageSize);
//            for (unsigned int z = 0; z < imageSize.z; ++z) {
//              for (unsigned int y = 0; y < imageSize.y; ++y) {
//                for (unsigned int x = 0; x < imageSize.x; ++x) {
//                  by(x,y,z) = vf(x,y,z).normL2();
//                }
//              }
//            }
//            //Array3DUtils::rescaleElements(by, 0.0F, 1.0F);
//            std::ostringstream tmpnameby;
//            tmpnameby << "BodyForceY" << iter;
//            Array3DIO::writeMETAVolume(by, tmpnameby.str().c_str());      
          //
	}

      //
      // compute the velocity field
      //
      _time("computing velocity field...");
      _computeVelocityField(vf, imageSize, parameters, lut,
			    fftwForwardPlan, fftwBackwardPlan);
      _stopTime();

      //
      // debug: write velocity for display
      //
//       Array3D<CoordinateType> vy(imageSize);
//       for (unsigned int z = 0; z < imageSize.z; ++z) {
//         for (unsigned int y = 0; y < imageSize.y; ++y) {
//           for (unsigned int x = 0; x < imageSize.x; ++x) {
//             vy(x,y,z) = vf(x,y,z).normL2();
//           }
//         }
//       }
//       //Array3DUtils::rescaleElements(vy, 0.0F, 1.0F);
//       std::ostringstream tmpnamevy;
//       tmpnamevy << "VelocityY" << iter;
//       Array3DIO::writeMETAVolume(vy, tmpnamevy.str().c_str());      
      //

      //
      // compute delta (currently only once at the beginning)
      //
      if (!haveDelta)
	{
	  _time("finding delta...");
	  delta = _computeDelta(vf, imageSize, parameters.maxPerturbation);
	  _stopTime();

	  std::ostringstream deltaString;
	  deltaString << "Computed delta: " << delta << std::endl;
	  _report(deltaString, FW_OUTPUT_TYPE_VERBOSE);
	  haveDelta = true;

          if (delta > 100)
            {
              std::cerr << "Images are too similar to process." << std::endl;
              break;
            }
	}

      //
      // update current h and hinv fields
      //
      _time("updating h fields...");

	  _updateHField(*h, vf, delta);

      _stopTime();

      //
      // output iter results
      //
      float rmsIntensityError = sqrt(intensitySquaredError / 
			    (imageSize.x * imageSize.y * imageSize.z));

      float rmsDVFError = sqrt(dvfSquaredError/ 
			    (imageSize.x * imageSize.y * imageSize.z));
    std::ostringstream oss;
    oss << "[" << iter <<  "] "
      << totalTimer.getSeconds()<< " "
      << std::setprecision(5)
      << delta << " "
      << intensitySquaredError<< " "
      << dvfSquaredError<< " "
      << rmsIntensityError<<" "
      << rmsDVFError<<std::endl;
  _report(oss, FW_OUTPUT_TYPE_STANDARD);

  //
  // flag increases in squaredError
  //
  if (iter != 0 && (intensitySquaredError+dvfSquaredError)> lastSquaredError)
    {
      _report("  <<---<<---<<---<<---<<", FW_OUTPUT_TYPE_STANDARD);
    }
     
      //adaptively change the stepsize
	  if ( iter != 0 && intensitySquaredError >  lastSquaredError){
		  parameters.maxPerturbation =parameters.maxPerturbation/2;
		  std::cerr <<"---->>----reduce maxPerturbation to half:" <<parameters.maxPerturbation<< std::endl;
		  if (parameters.maxPerturbation< 0.01)
		  {
			  parameters.numIterations=0;
			  std::cerr << "---->>----stopped---->>----" << std::endl;
		  }
	  }
      lastSquaredError = dvfSquaredError+intensitySquaredError;


      //
      // optionally write out volumes and/or slices for movies, etc.
      //
      _time("Writing per iteration data (if requested)...");
      //_writePerIterationData(iter, 0, 1, rmsError, def, fixed, *h, NULL);
      _stopTime();
    }

  _closeSpreadsheetOfstream();

  //
  // clean up memory
  //
  fftwf_destroy_plan(fftwForwardPlan);
  fftwf_destroy_plan(fftwBackwardPlan);

  //
  // report the total time
  //
  totalTimer.stop();
  std::ostringstream totalTimeString;
  totalTimeString << "Total Time: " << totalTimer.getSeconds() 
		  << " (sec)" << std::endl;
  _report(totalTimeString, FW_OUTPUT_TYPE_STANDARD);
}



void 
ConstrainedFluidWarp
::_generateBodyForce(const ImageType& fixed,
                     const ImageType& def,
                     VectorFieldType  * h_constrain,
                     VectorFieldType  * h,
                     VectorFieldType& gradToBodyForce,
                     double& intensitySquaredError,
                     double& dvfSquaredError)
{
  Vector3D<unsigned int> size = fixed.getSize();
  double di;
  intensitySquaredError = 0;
  dvfSquaredError= 0;

	// put xyz channels of the hfield into three images
	ImageType ** hFieldImages = new ImageType*[3];
	hFieldImages[0] = new ImageType(size.x, size.y, size.z);
	hFieldImages[1] = new ImageType(size.x, size.y, size.z);
	hFieldImages[2] = new ImageType(size.x, size.y, size.z);
	for (unsigned int z = 0; z < size.z; ++z) {
		for (unsigned int y = 0; y < size.y; ++y) {
			for (unsigned int x = 0; x < size.x; ++x) {
				hFieldImages[0]->set(x,y,z,h->get(x, y, z).x);
				hFieldImages[1]->set(x,y,z,h->get(x, y, z).y);
				hFieldImages[2]->set(x,y,z,h->get(x, y, z).z);
			}
		}
	}

  Vector3D<float> zerov;
  zerov[0]=0.0;
  zerov[1]=0.0;
  zerov[2]=0.0;
	VectorFieldType * *  tmpScratchVectorFieldPointer = new VectorFieldType*[3];
  for (int i =0; i<3; i++){
		tmpScratchVectorFieldPointer[i] = new VectorFieldType(size.x,size.y, size.z);
		tmpScratchVectorFieldPointer[i]->fill(zerov);  //store  the gradient
		Array3DUtils::computeGradient(*hFieldImages[i],*tmpScratchVectorFieldPointer[i]);
	}

//combine intensity and hfield forces together
  for (unsigned int z = 0; z < size.z; ++z) {
    for (unsigned int y = 0; y < size.y; ++y) {
      for (unsigned int x = 0; x < size.x; ++x) {
        di = 0;
        double dc[3];
        dc[0] = 0.0;
        dc[1] = 0.0;
        dc[2] = 0.0;
        di = fixed(x, y, z) - def(x, y, z);
        intensitySquaredError += di * di;
        gradToBodyForce(x, y, z) *= di;
        for (int i=0;i<3;i++)
           {
            dc[i] = h_constrain->get(x,y,z)[i] - h->get(x,y,z)[i];
            dvfSquaredError += dc[i]*dc[i];
            gradToBodyForce(x, y, z) += DVFSigma[i] * (*tmpScratchVectorFieldPointer[i])(x,y,z) * dc[i];
           }
      }
    }
  }
}



