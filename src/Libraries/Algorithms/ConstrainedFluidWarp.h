#ifndef ConstrainedFluidWarp_h
#define ConstrainedFluidWarp_h


#include "FluidWarp.h"

class ConstrainedFluidWarp : public FluidWarp {

public:
  ConstrainedFluidWarp();
  ~ConstrainedFluidWarp();
  void SetDVFSigma(double* sigma){
      DVFSigma[0]=sigma[0];
      DVFSigma[1]=sigma[1];
      DVFSigma[2]=sigma[2];
   }
  double * GetDVFSigma(){return DVFSigma;}


void computeHFieldAsymmetric(const ImageType& fixed,
			  const ImageType& moving,
      VectorFieldType * h_constrain,
			Parameters& parameters,
			  VectorFieldType * h);

protected:
void _generateBodyForce(const ImageType& fixed,
                     const ImageType& def,
                     VectorFieldType * h_constrain,
                     VectorFieldType * h,
                     VectorFieldType& gradToBodyForce,
                     double& intensitySquaredError,double& dvfSquaredError);
double DVFSigma[3];
private:

};

#endif
