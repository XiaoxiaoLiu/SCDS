
#ifndef PCA_STATS_H
#define PCA_STATS_H
#include "vnl/vnl_matrix_fixed.h"
#include "vnl/vnl_math.h"


typedef vnl_vector< double>				 VectorType;
typedef vnl_matrix< double>				 MatrixType;


//PCA stastitsics
class PCAStats{

public:
	PCAStats();
   ~PCAStats();
	
    bool readStats(char * FileName);
	bool writeStats(char * FileName);
	VectorType  GetMean(){return mean;}
	int GetNumPCs(){return numPCs;}
	int GetDim(){return dim;}
	MatrixType  GetPCs(){return PCs;}
	VectorType  GetSigmas(){return sigmas;}

	VectorType  GenObjVec(VectorType position);

protected:
	int numPCs;
	int dim;
	MatrixType  PCs;
	VectorType  mean;
	VectorType  sigmas; // standard deviation for each PC
	

};
#endif