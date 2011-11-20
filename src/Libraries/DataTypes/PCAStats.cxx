#include "PCAStats.h"
#include <iostream>
#include <fstream>

PCAStats:: PCAStats(){

     numPCs = 0;
	 dim = 0;
	 PCs = NULL;
     mean = NULL;
	 sigmas =NULL;  
}


PCAStats:: ~PCAStats(){

}

VectorType PCAStats:: GenObjVec(VectorType position){

   
	VectorType tVec(dim);
	tVec.fill(0);
	tVec = mean;

	// generate the locations using pcastats
	for (int i = 0; i < position.size(); i++){

		for (int j= 0; j<dim; j++){

			tVec[j] = tVec[j] + position[i] * sigmas[i] * PCs[j][i];
		}
	}
   return tVec;

}



bool  PCAStats:: readStats(char * fileName){

	FILE * pFile = fopen(fileName,"r");
	if (pFile == NULL) {std::cerr <<"File error"<< std::endl;exit(1);}


	if (! fscanf(pFile,"Dim = %d\n",&dim)) //Dim = 1000\n
	{ std::cerr<<"read error"<< std::endl;exit(1);} 


	//read mean
	char tmp[1];
	fscanf(pFile, "Mean = "); //Mean =\n


	mean = VectorType(dim);

	for(int i = 0 ; i < dim; i++) 
	{
		if(! fscanf(pFile, "%lf",& mean[i] ))
		{ 
			std::cerr<<"read error"<< std::endl;exit(1);

		} 
	}
	fscanf(pFile, "\n");


	//read PCs
	if (! fscanf( pFile,"numPCs = %d\n",&numPCs)) //numPCs = 3\n
	{
		std::cerr<<"read error"<< std::endl;exit(1);
	} 


	//read sigmas
	fscanf(pFile, "sigmas = "); 


	sigmas = VectorType(numPCs);

	for(int i = 0 ; i < numPCs; i++) 
	{
		if(! fscanf(pFile, "%lf",& sigmas[i] ))
		{ 
			std::cerr<<"read error"<< std::endl;exit(1);

		} 
	}
	fscanf(pFile, "\n");



	PCs = MatrixType(dim ,numPCs);//column matrix

	for (int j = 0; j <numPCs; j++)
	{
		fscanf(pFile,"PC[%c] = ", tmp); 

		for (int i =0; i < dim; i++)
		{
			fscanf(pFile, "%lf", &PCs(i, j) );
		}
		fscanf(pFile, "\n");
	}
return true;

}


bool PCAStats::writeStats(char * FileName){
//TODO

	return true;

}
