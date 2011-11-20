
#include "DeformWithPrior.h"
#include <iostream>
#include <fstream>
#include <exception>

#include "argio.h"
#include "itkImageFileReader.h"
#include "itkMesh.h"





#ifdef _MSC_VER
#pragma warning ( disable : 4786 )
#pragma warning ( disable : 4996 )
#endif



#define appout std::cerr
const std::string PROGRAM_NAME = "pDeform";


// reads in the VTK mesh file 
MeshType::Pointer readMesh(const char * fileName){

	MeshType::Pointer m_mesh = MeshType::New();

	std::ifstream infile ;
	std::string vtkFileName (fileName);
	char line[200] ;

	int index, v1, v2, v3 ;
	int nPts, nTris ;
	int i ;


	infile.open ( vtkFileName.c_str() ) ;


	// skip header - all we need here is the number of vertices
	for (int i = 0 ; i < 5 ; i++)
		infile.getline (line, 200 ) ;
	nPts = atoi(strchr(line,' ')) ;

	MeshType::PointType p;
	for ( i = 0 ; i < nPts ; i++ )
	{
		infile >> 	p[0] >> 	p[1] >> 	p[2] ;
		m_mesh->SetPoint(i, p);

	}

	// skip the next couple of lines - all we need is the number of triangles
	infile.getline ( line, 200 ) ;
	while ( strncmp (line, "POLYGONS", 8) )
		infile.getline ( line, 200 ) ;
	nTris = atoi(strchr(line,' ')) ;


	CellType ::CellAutoPointer tri;

	for ( i = 0 ; i < nTris ; i++)
	{
		infile >> index >> v1 >> v2 >> v3 ;
		assert ( index == 3 ) ;

		tri.TakeOwnership(new TriangleType );

		tri->SetPointId(0,v1);
		tri->SetPointId(1,v2);
		tri->SetPointId(2,v3);
		m_mesh->SetCell(i,tri);
	}

	//read normals?



	infile.close () ;


	return m_mesh;


}


bool writeMesh(MeshType::Pointer  m_mesh,MeshType ::Pointer normals, char * fileName){
	//write the resulting object into a mesh file

	std::ofstream outfile ;
	int nPts, nTris ;

	std::string outfileName(fileName) ;

	outfile.open ( outfileName.c_str() ) ;

	nPts = m_mesh->GetNumberOfPoints();
	nTris = m_mesh->GetNumberOfCells();


	// output header
	outfile << "# vtk DataFile Version 1.0" << std::endl << "vtk output" << std::endl << "ASCII" <<std:: endl << "DATASET POLYDATA" ;
	outfile << std::endl << "POINTS " << nPts << " double" << std::endl ;



	PointsContainerIterator pointIterator = m_mesh->GetPoints()->Begin();
	while (pointIterator!=m_mesh->GetPoints()->End())
	{
		PointType p = pointIterator.Value();
		outfile << p[0]<< " " <<p[1]<< " " << p[2] << std::endl;
		++pointIterator;
	}

	outfile<<std::endl << "POLYGONS " << nTris << " " << 4*nTris << std::endl ;



	CellsContainerIterator cellIterator = m_mesh->GetCells()->Begin();
	while (cellIterator!= m_mesh->GetCells()->End())
	{
		CellType * tri = cellIterator.Value();

		typedef CellType::PointIdIterator PointIdIterator;
		PointIdIterator  pointIdIterator=tri->PointIdsBegin();
		outfile<< "3 " ;
		while(pointIdIterator != tri->PointIdsEnd())
		{
			outfile << *pointIdIterator<< " ";
			++pointIdIterator;
		}

		outfile<<std:: endl ;
		++cellIterator;
	}


	//write normals

	if (normals->GetNumberOfPoints()>0){

		outfile << std::endl << "POINT_DATA " << nPts ;
		outfile << std::endl << "NORMALS normals double"<<std::endl;


		PointsContainerIterator normalIterator = normals->GetPoints()->Begin();
		while (normalIterator!=normals->GetPoints()->End())
		{
			PointType n = normalIterator.Value();
			outfile << n[0]<< " " <<n[1]<< " " << n[2] << std::endl;
			++normalIterator;
		}

	}

	outfile.close () ;
	return true;

}


bool writePts(MeshType::Pointer  m_Pts,char * fileName){
	//write the resulting object into a  lpts file

//TODO: check the file type

	std::ofstream outfile ;
	int nPts ;

	std::string outfileName(fileName) ;

	outfile.open ( outfileName.c_str() ) ;

	nPts = m_Pts->GetNumberOfPoints();
	

	PointsContainerIterator pointIterator = m_Pts->GetPoints()->Begin();
	while (pointIterator!=m_Pts->GetPoints()->End())
	{
		PointType p = pointIterator.Value();
		outfile << p[0]<< " " <<p[1]<< " " << p[2] << std::endl;
		++pointIterator;
	}

	outfile.close () ;
	return true;

}



void printUsage()
{
	appout << "Usage: " 
		<< PROGRAM_NAME << "-s ShapeStatsFileName -m MeanMeshFileName  -i ImageFileName	 -o OutModelFileName [Options]...  "
		<< std::endl << std::endl;
	appout << "Options:" 
		<< std::endl
		<< "  -h, --help   show this help"
		<< std::endl
		<< "  -v, --verbose  verbose output"
		<< std::endl
		<< "  -pcaStats, shape statistics file name"
		<< std::endl
		<< "  -refMesh, mesh file name [if mesh model is used]"
		<< std::endl
		<< "  -image, image file name"
		<< std::endl
		<< "  -output, output model file name"
		<< std::endl
		<<"   -numIters, maximum of iterations: default 100"
		<< std::endl
		<<"   -weightShapePrior, weight of shape prior: default 100.0 [weight of the image match is 1.0]"
		<< std::endl
		<<"   -numPCs, number of PCs: default 3"
		<< std::endl
		<<"   -stepSizeD, step size for calculating derivative during optimization: default 0.1"
		<< std::endl
		<<"   -gaussSigma, gaussian sigma for calculating image gradient feature: default 7";
	appout << std::endl << std::endl;
	appout << "Xiaoxiao Liu @ UNC-ChapelHill (sharonxx@cs.unc.edu)" << std::endl;
}






int main(int argc, const char *argv[])
{
	//Input Paramters
	char * ImageFileName =  ipGetStringArgument(argv, "-image", NULL);  
	char * PCAStatsFileName = ipGetStringArgument(argv, "-pcaStats", NULL);  
	char * refMeshFileName = ipGetStringArgument(argv, "-refMesh", NULL); 
	char * OutputModelFileName = ipGetStringArgument(argv, "-output", NULL);  



	// make sure the key arguments are valid
	if (  ipExistsArgument(argv, "--help") ||ipExistsArgument(argv, "-h") ||!ImageFileName||!PCAStatsFileName ||!OutputModelFileName)
	{  
		printUsage(); 
		return 0;

	}

    bool VERBOSE = false;
    if (ipExistsArgument(argv, "-v"))
		VERBOSE = true;


	// default settings
	int max_iters = ipGetIntArgument(argv, "-numIters", 100);  
	double WEIGHT_shapePrior = ipGetDoubleArgument(argv, "-weightShapePrior", 100); 
	int NUMofPCs = ipGetIntArgument(argv, "-numPCs", 3);  
	float StepSize = ipGetFloatArgument(argv, "-stepSizeD", 0.1); 
	float SigmaGaussian = ipGetFloatArgument(argv, "-gaussSigma", 10); 


	//read PCA statistics
	
	PCAStats * pcaStats = new PCAStats();
    pcaStats->readStats(PCAStatsFileName);

	std::cout<<"PCA statstics loaded:"<<PCAStatsFileName<<std::endl;



	
	//read image
	typedef itk::ImageFileReader< ImageType > ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(ImageFileName);
	reader->Update();
	ImageType::Pointer image = reader->GetOutput();
	std::cout<<"target image loaded:"<<ImageFileName<<std::endl;


	//setup DeformableMeshSegWithPirorclass
	DeformWithPrior  * deformer =  new DeformWithPrior();

    //read reference mesh model, if using mesh models for shape
    int shapeTYPE;
	if (refMeshFileName){
	   shapeTYPE = MESH;

	   //load refecemesh, to get the trianglations
	   MeshType::Pointer refMesh = readMesh(refMeshFileName);
	   std::cout<<"reference mesh model loaded :"<<refMeshFileName<<std::endl;
	   deformer->initialize( pcaStats,refMesh, image,SigmaGaussian);
	}
	else{
		shapeTYPE = PTS;
		//mean pts contained in PCA
		deformer->initialize( pcaStats, image,SigmaGaussian);
	}



	//params setting
	deformer->setMaxInterations(max_iters);

	deformer->setNumberOfUsedPCs(NUMofPCs);

	deformer->getCostFunction()->setShapePriorWeight(WEIGHT_shapePrior);

	deformer->getCostFunction()->setImageMatchWeight(1.0);

	deformer->getCostFunction()->setStepSizeForDerivative(StepSize);

	deformer->getCostFunction()->setVerbose(VERBOSE);


	//run program: deformation starts using conjugate gradient methods
	std::cout<<"deformation starts..."<<std::endl;

	deformer->Deform();

	std::cout<<"deformation done."<<std::endl;


	//Output: generate and write the final deformed mesh 
	MeshType::Pointer targetShape = deformer->getTargetObject();
	MeshType::Pointer normals = deformer->getCostFunction()->getNormals();



	if (shapeTYPE == PTS)
	{
		writePts(targetShape, OutputModelFileName);
	}
	if (shapeTYPE ==MESH)
	{
		writeMesh(targetShape,normals, OutputModelFileName);
	}
	


	//clean up memory
	return 0;

}

