#include "argio.hh"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <algorithm>

#include "itkMeshTovtkPolyData.h"
#include "vtkPolyDataToitkMesh.h"
//#include "itkMesh3DProcrustesAlignFilter.h"

#include "vtkPolyData.h"

#include "vtkPolyDataNormals.h"
#include "vtkPolyDataMapper.h"
#include "vtkDataSet.h"
#include "vtkPointData.h"

#include <itkDefaultDynamicMeshTraits.h>
#include <itkMetaMeshConverter.h>
#include <itkVariableLengthVector.h>
#include <itkMetaArrayReader.h>
#include <itkMetaArrayWriter.h>
#include <itkSpatialObjectWriter.h>
#include <itkSpatialObjectReader.h>
#include <itkSceneSpatialObject.h>

#include <itkLinearInterpolateImageFunction.h>
#include <itkContinuousIndex.h>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include <itkImageFileWriter.h>
#include <itkDanielssonDistanceMapImageFilter.h>
#include "itkBinaryThresholdImageFilter.h"

#include "vtkPolyDataReader.h"
#include "vtkCurvatures.h"
#include "vtkDoubleArray.h"

// takes in two itk mesh files (.meta) and produces a .txt file containing 
// the difference between the two meshes

// itk typedefs 
//
// typedef vtkPolyDataReader
typedef itk::DefaultDynamicMeshTraits < float, 3, 3 ,float, float> MeshTraitsType ; 
typedef itk::Mesh < float, 3, MeshTraitsType > MeshType ;
typedef itk::MeshSpatialObject < MeshType > MeshSOType ;
typedef itk::MetaMeshConverter < 3, float, MeshTraitsType > MeshConverterType ;
typedef MeshTraitsType::PointType PointType;

typedef MeshTraitsType::CellType CellType;

typedef itk::SpatialObjectWriter<3,float,MeshTraitsType> MeshWriterType;
typedef itk::SceneSpatialObject<3> SceneSpatialObjectType;  

typedef itk::DefaultDynamicMeshTraits<double, 3, 3,double,double> TriangleMeshTraits;
typedef itk::Mesh<double,3, TriangleMeshTraits> TriangleMeshType;
typedef TriangleMeshTraits::PointType PointTriangleType;
typedef itk::MetaMeshConverter < 3, double,TriangleMeshTraits > TriangleMeshConverterType ;
typedef itk::MeshSpatialObject < TriangleMeshType > TriangleMeshSOType ;  

typedef float ValueType;
typedef itk::VariableLengthVector<ValueType>       MeasurementVectorType;  
typedef itk::MetaArrayReader                       VectorReaderType;
typedef itk::MetaArrayWriter                       VectorWriterType;

//typedef itk::Mesh3DProcrustesAlignFilter<MeshType, MeshType> ProcrustesFilterType;
typedef itk::SpatialObjectReader<3,float,MeshTraitsType> ReaderType;

typedef MeshType::CellType CellType;
typedef itk::LineCell< CellType > LineType;
typedef itk::TriangleCell<CellType> TriangleType;
typedef MeshType::CellsContainer::ConstIterator CellIterator;
typedef  MeshType::Pointer MeshPointer;

typedef float PixelType;
typedef itk::Image < PixelType, 3 > ImageType;
typedef itk::LinearInterpolateImageFunction < ImageType,double > InterpolatorType;
typedef itk::ImageFileReader< ImageType > ImageReaderType;
typedef itk::ImageFileWriter< ImageType > ImageWriterType;

typedef itk::Image<short, 3 > ImageTypeShort;
typedef itk::BinaryThresholdImageFilter<ImageType,ImageTypeShort> Threshold;
typedef itk::DanielssonDistanceMapImageFilter<ImageTypeShort, ImageTypeShort> Distance;


int main(int argc, const char **argv)
{
	// make sure the arguments are valid
	if ( argc < 3 || ipExistsArgument(argv, "-usage") || ipExistsArgument(argv, "-help") )
	{
		std::cout <<"Usage:" << std::endl ;
		std::cout << argv[0] << " inputmesh/inputarray OutputFileName [options]" << std::endl ;
		std::cout <<"     -subtract <meshfile>      Subtract mesh from inputmesh, write a KWMeshVisu readable text file" << std::endl ;
		std::cout <<"     -magnitude                Magnitude of the input metaArray file (mvh/mva) and writes a KWMeshVisu readable file" << std::endl ;
		std::cout <<"     -scaleMVA <double>        Scales the input metaArray file (mvh/mva) and writes a KWMeshVisu readable file" << std::endl ;
		std::cout <<"     -scaleMesh <double>       Scales the input mesh file" << std::endl ;
		std::cout <<"     -avgMesh <Meshfile1> <Meshfile2> ..."  << std::endl;    
		std::cout <<"          Compute the average mesh from inputmesh file1, file2..." << std::endl;
		std::cout <<"     -ave <Vectorfile1> <Vectorfile2> ..."  << std::endl;
		std::cout <<"          Compute the average vector field from file1, file2... generated with -substract" << std::endl;
		std::cout <<"     -normave  <Vectorfile1> <Vectorfile2> ..." << std::endl;
		std::cout <<"          Works as the \"-ave\" option, but the average vector are projected on the normal at each point " << std::endl;
		std::cout <<"     -InvVect <VectorFile>     Invert all the vectors created with the -substract option and write a KWMeshVisu readable file" << std::endl;
		std::cout <<"     -magdir <VectorFile>      Compute the signed magnitude of each of the vector from the vector field.(+ if in the normal direction, - otherwise)" << std::endl;
		std::cout <<"     -applyVec <VectorFile> Deforme the mesh according to the vector field specified as input" << std::endl;
		std::cout <<"     -meshValues                Find the points and cells in a mesh. The outputfile is a textfile with the values" << std::endl;
		std::cout <<"      -avgGaussMesh <Meshfile1> <Meshfile2> ... -gaussMeshPara <mean>,<stdev>,<val1>,<val2>,... "  << std::endl; 
		std::cout <<"          Compute the gaussian average for mesh files." << std::endl;
		std::cout <<"          The first parameter is the average, then the standard deviation of the Gaussian model and the rest are the values associated with the files" << std::endl;
		std::cout <<"     -avgGaussKWM <txtfile1> <txtfile2>... -gaussKWMPara <mean>,<stdev>,<val1>,<val2>,... "<< std::endl;
		std::cout <<"          Compute the gaussian average for KWMeshVisu files." << std::endl;
		std::cout <<"          The first parameter is the average, then the standard deviation of the Gaussian model and the rest are the values associated with the files" << std::endl;
		std::cout <<"     -alignMesh <Meshfile1> <Meshfile2> ... Align all of the meshes"  << std::endl; 
		std::cout <<"     -BadTriangle <thresh value> [-correctMesh correctFilename] "  << std::endl; 
		std::cout <<"          Find the bad triangles in a Mesh. The <thresh value> is the value of the threshFactor to calculate the standard deviation for the bad triangles. The output is a KWMeshVisu text file with the values of the average of the triangles of the mesh "<< std::endl;
		std::cout <<"          To have a new Mesh with the correct triangles -correctMesh "<<  std::endl; 
		std::cout <<"     -extraction extractFilename [-extractClosest] "<< std::endl; 
		std::cout<< "          To extract an attribute.The Input is the Mesh, the extractFilename is the attribute image and the Output is a KWMeshVisu text file with the attribute extraction"<< std::endl; 
		std::cout <<"     -value <file1> <file2>... "<< std::endl;
		std::cout <<"          Extract the 5th column from a textfile and write a KWMeshVisu file with the values obtained"<<  std::endl;
		std::cout <<"     -subKWM <textname>...  Difference between 2 KWMeshVisu files"<< std::endl; 
		std::cout <<"     -MaxColor <textfile>...  Compare each points in every files, find a max for every points, keep 5% near the max, the other values will be 0"<< std::endl; 
		std::cout <<"     -dist_absolute <textfile>,<textfile>...  -result_relative <textfile>,<textfile>... Absolute distance map between  KWMeshVisu files"<< std::endl; 
		std::cout <<"     -dist_relative <textfile>,<textfile>...  -result_absolute <textfile>,<textfile>...Relative distance map between KWMeshVisu files (values between -1 & 1)"<< std::endl; 
		std::cout <<"     -label <textfile>... Separate every labels, find the mean..."<< std::endl; 
		std::cout <<"     -color -val <number_of_label>,<value_label>... -oldval <number_of_old_label>,<old_value_label>..."<< std::endl; 
		std::cout <<"          To change the value of labels to see the evolution with KWMeshVisu. "<<  std::endl; 
		std::cout <<"          Value_label is when the label grow up. "<<  std::endl; 
		std::cout <<"          Old_value_label is for the label wich has already grown up."<<  std::endl; 
		std::cout <<"     -first <textfile>... Convert a column file into a line file with a comma between every values"<< std::endl; 
		std::cout <<"     -v                        Verbose output" << std::endl;
		//cchou MC2Origin       
		std::cout <<"	-MC2Origin		Translate the Center of Mass to the Origin"<< std::endl;     
		//xxLiu
		std::cout <<" -gaussianCurvature <inputFile >  <outputFile>   : extract the gaussian curvatures into an attribut file "<<std::endl;	
		std::cout <<" -meanCurvature <inputFile >  <outputFile>   : extract the mean curvatures into an attribut file "<<std::endl;
		return 0 ;
	}

	// get the arguments  
	char * inputFilename = strdup(argv[1]);  
	char * outputFilename = strdup(argv[2]); 
	std::string outputFilename2;
	int DotLoc;

	bool colorOn=  ipExistsArgument(argv, "-color");
	bool valOn=  ipExistsArgument(argv, "-val");
	float *num = new float [30];
	float *old_num = new float [30];
	if(valOn)
	{
		char * tmp_str = ipGetStringArgument(argv, "-val", NULL);
	}
	bool oldvalOn=  ipExistsArgument(argv, "-oldval");
	if (oldvalOn)
	{
		char * tmp_str = ipGetStringArgument(argv, "-oldval", NULL);
	}

	const int maxNumFiles = 1000;
	int nbfile = 0;
	char * files[maxNumFiles];

	bool MaxColorOn= ipExistsArgument(argv, "-MaxColor");  
	if (MaxColorOn) 
		nbfile = ipGetStringMultipArgument(argv, "-MaxColor", files, maxNumFiles); 
	std::vector<std::string> MaxFiles;
	if(MaxColorOn) {
		for(int i = 0 ; i < nbfile ; i++) 
			MaxFiles.push_back(files[i]);
	}

	bool distAbsOn= ipExistsArgument(argv, "-dist_absolute");  
	if (distAbsOn) 
		nbfile = ipGetStringMultipArgument(argv, "-dist_absolute", files, maxNumFiles); 
	std::vector<std::string> distFiles;
	if(distAbsOn) {
		for(int i = 0 ; i < nbfile ; i++) 
			distFiles.push_back(files[i]);
	}
	bool resultAbsOn= ipExistsArgument(argv, "-result_absolute");  
	if (distAbsOn) 
		nbfile = ipGetStringMultipArgument(argv, "-result_absolute", files, maxNumFiles); 
	std::vector<std::string> resultFiles;
	if(resultAbsOn) {
		for(int i = 0 ; i < nbfile ; i++) 
			resultFiles.push_back(files[i]);
	}

	bool distRelOn= ipExistsArgument(argv, "-dist_relative");  
	if (distRelOn) 
		nbfile = ipGetStringMultipArgument(argv, "-dist_relative", files, maxNumFiles); 
	std::vector<std::string> distReFiles;
	if(distRelOn) {
		for(int i = 0 ; i < nbfile ; i++) 
			distReFiles.push_back(files[i]);
	}
	bool resultRelOn= ipExistsArgument(argv, "-result_relative");  
	if (distRelOn) 
		nbfile = ipGetStringMultipArgument(argv, "-result_relative", files, maxNumFiles); 
	std::vector<std::string> restFiles;
	if(resultRelOn) {
		for(int i = 0 ; i < nbfile ; i++) 
			restFiles.push_back(files[i]);
	}

	bool subtractOn = ipExistsArgument(argv, "-subtract");
	char * subtractFile = ipGetStringArgument(argv, "-subtract", NULL); 
	bool magnitudeOn = ipExistsArgument(argv, "-magnitude");
	bool scaleMVAOn = ipExistsArgument(argv, "-scaleMVA");
	bool scaleOn = ipExistsArgument(argv, "-scaleMesh");
	double scaleFactor = 1.0;
	if (scaleMVAOn)
		scaleFactor = ipGetDoubleArgument(argv,"-scaleMVA",1.0);
	else 
		scaleFactor = ipGetDoubleArgument(argv,"-scaleMesh",1.0);
	bool aveOn = ipExistsArgument(argv, "-ave");
	bool normaveOn = ipExistsArgument(argv, "-normave");
	bool averageOn = aveOn || normaveOn;
	if (aveOn) 
		nbfile = ipGetStringMultipArgument(argv, "-ave", files, maxNumFiles); 
	if (normaveOn)
		nbfile = ipGetStringMultipArgument(argv, "-normave", files, maxNumFiles);
	std::vector<std::string> AveFiles;
	if(averageOn) {
		for(int i = 1 ; i < nbfile ; i++) 
			AveFiles.push_back(files[i]);
	}

	bool meshValueOn =  ipExistsArgument(argv, "-meshValues");
	bool alignMeshOn =  ipExistsArgument(argv, "-alignMesh");
	if (alignMeshOn)
		nbfile = ipGetStringMultipArgument(argv, "-alignMesh", files, maxNumFiles); 
	std::vector<std::string> AlignMeshFiles;
	std::vector<std::string> AlignMeshOutputFiles;
	if(alignMeshOn) {
		for(int i = 0 ; i < nbfile ; i++) 
			AlignMeshFiles.push_back(files[i]);
		AlignMeshOutputFiles.assign(AlignMeshFiles.begin(), AlignMeshFiles.end());
		std::vector<std::string>::iterator inputFile= AlignMeshOutputFiles.begin();	
		AlignMeshOutputFiles.insert(inputFile,inputFilename );
		for( int i = 0; i < nbfile+1; i++ ) {
			AlignMeshOutputFiles[i].erase( AlignMeshOutputFiles[i].size()-5);
			AlignMeshOutputFiles[i].insert(AlignMeshOutputFiles[i].size(),"_align.meta");
		}
	}

	bool valueOn = ipExistsArgument(argv, "-value");
	///////////////uncomment if you want to select with the first column (like only one label)
	//double Index = ipGetDoubleArgument(argv, "-value",5.0);

	bool labelOn = ipExistsArgument(argv, "-label");
	char * labelFilename = ipGetStringArgument(argv, "-label",NULL);

	bool subKWMOn = ipExistsArgument(argv, "-subKWM");
	char * subKWMFilename = ipGetStringArgument(argv, "-subKWM",NULL);

	bool badtriangleMeshOn = ipExistsArgument(argv, "-BadTriangle");
	double threshFactor = ipGetDoubleArgument(argv,"-BadTriangle",5.0);
	char * correctFilename = ipGetStringArgument(argv, "-BadTriangle",NULL);
	bool correctbadtriangleMeshOn = ipExistsArgument(argv, "-correctMesh");

	bool extractionOn = ipExistsArgument(argv, "-extraction");
	char * extractFilename = ipGetStringArgument(argv, "-extraction", NULL);
	bool extractClosestOn= ipExistsArgument(argv, "-extractClosest");

	bool avgMeshOn = ipExistsArgument(argv, "-avgMesh");
	if (avgMeshOn) 
		nbfile = ipGetStringMultipArgument(argv, "-avgMesh", files, maxNumFiles); 
	std::vector<std::string> AvgMeshFiles;
	if(avgMeshOn) {
		for(int i = 0 ; i < nbfile ; i++) 
			AvgMeshFiles.push_back(files[i]);
	}

	bool firstOn = ipExistsArgument(argv, "-first");  
	if (firstOn) {
		nbfile = ipGetStringMultipArgument(argv, "-first", files, maxNumFiles);
	}
	std::vector<std::string> firstfiles;
	if(firstOn) {
		for(int i = 0 ; i < nbfile ; i++) 
			firstfiles.push_back(files[i]);
	}

	bool avgGaussMeshOn = ipExistsArgument(argv, "-avgGaussMesh");  
	if (avgGaussMeshOn) {
		nbfile = ipGetStringMultipArgument(argv, "-avgGaussMesh", files, maxNumFiles);
	}
	std::vector<std::string> AvgGaussMeshFiles;
	if(avgGaussMeshOn) {
		for(int i = 0 ; i < nbfile ; i++) 
			AvgGaussMeshFiles.push_back(files[i]);
	}
	bool gaussMeshParaOn =  ipExistsArgument(argv, "-gaussMeshPara");
	int numparams = nbfile+1+2;
	float *GaussParams = new float [numparams];
	float *Age = new float [nbfile+1];
	if(gaussMeshParaOn) {
		char * tmp_str = ipGetStringArgument(argv, "-gaussMeshPara", NULL);
		if (tmp_str) {
			int nbage = ipExtractFloatTokens(GaussParams, tmp_str, numparams);
			if (numparams != nbage) {
				cerr << "we need" << numparams<< " comma separated entries" << endl;
				exit(1);
			} 
			else {
				nbage = nbage - 2;
				for(int j = 0 ; j < nbage ; j++) 
					Age[j] = GaussParams[j+2];
			}
		}
		else {
			std::cerr << " error: need gaussian parameters " << std::endl;
			exit(-1);
		}
	}

	bool avgGaussKWMOn =  ipExistsArgument(argv, "-avgGaussKWM");  
	if(avgGaussKWMOn){
		nbfile = ipGetStringMultipArgument(argv, "-avgGaussKWM", files, maxNumFiles);
	}
	std::vector<std::string> avgGaussKWMFiles;
	if(avgGaussKWMOn) {
		for(int i = 0 ; i < nbfile ; i++) 
			avgGaussKWMFiles.push_back(files[i]);
	}

	bool gaussKWMParaOn =  ipExistsArgument(argv, "-gaussKWMPara");
	int numpara = nbfile+1+2;
	float *GaussKWMParams = new float [numpara];
	float *AgeKWM = new float [nbfile+1];
	if(gaussKWMParaOn) {
		char * tmp_str = ipGetStringArgument(argv, "-gaussKWMPara", NULL);
		if (tmp_str) {
			int nbage = ipExtractFloatTokens(GaussKWMParams, tmp_str, numpara);
			if (numpara != nbage) {
				cerr << "we need" << numpara<< " comma separated entries" << endl;
				exit(1);
			} 
			else {
				nbage = nbage - 2;
				for(int j = 0 ; j < nbage ; j++) 
					AgeKWM[j] = GaussKWMParams[j+2];
			}
		}
		else {
			std::cerr << " error: need gaussian parameters " << std::endl;
			exit(-1);
		}
	}

	char *VectFile;
	bool invvectOn = ipExistsArgument(argv, "-InvVect");
	VectFile = ipGetStringArgument(argv, "-InvVect", NULL); 

	bool magdirOn = ipExistsArgument(argv, "-magdir");
	VectFile = ipGetStringArgument(argv, "-magdir", NULL); 

	bool applyvecOn =  ipExistsArgument(argv, "-applyVec");
	VectFile = ipGetStringArgument(argv, "-applyVec", NULL); 

	bool debug = ipExistsArgument(argv, "-v");

	//cchou: MC2OriginOn
	bool MC2OriginOn = ipExistsArgument(argv, "-MC2Origin");

	bool GetGaussianCurvaturesOn = ipExistsArgument(argv, "-gaussianCurvature");

	bool GetMeanCurvaturesOn = ipExistsArgument(argv, "-meanCurvature");
	if (subtractOn) {
		// read in the input files
		MeshConverterType * converter = new MeshConverterType () ;

		if (debug)   std::cout << "Reading input mesh " << inputFilename << std::endl;
		MeshSOType::Pointer inputMeshSO = converter->ReadMeta (inputFilename) ;
		MeshType::Pointer inputMesh = inputMeshSO->GetMesh() ;

		if (debug) std::cout << "subtracting Mesh " << std::endl;
		MeshSOType::Pointer sutractMeshSO = converter->ReadMeta (subtractFile) ;
		MeshType::Pointer subtractMesh = sutractMeshSO->GetMesh() ;

		// make sure the two meshes have the same number of verts
		if (inputMesh->GetNumberOfPoints() != subtractMesh->GetNumberOfPoints()) {
			std::cout << "Meshes do not have same number of vertices: " << inputMesh->GetNumberOfPoints() <<
				" vs " << subtractMesh->GetNumberOfPoints() << std::endl;
			exit(-1);
		}

		std::ofstream outfile ;
		outfile.open ( outputFilename ) ;

		PointType *point1 = new PointType ;
		PointType *point2 = new PointType ;
		itk::Vector<float, 3> diff ;

		// print the header
		outfile << "NUMBER_OF_POINTS=" << inputMesh->GetNumberOfPoints() << std::endl ;
		outfile << "DIMENSION=" << 3 << std::endl ;
		outfile << "TYPE=Vector" << std::endl ;

		for ( unsigned int i = 0 ; i < inputMesh->GetNumberOfPoints()  ; i++ )
		{
			if ( inputMesh->GetPoint ( i, point1 ) && subtractMesh->GetPoint (i, point2) )
			{
				diff = *point1 - *point2 ;
				outfile << diff[0] << " " << diff[1] << " " << diff[2] << std::endl ;
			}
			else
				return 0 ;
		}
		outfile.close ();

		outputFilename2 = std::string(outputFilename);
		DotLoc = outputFilename2.find(".txt",0);
		if(DotLoc != -1)
			outputFilename2.insert(DotLoc,".mag");
		else
			outputFilename2 = std::string(outputFilename) + std::string(".mag");

		outfile.open ( outputFilename2.c_str() ) ;
		//  outputFilename => Filename
		// print the header
		outfile << "NUMBER_OF_POINTS=" << inputMesh->GetNumberOfPoints() << std::endl ;
		outfile << "DIMENSION=" << 1 << std::endl ;
		outfile << "TYPE=Scalar" << std::endl ;

		for ( unsigned int i = 0 ; i < inputMesh->GetNumberOfPoints()  ; i++ )
		{
			if ( inputMesh->GetPoint ( i, point1 ) && subtractMesh->GetPoint (i, point2) )
			{
				diff = *point1 - *point2 ;
				outfile << sqrt(diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]) << std::endl ;
			}
			else
				return 0 ;
		}
		outfile.close () ;

		delete ( point1 ) ;
		delete ( point2 ) ;
	} else if (scaleOn) {

		MeshConverterType * converter = new MeshConverterType () ;
		if (debug)   std::cout << "Reading input mesh " << inputFilename << std::endl;
		MeshSOType::Pointer inputMeshSO = converter->ReadMeta (inputFilename) ;
		MeshType::Pointer inputMesh = inputMeshSO->GetMesh() ;
		MeshType::PointsContainerPointer points = inputMesh->GetPoints();

		MeshType::PointsContainerPointer pointsTmp = MeshType::PointsContainer::New();
		for (unsigned int pointID = 0; pointID < points->Size();pointID++) {
			PointType curPoint =  points->GetElement(pointID);
			PointType vert;
			for (unsigned int dim = 0; dim < 3; dim++) 
				vert[dim] = curPoint[dim] * scaleFactor;
			pointsTmp->InsertElement(pointID, vert);
		}

		inputMesh->SetPoints(pointsTmp); 
		MeshWriterType::Pointer writer = MeshWriterType::New();
		writer->SetInput(inputMeshSO);
		writer->SetFileName(outputFilename);
		writer->Update();

	} else if (magnitudeOn) {
		if (debug) std::cout << "Magnitude of Vector field" << std::endl;
		MeasurementVectorType mv;
		// load the metaArray
		VectorReaderType::Pointer vectorReader = VectorReaderType::New();
		vectorReader->SetFileName(inputFilename);
		vectorReader->Update();
		vectorReader->GetOutput(MET_FLOAT, &mv, true);

		std::ofstream outfile ;
		outfile.open ( outputFilename ) ;

		// print the header
		outfile << "NUMBER_OF_POINTS = " << mv.GetSize()/3 << std::endl ;
		outfile << "DIMENSION = " << 1 << std::endl ;
		outfile << "TYPE=Scalar" << std::endl ;

		for ( unsigned int i = 0 ; i < mv.GetSize()  ; i += 3)
		{
			double length = sqrt(mv[i] * mv[i] + mv[i+1] * mv[i+1] + mv[i+2] * mv[i+2]) ;
			if (scaleMVAOn) {
				outfile << length * scaleFactor  << std::endl ;
			} else {
				outfile << length << std::endl ;
			}
		}

	} else if (scaleMVAOn) {
		if (debug) std::cout << "Scaling of Vector field" << std::endl;
		MeasurementVectorType mv;
		// load the metaArray
		VectorReaderType::Pointer vectorReader = VectorReaderType::New();
		vectorReader->SetFileName(inputFilename);
		vectorReader->Update();
		vectorReader->GetOutput(MET_FLOAT, &mv, true);

		std::ofstream outfile ;
		outfile.open ( outputFilename ) ;

		// print the header
		outfile << "NUMBER_OF_POINTS = " << mv.GetSize()/3 << std::endl ;
		outfile << "DIMENSION = " << 3 << std::endl ;
		outfile << "TYPE = Vector" << std::endl ;

		for ( unsigned int i = 0 ; i < mv.GetSize()  ; i += 3)
		{
			outfile << mv[i] * scaleFactor << " " << mv[i+1] * scaleFactor << " " 
				<< mv[i+2] * scaleFactor << std::endl ;
		}

	}if (meshValueOn){

		int num =0;
		//read input mesh
		MeshConverterType * converter = new MeshConverterType () ;
		MeshSOType::Pointer SOMesh = converter->ReadMeta (inputFilename) ;
		MeshType::Pointer mesh = SOMesh->GetMesh() ;
		MeshType::PointsContainerPointer Points = mesh->GetPoints();
		int numberOfPoints = Points->Size();

		//open output file
		std::ofstream outfileNor;
		outfileNor.open ( outputFilename ) ;
		outfileNor << "Number of  Points = "<<numberOfPoints <<std::endl;

		//find points
		for (unsigned int pointID = 0; pointID < Points->Size();pointID++){
			PointType curPoint =  Points->GetElement(pointID);
			outfileNor << "Points "<< pointID<< " " << curPoint<< std::endl;
		}

		//find number&points of cells
		typedef CellType::PointIdIterator PointIdIterator;
		CellIterator cellIterator = mesh->GetCells()->Begin();
		CellIterator cellEnd      = mesh->GetCells()->End();
		PointType curPoint;
		while( cellIterator != cellEnd )
		{
			CellType * cell = cellIterator.Value();
			TriangleType * line = dynamic_cast<TriangleType *> (cell);
			LineType::PointIdIterator pit = line->PointIdsBegin();
			int pIndex1,pIndex2,pIndex3;

			pIndex1= *pit;
			++pit;
			pIndex2= *pit;
			++pit;
			pIndex3= *pit;

			outfileNor << "Cells  "<<num<<" "<<pIndex1<< " "<< pIndex2 <<" "<<pIndex3<< std::endl;

			++cellIterator;
			++num;
		}

		//close outputFile
		outfileNor.close();

	} else if(extractionOn){

		PointType point; 
		PointType pixel;
		InterpolatorType::ContinuousIndexType point_index;   

		//read input 
		MeshConverterType * converter = new MeshConverterType();
		MeshSOType::Pointer SOMesh = converter->ReadMeta (inputFilename);
		MeshType::Pointer mesh = SOMesh->GetMesh();
		MeshType::PointsContainerPointer points = mesh->GetPoints();
		int numPoints = points->Size();
		double **pointIndex;
		pointIndex = new double* [numPoints];
		for (int c = 0; c < numPoints; c++)
		{    
			pointIndex[c]= new double [3];
			for (int d = 0; d < 3 ; d++)
				pointIndex[c][d]=0.0;
		}

		double extractValue[numPoints];
		for (int i = 0; i < numPoints ; i++) 
			extractValue[i]=0.0;

		//load attribute image
		ImageReaderType::Pointer reader = ImageReaderType::New();
		reader->SetFileName(extractFilename);
		reader->Update();

		ImageType::Pointer image = reader->GetOutput();
		InterpolatorType::Pointer   interpolator  = InterpolatorType::New();
		interpolator->SetInputImage(image);
		Distance::VectorImageType::Pointer closestVectorImage;

		if( extractClosestOn) {
			Threshold::Pointer threshold = Threshold::New();
			threshold->SetInput(image);
			threshold->SetLowerThreshold (1);
			threshold->SetUpperThreshold (32000);

			Distance::Pointer distance = Distance::New();
			distance->UseImageSpacingOn();
			distance->SquaredDistanceOn();
			distance->SetInput(threshold->GetOutput());
			distance->Update();
			closestVectorImage = distance->GetVectorDistanceMap();
		}

		// for all points in the mesh, extract image attribute
		for (int index = 0; index < numPoints ; index++) {

			//Origin & Spacing
			ImageType::SpacingType  spacing = image->GetSpacing();
			ImageType::PointType origin = image->GetOrigin();

			PointType curPoint =  points->GetElement(index);
			for (unsigned int dim = 0; dim < 3; dim++) {
				pixel[dim] = (curPoint[dim]- origin[dim]) / spacing[dim];
				point_index[dim]=pixel[dim];
				pointIndex[index][dim]=point_index[dim];
			}

			extractValue[index] = interpolator->EvaluateAtContinuousIndex(point_index);

			if(extractValue[index]==0 && extractClosestOn) {
				Distance::VectorImageType::PixelType vectorPixel;
				Distance::VectorImageType::IndexType vectorIndex;
				for (int indexDim = 0; indexDim < 3; indexDim++) {
					vectorIndex[indexDim] =  (long int) round( point_index[indexDim] );
				}

				vectorPixel = closestVectorImage->GetPixel(vectorIndex);

				for (unsigned int dim = 0; dim < 3; dim++) {
					point_index[dim]= vectorIndex[dim] + vectorPixel[dim];
					pointIndex[index][dim]=point_index[dim];
				}

				extractValue[index]=interpolator->EvaluateAtContinuousIndex(point_index);
				if (extractValue[index] == 0 ) std::cout << "error at extracting closest point" << std::endl;
			}
		}

		//write ouput KWMeshVisu text file, attribute extraction
		std::ofstream outfileNor;
		outfileNor.open ( outputFilename ) ;
		outfileNor << "NUMBER_OF_POINTS=" << numPoints << std::endl ;
		outfileNor << "DIMENSION=" << 1 << std::endl ;
		outfileNor << "TYPE=Scalar" << std::endl ;
		for (int i = 0; i < numPoints; i++)
			outfileNor  << extractValue[i] << std::endl;
		//close outputFile
		outfileNor.close();

		delete[] pointIndex;

	}else if(avgMeshOn) {

		MeshConverterType * converter = new MeshConverterType () ;

		//read input
		if (debug)   std::cout << "Reading  mesh " << inputFilename << std::endl;
		MeshSOType::Pointer SOMesh = converter->ReadMeta (inputFilename) ;
		MeshType::Pointer surfaceMesh = SOMesh->GetMesh() ;
		MeshType::PointsContainerPointer avgPoints = surfaceMesh->GetPoints();

		int numMeshes = AvgMeshFiles.size();

		// for all meshes
		for (int index = 0; index < numMeshes; index++) {
			try
			{
				if (debug)   std::cout << "Reading  mesh " << AvgMeshFiles[index] << std::endl;
				MeshSOType::Pointer inputMeshSO = converter->ReadMeta (AvgMeshFiles[index].c_str()) ;
				MeshType::Pointer inputMesh = inputMeshSO->GetMesh() ;
				MeshType::PointsContainerPointer points = inputMesh->GetPoints();
				MeshType::PointsContainerPointer pointsTmp = MeshType::PointsContainer::New();

				for (unsigned int pointID = 0; pointID < points->Size();pointID++) {
					PointType curPoint =  points->GetElement(pointID);
					PointType vert;
					PointType curPointAvg =  avgPoints->GetElement(pointID);
					for (unsigned int dim = 0; dim < 3; dim++) 
						vert[dim] = curPoint[dim] + curPointAvg[dim];
					pointsTmp->InsertElement(pointID, vert);
				}

				avgPoints = pointsTmp;
			}
			catch(itk::ExceptionObject ex)
			{
				std::cout<< "Error reading meshfile:  "<< AvgMeshFiles[index] << std::endl << "ITK error: " << ex.GetDescription()<< std::endl;
				exit(-3);
			}
		}

		MeshType::PointsContainerPointer pointsTmp = MeshType::PointsContainer::New();
		for (unsigned int pointID = 0; pointID < avgPoints->Size();pointID++) {
			PointType curPoint =  avgPoints->GetElement(pointID);
			PointType vert;
			for (unsigned int dim = 0; dim < 3; dim++) 
				vert[dim] = curPoint[dim] / (numMeshes + 1);
			pointsTmp->InsertElement(pointID, vert);
		}
		avgPoints = pointsTmp;

		surfaceMesh->SetPoints(avgPoints); 
		SOMesh->SetMesh(surfaceMesh);
		MeshWriterType::Pointer writer = MeshWriterType::New();
		writer->SetInput(SOMesh);
		writer->SetFileName(outputFilename);
		writer->Update();

	} /* else if(alignMeshOn) {

	     MeshConverterType * converter = new MeshConverterType () ;

	//read input
	MeshSOType::Pointer SOMesh = converter->ReadMeta (inputFilename) ;
	MeshType::Pointer surfaceMesh = SOMesh->GetMesh() ;
	MeshType::PointsContainerPointer Points = surfaceMesh->GetPoints();
	MeshType::Pointer newPoints = MeshType::New();

	int numMeshes = AlignMeshFiles.size();

	//align every mesh on the first one
	ProcrustesFilterType::Pointer procrustesFilter = ProcrustesFilterType::New();
	procrustesFilter->SetNumberOfInputs(numMeshes+1);
	procrustesFilter->SetUseInitialAverageOn();
	procrustesFilter->SetUseNormalizationOff();
	procrustesFilter->SetUseScalingOff();
	procrustesFilter->SetUseSingleIterationOn();

	procrustesFilter->SetInput(0,surfaceMesh);

	// read meshes and use Mesh3DProcrustesAlignFilter to align
	for (int index = 0; index < numMeshes; index++) {
	MeshSOType::Pointer inputMeshSO = converter->ReadMeta(AlignMeshFiles[index].c_str()) ;
	MeshType::Pointer inputMesh = inputMeshSO->GetMesh() ;
	procrustesFilter->SetInput(index+1, inputMesh);
	}

	procrustesFilter->Update();

	//create the align MeshFile for the input
	MeshWriterType::Pointer writer = MeshWriterType::New();
	MeshType::Pointer  RegisteredMesh = procrustesFilter->GetOutput(0);
	SOMesh->SetMesh(RegisteredMesh);
	writer->SetInput(SOMesh);
	writer->SetFileName(AlignMeshOutputFiles[0].c_str());
	writer->Update();

	//create align MeshFiles
	for (int index = 0; index < numMeshes; index++) {
	MeshType::Pointer  RegisteredMesh = procrustesFilter->GetOutput(index+1); 
	MeshSOType::Pointer inputMeshSO = converter->ReadMeta (AlignMeshFiles[index].c_str()) ;
	inputMeshSO->SetMesh(RegisteredMesh);
	writer->SetInput(inputMeshSO);
	writer->SetFileName(AlignMeshOutputFiles[index+1].c_str()); 
	writer->Update();
	}

	}*/else if(badtriangleMeshOn){

		//read input Mesh
		MeshConverterType * converter = new MeshConverterType();
		MeshSOType::Pointer SOMesh = converter->ReadMeta (inputFilename);
		MeshType::Pointer mesh = SOMesh->GetMesh();
		MeshType::PointsContainerPointer points = mesh->GetPoints();
		int numPoints = points->Size();
		double *avgTriangleArea = new double [numPoints];
		int *numTriangles = new int [numPoints];
		int *neighboor = new int [numPoints];
		for (int i = 0; i < numPoints; i++) {
			avgTriangleArea[i] = 0.0;
			numTriangles[i] = 0;
			neighboor[i] = 0;
		}
		int num=0;
		float std=0;
		double sum=0;
		float avgTriangle=0;
		double max=0, min=0, dif=0;
		double a=0, b=0, c=0;
		bool meshChanged =true;
		int z=0;

		while (meshChanged) {
			meshChanged = false;
			MeshType::PointType point;
			for (int dim = 0; dim < 3; dim++)
				point[dim] = 0;

			//z!=0 load file (correctFilename) wich is already correct to improve the correction
			if(z!=0){
				MeshConverterType * converter = new MeshConverterType();
				MeshSOType::Pointer SOMesh = converter->ReadMeta (correctFilename);
				MeshType::Pointer mesh = SOMesh->GetMesh();
				MeshType::PointsContainerPointer points = mesh->GetPoints();
			}

			typedef CellType::PointIdIterator PointIdIterator;
			CellIterator cellIterator = mesh->GetCells()->Begin();
			CellIterator cellEnd      = mesh->GetCells()->End();
			PointType curPoint;

			//calculate average triangle for all points, max, min
			while( cellIterator != cellEnd )
			{
				CellType * cell = cellIterator.Value();
				TriangleType * line = dynamic_cast<TriangleType *> (cell);
				LineType::PointIdIterator pit = line->PointIdsBegin();
				MeshType::PointType p1, p2, p3;
				int pIndex1,pIndex2,pIndex3;

				for (int dim = 0; dim < 3; dim++) {
					pIndex1=*pit;
					PointType curPoint =  points->GetElement(pIndex1);
					p1[dim]=curPoint[dim];
				}

				++pit;

				for (int dim = 0; dim < 3; dim++) {
					pIndex2=*pit;
					PointType curPoint =  points->GetElement(pIndex2);
					p2[dim]=curPoint[dim];
				}

				++pit;

				for (int dim = 0; dim < 3; dim++) {
					pIndex3=*pit;
					PointType curPoint =  points->GetElement(pIndex3);
					p3[dim]=curPoint[dim];
				}

				a = (p2[0]-p1[0])*(p2[0]-p1[0])+(p2[1]-p1[1])*(p2[1]-p1[1])+(p2[2]-p1[2])*(p2[2]-p1[2]);
				c = (p3[0]-p1[0])*(p3[0]-p1[0])+(p3[1]-p1[1])*(p3[1]-p1[1])+(p3[2]-p1[2])*(p3[2]-p1[2]);
				b = (p3[0]-p2[0])*(p3[0]-p2[0])+(p3[1]-p2[1])*(p3[1]-p2[1])+(p3[2]-p2[2])*(p3[2]-p2[2]);

				double triangleArea =sqrt(2*b*c+2*c*a+2*a*b-a*a-b*b-c*c)*0.25;

				avgTriangleArea[pIndex1] += triangleArea;
				avgTriangleArea[pIndex2] += triangleArea;
				avgTriangleArea[pIndex3] += triangleArea;
				numTriangles[pIndex1] ++;
				numTriangles[pIndex2] ++;
				numTriangles[pIndex3] ++;

				++cellIterator;
				++num;
			}
			//calculate the average of the triangles
			for (int i = 0; i < numPoints; i++) {
				avgTriangleArea[i] /= numTriangles[i];     
			}
			//find max, min of all of the triangle area and calculate the total average triangle
			max=avgTriangleArea[0];
			min=avgTriangleArea[0];
			for(int i=0;i<numPoints;i++) {
				if( avgTriangleArea[i]>max){
					max=avgTriangleArea[i];
				}
				if( avgTriangleArea[i]<min){
					min=avgTriangleArea[i];
				}
				avgTriangle= avgTriangle + avgTriangleArea[i];
			}
			avgTriangle=avgTriangle/numPoints;
			//calculate the standard deviation
			for(int i=0;i<numPoints;i++) {
				dif= avgTriangleArea[i] - avgTriangle;
				sum = sum + dif*dif;
			}
			std=sqrt(sum/(numPoints-1)) ;

			double areaThreshold = avgTriangle + threshFactor*std;
			int n=0;
			int numNeighboor=0;
			//correct points if it's necessary
			for(int i=0;i<numPoints;i++) {
				//if point i is a bad point
				if( avgTriangleArea[i]  >= areaThreshold ){
					for (int dim = 0; dim < 3; dim++) 
						point[dim] = 0;
					n=0;
					numNeighboor=0;
					for (int c = 0; c < numPoints; c++) 
						neighboor[c] = 0;
					PointType curPoint =  points->GetElement(i);
					CellIterator cellIterator = mesh->GetCells()->Begin();
					CellIterator cellEnd      = mesh->GetCells()->End();
					bool p2neighboor=0,p3neighboor=0;
					//find the neighboors and calculate the new point with the neighboors
					while( cellIterator != cellEnd )
					{
						CellType * cell = cellIterator.Value();
						TriangleType * line = dynamic_cast<TriangleType *> (cell);
						LineType::PointIdIterator pit = line->PointIdsBegin();
						MeshType::PointType p1, p2, p3,p;
						int pIndex1,pIndex2,pIndex3;
						//find the neighboors
						for (int dim = 0; dim < 3; dim++) {
							pIndex1=*pit;
							PointType curPoint =  points->GetElement(pIndex1);
							p1[dim]=curPoint[dim];
						}
						++pit;
						for (int dim = 0; dim < 3; dim++) {
							pIndex2=*pit;
							PointType curPoint =  points->GetElement(pIndex2);
							p2[dim]=curPoint[dim];
						}
						++pit;
						for (int dim = 0; dim < 3; dim++) {
							pIndex3=*pit;
							PointType curPoint =  points->GetElement(pIndex3);
							p3[dim]=curPoint[dim];
						}
						//if pIndex1 is the bad point
						if(pIndex1==i){
							//know if pIndex2 and pIndex3 are already neighboors
							for(int b=0; b< n;b++){
								if(pIndex2 == neighboor[b])
									p2neighboor=1;
								if(pIndex3 == neighboor[b])
									p3neighboor=1;
							}
							//nothing change
							if(p2neighboor && p3neighboor){
								for (int dim = 0; dim < 3; dim++) {
									point[dim] = point[dim];
								}
							}
							//if not calculate the new sum of neighboors, and the new num of neighboors
							if(!p3neighboor) {
								for (int dim = 0; dim < 3; dim++) {
									point[dim] += p3[dim];
								}
								++numNeighboor;
							}
							if(!p2neighboor){
								for (int dim = 0; dim < 3; dim++) {
									point[dim] += p2[dim];
								}
								++numNeighboor;
							}
							neighboor[n]=pIndex2;
							++n;
							neighboor[n]=pIndex3;
							++n;
						} else if(pIndex2==i){
							//know if pIndex1 and pIndex3 are already neighboors
							for(int b=0; b< n;b++){
								if(pIndex1 == neighboor[b])
									p2neighboor=1;
								if(pIndex3 == neighboor[b])
									p3neighboor=1;
							}
							//nothing change
							if(p2neighboor && p3neighboor){
								for (int dim = 0; dim < 3; dim++) {
									point[dim] = point[dim];
								}
							}
							//if not calculate the new sum of neighboors, and the num of neighboors
							if(!p3neighboor) {
								for (int dim = 0; dim < 3; dim++) {
									point[dim] += p3[dim];
								}
								++numNeighboor;
							}
							if(!p2neighboor){
								for (int dim = 0; dim < 3; dim++) {
									point[dim] += p1[dim];
								}
								++numNeighboor;
							}
							neighboor[n]=pIndex1;
							++n;
							neighboor[n]=pIndex3;
							++n;
						} else if(pIndex3==i){
							//know if pIndex1 and pIndex3 are already neighboors
							for(int b=0; b< n;b++){
								if(pIndex1 == neighboor[b])
									p2neighboor=1;
								if(pIndex2 == neighboor[b])
									p3neighboor=1;
							}
							//nothing change
							if(p2neighboor && p3neighboor){
								for (int dim = 0; dim < 3; dim++) {
									point[dim] = point[dim];
								}
							}
							//if not calculate the new sum of neighboors, and the num of neighboors
							if(!p3neighboor) {
								for (int dim = 0; dim < 3; dim++) {
									point[dim] += p2[dim];
								}
								++numNeighboor;
							}
							if(!p2neighboor){
								for (int dim = 0; dim < 3; dim++) {
									point[dim] += p1[dim];
								}
								++numNeighboor;
							}
							neighboor[n]=pIndex1;
							++n;
							neighboor[n]=pIndex2;
							++n;
						}
						++cellIterator;
						++num;
					}
					//calculate the new point : the average of the neighboors
					for (int dim = 0; dim < 3; dim++) {
						point[dim] /= numNeighboor;
					}

					mesh->SetPoint(i,point);

					//write a new Mesh with correct triangles
					if(correctbadtriangleMeshOn){
						SOMesh->SetMesh(mesh);
						MeshWriterType::Pointer writer = MeshWriterType::New();
						writer->SetInput(SOMesh);
						writer->SetFileName(correctFilename);
						writer->Update();
					}
					//know if the mesh has changed
					if(curPoint[0]!=point[0] || curPoint[1]!=point[1] || curPoint[2]!=point[2])
						meshChanged =true;     
				}
			}
			++z;
		}

		//write the correct Mesh with all correct triangles
		if(correctbadtriangleMeshOn){
			SOMesh->SetMesh(mesh);
			MeshWriterType::Pointer writer = MeshWriterType::New();
			writer->SetInput(SOMesh);
			writer->SetFileName(correctFilename);
			writer->Update();
		}

		//write the correct values of the average of triangles in a text file
		if(!correctbadtriangleMeshOn){
			//open outputFile
			std::ofstream outfileNor;
			outfileNor.open (outputFilename ) ;
			outfileNor << "NUMBER_OF_POINTS=" << numPoints << std::endl ;
			outfileNor << "DIMENSION=" << 1 << std::endl ;
			outfileNor << "TYPE=Scalar" << std::endl ;
			for (int i = 0; i < numPoints; i++)
				outfileNor << avgTriangleArea[i] << std::endl;
			//close outputFile
			outfileNor.close();
		}

		delete[] avgTriangleArea;
		delete[] numTriangles;
		delete[] neighboor;

	} else if(valueOn){

		char line[1000];

		//open inputFile
		std::ifstream input, input2;
		input2.open ( inputFilename ); 
		input2.seekg(0,std::ios::beg);
		//calculate number of lines of the inputFilename
		int numberOfLines = 0;
		while ((input2.getline(line, 1000)) != NULL)
			numberOfLines++;
		input2.close () ;
		//numberOfLines-3 corresponds to the number of points (without the three first lines of KWMeshVisu textfile)
		numberOfLines=numberOfLines-3;

		//open inputFile
		input.open ( inputFilename ); 
		input.seekg(0,std::ios::beg);

		//create output (KWMeshVisu textFile)
		std::ofstream output;
		output.open ( outputFilename ) ;
		output << "NUMBER_OF_POINTS= " << numberOfLines << std::endl;
		output << "DIMENSION=1" << std::endl;
		output << "TYPE=Scalar" << std::endl;

		//extract column 5 and write the values
		while ((input.getline(line, 1000)) != NULL)
		{
			float dummy, extractValue, point;
			sscanf(line, " %f %f %f %f %f ", &point, &dummy, &dummy, &dummy, &extractValue);
			//////////uncomment if you want to select with the first column (like only one label)
			// if(point == Index)
			// 	output <<extractValue << std::endl;
		}

		//close input&output
		input.close () ;
		output.close();

	}else if(labelOn) {

		char line[1000],line1[1000],line2[1000];
		char * valuePtr;
		int nPts;

		std::ifstream labelFile;
		labelFile.open ( labelFilename ); 
		labelFile.seekg(0,std::ios::beg);
		bool found=false;
		while ( !found && !labelFile.eof())
		{ labelFile.getline ( line, 1000 ) ;
			if (line[0] != '#' && strstr ( line, "NUMBER_OF_POINTS" )) found = true;
		}
		valuePtr=strchr(line, '=');
		if (!valuePtr) return 0;
		valuePtr++;
		sscanf(valuePtr, " %d ", &nPts);
		labelFile.close();

		std::ofstream output;
		output.open ( outputFilename ) ;

		std::ifstream inputFile;
		inputFile.open ( inputFilename ); 
		inputFile.seekg(0,std::ios::beg);

		std::ifstream labelFile1;
		labelFile1.open ( labelFilename ); 
		labelFile1.seekg(0,std::ios::beg);
		float label[nPts];
		float thick[nPts];
		int i=0;

		found=false;
		while ( !found && !labelFile1.eof())
		{
			labelFile1.getline(line1, 1000);
			if ( strstr ( line1, "TYPE=Scalar" )) found= true;
		}

		while ((labelFile1.getline(line1, 1000)) != NULL)
		{
			float Value;
			sscanf(line1, " %f ", &Value);
			label[i]=Value;
			i++;
		}

		found=false;
		while ( !found && !inputFile.eof())
		{
			inputFile.getline(line2, 1000);
			if ( strstr ( line2, "TYPE=Scalar" )) found= true;
		}

		i=0;
		while ((inputFile.getline(line2, 1000)) != NULL)
		{
			float Value1;
			sscanf(line2, " %f ", &Value1);
			thick[i]=Value1;
			i++;
		}

		///////separeted all labels in 1file
		//     int nu=0;
		//     float label_new[nPts];
		//     int d=0;
		//     float num[nPts];
		//     for (int c = 0; c < nPts; c++) 
		//       {
		// 	for(int a =0; a<29;a++)
		// 	  {
		// 	    if(label[c]==a)
		// 	      {
		// 		num[a]++;
		// 	      }
		// 	  }
		//       }

		//     for(int b =0; b<29;b++)
		//       {
		// 	if(num[b]!=0)
		// 	  {
		// 	    label_new[d]=b;
		// 	    d++;
		// 	  }
		//       }

		//       ///////num labels
		//       //     for(int g =0; g<25;g++)
		//       //       {
		//       // 	std::cout <<g <<" " << label_new[g] <<std::endl;
		//       //       }

		//     for(int a =0; a<25;a++)
		//       {
		// 	nu=0;
		// 	output << label_new[a] <<std::endl;
		// 	for (int c = 0; c < nPts; c++) 
		// 	  {
		// 	    if(label[c]==label_new[a])
		// 	      {
		// 		//std::cout << thick[c] <<",";
		// 		output << thick[c] <<",";
		// 		nu++;
		// 	      }
		// 	  }
		// 	output << "" <<std::endl;
		// 	output << nu <<std::endl;
		// 	output << "" <<std::endl;
		//       }



		////// average all labels :
		//     float thicknew[nPts];
		//     float num[nPts];
		//     for (int c = 0; c < nPts; c++) 
		//       {
		// 	for(int a =0; a<29;a++)
		// 	  {
		// 	    if(label[c]==a)
		// 	      {
		// 		thicknew[a]=thick[c]+thicknew[a];
		// 		num[a]++;
		// 	      }
		// 	  }
		//       }
		//     for(int b =0; b<29;b++)
		//       {
		// 	if(num[b]!=0)

		// 	  {
		// 	    thicknew[b]=thicknew[b]/num[b];
		// 	    output <<thicknew[b] <<std::endl;
		// 	  }
		//       }


		//////mean for one label : label == 11
		output <<"NUMBER_OF_POINTS="<<1<<std::endl;
		output <<"DIMENSION=1 "<<std::endl;
		output <<"TYPE=Scalar "<<std::endl;

		float thicknew=0;
		int num=0;
		for (int c = 0; c < nPts; c++) 
		{
			if(label[c]==11)	
			{
				thicknew=thick[c]+thicknew;
				num++;
				//every values of one label : label == 11
				//output <<thick[c] << std::endl;
			}
		}
		thicknew=thicknew/num;
		output <<thicknew <<std::endl;

		inputFile.close();
		labelFile1.close();
		output.close ();

	}else if(colorOn){
		char line[1000];
		bool found=false;

		std::ifstream input;
		input.open ( inputFilename ); 
		input.seekg(0,std::ios::beg);

		std::ofstream output;
		output.open ( outputFilename );

		while ( !found && !input.eof())
		{
			input.getline(line, 1000);
			output <<line<<std::endl;
			if ( strstr ( line, "LOOKUP_TABLE default" )) found = true;
		}

		float *val = new float [30];
		int i=1;
		int max=0;
		for(int y=1; y<=num[0];y++)
		{
			val[i]=num[y];
			i++;
		}

		val[i]=num[1];
		i++;
		val[i]=num[1];
		i++;
		val[i]=num[1];
		i++;

		for( int z=1;z<=old_num[0];z++)
		{
			val[i]=old_num[z];
			i++;
			max=i;
		}

		while ((input.getline(line, 1000)) != NULL )
		{
			float a,b,c,d,e,f,g,h,i;
			bool Afound=false,Bfound=false,Cfound=false,Dfound=false,Efound=false,Ffound=false,Gfound=false,Hfound=false,Ifound=false;
			sscanf(line, " %f %f %f %f %f %f %f %f %f", &a, &b, &c, &d, &e, &f, &g, &h, &i );
			for(int q=1;q<max;q++){
				if(a==val[q] && Afound==false) {
					a=q;
					Afound=true;}
					if(b==val[q] && Bfound==false) {
						b=q;
						Bfound=true; }
						if(c==val[q] && Cfound==false) {
							c=q;
							Cfound=true; }
							if(d==val[q] && Dfound==false) {
								d=q;
								Dfound=true; }
								if(e==val[q] && Efound==false) {
									e=q;
									Efound=true; }
									if(f==val[q] && Ffound==false) {
										f=q;
										Ffound=true; }
										if(g==val[q] && Gfound==false) {
											g=q;
											Gfound=true; }
											if(h==val[q] && Hfound==false) {
												h=q;
												Hfound=true; }
												if(i==val[q] && Ifound==false) {
													i=q;
													Ifound=true; }
			}

			if(Afound==false)
				a=0;
			if(Bfound==false)
				b=0;
			if(Cfound==false)
				c=0;
			if(Dfound==false)
				d=0;
			if(Efound==false)
				e=0;
			if(Ffound==false)
				f=0;
			if(Gfound==false)
				g=0;
			if(Hfound==false)
				h=0;
			if(Ifound==false)
				i=0;

			output <<a<<" "<<b<<" "<<c<<" "<<d<<" "<<e<<" "<<f<<" "<<g<<" "<<h<<" "<<i<<std::endl;
		}

		input.close();
		output.close ();

	} else if(MaxColorOn){

		char line[1000], line2[1000];
		bool found=false;

		std::ifstream in;
		in.open ( inputFilename ); 
		in.seekg(0,std::ios::beg);
		int numberOfLines = 0;
		//find number of lines
		while ((in.getline(line, 1000)) != NULL)
			numberOfLines++;
		in.close () ;
		//numberOfLines-3 corresponds to the number of points (without the three first lines of KWMeshVisu textfile)
		numberOfLines=numberOfLines-3;

		std::ifstream input;
		input.open ( inputFilename ); 
		input.seekg(0,std::ios::beg);
		//pass the three first lines of KWMeshVisu textfile
		while ( !found && !input.eof())
		{
			input.getline(line, 1000);
			if ( strstr ( line, "TYPE=Scalar" )) found= true;
		}
		//numFiles corresponds to number of files -1 (without the inputfile)
		int numFiles = MaxFiles.size();
		int index=0;
		int val=0;

		std::ifstream attrFile[numFiles];

		double **value;
		value = new double* [numberOfLines];
		for (int c = 0; c <numberOfLines ; c++)
		{    
			value[c]= new double [numFiles +1];
			for (int d = 0; d < numFiles +1 ; d++)
				value[c][d]=0.0;
		}

		std::ofstream output,output2;

		for (int num = 0; num < numFiles; num++) {
			attrFile[num].open (MaxFiles[num].c_str()) ; 
			attrFile[num].seekg(0,std::ios::beg);
			for(val=0;val<2 ;val++)
			{
				attrFile[num].getline(line2, 1000);
			}
		}

		while ((input.getline(line, 1000)) != NULL )
		{
			float point;
			sscanf(line, " %f ", &point);
			value[index][0]=point;

			for (int num = 0; num < numFiles; num++) {
				float point2;
				attrFile[num].getline(line2, 1000);
				sscanf(line2, " %f ", &point2);
				value[index][num+1]=point2;
			}
			++index;
		}

		for (int num = 0; num < numFiles; num++) {
			attrFile[num].close();
		}
		input.close();

		for(int a=0; a <numberOfLines ; a++)
		{
			float max=0.0,maxI=0.0;

			for (int b = 0; b < numFiles +1 ; b++)
			{
				if( value[a][b]>max)
					max=value[a][b];
			}

			maxI=max-max*5/100;

			for (int e = 0; e < numFiles +1 ; e++)
			{
				if( value[a][e]<=maxI)
					value[a][e]=0.0;
			}
		}

		//write in the first file -> inputfile
		output.open(inputFilename);
		output << "NUMBER_OF_POINTS="<<numberOfLines<< std::endl;
		output << "DIMENSION=1" << std::endl;
		output << "TYPE=Scalar" << std::endl;
		for(int z=0;z<numberOfLines;z++)
			output << value[z][0]<< std::endl;
		output.close();

		//write in the other files
		for(int x = 0; x < numFiles +1 ; x++)
		{
			output2.open (MaxFiles[x].c_str()) ; 
			output2 << "NUMBER_OF_POINTS="<<numberOfLines<< std::endl;
			output2 << "DIMENSION=1" << std::endl;
			output2 << "TYPE=Scalar" << std::endl;
			for(int y=0; y <numberOfLines ; y++)
			{
				output2 << value[y][x+1]<< std::endl;
			}
			output2.close();
		}

		delete[] value;

	}  else if(distAbsOn){

		char line[1000], line2[1000];
		bool found=false;

		std::ifstream in;
		in.open ( inputFilename ); 
		in.seekg(0,std::ios::beg);
		int numberOfLines = 0;
		//find number of lines
		while ((in.getline(line, 1000)) != NULL)
			numberOfLines++;
		in.close () ;
		//numberOfLines-3 corresponds to the number of points (without the three first lines of KWMeshVisu textfile)
		numberOfLines=numberOfLines-3;

		std::ifstream input;
		input.open ( inputFilename ); 
		input.seekg(0,std::ios::beg);
		//pass the three first lines of KWMeshVisu textfile

		while ( !found && !input.eof())
		{
			input.getline(line, 1000);
			if ( strstr ( line, "TYPE=Scalar" )) found= true;
		}

		int index=0;
		int val=0;
		int numFiles=distFiles.size();
		double **value;
		value = new double* [numberOfLines];
		for (int c = 0; c <numberOfLines ; c++)
		{    
			value[c]= new double [numFiles +1];
			for (int d = 0; d < numFiles +1 ; d++)
				value[c][d]=0.0;
		}
		double **value_new;
		value_new = new double* [numberOfLines];
		for (int u = 0; u <numberOfLines ; u++)
		{    
			value_new[u]= new double [numFiles +1];
			for (int k = 0; k < numFiles +1 ; k++)
				value_new[u][k]=0.0;
		}
		std::ifstream attrFile[numFiles];
		for (int num = 0; num < numFiles; num++) {
			attrFile[num].open (distFiles[num].c_str()) ; 
			attrFile[num].seekg(0,std::ios::beg);
			for(val=0;val<2 ;val++)
			{
				attrFile[num].getline(line2, 1000);
			}
		}

		while ((input.getline(line, 1000)) != NULL )
		{
			float point;
			sscanf(line, " %f ", &point);
			value[index][0]=point;
			for (int num = 0; num < numFiles; num++) {
				float point2;
				attrFile[num].getline(line2, 1000);
				sscanf(line2, " %f ", &point2);
				value[index][num+1]=point2;
			}
			++index;
		}

		for (int num = 0; num < numFiles; num++) {
			attrFile[num].close();
		}
		input.close();

		for(int a=0; a <numberOfLines ; a++)  {

			for (int b = 1; b < numFiles; b++)
				value_new[a][b]=value[a][b+1]-value[a][b-1];
		}

		//writte output file
		std::ofstream output;
		for(int x = 0; x < numFiles  ; x++)
		{	
			output.open (resultFiles[x].c_str()) ; 
			output << "NUMBER_OF_POINTS="<<numberOfLines<< std::endl;
			output << "DIMENSION=1" << std::endl;
			output << "TYPE=Scalar" << std::endl;

			for(int p=0; p < numberOfLines ; p++)
				output << value_new[p][x+1]<< std::endl;

			output.close();
		}

		delete[] value;
		delete[] value_new;

	} else if(distRelOn){

		char line[1000], line2[1000];
		bool found=false;

		std::ifstream in;
		in.open ( inputFilename ); 
		in.seekg(0,std::ios::beg);
		int numberOfLines = 0;
		//find number of lines
		while ((in.getline(line, 1000)) != NULL)
			numberOfLines++;
		in.close () ;
		//numberOfLines-3 corresponds to the number of points (without the three first lines of KWMeshVisu textfile)
		numberOfLines=numberOfLines-3;

		std::ifstream input;
		input.open ( inputFilename ); 
		input.seekg(0,std::ios::beg);
		//pass the three first lines of KWMeshVisu textfile

		while ( !found && !input.eof())
		{
			input.getline(line, 1000);
			if ( strstr ( line, "TYPE=Scalar" )) found= true;
		}

		int index=0;
		int val=0;
		int numFiles=distReFiles.size();
		double **value;
		value = new double* [numberOfLines];
		for (int c = 0; c <numberOfLines ; c++)
		{    
			value[c]= new double [numFiles +1];
			for (int d = 0; d < numFiles +1 ; d++)
				value[c][d]=0.0;
		}
		double **value_new;
		value_new = new double* [numberOfLines];
		for (int u = 0; u <numberOfLines ; u++)
		{    
			value_new[u]= new double [numFiles +1];
			for (int k = 0; k < numFiles +1 ; k++)
				value_new[u][k]=0.0;
		}
		std::ifstream attrFile[numFiles];
		for (int num = 0; num < numFiles; num++) {
			attrFile[num].open (distReFiles[num].c_str()) ; 
			attrFile[num].seekg(0,std::ios::beg);
			for(val=0;val<2 ;val++)
			{
				attrFile[num].getline(line2, 1000);
			}
		}

		while ((input.getline(line, 1000)) != NULL )
		{
			float point;
			sscanf(line, " %f ", &point);
			value[index][0]=point;
			for (int num = 0; num < numFiles; num++) {
				float point2;
				attrFile[num].getline(line2, 1000);
				sscanf(line2, " %f ", &point2);
				value[index][num+1]=point2;
			}
			++index;
		}

		for (int num = 0; num < numFiles; num++) {
			attrFile[num].close();
		}
		input.close();

		for(int a=0; a <numberOfLines ; a++)  {

			for (int b = 1; b < numFiles; b++)
				value_new[a][b]=value[a][b+1]-value[a][b-1];

			float max=0.0;
			float min=10.0;
			for (int h = 1; h < numFiles  ; h++)
			{
				if( value_new[a][h]>max){
					max=value_new[a][h];
				}
				if(value[a][h]<min){
					min=value_new[a][h];
				}
			}

			for (int g = 1; g < numFiles  ; g++)
			{
				value_new[a][g]=(value_new[a][g]-min)/(max-min);
			}
		}

		//writte output file
		std::ofstream output;
		for(int x = 0; x < numFiles  ; x++)
		{	
			output.open (restFiles[x].c_str()) ; 
			output << "NUMBER_OF_POINTS="<<numberOfLines<< std::endl;
			output << "DIMENSION=1" << std::endl;
			output << "TYPE=Scalar" << std::endl;

			for(int p=0; p < numberOfLines ; p++)
				output << value_new[p][x+1]<< std::endl;

			output.close();
		}

		delete[] value;
		delete[] value_new;

	} else if(firstOn){

		char line[1000],line2[1000];
		int numfile = firstfiles.size();

		std::ifstream input;
		input.open ( inputFilename ); 
		input.seekg(0,std::ios::beg);

		std::ofstream output;
		output.open ( outputFilename ) ;
		output << "NUMBER_OF_POINTS=" <<numfile+2<< std::endl;
		output << "DIMENSION=1" << std::endl;
		output << "TYPE=Scalar" << std::endl;
		const int numEntries = 3;
		int counter = 0;

		while ( counter < numEntries && !input.eof())
		{ input.getline ( line, 1000 ) ;
			if ((line[0] != '#')) counter++;
		}

		while ((input.getline(line, 1000)) != NULL)
		{
			float point;
			sscanf(line, " %f ", &point);
			output <<point << ",";
		}
		input.close () ;

		for (int index = 0; index < numfile; index++) {

			std::ifstream input2;
			counter=0;
			input2.open (firstfiles[index].c_str()) ; 
			input2.seekg(0,std::ios::beg);

			while ( counter < numEntries && !input2.eof())
			{ input2.getline ( line2, 1000 ) ;
				if ((line2[0] != '#')) counter++;
			}

			while ((input2.getline(line2, 1000)) != NULL)
			{
				float point2;
				sscanf(line2, " %f ", &point2);
				output <<point2 <<",";
			}
			input2.close();
		}
		output <<" "<<std::endl;
		output.close();

	}else if(subKWMOn){

		char line1[1000],line2[1000];
		float value1=0.0;
		float value2=0.0;
		float diff=0.0;
		bool found;
		char * valuePtr;
		int nPts1, nPts2, nDim;
		char typeString[1001], line[1001];

		//open input
		std::ifstream input1,input2 ;
		input1.open ( inputFilename ) ; 
		input1.seekg(0,std::ios::beg);
		found = false ;
		//find number of points
		while ( !found && !input1.eof())
		{ input1.getline ( line, 1000 ) ;
			if (line[0] != '#' && strstr ( line, "NUMBER_OF_POINTS" )) found = true;
		}
		valuePtr=strchr(line, '=');
		if (!valuePtr) return 0;
		valuePtr++;
		sscanf(valuePtr, " %d ", &nPts1);

		input1.seekg(0,std::ios::beg);
		found = false ;
		while ( !found && !input1.eof())
		{ input1.getline ( line, 1000 ) ;
			if (line[0] != '#' && strstr ( line, "DIMENSION" )) found = true;
		}
		valuePtr=strchr(line, '=');
		if (!valuePtr) return 0;
		valuePtr++;
		sscanf(valuePtr, " %d ", &nDim);

		input1.seekg(0,std::ios::beg);
		found = false ;
		while ( !found && !input1.eof())
		{ input1.getline ( line, 1000 ) ;
			if (line[0] != '#' && strstr ( line, "TYPE" )) found = true;
		}
		valuePtr=strchr(line, '=');
		if (!valuePtr) return 0;
		valuePtr++;
		sscanf(valuePtr, " %s ", typeString);
		assert ( nDim == 1 ) ;
		assert ( strcmp ( typeString, "Scalar" ) == 0 ) ;

		//open subKWMFile
		input2.open ( subKWMFilename ) ; 
		input2.seekg(0,std::ios::beg);
		found = false ;
		//find the number of points
		while ( !found && !input2.eof())
		{ input2.getline ( line, 1000 ) ;
			if (line[0] != '#' && strstr ( line, "NUMBER_OF_POINTS" )) found = true;
		}
		valuePtr=strchr(line, '=');
		if (!valuePtr) return 0;
		valuePtr++;
		sscanf(valuePtr, " %d ", &nPts2);

		input2.seekg(0,std::ios::beg);
		found = false ;
		while ( !found && !input2.eof())
		{ input2.getline ( line, 1000 ) ;
			if (line[0] != '#' && strstr ( line, "DIMENSION" )) found = true;
		}
		valuePtr=strchr(line, '=');
		if (!valuePtr) return 0;
		valuePtr++;
		sscanf(valuePtr, " %d ", &nDim);

		input2.seekg(0,std::ios::beg);
		found = false ;
		while ( !found && !input2.eof())
		{ input2.getline ( line, 1000 ) ;
			if (line[0] != '#' && strstr ( line, "TYPE" )) found = true;
		}
		valuePtr=strchr(line, '=');
		if (!valuePtr) return 0;
		valuePtr++;
		sscanf(valuePtr, " %s ", typeString);
		//number of points of the first file = number of points of the second file
		assert ( nPts2 == nPts1 ) ;
		assert ( nDim == 1 ) ;
		assert ( strcmp ( typeString, "Scalar" ) == 0 ) ;

		//create ouput file
		std::ofstream output;
		output.open ( outputFilename ) ;
		output << "NUMBER_OF_POINTS=" << nPts1 << std::endl ;
		output << "DIMENSION=" << 1 << std::endl ;
		output << "TYPE=Scalar" << std::endl ;
		//find the values of the two files, make the difference and write them into the ouputFile
		while(((input1.getline(line1, 1000)) != NULL) && (input2.getline(line2, 1000) != NULL)) {	
			value1=atof(line1);
			value2=atof(line2);
			diff=value1-value2;
			output  << diff << std::endl;
		}

		//close inputs and ouput
		input1.close () ;
		input2.close();
		output.close();

	}else if(avgGaussKWMOn){

		float sum =0;
		float gaussKWM =0;
		int nPts ;     
		int nPtsFile, nDim ;
		bool found;
		char * valuePtr;
		char typeString[1000], line[1000];
		double Gaussvalue;

		//open file(input)
		std::ifstream input ;
		input.open ( inputFilename ) ; 
		input.seekg(0,std::ios::beg);
		found = false ;
		//find the number of points
		while ( !found && !input.eof())
		{ input.getline ( line, 1000 ) ;
			if (line[0] != '#' && strstr ( line, "NUMBER_OF_POINTS" )) found = true;
		}
		valuePtr=strchr(line, '=');
		if (!valuePtr) return 0;
		valuePtr++;
		sscanf(valuePtr, " %d ", &nPtsFile);
		nPts=atoi(valuePtr);
		input.seekg(0,std::ios::beg);
		found = false ;
		while ( !found && !input.eof())
		{ input.getline ( line, 1000 ) ;
			if (line[0] != '#' && strstr ( line, "DIMENSION" )) found = true;
		}
		valuePtr=strchr(line, '=');
		if (!valuePtr) return 0;
		valuePtr++;
		sscanf(valuePtr, " %d ", &nDim);
		input.seekg(0,std::ios::beg);
		found = false ;
		while ( !found && !input.eof())
		{ input.getline ( line, 1000 ) ;
			if (line[0] != '#' && strstr ( line, "TYPE" )) found = true;
		}
		valuePtr=strchr(line, '=');
		if (!valuePtr) return 0;
		valuePtr++;
		sscanf(valuePtr, " %s ", typeString);
		assert ( nPtsFile == nPts ) ;
		assert ( nDim == 1 ) ;
		assert ( strcmp ( typeString, "Scalar" ) == 0 ) ;


		float *avgGaussKWM= new float [nPts];
		float val;

		input.seekg(0,std::ios::beg);
		const int numEntries = 3;
		int counter = 0;
		while ( counter < numEntries && !input.eof())
		{ input.getline ( line, 1000 ) ;
			if ((line[0] != '#')) counter++;
		}
		//calculate gauss coefficient
		gaussKWM = (1/(GaussKWMParams[1]*sqrt(2*M_PI))*exp(-(AgeKWM[0]-GaussKWMParams[0])*(AgeKWM[0]-GaussKWMParams[0])/(GaussKWMParams[1]*GaussKWMParams[1]*2)))/0.0797885;
		sum=gaussKWM;

		//find values and calculate the gauss values
		for (int i = 0 ; i < nPts ; i++ )
		{
			input.getline(line, 1000);
			val = atof(line);
			avgGaussKWM[i]=gaussKWM*val;
		}
		//close file
		input.close () ;

		int numFiles = avgGaussKWMFiles.size();
		float newval;

		//for all files
		for (int index = 0; index < numFiles; index++) {

			//open file
			input.open (avgGaussKWMFiles[index].c_str()) ; 
			input.seekg(0,std::ios::beg);
			found = false ;
			//find the number of points
			while ( !found && !input.eof())
			{ input.getline ( line, 1000 ) ;
				if (line[0] != '#' && strstr ( line, "NUMBER_OF_POINTS" )) found = true;
			}
			valuePtr=strchr(line, '=');
			if (!valuePtr) return 0;
			valuePtr++;
			sscanf(valuePtr, " %d ", &nPtsFile);
			nPts=atoi(valuePtr);
			input.seekg(0,std::ios::beg);
			found = false ;
			while ( !found && !input.eof())
			{ input.getline ( line, 1000 ) ;
				if (line[0] != '#' && strstr ( line, "DIMENSION" )) found = true;
			}
			valuePtr=strchr(line, '=');
			if (!valuePtr) return 0;
			valuePtr++;
			sscanf(valuePtr, " %d ", &nDim);
			input.seekg(0,std::ios::beg);
			found = false ;
			while ( !found && !input.eof())
			{ input.getline ( line, 1000 ) ;
				if (line[0] != '#' && strstr ( line, "TYPE" )) found = true;
			}
			valuePtr=strchr(line, '=');
			if (!valuePtr) return 0;
			valuePtr++;
			sscanf(valuePtr, " %s ", typeString);
			assert ( nPtsFile == nPts ) ;
			assert ( nDim == 1 ) ;
			assert ( strcmp ( typeString, "Scalar" ) == 0 ) ;
			input.seekg(0,std::ios::beg);
			const int numEntries = 3;
			int counter = 0;
			while ( counter < numEntries && !input.eof())
			{ input.getline ( line, 1000 ) ;
				if ((line[0] != '#')) counter++;
			}

			//calculate gauss coefficient
			gaussKWM = (1/(GaussKWMParams[1]*sqrt(2*M_PI))*exp(-(AgeKWM[index+1]-GaussKWMParams[0])*(AgeKWM[index+1]-GaussKWMParams[0])/(GaussKWMParams[1]*GaussKWMParams[1]*2)))/0.0797885;
			sum=sum+gaussKWM;

			//find values and calculate the gauss values
			for (int i = 0 ; i < nPts ; i++ )
			{
				input.getline(line, 1000);
				newval = atof(line);
				Gaussvalue = gaussKWM*newval;
				avgGaussKWM[i]=avgGaussKWM[i]+Gaussvalue;
			}

			//close file
			input.close () ;

		}
		//calculate average (divide the values by sum)
		for(int i =0; i < nPts;i++)
			avgGaussKWM[i]=avgGaussKWM[i]/sum;
		//write KWMeshVisu textFile with new values  
		std::ofstream output;
		output.open ( outputFilename ) ;
		output << "NUMBER_OF_POINTS=" << nPts << std::endl ;
		output << "DIMENSION=" << 1 << std::endl ;
		output << "TYPE=Scalar" << std::endl ;
		for (int i = 0; i < nPts; i++) 
			output  << avgGaussKWM[i] << std::endl;
		//close output
		output.close();

		delete[] avgGaussKWM;

	} else if(avgGaussMeshOn) {

		float sum =0;
		float gauss=0;

		MeshConverterType * converter = new MeshConverterType () ;
		//read input
		MeshSOType::Pointer SOMesh = converter->ReadMeta (inputFilename) ;
		MeshType::Pointer surfaceMesh = SOMesh->GetMesh() ;
		MeshType::PointsContainerPointer avgPoints = surfaceMesh->GetPoints();
		int numberOfPoints = avgPoints->Size();
		{
			MeshType::PointsContainerPointer pointsTmp = MeshType::PointsContainer::New();
			//calculate the gauss coefficient with mean,std,age...
			gauss = (1/(GaussParams[1]*sqrt(2*M_PI))*exp(-(Age[0]-GaussParams[0])*(Age[0]-GaussParams[0])/(GaussParams[1]*GaussParams[1]*2)))/0.0797885;
			//calculate the average values
			for (unsigned int pointID = 0; pointID < avgPoints->Size();pointID++) {
				PointType vert;
				PointType curPointAvg =  avgPoints->GetElement(pointID);
				for (unsigned int dim = 0; dim < 3; dim++) 
					vert[dim] = gauss*curPointAvg[dim];
				pointsTmp->InsertElement(pointID, vert);
			}
			sum=gauss;
			avgPoints = pointsTmp;
		}

		//for all meshes
		int numMeshes = AvgGaussMeshFiles.size();
		for (int index = 0; index < numMeshes; index++) {
			try
			{
				MeshSOType::Pointer inputMeshSO = converter->ReadMeta (AvgGaussMeshFiles[index].c_str()) ;
				MeshType::Pointer inputMesh = inputMeshSO->GetMesh() ;
				MeshType::PointsContainerPointer points = inputMesh->GetPoints();
				MeshType::PointsContainerPointer pointsTmp = MeshType::PointsContainer::New();
				//calculate the gauss coefficient with mean,std,age...
				gauss = (1/(GaussParams[1]*sqrt(2*M_PI))*exp(-(Age[index+1]-GaussParams[0])*(Age[index+1]-GaussParams[0])/(GaussParams[1]*GaussParams[1]*2)))/0.0797885;
				sum=sum+gauss;

				if (points->Size() != (unsigned int) numberOfPoints) {
					std::cout << "Number of Points do not agree, expected : " << numberOfPoints << " , got " << points->Size() << " File: " << AvgGaussMeshFiles[index] <<  std::endl;
					exit (-1);
				}
				//calculate the average values
				for (unsigned int pointID = 0; pointID < points->Size();pointID++) {
					PointType curPoint =  points->GetElement(pointID);
					PointType vert;
					PointType curPointAvg =  avgPoints->GetElement(pointID);
					for (unsigned int dim = 0; dim < 3; dim++) 
						vert[dim] = gauss*curPoint[dim] + curPointAvg[dim];
					pointsTmp->InsertElement(pointID, vert);
				}
				avgPoints = pointsTmp;
			}
			catch(itk::ExceptionObject ex)
			{
				std::cout<< "Error reading meshfile:  "<< AvgGaussMeshFiles[index] << std::endl << "ITK error: " << ex.GetDescription()<< std::endl;
				exit(-3);
			}
		}
		//calculate the average (values divide by the sum)	    
		MeshType::PointsContainerPointer pointsTmp = MeshType::PointsContainer::New();
		for (unsigned int pointID = 0; pointID < avgPoints->Size();pointID++) {
			PointType curPoint =  avgPoints->GetElement(pointID);
			PointType vert;
			for (unsigned int dim = 0; dim < 3; dim++) 
				vert[dim] = curPoint[dim]/sum;
			pointsTmp->InsertElement(pointID, vert);
		}
		avgPoints = pointsTmp;

		//create new Mesh with the new values
		surfaceMesh->SetPoints(avgPoints); 
		SOMesh->SetMesh(surfaceMesh);
		MeshWriterType::Pointer writer = MeshWriterType::New();
		writer->SetInput(SOMesh);
		writer->SetFileName(outputFilename);
		writer->Update();
		delete []Age;

	} else if(averageOn) {

		int NbAveFile = AveFiles.size();
		if(NbAveFile == 0)
		{
			std::cout << "-ave option has to have at least one file as input" << std::endl;
			return 0;
		}

		// read in the input files
		MeshConverterType * converter = new MeshConverterType () ;
		if (debug)   std::cout << "Reading input mesh " << inputFilename << std::endl;
		MeshSOType::Pointer inputMeshSO = converter->ReadMeta (inputFilename);
		MeshType::Pointer inputMesh = inputMeshSO->GetMesh();
		int NbOfPoint = inputMesh->GetNumberOfPoints();

		//Convert the meta data into the VTK Poly Data format
		if (debug)   std::cout << "Computing normals" << std::endl;
		TriangleMeshConverterType * Triangleconverter = new TriangleMeshConverterType () ;
		TriangleMeshSOType::Pointer inputTriangleMeshSO = Triangleconverter->ReadMeta (inputFilename);
		TriangleMeshType::Pointer inputTriangleMesh = inputTriangleMeshSO->GetMesh();

		itkMeshTovtkPolyData *convertMeshToVTK = new itkMeshTovtkPolyData();
		convertMeshToVTK->SetInput(inputTriangleMesh);

		vtkPolyData * vtkMesh = convertMeshToVTK->GetOutput();
		vtkPolyDataNormals *MeshNormals = vtkPolyDataNormals::New();

		MeshNormals->SetComputePointNormals(1);
		MeshNormals->SetComputeCellNormals(0);
		MeshNormals->SetSplitting(0);
		MeshNormals->SetInput(vtkMesh);
		MeshNormals->Update();
		vtkPolyData * vtkMeshNormals = MeshNormals->GetOutput();
		vtkMeshNormals->Update();

		vtkPointData * NormalPoints = vtkMeshNormals->GetPointData();
		vtkDataArray * ArrayNormal = NormalPoints->GetNormals();

		/*    double *n = new double[3];
		      for(int NbNorm = 0 ; NbNorm < vtkMeshNormals->GetNumberOfPoints() ; NbNorm++)
		      {
		      n = ArrayNormal->GetTuple3(NbNorm);
		      std::cout << "Coord[" << NbNorm << "]= " << vtkMesh->GetPoint(NbNorm)[0] << " | " << vtkMesh->GetPoint(NbNorm)[1] << " | " << vtkMesh->GetPoint(NbNorm)[2] << std::endl;
		      std::cout << "Normal[" << NbNorm << "]= " << n[0]  << " | " << n[1] << " | " << n[2] << std::endl;
		      }
		      */
		char output[512];
		float *AverageX = new float[NbOfPoint];
		float *AverageY = new float[NbOfPoint];
		float *AverageZ = new float[NbOfPoint];
		if(debug) std::cout << "Vector initialization" << std::endl;
		for(int i = 0 ; i < NbOfPoint ; i++)
			AverageX[i] = AverageY[i] = AverageZ[i] = 0;

		for(int FileId = 0 ; FileId < (int)(AveFiles.size()) ; FileId++)
		{
			if(debug) std::cout << "Reading file number " << FileId << std::endl;
			//Open the file
			std::ifstream avefile ;
			avefile.open(AveFiles[FileId].c_str());
			//Get the number of point in the file
			avefile.getline(output,32);
			std::string Line = output;	
			int loc = Line.find("=",1);
			std::string NbPoint = Line.erase(0,loc+1);
			int Val = atoi(NbPoint.c_str());
			if(Val != NbOfPoint)
			{
				std::cout << "The file must have the same number of point as the meta file" << std::endl;
				return 0;
			}else {

				//Read the two next lines
				avefile.getline(output,64);
				avefile.getline(output,64);
				//Then get the numbers
				for(int Line = 0 ; Line < NbOfPoint ; Line++)
				{
					float X,Y,Z;
					avefile.getline(output,64);	  
					sscanf(output,"%f %f %f",&X,&Y,&Z);
					AverageX[Line] += X;
					AverageY[Line] += Y;
					AverageZ[Line] += Z;
				}
			}
		}

		if(debug) std::cout << "Writing outputs" << std::endl;
		//Computing the average of each vector and writing it down in the file along with the norm values
		std::ofstream outfileVec ;
		outfileVec.open ( outputFilename ) ;
		// print the header
		outfileVec << "NUMBER_OF_POINTS = " << NbOfPoint << std::endl ;
		outfileVec << "DIMENSION = " << 3 << std::endl ;
		outfileVec << "TYPE = Vector" << std::endl ;

		std::ofstream outfileNor;
		outfileNor.open ( outputFilename2.c_str() ) ;  
		// print the header
		outfileNor << "NUMBER_OF_POINTS=" << NbOfPoint << std::endl ;
		outfileNor << "DIMENSION=" << 1 << std::endl ;
		outfileNor << "TYPE=Scalar" << std::endl ;

		double *Normal = new double[3];

		for(int Line = 0 ; Line < NbOfPoint ; Line++)
		{
			float MeanX = AverageX[Line] / NbAveFile;
			float MeanY = AverageY[Line] / NbAveFile;
			float MeanZ = AverageZ[Line] / NbAveFile;
			//if the mean vector is projected on the normal at each point
			if(normaveOn){
				Normal = ArrayNormal->GetTuple3(Line);
				float DotProd = (float)Normal[0] * MeanX + (float)Normal[1] * MeanY + (float)Normal[2] * MeanZ;
				float MeanXProj = (float)Normal[0] * DotProd;
				float MeanYProj = (float)Normal[1] * DotProd;
				float MeanZProj = (float)Normal[2] * DotProd;
				outfileVec << MeanXProj << " " << MeanYProj << " " << MeanZProj << std::endl ;
				float Norme = sqrt(MeanXProj * MeanXProj + MeanYProj * MeanYProj + MeanZProj * MeanZProj);
				outfileNor << Norme << std::endl;
			}else{       
				outfileVec << MeanX << " " << MeanY << " " << MeanZ << std::endl ;
				float Norme = sqrt(MeanX * MeanX + MeanY * MeanY + MeanZ * MeanZ);
				outfileNor << Norme << std::endl;
			}
		}

		outfileVec.close();
		outfileNor.close();

	} else if(invvectOn) {

		std::ifstream vectorFile ;
		vectorFile.open(VectFile);
		std::ofstream outfileVec ;
		outfileVec.open ( outputFilename ) ;
		//Get the number of point in the file
		char output[64];
		vectorFile.getline(output,32);
		outfileVec << output << std::endl;
		std::string Line = output;	
		int loc = Line.find("=",1);
		std::string NbPoint = Line.erase(0,loc+1);
		int Val = atoi(NbPoint.c_str());
		vectorFile.getline(output,32);
		outfileVec << output << std::endl;
		vectorFile.getline(output,32);
		outfileVec << output << std::endl;

		float X,Y,Z;
		for(int Line = 0 ; Line < Val ; Line++)
		{
			vectorFile.getline(output,64);	  
			sscanf(output,"%f %f %f",&X,&Y,&Z);
			X = -1 * X;
			Y = -1 * Y;
			Z = -1 * Z;
			outfileVec << X << " " << Y << " " << Z << std::endl;
		}
		outfileVec.close();

	} else if (magdirOn) {

		//Read the meta file
		MeshConverterType * converter = new MeshConverterType () ;
		if (debug)   std::cout << "Reading input mesh " << inputFilename << std::endl;
		MeshSOType::Pointer inputMeshSO = converter->ReadMeta (inputFilename);
		MeshType::Pointer inputMesh = inputMeshSO->GetMesh();
		int NbOfPoint = inputMesh->GetNumberOfPoints();

		//Convert the meta data into the VTK Poly Data format
		if (debug)   std::cout << "Computing normals" << std::endl;
		TriangleMeshConverterType * Triangleconverter = new TriangleMeshConverterType () ;
		TriangleMeshSOType::Pointer inputTriangleMeshSO = Triangleconverter->ReadMeta (inputFilename);
		TriangleMeshType::Pointer inputTriangleMesh = inputTriangleMeshSO->GetMesh();

		itkMeshTovtkPolyData *convertMeshToVTK = new itkMeshTovtkPolyData();
		convertMeshToVTK->SetInput(inputTriangleMesh);

		//compute the normal
		vtkPolyData * vtkMesh = convertMeshToVTK->GetOutput();
		vtkPolyDataNormals *MeshNormals = vtkPolyDataNormals::New();

		MeshNormals->SetComputePointNormals(1);
		MeshNormals->SetComputeCellNormals(0);
		MeshNormals->SetSplitting(0);
		MeshNormals->SetInput(vtkMesh);
		MeshNormals->Update();
		vtkPolyData * vtkMeshNormals = MeshNormals->GetOutput();
		vtkMeshNormals->Update();

		vtkPointData * NormalPoints = vtkMeshNormals->GetPointData();
		vtkDataArray * ArrayNormal = NormalPoints->GetNormals();
		double * MagList = new double[NbOfPoint];

		//Read the vector field file
		char out[64];
		if(debug) std::cout << "Reading vector field " << VectFile << std::endl;
		//Open the file
		std::ifstream avefile ;
		avefile.open(VectFile);
		//Get the number of point in the file
		avefile.getline(out,32);
		std::string Line = out;	
		int loc = Line.find("=",1);
		std::string NbPoint = Line.erase(0,loc+1);
		int Val = atoi(NbPoint.c_str());
		if(Val != NbOfPoint)
		{
			std::cout << "The file must have the same number of point as the meta file" << std::endl;
			return 0;
		}else {
			double *Normal = new double[3];

			//Read the two next lines
			avefile.getline(out,64);
			avefile.getline(out,64);
			//Then get the numbers
			for(int Line = 0 ; Line < NbOfPoint ; Line++)
			{
				float X,Y,Z;
				avefile.getline(out,64);	  
				sscanf(out,"%f %f %f",&X,&Y,&Z);
				double magnitude = sqrt((double)X*(double)X + (double)Y*(double)Y + (double)Z*(double)Z);
				MagList[Line] = magnitude;
				Normal = ArrayNormal->GetTuple3(Line);
				double DotProduct = X*Normal[0] + Y*Normal[1] + Z*Normal[2];
				if(DotProduct < 0) MagList[Line] = -1 * MagList[Line];
			}
		}
		if(debug) std::cout << "Write the magnitude file" << std::endl;
		std::ofstream outfileNor;
		outfileNor.open ( outputFilename ) ;  
		// print the header
		outfileNor << "NUMBER_OF_POINTS=" << NbOfPoint << std::endl ;
		outfileNor << "DIMENSION=" << 1 << std::endl ;
		outfileNor << "TYPE=Scalar" << std::endl ;
		for(int OutLine = 0 ; OutLine < NbOfPoint ; OutLine++)
			outfileNor << MagList[OutLine] << std::endl;
		outfileNor.close();

	} else if(applyvecOn){

		MeshConverterType * converter = new MeshConverterType () ;  
		if (debug)   std::cout << "Reading input mesh " << inputFilename << std::endl;
		MeshSOType::Pointer inputMeshSO = converter->ReadMeta (inputFilename) ;
		MeshType::Pointer inputMesh = inputMeshSO->GetMesh() ;
		TriangleMeshConverterType * Triangleconverter = new TriangleMeshConverterType () ;
		TriangleMeshSOType::Pointer outputTriangleMeshSO = Triangleconverter->ReadMeta (inputFilename);
		TriangleMeshType::Pointer outputTriangleMesh = outputTriangleMeshSO->GetMesh();

		char output[64];
		std::ifstream vectorFile ;
		vectorFile.open(VectFile);
		if (debug)   std::cout << "Reading the vector file " << VectFile << std::endl;
		vectorFile.getline(output,64);
		std::string Line = output;	
		int loc = Line.find("=",1);
		std::string NbPoint = Line.erase(0,loc+1);
		unsigned int Val = atoi(NbPoint.c_str());

		// make sure the two meshes have the same number of verts
		if (inputMesh->GetNumberOfPoints() != Val) {
			std::cout << "The mesh and the vector field has to have the same number of point" << std::endl;
			std::cout << "The mesh has " << inputMesh->GetNumberOfPoints() << " points and the vector field has " << Val << " points." << std::endl;
			exit(-1);
		}
		//Skipping the second and third lines of the file
		vectorFile.getline(output,32);
		vectorFile.getline(output,32);

		float X,Y,Z;
		std::vector<float> Xval;
		std::vector<float> Yval;
		std::vector<float> Zval;

		for(unsigned int Line = 0 ; Line < Val ; Line++)
		{
			vectorFile.getline(output,64);	  
			sscanf(output,"%f %f %f",&X,&Y,&Z);
			Xval.push_back(X);
			Yval.push_back(Y);
			Zval.push_back(Z);
		}

		PointType *point1 = new PointType ;
		PointTriangleType *Tripoint = new PointTriangleType;

		for ( unsigned int i = 0 ; i < inputMesh->GetNumberOfPoints()  ; i++ )
		{
			if ( inputMesh->GetPoint ( i, point1 ) && &Xval[i] != NULL )
			{	    
				if(debug) std::cout << "Original Point[" << i << "]= " << point1->GetElement(0) << " | " << point1->GetElement(1) << " | " << point1->GetElement(2) << std::endl;
				Tripoint->SetElement(0,point1->GetElement(0) + Xval[i]);
				Tripoint->SetElement(1,point1->GetElement(1) + Yval[i]);
				Tripoint->SetElement(2,point1->GetElement(2) + Zval[i]);
				if(debug) std::cout << "New Point[" << i << "]= " << Tripoint->GetElement(0) << " | " << Tripoint->GetElement(1) << " | " << Tripoint->GetElement(2) << std::endl;
				outputTriangleMesh->SetPoint(i,*Tripoint);
			}
			else
				return 0 ;
		}

		//Convert the outputTriangleMesh to a vtk poly data.
		itkMeshTovtkPolyData *convertMeshToVTK = new itkMeshTovtkPolyData();
		convertMeshToVTK->SetInput(outputTriangleMesh);
		vtkPolyData * vtkPolyData = convertMeshToVTK->GetOutput();

		//Convert the vtk poly data to itk mesh data structure
		vtkPolyDataToitkMesh *vtkItkConverter = new vtkPolyDataToitkMesh () ;
		vtkItkConverter->SetInput(vtkPolyData);

		// Convert the itk mesh data in Spatial Object mesh
		// write out the itk spatial object meta mesh file
		TriangleMeshSOType::Pointer meshSO = TriangleMeshSOType::New();
		meshSO->SetMesh(vtkItkConverter->GetOutput());
		Triangleconverter->WriteMeta(meshSO,outputFilename);

	} else if(MC2OriginOn) //cchou
	{
		MeshConverterType * converter = new MeshConverterType () ;
		MeshSOType::Pointer SOMesh = converter->ReadMeta(inputFilename);
		MeshType::Pointer surfaceMesh = SOMesh->GetMesh () ;
		MeshType::PointsContainerPointer points = surfaceMesh->GetPoints();

		//Sum up Original Points	
		double sum[3];
		for(unsigned int pointID = 0; pointID < points->Size(); pointID++)
		{
			PointType curPoint = points->GetElement(pointID);
			for(unsigned int dim = 0; dim < 3; dim++)
			{
				sum[dim] += curPoint[dim];
			}	
		}	
		//Calculate the Center of Mass
		double MC[3];
		for(unsigned int dim=0; dim < 3; dim++)
		{
			MC[dim] = (double) sum[dim] / (points->Size()+1) ;
		}
		//Create a New Point Set "newpts" with the Shifted Values
		MeshType::PointsContainerPointer newpts = MeshType::PointsContainer::New();
		for(unsigned int pointID = 0; pointID < points->Size(); pointID++)
		{
			PointType curPoint = points->GetElement(pointID);
			PointType sftPoint;
			for(unsigned int dim=0; dim < 3; dim++)
			{
				sftPoint[dim] = curPoint[dim] - MC[dim];
			}
			newpts->InsertElement(pointID,sftPoint);
		}
		//Set the Shifted Points back to "points"
		points = newpts;
		surfaceMesh->SetPoints(points);
		SOMesh->SetMesh(surfaceMesh);
		MeshWriterType::Pointer writer = MeshWriterType::New();
		writer->SetInput(SOMesh);
		writer->SetFileName(outputFilename);
		writer->Update();

	}else if(GetGaussianCurvaturesOn || GetMeanCurvaturesOn) 
	{
		std::cout << "Arguments " << " " << argv[0] << " " << argv[1] << " " << argv[2] << " " << argv[3] << std::endl;

		// Read PolyData info
		vtkPolyDataReader *polyIn = vtkPolyDataReader::New();
		polyIn->SetFileName(argv[2]); 
		polyIn->Update();
		vtkPolyData* polydata = polyIn->GetOutput();

	
		//Calculate curvatures
		vtkCurvatures* curve = vtkCurvatures::New();
		curve->SetInput(polydata);

		if( GetGaussianCurvaturesOn){
			curve->SetCurvatureTypeToGaussian();}
		if (GetMeanCurvaturesOn){
			curve->SetCurvatureTypeToMean();
		}

		//curve->SetCurvatureTypeToMaximum ();
		//curve->SetCurvatureTypeToMinimum ();


		//Writing curvatures out
		std::ofstream output ( argv[3] ) ;
		curve->Update();
		
		vtkPolyData* polydataCurv = curve->GetOutput();
	

		unsigned int nPoints = polydataCurv->GetNumberOfPoints();
		//std::cout << polydataCurv->GetPointData()->GetScalars() << std::endl;
		
		/*vtkDataArray *Array;
		if(  GetGaussianCurvaturesOn){
			Array = polydataCurv->GetPointData()->GetArray("Gauss_Curvature");
		}
		if (GetMeanCurvaturesOn){
			Array =polydataCurv->GetPointData()->GetArray("Mean_Curvature");
		}
		*/
		vtkDoubleArray *Array = vtkDoubleArray::SafeDownCast(polydataCurv->GetPointData()->GetScalars());
		

		output << "NUMBER_OF_POINTS=" << nPoints << std::endl ;
		output << "DIMENSION=1" << std::endl << "TYPE=Scalar" << std::endl;

		if(Array)
		{
			for(unsigned int i = 0; i < nPoints; i++)
			{
				//std::cout << "Got array." << std::endl;
				double curv;
				//debug
				curv= Array->GetValue(i);
				//std::cout << "Curvature1: " << curv << std::endl;
				//Array->GetTuple(i, &curv);
				output << float(curv*1000) << std::endl ;
				std::cout << "Curvature: " << curv << std::endl;
			}
		}
		else
			std::cout << "Got no array? Error!" << std::endl;


		/*output << "NUMBER_OF_POINTS=" << nPoints << std::endl ;
		  output << "DIMENSION=1" << std::endl << "TYPE=Scalar" << std::endl ;

		  for ( int i = 0 ; i < nPoints ; i++ )
		  {
		//polydata >> point[0] >> point[1] >> point[2] ;
		//output << point[channel] << std::endl ;
		} */

		output.close ();


		/*vtkPolyDataWriter *SurfaceWriter = vtkPolyDataWriter::New();
		  SurfaceWriter->SetInput(curve->GetOutput());
		  SurfaceWriter->SetFileName(argv[3]);
		  SurfaceWriter->Update();
		  std::cout << "Writing new mesh " << argv[3] << std::endl;*/

		std::cout << "Finished!!" << std::endl;

	}


	else{
		std::cout << "No operation to do -> exiting" << std::endl;
		exit(-1);
	}

	return 0 ; 
}
