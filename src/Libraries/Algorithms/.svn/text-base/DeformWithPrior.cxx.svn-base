
#include "DeformWithPrior.h"
#include "itkMesh.h"
#include "itkNumericTraits.h"
#include <itkConjugateGradientOptimizer.h>
#include <vnl/vnl_math.h>
#include "vnl/vnl_matrix_fixed.h"
#include "vnl/vnl_math.h"
#include <cstdlib>
#include "itkGradientRecursiveGaussianImageFilter.h"
#include "itkGradientImageFilter.h"
#include "itkGradientToMagnitudeImageFilter.h"
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkMeshToMeshFilter.h"

#define MAX_OBJ_VAL 1000000


void DeformationCostFunction
::initialize(PCAStats * pcaStats,   MeshType::Pointer meanObj, 
			 const  ImageType::Pointer image, const int NumberOfUsedPCs, const float sigmaGaussian)
{  
	if (meanObj->GetCells() ) 
	{
		shapeTYPE = MESH;
	}else
	{
		shapeTYPE= PTS;
	}

	m_image = image;
	m_pcaStats = pcaStats;
	m_meanObj = meanObj;
	m_SigmaGaussian = sigmaGaussian;

	m_ShapePriorWeight = 1;
	m_ImageMatchWeight = 1;
	m_StepSize = 0.1;
	VERBOSE = false;

	m_NumberOfUsedPCs = NumberOfUsedPCs;
	m_NumberOfNodes = meanObj->GetNumberOfPoints();

	//initiallize  forces 
	m_forces = MeshType::New();
	PointsContainer::Pointer myForces = m_forces->GetPoints();
	myForces->Reserve(m_NumberOfNodes);
	PointsContainerIterator  forces = myForces->Begin();



	PointType d;
	d.Fill(0.0);

	PointsContainerIterator points = meanObj->GetPoints()->Begin();

	while(points != meanObj->GetPoints()->End())
	{  
		forces.Value()=d;
		++points;
		++forces;
	}

	if (shapeTYPE== MESH)
	{
		// initialize normals for mesh object 
		m_normals = MeshType::New();
		PointsContainer::Pointer myNormals = m_normals->GetPoints();
		myNormals->Reserve(m_NumberOfNodes);
		PointsContainerIterator normals = myNormals->Begin();

		PointsContainerIterator points = meanObj->GetPoints()->Begin();
		while(points != meanObj->GetPoints()->End())
		{  
			normals.Value()=d;
			++points;

			++normals;
		}
		m_normals = ComputeNormals(meanObj);//m_normals are computed from meanMesh
	}
	else{
		m_normals=NULL;
	}

	//clculate gradient image
	ComputeGradientImage(m_image);//m_gradientImageImage
}






//----------------------- SCDSCostFunction----------------------
double DeformationCostFunction::GetValue( const ParametersType & position ) const
{ //energy function value, for debugging and printing

	double objFuncVal =0;
	double imageMatchFuncVal = 0;
	double shapePriorFuncVal = MahalDistance( position);
	//std::cout<<"mahal dis = "<<shapePriorFuncVal<<std::endl;

	if (shapePriorFuncVal/position.size()>5)//too far off mean
	{
		imageMatchFuncVal = MAX_OBJ_VAL;
	}else
	{
		MeshType::Pointer tObject = GenerateObject(position);//update normals, locations

		imageMatchFuncVal = ImageMatchValue(tObject);//project normals onto second-order gradients
		//std::cout<<"imageMatchFuncVal = "<<imageMatchFuncVal<<std::endl;
	}

	objFuncVal = m_ShapePriorWeight * shapePriorFuncVal + m_ImageMatchWeight * imageMatchFuncVal;

	if (VERBOSE){
		std::cout<<"objFuncVal = "<<objFuncVal<<" "<<position <<"  : mahal dis-- "<< m_ShapePriorWeight *shapePriorFuncVal<<"  imageMatchFuncVal- "<<imageMatchFuncVal<<std::endl;
	}

	return objFuncVal;

}



MeshType::Pointer  DeformationCostFunction::
GenerateObject(ParametersType position) const
{
	MeshType::Pointer tObj = MeshType::New();

	VectorType tVec = m_pcaStats->GenObjVec(*(dynamic_cast<VectorType *> (&position)));




	PointsContainer::Pointer myLocations = tObj->GetPoints();
	myLocations->Reserve(m_NumberOfNodes);
	PointsContainerIterator locations = myLocations->Begin();

	int k = 0;
	while (locations!=myLocations->End()){
		locations.Value()[0]=tVec[k++];
		locations.Value()[1]=tVec[k++];
		locations.Value()[2]=tVec[k++];
		++locations;
	}


	if (shapeTYPE == MESH){
		tObj->SetCells(m_meanObj->GetCells());
	}

	return tObj;

}


void DeformationCostFunction::GetDerivative( const ParametersType & position, 
											DerivativeType & derivative ) const
{//derivatives used for the search

	derivative = DerivativeType(m_NumberOfUsedPCs);

	for(int i= 0; i< m_NumberOfUsedPCs; i++){

		ParametersType newP = position;
		newP[i] = position[i]  + 0.5 * m_StepSize;
		double upper = GetValue(newP);

		newP[i] = position[i] - 0.5 * m_StepSize;
		double lower = GetValue(newP);
		derivative[i] = (upper - lower )/m_StepSize;

	}
	if (VERBOSE)
	{
		std::cout<< "der="<<derivative<< std::endl;
	}

}



MeshType::Pointer DeformationCostFunction::GradientFit(MeshType::Pointer tObj)const 
{

	VectorType vec_nor(3);
	VectorType vec_loc(3);


	MeshType::Pointer tForces = MeshType::New();
	PointsContainer::Pointer myForces = tForces->GetPoints();
	myForces->Reserve( m_NumberOfNodes);
	PointsContainerIterator  forces = myForces->Begin();

	PointsContainer::Pointer   myLocations = tObj->GetPoints();
	PointsContainerIterator  locations = myLocations->Begin();

	MeshType::Pointer tNormals = ComputeNormals(tObj);
	PointsContainer::Pointer   myNormals = tNormals->GetPoints();
	PointsContainerIterator   normals = myNormals->Begin();


	// New gradient fit method testing. 
	long int count = 0;

	while( forces != myForces->End() ) {



		GradientImageType::PointType point;
		point[0]= locations.Value()[0];
		point[1]= locations.Value()[1];
		point[2]= locations.Value()[2];


		GradientImageType::IndexType pixelIndex;

		bool  isInside = m_gradientImage->TransformPhysicalPointToIndex(point, pixelIndex);


		GradientType vec_for;
		if ( !isInside) {    
			// if the model points goes outside of the image , find the nearest boundaoy point
			GradientImageType::SizeType size = m_gradientImage->GetBufferedRegion().GetSize();

			float deltaWeight = 0;
			GradientImageType::IndexType delta;
			delta[0]=0;delta[1] = 0; delta[2]=0;
			for( int i=0; i<3; i++)
			{
				if (pixelIndex[i] > long(size[i]-1)){
					delta[i] = pixelIndex[i] - size[i]+1;
					pixelIndex[i] = long(size[i])-1;}

				if (pixelIndex[i] < 0){
					delta[i] = - pixelIndex[i];
					pixelIndex[i] = 0;

				}
				deltaWeight = deltaWeight +  delta[i] * delta[i];// quadratic penalty, further distance:bigger penalty

			}

			count ++;

			vec_for = (1 + deltaWeight) *  m_gradientImage->GetPixel(pixelIndex);



		}else{
			vec_for = m_gradientImage->GetPixel(pixelIndex);
		}


		/*if ( shapeTYPE == MESH) {


			vec_nor[0] = normals.Value()[0];
			vec_nor[1] = normals.Value()[1];
			vec_nor[2] = normals.Value()[2];

			//project the forces onto the normals
			double mag = vec_for[0] * vec_nor[0] + vec_for[1] * vec_nor[1] + vec_for[2] * vec_nor[2];


			forces.Value()[0] = mag * vec_nor[0];
			forces.Value()[1] = mag * vec_nor[1]; 
			forces.Value()[2] = mag * vec_nor[2]; 
			++normals;
		}
		else{//for PTs*/
			//only use the gradient magnitude
			forces.Value()[0] = vec_for[0];
			forces.Value()[1] = vec_for[1]; 
			forces.Value()[2] = vec_for[2]; 
		//}
		++forces;
		++locations;

	}

	if (count>0 && VERBOSE){
		std::cout<<count<< " pts outside of the image"<<std::endl;
		/*if (count>100)
		{
		forces = myForces->Begin();
		while( forces != myForces->End() ) {

		forces.Value()[0] = 1000;
		forces.Value()[1] = 1000; 
		forces.Value()[2] = 1000; 
		++forces;
		}
		}
		*/

	}

	return tForces;

}




/* Compute normals. */
MeshType::Pointer DeformationCostFunction::ComputeNormals(MeshType::Pointer tObj) const
{
	const unsigned long *tp;
	PointType v1, v2, v3, v4, d;
	v1.Fill(0.);
	v2.Fill(0.);
	v3.Fill(0.);
	d.Fill(0.);

	double coa, cob, coc ;
	double absvec ;

	CellsContainer::Pointer  myCells = m_meanObj->GetCells();// tObj only contains locations
	CellsContainerIterator   cells = myCells->Begin();


	MeshType::Pointer Normals = MeshType::New();
	PointsContainer::Pointer myNormals = Normals->GetPoints();
	myNormals->Reserve(m_NumberOfNodes);
	PointsContainerIterator normals = myNormals->Begin();


	while( normals != myNormals->End() ) {
		normals.Value() = d;
		++ normals;
	}

	while ( cells != myCells->End() ) 
	{
		tp = cells.Value()->GetPointIds();
		++cells;

		tObj->GetPoint (tp[0], &v1);
		tObj->GetPoint (tp[1], &v2);
		tObj->GetPoint (tp[2], &v3);

		coa = -(v1[1]*(v2[2]-v3[2]) + 
			v2[1]*(v3[2]-v1[2]) +
			v3[1]*(v1[2]-v2[2])) ;
		cob = -(v1[2] * (v2[0]-v3[0]) +
			v2[2]*(v3[0]-v1[0]) +
			v3[2]*(v1[0]-v2[0])) ;
		coc = -(v1[0] * (v2[1]-v3[1]) +
			v2[0]*(v3[1]-v1[1]) +
			v3[0]*(v1[1]-v2[1])) ;

		absvec = -sqrt ((double) ((coa*coa) + (cob*cob) + (coc*coc))) ;
		if (absvec==0)
		{//hack
			coa =0; cob=0;coc=1;absvec=1;

		}
		v4[0] = coa/absvec;
		v4[1] = cob/absvec;
		v4[2] = coc/absvec;
		Normals->GetPoint (tp[0], &v1);
		Normals->GetPoint (tp[1], &v2);
		Normals->GetPoint (tp[2], &v3);

		v1[0] += v4[0];
		v1[1] += v4[1];
		v1[2] += v4[2];

		v2[0] += v4[0];
		v2[1] += v4[1];
		v2[2] += v4[2];

		v3[0] += v4[0];
		v3[1] += v4[1];
		v3[2] += v4[2];

		Normals->SetPoint (tp[0], v1);
		Normals->SetPoint (tp[1], v2);
		Normals->SetPoint (tp[2], v3);

	}

	normals = myNormals->Begin();
	while( normals != myNormals->End() ) 
	{
		v1 = normals.Value();

		absvec = sqrt((double) ((v1[0]*v1[0]) + (v1[1]*v1[1]) + 
			(v1[2]*v1[2])));

		if (absvec==0)
		{//hack
			v1[0] =0; v1[1]=0; v1[2]=1;absvec=1;
		}

		v1[0] = v1[0]/absvec;
		v1[1] = v1[1]/absvec;
		v1[2] = v1[2]/absvec;

		normals.Value() = v1;
		++normals;
	}

	return Normals;
}


void DeformationCostFunction::ComputeGradientImage(ImageType::Pointer image){

	//only for debug purpose
	bool USE_GENERATED_GRAD_IMAGE = false;

	if (USE_GENERATED_GRAD_IMAGE)
	{
		std::cout<<"read gradient image.."<<std::endl;
		typedef itk::ImageFileReader <GradientImageType> ImageReaderType;
		ImageReaderType ::Pointer reader = ImageReaderType::New();
		reader->SetFileName("gradientImage_2order.mhd");
		try {
			reader->Update();
		}
		catch(itk::ExceptionObject &err){
			std::cerr<<err<<std::endl;
		}

		m_gradientImage = reader->GetOutput();

	}else{ 
		//generate the second order gradient images
		//first: guassian smooth the image
		//second: get teh magnitude of the first
		//third: get the gradient of the second
		std::cout<<"generate second derivative image..."<<std::endl;

		typedef itk::GradientRecursiveGaussianImageFilter <ImageType,GradientImageType>  grgFilterType;
		typedef itk::GradientImageFilter  <ImageType, float, float>  gFilterType;
		typedef itk::GradientToMagnitudeImageFilter <GradientImageType, ImageType> g2mFilterType;

		grgFilterType ::Pointer grgfilter = grgFilterType::New();
		gFilterType :: Pointer gfilter = gFilterType::New();
		g2mFilterType:: Pointer g2mfilter = g2mFilterType::New();

		grgfilter->SetInput(image);
		grgfilter->SetSigma(m_SigmaGaussian);

		g2mfilter->SetInput(grgfilter->GetOutput());

		gfilter->SetInput(g2mfilter->GetOutput());
		gfilter->Update();
		m_gradientImage = gfilter->GetOutput();
		std::cout<<"derivative image done"<<std::endl;

//# image window  0-256


		/*typedef itk::ImageFileWriter <ImageType> WriterType;
		WriterType ::Pointer writer = WriterType::New();
		writer->SetFileName("gradientMagImage.mhd");
		writer->SetInput(g2mfilter->GetOutput());
		try {
			writer->Update();
		}
		catch(itk::ExceptionObject &err){
			std::cerr<<err<<std::endl;
		}
		typedef itk::ImageFileWriter <GradientImageType> WriterType2;
		WriterType2 ::Pointer writer2 = WriterType2::New();
		writer2->SetFileName("gradientImage_2order.mhd");
		writer2->SetInput(m_gradientImage);
		try {
			writer2->Update();
		}
		catch(itk::ExceptionObject &err){
			std::cerr<<err<<std::endl;
		}*/
	}
}


double DeformationCostFunction::MahalDistance(ParametersType position) const
{
	int size = position.GetSize();
	double dis = 0;
	VectorType var(size);
	for( int i=0; i<size; i++){
		var[i] = position[i]; //* m_pcaStats.sigmas[i];
		dis = dis + var[i]*var[i];

	}
	return sqrt(dis);


}


double DeformationCostFunction::ImageMatchValue(MeshType::Pointer tObj)const
{
	//how well the image featuer matches with the current model
	//m_gradientImage :second order gradient derivative

	double val = 0;
	double absvec; 
	PointType v;

	MeshType::Pointer  tForces = GradientFit(tObj);//project surface normal onto g2 to get m_froces

	PointsContainer::Pointer myForces = tForces->GetPoints();

	PointsContainerIterator forces = myForces->Begin();

	//sum over the forces to get val
	while(forces != myForces->End())
	{
		v = forces.Value();

		absvec = sqrt((double) ((v[0]*v[0]) + (v[1]*v[1]) + (v[2]*v[2])));
		val = val + absvec;
		++forces;
	}

	val = val/tForces->GetNumberOfPoints() * 1000; //normalize, independe of the points set size

	return val;
}


//////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////



DeformWithPrior
::DeformWithPrior()
{
	m_deformationOptimizer = OptimizerType::New();
	m_MaxIterations = 100;
	m_image = NULL; 
}


DeformWithPrior
::~DeformWithPrior()
{//clean up memory

	if (&m_pcaStats!=NULL)
		delete (&m_pcaStats);
}




void DeformWithPrior
::initialize( PCAStats * pcaStats, const  MeshType::Pointer refMesh, const ImageType::Pointer image, float sigmaGuassian){
	m_image = image;
	m_pcaStats = pcaStats;
	m_meanObj = refMesh;
	m_targetObject = refMesh;
	m_NumberOfNodes = refMesh->GetNumberOfPoints();
	m_NumberOfUsedPCs = pcaStats->GetNumPCs();
	m_costFunction = DeformationCostFunction::New();
	m_costFunction->initialize(m_pcaStats,m_meanObj, m_image,m_NumberOfUsedPCs,sigmaGuassian);
	shapeTYPE= MESH;

}

void DeformWithPrior::
initialize(PCAStats *  pcaStats,  const  ImageType::Pointer image, const float sigmaGaussian){

	m_image = image;
	m_pcaStats = pcaStats;
	m_NumberOfNodes = pcaStats->GetDim()/3;

	//put the mean pts into the data structure of MeshType, leave the  cells empty
	m_meanObj = MeshType::New();
	PointsContainer::Pointer myPts = m_meanObj->GetPoints();
	myPts->Reserve(m_NumberOfNodes);
	PointsContainerIterator pts = myPts->Begin();

	int i=0;
	VectorType smean= pcaStats->GetMean();
	while( pts != myPts->End() ) {
		pts.Value()[0] = smean[i++];
		pts.Value()[1] = smean[i++];
		pts.Value()[2] = smean[i++];
		++ pts;
	}


	m_targetObject = m_meanObj;

	m_NumberOfUsedPCs = pcaStats->GetNumPCs();
	m_costFunction = DeformationCostFunction::New();
	m_costFunction->initialize(m_pcaStats,m_meanObj, m_image,m_NumberOfUsedPCs,sigmaGaussian);

	shapeTYPE= PTS;

}


void DeformWithPrior
::Deform()
{

	m_deformationOptimizer -> SetCostFunction(m_costFunction.GetPointer());//GetPointer()?


	vnlOptimizerType * vnlOptimizer = m_deformationOptimizer->GetOptimizer();

	const double F_Tolerance      = 1e-2;  // Function value tolerance
	const double G_Tolerance      = 1e-2;  // Gradient magnitude tolerance 
	const double X_Tolerance      = 1e-2;  // Search space tolerance
	const double Epsilon_Function = 1e-5;//1e-10; // Step
	const int    Max_Iterations   =  m_MaxIterations; // Maximum number of iterations

	vnlOptimizer->set_f_tolerance( F_Tolerance );
	vnlOptimizer->set_g_tolerance( G_Tolerance );
	vnlOptimizer->set_x_tolerance( X_Tolerance ); 
	vnlOptimizer->set_epsilon_function( Epsilon_Function );
	vnlOptimizer->set_max_function_evals( Max_Iterations );

	vnlOptimizer->set_check_derivatives(1);


	OptimizerType::ParametersType initialValue(m_NumberOfUsedPCs);       // constructor requires vector size

	initialValue.fill(0);

	m_deformationOptimizer->SetInitialPosition( initialValue );


	CommandIterationUpdateConjugateGradient::Pointer observer = 
		CommandIterationUpdateConjugateGradient::New();
	m_deformationOptimizer->AddObserver( itk::IterationEvent(), observer );
	//m_deformationOptimizer->AddObserver( itk::FunctionEvaluationIterationEvent(), observer );


	try 
	{
		m_deformationOptimizer->StartOptimization();
	}
	catch( itk::ExceptionObject & e )
	{
		std::cout << "Exception thrown ! " << std::endl;
		std::cout << "An error ocurred during Optimization" << std::endl;
		std::cout << "Location    = " << e.GetLocation()    << std::endl;
		std::cout << "Description = " << e.GetDescription() << std::endl;
		exit(-1);
	}


	std::cout << "Number of iters = " << m_deformationOptimizer->GetCurrentIteration()  << std::endl;
	std::cout << "Number of evals = " << vnlOptimizer->get_num_evaluations() << std::endl;    

	std::cout << "Report from vnl optimizer: " << std::endl;

	vnlOptimizer->diagnose_outcome( std::cout );

	std::cout << std::endl;

	//
	// check results to see if it is within range
	//

	OptimizerType::ParametersType finalPosition;
	finalPosition = m_deformationOptimizer->GetCurrentPosition();


	m_targetObject = m_costFunction->GenerateObject(finalPosition);



	std::cout << "Solution        = (";
	std::cout << finalPosition <<")"<< std::endl;  



	// Get the final value of the optimizer
	std::cout << "energy function : ";
	OptimizerType::MeasureType finalValue = m_deformationOptimizer->GetValue();

	std::cout << finalValue << std::endl;




}

