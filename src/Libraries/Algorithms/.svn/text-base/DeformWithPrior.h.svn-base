#ifndef DEFORM_WITH_PRIOR_H
#define DEFORM_WITH_PRIOR_H


#include "itkTriangleCell.h"
#include "itkImage.h"
#include "itkImageRegionIterator.h"
#include "itkCovariantVector.h"
#include "itkCommand.h"
#include "itkConjugateGradientOptimizer.h"
#include "itkMeshSpatialObject.h"
#include "itkDefaultDynamicMeshTraits.h"
#include "itkPointSet.h"
#include "PCAStats.h"
//#include "itkVector.h"

#define MESH 1
#define PTS 2


typedef itk::CovariantVector<float, 3> GradientType;
typedef itk::Image<GradientType, 3>   GradientImageType;


typedef float								PixelType;
typedef itk::Image<PixelType, 3>			ImageType;


typedef itk::ConjugateGradientOptimizer  OptimizerType;


//  itkMesh object, subclass of itkPointSet 
//both point set and triangle meshe represnetation are in MeshType )
typedef itk::DefaultDynamicMeshTraits<double, 3, 3> MeshTraits;
typedef itk::Mesh<double,3, MeshTraits>				MeshType;
typedef MeshType::CellType							CellType;

typedef itk::TriangleCell< CellType >				TriangleType;

typedef MeshType::CellsContainer          CellsContainer;
typedef MeshType::CellsContainer::Iterator			CellsContainerIterator;

typedef MeshType::PointType		PointType;
typedef MeshType::PointsContainer    PointsContainer;
typedef MeshType::PointsContainer::Iterator   PointsContainerIterator;

typedef  OptimizerType::InternalOptimizerType  vnlOptimizerType;



//cost function for itk conjugate gradient filter
class DeformationCostFunction: public itk::SingleValuedCostFunction{

	public:
		typedef DeformationCostFunction           Self;
		typedef itk::SingleValuedCostFunction     Superclass;
		typedef itk::SmartPointer<Self>           Pointer;
		typedef itk::SmartPointer<const Self>     ConstPointer;
		itkNewMacro( Self );
		itkTypeMacro( conjugateCostFunction, SingleValuedCostFunction );

		typedef Superclass::ParametersType        ParametersType;
		typedef Superclass::DerivativeType        DerivativeType;

		typedef double MeasureType ;


		DeformationCostFunction(){}

		void initialize( PCAStats * pcaStats,  MeshType::Pointer meanObj, 
			const  ImageType::Pointer image, const int NumberOfSsedPCs, const float sigmaGaussian);

		//energy function value, for debugging and printing
		double GetValue( const ParametersType & position )const;


		//derivatives used for the search
		void GetDerivative( const ParametersType & position, DerivativeType & derivative )const ;

	
		unsigned int GetNumberOfParameters(void) const{return m_NumberOfUsedPCs;}

		MeshType::Pointer GenerateObject(ParametersType position)const ;

		MeshType::Pointer getNormals()const{return m_normals;}


		void setShapePriorWeight(float w) { m_ShapePriorWeight = w;}
		float getShapePriorWeight()const { return m_ShapePriorWeight;}

		void setImageMatchWeight(float w){ m_ImageMatchWeight = w;}
		float getImageMatchWeight()const { return m_ImageMatchWeight;}

		void setStepSizeForDerivative (float w) { m_StepSize = w;}
		float getStepSizeForDerivative()const { return m_StepSize;}

		void setNumberOfUsedPCs(int w) { m_NumberOfUsedPCs = w;}
		float getNumberOfUsedPCs()const { return m_NumberOfUsedPCs;}

		//void setGaussianSigma(float w){ m_SigmaGaussian =w;}
		float getGaussianSigma() const { return m_SigmaGaussian;}

		void setVerbose(bool w){VERBOSE = w;}

         int shapeTYPE;

	private:
		

		MeshType::Pointer ComputeNormals(MeshType::Pointer tObj) const;//compute mesh surface normals

		double MahalDistance(ParametersType position) const;// compute model penalty

		double ImageMatchValue(MeshType::Pointer tObj) const;

		MeshType::Pointer GradientFit(MeshType::Pointer tObj)const; //compute exteral force (  LoG operator, doc product with normals )

		void ComputeGradientImage(ImageType::Pointer image);//second-order gradient

		
		MeshType::Pointer m_meanObj;

        PCAStats * m_pcaStats;

		ImageType::Pointer m_image;

		GradientImageType::Pointer m_gradientImage;

		MeshType::Pointer m_normals;

		MeshType::Pointer m_forces;


		int       m_NumberOfNodes;
		int       m_NumberOfCells;
		int       m_ImageWidth;      /** Image size */
		int       m_ImageHeight;
		int       m_ImageDepth;
		int       m_StepThreshold;
		int       m_NumberOfUsedPCs;
		float	m_ShapePriorWeight;
		float	m_ImageMatchWeight;
		float	m_StepSize;
		float	m_SigmaGaussian;

		int shapeType; // mesh or pts 
		bool VERBOSE;

};








//for checking the optimizer iterations

class CommandIterationUpdateConjugateGradient : public itk::Command 
{
	public:
		typedef  CommandIterationUpdateConjugateGradient   Self;
		typedef  itk::Command             Superclass;
		typedef itk::SmartPointer<Self>  Pointer;
		itkNewMacro( Self );
	protected:
		CommandIterationUpdateConjugateGradient() 
		{
			m_IterationNumber=0;
		}
	public:
		typedef itk::ConjugateGradientOptimizer   OptimizerType;
		typedef   const OptimizerType   *    OptimizerPointer;

		void Execute(itk::Object *caller, const itk::EventObject & event)
		{
			Execute( (const itk::Object *)caller, event);
		}

		void Execute(const itk::Object * object, const itk::EventObject & event)
		{
			OptimizerPointer optimizer = 
				dynamic_cast< OptimizerPointer >( object );
			/*if( m_FunctionEvent.CheckEvent( &event ) )
			  {
			  std::cout << m_IterationNumber++ << "   ";
			  std::cout << optimizer->GetCachedValue() << "   ";
			  std::cout << optimizer->GetCachedCurrentPosition() << std::endl;
			  }
			  else */if( m_GradientEvent.CheckEvent( &event ) )
			{
				std::cout << m_IterationNumber++ << "  ";
				std::cout << "energy =" <<optimizer->GetCachedValue() << std::endl;
				std::cout << "position = "<<optimizer->GetCachedCurrentPosition() << std::endl;
				//}
				std::cout << "Gradient = " << optimizer->GetCachedDerivative() << std::endl;
		}

}
private:
unsigned long m_IterationNumber;

itk::FunctionEvaluationIterationEvent m_FunctionEvent;
itk::GradientEvaluationIterationEvent m_GradientEvent;
};










class DeformWithPrior{

	public:

		DeformWithPrior();
		~DeformWithPrior();

		void initialize( PCAStats *  pcaStats,  const  ImageType::Pointer image, const float sigmaGaussian);
		void initialize( PCAStats *  pcaStats,  const MeshType::Pointer meanObj,  const  ImageType::Pointer image, const float sigmaGaussian);

		void Deform();


		DeformationCostFunction::Pointer getCostFunction(){return m_costFunction;}	



		void setNumberOfUsedPCs(int w) {m_NumberOfUsedPCs = w;m_costFunction->setNumberOfUsedPCs(w);}
		int getNumberOfUsedPCs()const{return m_NumberOfUsedPCs;}

		void setMaxInterations(unsigned int t){m_MaxIterations = t;}
		unsigned int getMaxInterations(){return m_MaxIterations;}

		MeshType::Pointer getTargetObject(){return m_targetObject;}

         int shapeTYPE;// MESH/PTS

	protected:

		MeshType::Pointer m_meanObj;
		MeshType::Pointer m_targetObject;
		OptimizerType::Pointer m_deformationOptimizer;
		PCAStats * m_pcaStats;
		ImageType::Pointer m_image;
		DeformationCostFunction::Pointer m_costFunction;

		unsigned int m_MaxIterations;         /** Number of iterations */
		int       m_NumberOfNodes;
		int       m_NumberOfCells;
		int       m_StepThreshold;
		int       m_NumberOfUsedPCs;


	private:


};




#endif
