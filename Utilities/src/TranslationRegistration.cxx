
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif




#include "itkImageRegistrationMethod.h"
#include "itkMeanSquaresImageToImageMetric.h"
#include "itkMattesMutualInformationImageToImageMetric.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkImage.h"

#include "itkRescaleIntensityImageFilter.h"
#include "itkCenteredRigid2DTransform.h"
#include "itkTranslationTransform.h"


#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkResampleImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"



#include "itkCommand.h"
class CommandIterationUpdate : public itk::Command 
{
	public:
		typedef  CommandIterationUpdate   Self;
		typedef  itk::Command             Superclass;
		typedef itk::SmartPointer<Self>  Pointer;
		itkNewMacro( Self );
	protected:
		CommandIterationUpdate() {};
	public:
		typedef itk::RegularStepGradientDescentOptimizer     OptimizerType;
		typedef   const OptimizerType   *    OptimizerPointer;

		void Execute(itk::Object *caller, const itk::EventObject & event)
		{
			Execute( (const itk::Object *)caller, event);
		}

		void Execute(const itk::Object * object, const itk::EventObject & event)
		{
			OptimizerPointer optimizer = 
				dynamic_cast< OptimizerPointer >( object );
			if( ! itk::IterationEvent().CheckEvent( &event ) )
			{
				return;
			}
			std::cout << optimizer->GetCurrentIteration() << "   ";
			std::cout << optimizer->GetValue() << "   ";
			std::cout << optimizer->GetCurrentPosition() << std::endl;
		}
};


int main( int argc, char *argv[] )
{
	if( argc < 4 )
	{
		std::cerr << "Missing Parameters " << std::endl;
		std::cerr << "Usage: " << argv[0];
		std::cerr << " fixedImageFile  movingImageFile ";
		std::cerr << " outputImagefile  [differenceAfterRegistration] ";
		std::cerr << " [differenceBeforeRegistration] ";
		std::cerr << " [initialStepLength]  "<< std::endl;
		return EXIT_FAILURE;
	}

	const    unsigned int    Dimension = 3;
	typedef  float    PixelType;

	typedef itk::Image< PixelType, Dimension >  FixedImageType;
	typedef itk::Image< PixelType, Dimension >  MovingImageType;



	typedef itk::TranslationTransform< double > TransformType;


	typedef itk::RegularStepGradientDescentOptimizer       OptimizerType;

	typedef itk:: LinearInterpolateImageFunction< 
		MovingImageType,
		double          >    InterpolatorType;
	typedef itk::ImageRegistrationMethod< 
		FixedImageType, 
		MovingImageType >    RegistrationType;

	OptimizerType::Pointer      optimizer     = OptimizerType::New();
	InterpolatorType::Pointer   interpolator  = InterpolatorType::New();
	RegistrationType::Pointer   registration  = RegistrationType::New();

	typedef itk::MattesMutualInformationImageToImageMetric< 
		FixedImageType, 
		MovingImageType >    MetricType;
	MetricType::Pointer         metric        = MetricType::New();
	registration->SetMetric(        metric        );

	//MetricType::Pointer         metric        = MetricType::New();

	// registration->SetMetric(        metric        );
	registration->SetOptimizer(     optimizer     );
	registration->SetInterpolator(  interpolator  );



	TransformType::Pointer  transform = TransformType::New();
	registration->SetTransform( transform );



	typedef itk::ImageFileReader< FixedImageType  > FixedImageReaderType;
	typedef itk::ImageFileReader< MovingImageType > MovingImageReaderType;

	FixedImageReaderType::Pointer  fixedImageReader  = FixedImageReaderType::New();
	MovingImageReaderType::Pointer movingImageReader = MovingImageReaderType::New();

	fixedImageReader->SetFileName(  argv[1] );
	movingImageReader->SetFileName( argv[2] );


	registration->SetFixedImage(    fixedImageReader->GetOutput()    );
	registration->SetMovingImage(   movingImageReader->GetOutput()   );
	fixedImageReader->Update();

	registration->SetFixedImageRegion( 
			fixedImageReader->GetOutput()->GetBufferedRegion() );


	fixedImageReader->Update();
	movingImageReader->Update();


	typedef FixedImageType::SpacingType    SpacingType;
	typedef FixedImageType::PointType      OriginType;
	typedef FixedImageType::RegionType     RegionType;
	typedef FixedImageType::SizeType       SizeType;


	FixedImageType::Pointer fixedImage = fixedImageReader->GetOutput();


	const SpacingType fixedSpacing = fixedImage->GetSpacing();
	const OriginType  fixedOrigin  = fixedImage->GetOrigin();
	const RegionType  fixedRegion  = fixedImage->GetLargestPossibleRegion(); 
	const SizeType    fixedSize    = fixedRegion.GetSize();


	MovingImageType::Pointer movingImage = movingImageReader->GetOutput();

	const SpacingType movingSpacing = movingImage->GetSpacing();
	const OriginType  movingOrigin  = movingImage->GetOrigin();
	const RegionType  movingRegion  = movingImage->GetLargestPossibleRegion();
	const SizeType    movingSize    = movingRegion.GetSize();



	transform->SetIdentity();

	registration->SetInitialTransformParameters( transform->GetParameters() );



	double initialStepLength = 5;


	if( argc > 6 )
	{
		initialStepLength = atof( argv[6] );
	}


	// optimizer->SetRelaxationFactor( 0.6 );
	optimizer->SetMaximumStepLength( initialStepLength ); 
	optimizer->SetMinimumStepLength( 0.01 );
	optimizer->SetNumberOfIterations( 200 );




	CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
	optimizer->AddObserver( itk::IterationEvent(), observer );

	try 
	{ 
		registration->StartRegistration(); 
	} 
	catch( itk::ExceptionObject & err ) 
	{ 
		std::cerr << "ExceptionObject caught !" << std::endl; 
		std::cerr << err << std::endl; 
		return EXIT_FAILURE;
	} 

	OptimizerType::ParametersType finalParameters = 
		registration->GetLastTransformParameters();


	const double finalTranslationX    = finalParameters[0];
	const double finalTranslationY    = finalParameters[1];
	const double finalTranslationZ    = finalParameters[2];


	const unsigned int numberOfIterations = optimizer->GetCurrentIteration();

	const double bestValue = optimizer->GetValue();





	std::cout << "Result = " << std::endl;

	std::cout << " Translation X = " << finalTranslationX  << std::endl;
	std::cout << " Translation Y = " << finalTranslationY  << std::endl;
	std::cout << " Translation Z = " << finalTranslationZ  << std::endl;
	std::cout << " Iterations  = " << numberOfIterations << std::endl;
	std::cout << " Metric value = " << bestValue          << std::endl;




	typedef itk::ResampleImageFilter< 
		MovingImageType, 
		FixedImageType >    ResampleFilterType;

	TransformType::Pointer finalTransform = TransformType::New();

	finalTransform->SetParameters( finalParameters );

	ResampleFilterType::Pointer resample = ResampleFilterType::New();

	resample->SetTransform( finalTransform );
	resample->SetInput( movingImageReader->GetOutput() );

	resample->SetSize(    fixedImage->GetLargestPossibleRegion().GetSize() );
	resample->SetOutputOrigin(  fixedImage->GetOrigin() );
	resample->SetOutputSpacing( fixedImage->GetSpacing() );
	resample->SetOutputDirection( fixedImage->GetDirection() );
	resample->SetDefaultPixelValue( 100 );

	typedef itk::ImageFileWriter< FixedImageType >  WriterFixedType;

	WriterFixedType::Pointer  writer =  WriterFixedType::New();

	writer->SetFileName( argv[3] );

	writer->SetInput( resample->GetOutput()   );

	try
	{
		writer->Update();
	}
	catch( itk::ExceptionObject & excp )
	{ 
		std::cerr << "ExceptionObject while writing the resampled image !" << std::endl; 
		std::cerr << excp << std::endl; 
		return EXIT_FAILURE;
	} 



	typedef itk::Image< float, Dimension > DifferenceImageType;

	typedef itk::SubtractImageFilter< 
		FixedImageType, 
		FixedImageType, 
		DifferenceImageType > DifferenceFilterType;

	DifferenceFilterType::Pointer difference = DifferenceFilterType::New();

	typedef  unsigned char  OutputPixelType;

	typedef itk::Image< OutputPixelType, Dimension > OutputImageType;

	typedef itk::RescaleIntensityImageFilter< 
		DifferenceImageType, 
		OutputImageType >   RescalerType;

	RescalerType::Pointer intensityRescaler = RescalerType::New();

	intensityRescaler->SetOutputMinimum(   0 );
	intensityRescaler->SetOutputMaximum( 255 );

	difference->SetInput1( fixedImageReader->GetOutput() );
	difference->SetInput2( resample->GetOutput() );

	resample->SetDefaultPixelValue( 1 );

	intensityRescaler->SetInput( difference->GetOutput() );  

	typedef itk::ImageFileWriter< OutputImageType >  WriterType;

	WriterType::Pointer      writer2 =  WriterType::New();

	writer2->SetInput( intensityRescaler->GetOutput() );


	try
	{
		// Compute the difference image between the 
		// fixed and moving image after registration.
		if( argc > 4 )
		{
			writer2->SetFileName( argv[4] );
			writer2->Update();
		}

		// Compute the difference image between the 
		/*   // fixed and resampled moving image after registration.
		     TransformType::Pointer identityTransform = TransformType::New();
		     identityTransform->SetIdentity();
		     resample->SetTransform( identityTransform );
		     if( argc > 5 )
		     {
		     writer2->SetFileName( argv[5] );
		     writer2->Update();
		     }
		     */  }
		catch( itk::ExceptionObject & excp )
		{
			std::cerr << "Error while writing difference images" << std::endl;
			std::cerr << excp << std::endl;
			return EXIT_FAILURE;
		}







		return EXIT_SUCCESS;
}


