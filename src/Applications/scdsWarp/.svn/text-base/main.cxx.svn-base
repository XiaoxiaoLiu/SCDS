#include <ctime>
#include "ApplicationUtils.h"
#include "Image.h"
#include "Vector3D.h"
#include "Array3D.h"
#include "Array3DUtils.h"
#include "ImageUtils.h"
#include "HField3DUtils.h"
#include "HField3DIO.h"
#include <string>
#include <sstream>
#include <vector>
#include <iostream>
#include <iomanip>
#include <exception>
#include "FluidWarp.h"
#include "FluidWarpParameters.h"
#include "SCDSAtlasBuilder.h"
#include "Timer.h"
#include <itkMultiThreader.h>
#include <itkTransformFileReader.h>
#include <itkTransformBase.h>
#include <itkAffineTransform.h>

#define appout std::cerr
const std::string PROGRAM_NAME = "scdsWarp";

typedef float VoxelType;

void initializeHField(Array3D< Vector3D<float> >* h,
		Array3D< Vector3D<float> >* hinv,
		itk::AffineTransform<float, 3>::Pointer affine,
		Vector3D<float> spacing,
		Vector3D<float> origin)
{
	Vector3D<unsigned int> size = h->getSize();
	itk::AffineTransform<float, 3>::Pointer invtransform = itk::AffineTransform<float, 3>::New();
	invtransform->SetCenter(affine->GetCenter());
	affine->GetInverse(invtransform);

	for (unsigned int z = 0; z < size.z; ++z) 
	{
		for (unsigned int y = 0; y < size.y; ++y) 
		{
			for (unsigned int x = 0; x < size.x; ++x) 
			{
				itk::Point<float, 3> p, tp, tinvp;
				p[0] = x * spacing[0] + origin[0];
				p[1] = y * spacing[1] + origin[1];
				p[2] = z * spacing[2] + origin[2];

				tp = affine->TransformPoint(p);
				tinvp = invtransform->TransformPoint(p);

				(*h)(x,y,z).set((tp[0] - origin[0]) /spacing[0],
						(tp[1] - origin[1]) /spacing[1],
						(tp[2] - origin[2]) /spacing[2]);

				if(hinv != NULL)
				{
					(*hinv)(x,y,z).set((tinvp[0] - origin[0]) /spacing[0],
							(tinvp[1] - origin[1]) /spacing[1],
							(tinvp[2] - origin[2]) /spacing[2]);
				}
			}
		}
	}

}

void printUsage()
{
	appout << "Usage: " 
		<< PROGRAM_NAME << " [OPTION]... transformation0   transformation2 .... "
		<< std::endl << std::endl;

	appout << "Options:" 
		<< std::endl
		<< "  -h, --help                         show this help"
		<< std::endl
		<< "  -v, --verbose                      verbose output"
		<< std::endl
		<< "  -I, --incompressible               compute incompressible flows"
		<< std::endl
		<< "   --numberOfImages=VAL              number of input images"   
		<< std::endl
		<< "  --initWithInputHField=BOOL          default true "
		<< std::endl
		<< "   --inputImageList image0 image1 ...   list of image file names"
		<< std::endl
		<< "   --inputHFieldList hfield0 hfield1..  list of hfield file names"
		<< std::endl
		<< "   --inputCompHFieldList hfield0 hfield1..  list of groud truth (for comparison) hfield file names"
		<< std::endl
		<< "   --optWeights <w1> <w2> ...     weights for body forces of image intensity and  prediced Hfields "
		<< std::endl
		<< "   --shapeScores <w1> <w2> ...       surrogate shapeScores for each image "
		<< std::endl
		<< "  -o, --outputImageFilenamePrefix=PREFIX          image filename prefix"
		<< std::endl
		<< "  -o, --outputDeformedImageFilenamePrefix=PREFIX  deformedimage filename prefix"
		<< std::endl
		<< "  -f, --outputHFieldFilenamePrefix=PREFIX    hfield filename prefix"
		<< std::endl
		<< "  -p, --outputHInvFieldFilenamePrefix=PREFIX inverse hfield filename prefix"
		<< std::endl
		<< "  -n, --intensityWindowMin=VAL       intensity window min value"
		<< std::endl
		<< "  -x, --intensityWindowMax=VAL       intensity window max value"
		<< std::endl
		<< "  -s, --scaleLevel=SCALE             add scale level (1,2,4,...)"
		<< std::endl << std::endl
		<< "  Algorithm parameters for the current scale level:"
		<< std::endl
		<< "      --alpha=VAL                    "
		<< std::endl
		<< "      --beta=VAL                     "
		<< std::endl
		<< "      --gamma=VAL                    "
		<< std::endl
		<< "      --maxPerturbation=VAL          "
		<< std::endl
		<< "  -i, --numberOfIterations=VAL       "
		<< std::endl
		<< std::endl
		<< "  LDMM parameters for the current scale level:"
		<< std::endl
		<< "      --numberOfTimeSteps=VAL        (NOT IMPLEMENTED!)"
		<< std::endl
		<< "      --epsilon=VAL                  (NOT IMPLEMENTED!)"
		<< std::endl
		<< "      --sigma=VAL                    (NOT IMPLEMENTED!)"
		<< std::endl << std::endl
		<< "  Global parameters:"
		<< std::endl
		<< "      --updateAverageAfterSubIterations=BOOL default is false"
		<< std::endl
		<< "      --deltaSelectionUseMean=BOOL   default uses individual deltas"
		<< std::endl
		<< "      --imageWeights <w1> <w2> ...   image weights for Frechet mean"
		<< std::endl << std::endl
		<< "  Miscellaneous parameters:          "
		<< std::endl
		<< "      --fftwMeasure=BOOL             default is true"
		<< std::endl
		<< "      --numberOfThreads=VAL          default is number of processors"
		<< std::endl
		<< "      --fftwPlan=FILENAME            FFTW plan (NOT IMPLEMENTED!)"
		<< std::endl << std::endl
		<< "  Extra output:" << std::endl 
		<< "      --extraOutputStride=VAL        write atlas data every VAL iterations"
		<< std::endl
		<< "      --writeVolume                  write atlas volume every Stride iterations"
		<< std::endl
		<< "      --writeXSlice=SLICE            write x slice of atlas volume every Stride iterations"
		<< std::endl
		<< "      --writeYSlice=SLICE            write y slice of atlas volume every Stride iterations"
		<< std::endl
		<< "      --writeZSlice=SLICE            write z slice of atlas volume every Stride iterations"
		<< std::endl;
	appout << std::endl << std::endl;

	appout << "Example: " << std::endl
		<< PROGRAM_NAME << "                                  \\"
		<< std::endl
		<< " --numberOfImages  2                              \\"
		<< std::endl
		<< " --inputImageList image0.mhd image1.mhd                 \\" 
		<< std::endl
		<< " --inputHFieldList hfield0.mhd hfield1.mhd              \\" 
		<< std::endl
		<< " --inputCompHFieldList gt_hfield0.mhd gt_hfield1.mhd              \\" 
		<< std::endl
		<< " --outputImageFilenamePrefix=averageImage_        \\"
		<< std::endl
		<< " --outputDeformedImageFilenamePrefix=deformedImage_ \\"
		<< std::endl
		<< " --outputHFieldFilenamePrefix=deformationField_   \\"
		<< std::endl
		<< " --scaleLevel=4 --numberOfIterations=100          \\"
		<< std::endl
		<< " --scaleLevel=2 --numberOfIterations=50           \\"
		<< std::endl
		<< " --scaleLevel=1 --numberOfIterations=25           \\"
		<< std::endl 
		<< std::endl
		<< "This downsamples twice, runs 100 iterations, upsamples to half "
		<< std::endl 
		<< "resolution, runs 50 iterations, and finally upsamples to full"
		<< std::endl
		<< "resolution and runs 25 iterations." << std::endl;
	appout << std::endl << std::endl;
	appout << " Xiaoxiao Liu (sharonxx@cs.unc.edu)" << std::endl;
	appout << " Brad Davis, Sarang Joshi (davisb@cs.unc.edu)" << std::endl;
	
}

int main(int argc, char **argv)
{
	Timer totalTimer;
	totalTimer.start();

	//
	// default parameters
	//
	enum FluidOutputMode {FluidSilent, FluidStandard, FluidVerbose};
	FluidOutputMode outputMode = FluidStandard;

	std::string outputImageFilenamePrefix = "";
	std::string outputDeformedImageFilenamePrefix = "";
	std::string outputHFieldFilenamePrefix = "";
	std::string outputHInvFieldFilenamePrefix = "";
	bool computeInverseHFields = false;
	std::vector<std::string> inputImageFilenames;
	std::vector<std::string> transformFilenames;
	std::vector<std::string> inputHFieldFilenames;
	bool COMPARE_HFIELD_ON=false;   
	std::vector<std::string> inputCompHFieldFilenames;

	std::vector<FluidWarpParameters> fluidParams;
	std::vector<float> scaleLevels;
	FluidWarpParameters defaultFluidParams;
	defaultFluidParams.alpha = 0.01;
	defaultFluidParams.beta  = 0.01;
	defaultFluidParams.gamma = 0.001;
	defaultFluidParams.maxPerturbation = 0.5;
	defaultFluidParams.numIterations = 50;
	defaultFluidParams.numBasis = 2000;
	defaultFluidParams.divergenceFree = false;

	bool useIntensityWindow = false;
	int extraOutputStride = -1;
	bool writeVolume = false;
	int writeXSlice = -1;
	int writeYSlice = -1;
	int writeZSlice = -1;
	float iwMin, iwMax;

	bool updateAverageAfterSubIterations = false;
	bool deltaSelectionUseMean = false;
	bool fftwMeasure = true;
	int  numberOfThreads = 
		itk::MultiThreader::GetGlobalDefaultNumberOfThreads();
	int  numberOfImages =0;
	double* imageWeights = NULL;
	double * optWeights =NULL;
	double * shapeScores = NULL;
	bool initWithInputHField = true;
	// parse command line
	//
	argv++; argc--;
	while (argc > 0)
	{
		std::string arg(argv[0]);

		if (arg.find("-h")==0 || arg.find("--help")==0)
		{
			// print help and exit
			printUsage();
			exit(0);
		}
		else if (arg.find("-v")==0 || arg.find("--verbose")==0)
		{
			// set output mode to verbose
			outputMode = FluidVerbose;
			ApplicationUtils::ParseOptionVoid(argc, argv);
		}
		else if (arg.find("-I")==0 || arg.find("--incompressible")==0)
		{
			// set output mode to incompressible
			defaultFluidParams.divergenceFree = true;
			ApplicationUtils::ParseOptionVoid(argc, argv);
		}	
		else if (arg.find("--numberOfImages")==0)
		{     
			int tmpImages = ApplicationUtils::ParseOptionInt(argc, argv);
			if (tmpImages > 0)
			{
				numberOfImages= tmpImages;
			}
		}
		else if(arg.find("--initWithInputHField")==0)
		{
			initWithInputHField = ApplicationUtils::ParseOptionBool(argc, argv);

		}
	    else if (arg.find("--inputImageList")==0)
		{     

			ApplicationUtils::ParseOptionVoid(argc, argv);

			//imageWeights = new double[numberOfImages];
			for (int i = 0; i < numberOfImages; ++i)
			{
				inputImageFilenames.push_back(argv[0]);
				--argc; ++argv;
			}

		}
		else if (arg.find("--inputHFieldList")==0)
		{     
			ApplicationUtils::ParseOptionVoid(argc, argv);
			for (int i = 0; i < numberOfImages; ++i)
			{
				inputHFieldFilenames.push_back(argv[0]);
				--argc; ++argv;
			}
		}else if (arg.find("--inputCompHFieldList")==0)
		{     
			ApplicationUtils::ParseOptionVoid(argc, argv);
			for (int i = 0; i < numberOfImages; ++i)
			{
				inputCompHFieldFilenames.push_back(argv[0]);
				--argc; ++argv;
			}
		}else if (arg.find("--optWeights")==0) 
		{  
			ApplicationUtils::ParseOptionVoid(argc, argv);
			optWeights = new double[4];
			for (int i = 0; i < 4; ++i)
			{
				optWeights[i] = atof(argv[0]);
				--argc; ++argv;
			}

		}else if (arg.find("--shapeScores")==0) 
		{  
			ApplicationUtils::ParseOptionVoid(argc, argv);
			shapeScores = new double[numberOfImages];
			for (int i = 0; i < numberOfImages; ++i)
			{
				shapeScores[i] = atof(argv[0]);
				--argc; ++argv;
			}

		}
		else if (arg.find("-s")==0 || arg.find("--scaleLevel")==0)
		{
			double newScaleLevel = 
				ApplicationUtils::ParseOptionDouble(argc, argv);
			scaleLevels.push_back(newScaleLevel);
			fluidParams.push_back(defaultFluidParams);
		}
		else if (arg.find("-o")==0 || arg.find("--outputImageFilenamePrefix")==0)
		{
			outputImageFilenamePrefix = 
				ApplicationUtils::ParseOptionString(argc, argv);
		}
		else if (arg.find("-d")==0 || arg.find("--outputDeformedImageFilenamePrefix")==0)
		{
			outputDeformedImageFilenamePrefix = 
				ApplicationUtils::ParseOptionString(argc, argv);
		}
		else if (arg.find("-f")==0 || arg.find("--outputHFieldFilenamePrefix")==0)
		{
			outputHFieldFilenamePrefix = 
				ApplicationUtils::ParseOptionString(argc, argv);
		}
		else if (arg.find("-p")==0 || 
				arg.find("--outputHInvFieldFilenamePrefix")==0)
		{
			outputHInvFieldFilenamePrefix = 
				ApplicationUtils::ParseOptionString(argc, argv);
			computeInverseHFields = true;
		}
		else if (arg.find("--extraOutputStride")==0)
		{
			extraOutputStride = ApplicationUtils::ParseOptionBool(argc, argv);
		}
		else if (arg.find("--writeVolume")==0)
		{
			writeVolume = ApplicationUtils::ParseOptionBool(argc, argv);
		}
		else if (arg.find("--writeXSlice")==0)
		{
			writeXSlice = ApplicationUtils::ParseOptionInt(argc, argv);
		}
		else if (arg.find("--writeYSlice")==0)
		{
			writeYSlice = ApplicationUtils::ParseOptionInt(argc, argv);
		}
		else if (arg.find("--writeZSlice")==0)
		{
			writeZSlice = ApplicationUtils::ParseOptionInt(argc, argv);
		}
		else if (arg.find("--alpha")==0)
		{
			double alpha = ApplicationUtils::ParseOptionDouble(argc, argv);

			if (fluidParams.size() == 0)
			{
				fluidParams.push_back(defaultFluidParams);
				scaleLevels.push_back(1.0F);
			}
			fluidParams.back().alpha = alpha;
		}
		else if (arg.find("--beta")==0)
		{
			double beta = ApplicationUtils::ParseOptionDouble(argc, argv);

			if (fluidParams.size() == 0)
			{
				fluidParams.push_back(defaultFluidParams);
				scaleLevels.push_back(1.0F);
			}
			fluidParams.back().beta = beta;
		}
		else if (arg.find("--gamma")==0)
		{
			double gamma = ApplicationUtils::ParseOptionDouble(argc, argv);

			if (fluidParams.size() == 0)
			{
				fluidParams.push_back(defaultFluidParams);
				scaleLevels.push_back(1.0F);
			}
			fluidParams.back().gamma = gamma;
		}
		else if (arg.find("--maxPerturbation")==0)
		{
			double maxPert = ApplicationUtils::ParseOptionDouble(argc, argv);

			if (fluidParams.size() == 0)
			{
				fluidParams.push_back(defaultFluidParams);
				scaleLevels.push_back(1.0F);
			}

			fluidParams.back().maxPerturbation = maxPert;
		}
		else if (arg.find("-i")==0 || arg.find("--numberOfIterations")==0)
		{
			int numIter = ApplicationUtils::ParseOptionInt(argc, argv);

			if (fluidParams.size() == 0)
			{
				fluidParams.push_back(defaultFluidParams);
				scaleLevels.push_back(1.0F);
			}
			fluidParams.back().numIterations = numIter;
		}
		else if (arg.find("--numberOfTimeSteps")==0)
		{
			appout << "numberOfTimeSteps is not currently implemented!" 
				<< std::endl;
			exit(0);
		}
		else if (arg.find("--epsilon")==0)
		{
			appout << "epsilon is not currently implemented!" 
				<< std::endl;
			exit(0);
		}
		else if (arg.find("--sigma")==0)
		{
			appout << "sigma is not currently implemented!" 
				<< std::endl;
			exit(0);
		}
		else if (arg.find("-n")==0 || arg.find("--intensityWindowMin")==0)
		{
			iwMin = ApplicationUtils::ParseOptionDouble(argc, argv);
			useIntensityWindow = true;
		}
		else if (arg.find("-x")==0 || arg.find("--intensityWindowMax")==0)
		{
			iwMax = ApplicationUtils::ParseOptionDouble(argc, argv);
			useIntensityWindow = true;
		}
		else if (arg.find("--updateAverageAfterSubIterations")==0)
		{
			updateAverageAfterSubIterations = 
				ApplicationUtils::ParseOptionBool(argc, argv);
		}
		else if (arg.find("--deltaSelectionUseMean")==0)
		{
			deltaSelectionUseMean = 
				ApplicationUtils::ParseOptionBool(argc, argv);
		}
		else if (arg.find("--fftwMeasure")==0)
		{
			fftwMeasure = ApplicationUtils::ParseOptionBool(argc, argv);
		}
		else if (arg.find("--numberOfThreads")==0)
		{
			int tmpThreads = ApplicationUtils::ParseOptionInt(argc, argv);
			if (tmpThreads > 0)
			{
				numberOfThreads = tmpThreads;
			}
		}
		else if (arg.find("--fftwPlan")==0)
		{
			appout << "The fftwPlan option is currently not implemented!" 
				<< std::endl;
			exit(0);
		}
		else if (arg.find("--imageWeights")==0)
		{
			ApplicationUtils::ParseOptionVoid(argc, argv);
			unsigned int numImages = inputImageFilenames.size();
			imageWeights = new double[numImages];
			for (int i = 0; i < numImages; ++i)
			{
				imageWeights[i] = atof(argv[0]);
				--argc; ++argv;
			}
		}

		else if (std::string(argv[0]).find(".txt") == 0)
		{
			transformFilenames.push_back(argv[0]);
			ApplicationUtils::ParseOptionVoid(argc, argv);
		}
		else{
			std::cerr<<"unrecogonized input parameters"<<std::endl;
			exit(0);
		}

	}

	unsigned int numImages = inputImageFilenames.size();
	if (numImages < 1)
	{
		printUsage();
		return 0;
	}

	if (fluidParams.size() == 0)
	{
		fluidParams.push_back(defaultFluidParams);
		scaleLevels.push_back(1.0F);
	}

	// Check # transforms == numimages if # transforms > 0
	bool initTransformFiles = !transformFilenames.empty();
	if(initTransformFiles)
	{
		if( transformFilenames.size() !=  inputImageFilenames.size())
		{
			std::cerr << "Affine initialization files provides, but # transforms (" << transformFilenames.size() << ") != # images (" << inputImageFilenames.size() << ")"  << std::endl;
			return EXIT_FAILURE;
		}
	}

	// Check # hfields == numimages if # hfields > 0
	bool initHFieldFiles = !inputHFieldFilenames.empty();
	if(initHFieldFiles)
	{
		if( inputHFieldFilenames.size() !=  inputImageFilenames.size())
		{
			std::cerr << "HField files provides, but # HFields (" << inputHFieldFilenames.size() << ") != # images (" << inputImageFilenames.size() << ")"  << std::endl;
			return EXIT_FAILURE;
		}
	}
	//
	// print parameters
	//
	using appout;
	using std::endl;
	appout << PROGRAM_NAME << " parameters..." << std::endl;
	appout << "Output Mode           : ";
	switch(outputMode)
	{
		case (FluidVerbose):
			appout << "Verbose" << std::endl;
			break;
		case (FluidStandard):
			appout << "Standard" << std::endl;
			break;
		case (FluidSilent):
			appout << "Silent" << std::endl;
			break;
		default:
			appout << "Unknown" << std::endl;
			printUsage();
			return(0);
	}      
	appout << "Num. Images           : " << numImages << std::endl;
	appout << "Weights               : ";
	if (imageWeights != NULL)
	{
		double imageWeightSum = 0;
		for(int i = 0; i < numImages; ++i)
		{
			appout << imageWeights[i] << ' ';
			imageWeightSum += imageWeights[i];
		}
		appout << std::endl;
		appout << "Sum of Weights        : " << imageWeightSum << std::endl;
	}
	else
	{
		//equal weights
		imageWeights= new double[numImages];
		for(int i = 0; i < numImages; ++i)
		{
			imageWeights[i] = 1.0/numImages;
		}

		appout << "(none)" << std::endl;




	}

	appout << "Output Image Prefix   : " 
		<< (outputImageFilenamePrefix.length() ? 
				outputImageFilenamePrefix : "(none)")
		<< std::endl;
	appout << "Output Deformed Image Prefix   : " 
		<< (outputDeformedImageFilenamePrefix.length() ? 
				outputDeformedImageFilenamePrefix : "(none)")
		<< std::endl;
	appout << "Output h Field Prefix : " 
		<< (outputHFieldFilenamePrefix.length() ? 
				outputHFieldFilenamePrefix : "(none)")
		<< std::endl;
	appout << "Output inverse h Field Prefix : " 
		<< (outputHInvFieldFilenamePrefix.length() ? 
				outputHInvFieldFilenamePrefix : "(none)")
		<< std::endl;
	appout << "Intensity Window      : ";
	if (useIntensityWindow)
	{
		appout << iwMin << "-" << iwMax << std::endl;
	}
	else
	{
		appout << "No Intensity Windowing" << std::endl;
	}

	appout << "FFTW Parameters       : "
		<< (fftwMeasure ? "Measure plan, " : "Don't measure plan, ")
		<< "Threads=" << numberOfThreads
		<< std::endl;
	appout << "Update atlas after every sub-iteration: "
		<< (updateAverageAfterSubIterations ? "true" : "false")
		<< std::endl;
	appout << "Delta selction use mean:                "
		<< (deltaSelectionUseMean ? "true" : "false")
		<< std::endl;

	for (int scale = 0; scale < (int) scaleLevels.size(); ++scale)
	{
		appout << "Scale                 : " 
			<< scaleLevels[scale] 
			<< std::endl;
		appout << "alpha                 : " 
			<< fluidParams[scale].alpha 
			<< std::endl;
		appout << "beta                  : " 
			<< fluidParams[scale].beta 
			<< std::endl;
		appout << "gamma                 : " 
			<< fluidParams[scale].gamma 
			<< std::endl;
		appout << "Max. Pert.            : " 
			<< fluidParams[scale].maxPerturbation 
			<< std::endl;
		appout << "Num. Iterations       : " 
			<< fluidParams[scale].numIterations 
			<< std::endl;
	}

	//
	// load images
	//
	appout << "Loading Images..." << std::endl;
	Image<VoxelType>** images = new Image<VoxelType>*[numImages];
	for (int i = 0; i < (int) numImages; ++i)
	{
		images[i] = new Image<VoxelType>;
		ApplicationUtils::LoadImageITK(inputImageFilenames[i].c_str(), *images[i]);
		appout << "   Loaded: " << inputImageFilenames[i] << std::endl;
		appout << "   Dimensions: " << images[i]->getSize() 
			<< std::endl;
		appout << "   Origin: " << images[i]->getOrigin() 
			<< std::endl;
		appout << "   Spacing: " << images[i]->getSpacing() 
			<< std::endl;
		VoxelType iMin, iMax;
		Array3DUtils::getMinMax(*images[i], iMin, iMax);
		appout << "   Intensity Range: " << iMin << "-" << iMax << std::endl;
		if (!useIntensityWindow)
		{
			iwMin = iMin;
			iwMax = iMax;
		}
		appout << "   Rescaling to [0,1]...";
		Array3DUtils::rescaleElements(*images[i], iwMin, iwMax, 0.0F, 1.0F);
		appout << "DONE" << std::endl;
		Array3DUtils::getMinMax(*images[i], iMin, iMax);
		appout << "   Intensity Range: " << iMin << "-" << iMax << std::endl;
	}  


	//
	//load hfields
	//
    Array3D< Vector3D<float> > ** inputHFields = NULL;
	if(inputHFieldFilenames.size()>0)
	{appout <<"Loading Hfields.." << std::endl;
	inputHFields = new Array3D< Vector3D<float> > * [numImages];
	for (int i = 0; i < (int) numImages; ++i)
	{
		Vector3D<unsigned int> hSize;
		Vector3D<double> hOrigin, hSpacing;
		inputHFields[i]=new Array3D<Vector3D<float> >(images[i]->getSize());
		std::cerr << "Loading hField: " << inputHFieldFilenames[i].c_str() << "...";
		HField3DIO::readMETAHeader(hSize,hOrigin,hSpacing,inputHFieldFilenames[i].c_str());
		HField3DIO::readMETA(*inputHFields[i], inputHFieldFilenames[i].c_str());
		appout  << "DONE" << std::endl;
		appout  << " Dimensions: " << hSize << std::endl;
		appout  << " Origin    : " << hOrigin << std::endl;
		appout  << " Spacing   : " << hSpacing << std::endl;

	}
	}
	//
	//load ground truth hfields
	//

	 Array3D< Vector3D<float> > ** inputCompHFields = NULL;
	if(!inputCompHFieldFilenames.empty())
	{
		COMPARE_HFIELD_ON=true;
			appout <<"Loading Comparing Hfields.." << std::endl;
	inputCompHFields = new Array3D< Vector3D<float> > * [numImages];
	for (int i = 0; i < (int) numImages; ++i)
	{
		Vector3D<unsigned int> hSize;
		Vector3D<double> hOrigin, hSpacing;
		inputCompHFields[i]=new Array3D<Vector3D<float> >(images[i]->getSize());
		std::cerr << "Loading hField: " << inputCompHFieldFilenames[i].c_str() << "...";
		HField3DIO::readMETAHeader(hSize,hOrigin,hSpacing,inputCompHFieldFilenames[i].c_str());
		HField3DIO::readMETA(*inputCompHFields[i], inputCompHFieldFilenames[i].c_str());
		appout  << "DONE" << std::endl;
		appout  << " Dimensions: " << hSize << std::endl;
		appout  << " Origin    : " << hOrigin << std::endl;
		appout  << " Spacing   : " << hSpacing << std::endl;

	}
	}

	//
	// build atlas at each scale level
	//
	Image<VoxelType>** scaledImages  = new Image<VoxelType>*[numImages];
	SCDSAtlasBuilder::VectorFieldType ** scaledInputHFields = new SCDSAtlasBuilder::VectorFieldType * [numImages];
	SCDSAtlasBuilder::VectorFieldType ** scaledInputCompHFields = new SCDSAtlasBuilder::VectorFieldType * [numImages];



	// write initial atlas image
	//
	// orignal size image
	// images and scaledImages are the same!!!
	// deformed/initilized scaledImages will be input to the atlasBuilder
	if (outputImageFilenamePrefix != "")
	{
		Image<VoxelType> iHat(*images[0]);



		//deform image if the intial hfield is not identity
		if(initWithInputHField == true){
			for(int imageIndex =0; imageIndex<numImages;imageIndex++)
			{
				scaledImages[imageIndex] = new Image<VoxelType>(*images[imageIndex]); //just for temporary usage
				HField3DUtils::apply(*images[imageIndex],*inputHFields[imageIndex],*scaledImages[imageIndex]);
			}



			Array3DUtils::weightedArithmeticMean(numImages, 
					(const Array3D<float>** const) scaledImages, 
					imageWeights,
					iHat);


			for (int i = 0; i < (int) numImages; ++i)
			{
				delete scaledImages[i];
			}  
		}else{

			Array3DUtils::weightedArithmeticMean(numImages, 
					(const Array3D<float>** const) images, 
					imageWeights,
					iHat);


		}



		std::ostringstream oss;
		oss << outputImageFilenamePrefix << "initial";
		appout << "Writing Average Image...";
		ImageUtils::writeMETA(iHat, oss.str().c_str());
		appout << "DONE" << std::endl;
	}







	//Image<VoxelType> iHat;
	SCDSAtlasBuilder::VectorFieldType ** h    = new   SCDSAtlasBuilder::VectorFieldType *[numImages];
	SCDSAtlasBuilder::VectorFieldType** hinv = 0;

	if (computeInverseHFields)
	{
		hinv = new   SCDSAtlasBuilder::VectorFieldType *[numImages];  
	}
	Vector3D<float> origin, spacing;
	origin = images[0]->getOrigin();
	for (int scale = 0; scale < (int) scaleLevels.size(); ++scale)
	{
		appout << "Scale: " << scaleLevels[scale] << std::endl;

		//
		// create images for this scale level
		//
		appout << "Downsampling Images...";
		int f = (int) scaleLevels[scale];
		if (f == 1) {
			appout << "Using Actual Images..." << std::endl;
		}
		else {
			appout << "sigma=" << f << "...";
		}
		for (int i = 0; i < (int) numImages; ++i)
		{
			if (f == 1) {//original image size
				scaledImages[i] = new Image<VoxelType>(*images[i]);
				scaledInputHFields[i] = new SCDSAtlasBuilder::VectorFieldType (*inputHFields[i]); 
				if (COMPARE_HFIELD_ON){
					scaledInputCompHFields[i] = new SCDSAtlasBuilder::VectorFieldType (*inputCompHFields[i]);
				}
			}
			else {
				scaledImages[i] = new Image<VoxelType>;
				ImageUtils::gaussianDownsample(*images[i],
						*scaledImages[i],
						Vector3D<int>(f, f, f),
						Vector3D<double>(f, f, f),
						Vector3D<int>(2*f, 2*f, 2*f));


				appout << "InputHFields Downsampling...";

				Vector3D<unsigned int> scaledImageSize = scaledImages[i]->getSize();
				scaledInputHFields[i] = new SCDSAtlasBuilder::VectorFieldType(scaledImageSize);

				HField3DUtils::resample(*inputHFields[i], *scaledInputHFields[i], scaledImageSize);
				if (COMPARE_HFIELD_ON){
				scaledInputCompHFields[i] = new SCDSAtlasBuilder::VectorFieldType(scaledImageSize);
				HField3DUtils::resample(*inputCompHFields[i], *scaledInputCompHFields[i], scaledImageSize);
				}

				if (hinv)
				{//todo:
					std::cerr<<"not implemented"<<std::endl;
				}

			}
		}
		Vector3D<unsigned int> scaledImageSize = scaledImages[0]->getSize();
		spacing = scaledImages[0]->getSpacing();
		appout << "DONE, size = " << scaledImageSize 
			<< ", spacing = " << spacing
			<< std::endl;



		//
		// create h fields for this scale level
		//
		appout << "Creating h-fields...";
		if (scale == 0)
		{
			// start with identity
			for (int i = 0; i < (int) numImages; ++i)
			{
				try 
				{
					h[i] = new Array3D<Vector3D<float> >(scaledImageSize);
					if (hinv)
					{
						hinv[i] = new Array3D<Vector3D<float> >(scaledImageSize);
					}

					// initalize deformation fields to that specified by the
					// affine transforms
					if(initTransformFiles)
					{

						itk::TransformFileReader::Pointer transformreader = itk::TransformFileReader::New();
						transformreader->SetFileName(transformFilenames[i].c_str());
						try 
						{
							transformreader->Update();
							itk::AffineTransform<float, 3>::Pointer transform = dynamic_cast<itk::AffineTransform<float, 3>* >((*transformreader->GetTransformList()->begin()).GetPointer());
							if(!transform)
							{
								std::cerr << "Could not read transform into an itk::AffineTransform<float, 3>" << std::endl;
								return EXIT_FAILURE;
							}

							if(hinv)
								initializeHField(h[i], hinv[i], transform, spacing, images[i]->getOrigin());
							else
								initializeHField(h[i], NULL, transform, spacing, images[i]->getOrigin());

						}
						catch (std::exception& e)      {
							appout << "Error reading transform: " << e.what() << std::endl;
							return EXIT_FAILURE;
						}
					}

					else
					{
						if (initWithInputHField)
						{ 
							//initialize with the input Hfields (predictions)

							*h[i] = *scaledInputHFields[i];
						}
						else{
							HField3DUtils::setToIdentity(*h[i]);
							std::cerr<<"start from identity hfields..."<<std::endl;
						}

						if (hinv)
						{//todo:
							std::cerr<<"not implemented"<<std::endl;

						}
					}


				} catch (std::exception& e)
				{
					appout << "Error creating h-field: " << e.what() << std::endl;
					exit(0);
				}
			}

		}
		else
		{
			appout << "Upsampling...";
			// upsample old xforms
			Array3D<Vector3D<float> > tmph(h[0]->getSize());

			for (int i = 0; i < (int) numImages; ++i)
			{
				tmph = *h[i];
				HField3DUtils::resample(tmph, *h[i], scaledImageSize);

				//smooth hfields
				//--------------------------------------------------------------------
				// put xyz channels of the hfield into three images

				Vector3D<unsigned int> size( scaledImageSize);
				Image<VoxelType> ** hFieldImages = new Image<VoxelType>*[3];
				hFieldImages[0] = new Image<VoxelType>(size);
				hFieldImages[1] = new Image<VoxelType>(size);
				hFieldImages[2] = new Image<VoxelType>(size);
				for (unsigned int z = 0; z < size.z; ++z) {
					for (unsigned int y = 0; y < size.y; ++y) {
						for (unsigned int x = 0; x < size.x; ++x) {
							hFieldImages[0]->set(x,y,z,  h[i]->get(x, y, z).x);
							hFieldImages[1]->set(x,y,z,  h[i]->get(x, y, z).y);
							hFieldImages[2]->set(x,y,z,  h[i]->get(x, y, z).z);
						}
					}
				}

				Image<VoxelType> *	tmpHFieldChanhnelImage = new Image<VoxelType> (*hFieldImages[0]);
				for (int chan=0; chan<3; chan++){
					ImageUtils::gaussianDownsample(*hFieldImages[chan],*tmpHFieldChanhnelImage,
							Vector3D<int>(1,1,1),//factor,original size 
							Vector3D<double>(2, 2, 2),//sigmas
							Vector3D<int>(5, 5, 5));//kernal size

					for (unsigned int z = 0; z < size.z; ++z) {
						for (unsigned int y = 0; y < size.y; ++y) {
							for (unsigned int x = 0; x < size.x; ++x) {
								h[i]->get(x,y,z)[chan]= tmpHFieldChanhnelImage->get(x,y,z);
							}
						}
					}
				}
				delete tmpHFieldChanhnelImage;

				//------------------------------------------------------------
				if (hinv)
				{  //to be implemented
					//tmph = *hinv[i];
					//HField3DUtils::resample(tmph, *hinv[i], scaledImageSize);          
				}
			}

			//use the same method to upsampling the inputH to debug 

		}
		appout << "DONE" << std::endl;

		//
		// create atlas at this scale level
		//
		appout << "Computing atlas at this scale level..." << std::endl;
		Image<VoxelType> iHat(*scaledImages[0]);
		SCDSAtlasBuilder atlasBuilder;

		ArithmeticMeanComputationStrategy<SCDSAtlasBuilder::VoxelType> meanCalculator;
		meanCalculator.SetNumberOfElements(numImages);
		//meanCalculator.SetWeightsToEqual();
		for (unsigned int i = 0; i < numImages; ++i)
		{
			meanCalculator.SetNthWeight(i, imageWeights[i]);
		}
		if( optWeights != NULL){
			atlasBuilder.SetSigmas(optWeights);
		}

		atlasBuilder.SetNumberOfInputImages(numImages);
		atlasBuilder.
			SetMeanComputationStrategy((SCDSAtlasBuilder::
						MeanComputationStrategyType*) 
					&meanCalculator);
		atlasBuilder.SetNumberOfThreads(numberOfThreads);
		atlasBuilder.SetFFTWNumberOfThreads(1);
		atlasBuilder.SetFFTWMeasure(fftwMeasure);
		atlasBuilder.SetLogOutputStream(std::cerr);
		atlasBuilder.
			SetUpdateAverageEverySubIteration(updateAverageAfterSubIterations);
		if (deltaSelectionUseMean)
		{
			atlasBuilder.SetDeltaSelectionToMean();
		}
		else
		{
			atlasBuilder.SetDeltaSelectionToIndividual();
		}
		if (computeInverseHFields)
			atlasBuilder.SetComputeInverseDeformationsOn();
		else
			atlasBuilder.SetComputeInverseDeformationsOff();
		atlasBuilder.SetFluidWarpParameters(fluidParams[scale]);
		atlasBuilder.SetAverageImage(&iHat);

		//the middle results output
		std::ostringstream oss;
		oss << outputDeformedImageFilenamePrefix << scale << "_";  
		atlasBuilder.SetOutputFileNamePrefix(oss.str().c_str());
		if (writeXSlice >=0 )atlasBuilder.SetWriteXSlice(writeXSlice/f);
		if (writeYSlice >=0 )atlasBuilder.SetWriteYSlice(writeYSlice/f);
		if (writeZSlice >=0 )atlasBuilder.SetWriteZSlice(writeZSlice/f);



		if( shapeScores !=NULL)
		{
			atlasBuilder.SetShapeScores(shapeScores);
		}
		if (COMPARE_HFIELD_ON)
				atlasBuilder.SetCompareHFieldOn();

		for (unsigned int imageIndex = 0; imageIndex < numImages; ++imageIndex)
		{
			atlasBuilder.SetNthInputImage(imageIndex, scaledImages[imageIndex]);
			atlasBuilder.SetNthDeformationField(imageIndex, h[imageIndex]);
			atlasBuilder.SetNthInputHField(imageIndex, scaledInputHFields[imageIndex]);
			if (COMPARE_HFIELD_ON){
				atlasBuilder.SetNthInputCompHField(imageIndex, scaledInputCompHFields[imageIndex]);
			}

			if (hinv)
			{   std::cout<<"error, not implemented hinv case" <<std::endl;
				atlasBuilder.SetNthDeformationFieldInverse(imageIndex, hinv[imageIndex]);
			}
		}

		atlasBuilder.GenerateAverage();

		if (outputDeformedImageFilenamePrefix != "")
		{
			for (unsigned int imageIndex = 0; imageIndex < numImages; ++imageIndex)
			{
				std::ostringstream oss;
				oss << outputDeformedImageFilenamePrefix << scale << "_"  
					<< std::setw(4) << std::setfill('0') << imageIndex;
				Image<VoxelType> deformedImage(*atlasBuilder.GetNthDeformedImage(imageIndex));
				deformedImage.setOrigin(scaledImages[imageIndex]->getOrigin());
				deformedImage.setSpacing(scaledImages[imageIndex]->getSpacing());
				ImageUtils::writeMETA(deformedImage, oss.str().c_str());
			}
		}

		appout << "DONE Computing Atlas." << std::endl;

		//
		// write atlas image
		//
		if (outputImageFilenamePrefix != "")
		{
			std::ostringstream oss;
			oss << outputImageFilenamePrefix << scale;
			appout << "Writing Average Image...";
			ImageUtils::writeMETA(iHat, oss.str().c_str());
			appout << "DONE" << std::endl;
		}

		//
		// delete scaled images
		//
		for (int i = 0; i < (int) numImages; ++i)
		{
			delete scaledImages[i];
			delete scaledInputHFields[i];
			 if (COMPARE_HFIELD_ON) {
				delete scaledInputCompHFields[i];}
		}          
	}

	//
	// write h fields at last scale level
	//
	if (outputHFieldFilenamePrefix != "") 
	{
		for (int i = 0; i < (int) numImages; ++i) 
		{
			std::stringstream ss;
			ss << outputHFieldFilenamePrefix << std::setw(4) << std::setfill('0') << i;
			appout << "Writing H Field " << outputHFieldFilenamePrefix << "...";
			HField3DIO::writeMETA(*h[i], origin, spacing, ss.str().c_str());
			appout << "DONE" << std::endl;      
		}
	}

	if (computeInverseHFields && outputHInvFieldFilenamePrefix != "") 
	{
		for (int i = 0; i < (int) numImages; ++i) 
		{
			std::stringstream ss;
			ss << outputHInvFieldFilenamePrefix << std::setw(4) << std::setfill('0') << i;
			appout << "Writing inverse H Field " << outputHInvFieldFilenamePrefix << "...";
			HField3DIO::writeMETA(*hinv[i], origin, spacing, ss.str().c_str());
			appout << "DONE" << std::endl;      
		}
	}

	//
	// clean up memory
	//
	for (int i = 0; i < (int) numImages;++i)
	{
		delete images[i];
		delete h[i];

		delete inputHFields[i];
   	    if (COMPARE_HFIELD_ON)
           delete inputCompHFields[i];
		if (hinv)
		{
			delete hinv[i];
		}
	}

	delete [] scaledImages;
	delete [] inputHFields;
	delete [] scaledInputHFields;
	if (COMPARE_HFIELD_ON)
		delete [] scaledInputCompHFields;
	delete [] images;
	delete [] h;
	if (hinv)
	{
		delete [] hinv;
	}

	appout << "Total Time: " << totalTimer.getTime() << std::endl;
	return 0;
}


