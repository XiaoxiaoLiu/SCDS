#include "Array3D.h"
#include "Array3DIO.h"
#include "Array3DUtils.h"
#include "HField3DUtils.h"
#include <iostream>
#include <string>
#include "Image.h"
#include "Surface.h"
#include "SurfaceUtils.h"
#include "HField3DUtils.h"
#include "HField3DIO.h"
#include "Anastruct.h"
#include "AnastructUtils.h"
#include "ImageUtils.h"
#include "math.h"

#include "ApplicationUtils.h"

typedef float VoxelType;

void printUsage()
{
  std::cerr << "Usage: txApply" << std::endl << std::endl
            << "Specify a transformation: " << std::endl
	    << "\t-h hFieldFileName          ---deformation field" << std::endl
	    << "\t-t tx ty tz                ---translation" << std::endl
            << std::endl
            << "Specify the direction to apply the transformation: " << std::endl
	    << "\t-f                         ---forwards (push, with arrows)"
            << std::endl
	    << "\t-b                         ---backwards (pull, against arrows)"
            << std::endl << std::endl
            << "Specify the object to transform: " << std::endl
	    << "\t-i imageFileName           ---image to be deformed"
            << std::endl
	    << "\t-s surfaceFileName         ---surface to be deformed"
            << std::endl
	    << "\t-a anastructFileName       ---anastruct to be deformed" 
            << std::endl << std::endl
            << "Specify the filename of the transformed object: " << std::endl
	    << "\t-o outFilePrefix                 ---filename prefix of deformed object" 
            << std::endl << std::endl
            << "Options: " << std::endl
	    << "\t-n ---use nearest neighbor interpolation; trilerp is default"
            << std::endl
	    << "\t-d srcOrigin srcSpacing dstOrigin dstSpacing" << std::endl
            << "           ---applies to vertices of surfaces and anastructs only" << std::endl
            << "             -transform must be applied to vertices in hfield index coordinates" << std::endl
            << "             -used to apply transform to a surface stored in world coordinates" << std::endl
            << "             -before transform, v=(v-srcOrigin)/srcSpacing" << std::endl
            << "             -after transform, v=v*dstSpacing+dstOrigin" << std::endl
            << "             -example: -d 0 0 0 1.5 1.5 1.5 0 0 0 1.5 1.5 1.5" << std::endl
	    << "\t-p alpha" << std::endl
            << "           ---for use with deformation fields only" << std::endl
            << "             -allows linear interpolation between the identity transform" << std::endl
            << "              and the given deformation field" << std::endl
            << "             -the interpolated deformation is applied to the object" << std::endl
            << "             -specify alpha between 0 (identity) and 1 (the full deformation)" << std::endl

            << std::endl;
}

int main(int argc, char **argv)
{
  //
  // parse first part of command line
  //
  std::string outputFilename = "";
  std::string hFieldFileName = "";
  std::string imageFileName  = "";
  std::string surfaceFileName = "";
  std::string anastructFileName = "";
  Vector3D<double> dstOrigin(0,0,0);  
  Vector3D<double> dstSpacing(1,1,1);
  Vector3D<double> srcOrigin(0,0,0);  
  Vector3D<double> srcSpacing(1,1,1);
  bool doDeformationInterpolation = false;
  double alpha; // deformation field interpolation parameter
  bool forward = false;
  bool backward = false;
  bool nearestNeighbor = false;
  bool doTranslation = false;
  Vector3D<double> translation;

  argv++; argc--;
  while (argc > 0)
    {
      if (std::string(argv[0]) == std::string("-o"))
	{
	  // outfile: make sure at least one argument follows
	  if (argc < 2)
	    {
	      printUsage();
	      return 0;
	    }
	  outputFilename = argv[1];
	  argv++; argc--;
	  argv++; argc--;
	}
      else if (std::string(argv[0]) == std::string("-i"))
	{
	  // image: make sure at least one argument follows
	  if (argc < 2)
	    {
	      printUsage();
	      return 0;
	    }
	  imageFileName = argv[1];
	  argv++; argc--;
	  argv++; argc--;
	}
      else if (std::string(argv[0]) == std::string("-d"))
	{
	  // src and dest origing and dimensions: make sure at least
	  // one argument follows
	  if (argc < 13)
	    {
	      printUsage();
	      return 0;
	    }
          srcOrigin.x = atof(argv[1]);
          srcOrigin.y = atof(argv[2]);
          srcOrigin.z = atof(argv[3]);
          srcSpacing.x = atof(argv[4]);
          srcSpacing.y = atof(argv[5]);
          srcSpacing.z = atof(argv[6]);

          dstOrigin.x = atof(argv[7]);
          dstOrigin.y = atof(argv[8]);
          dstOrigin.z = atof(argv[9]);
          dstSpacing.x = atof(argv[10]);
          dstSpacing.y = atof(argv[11]);
          dstSpacing.z = atof(argv[12]);

	  for (int i = 0; i < 13; ++i) 
            {
              argv++; argc--;
            }
	}
      else if (std::string(argv[0]) == std::string("-s"))
	{
	  // surface: make sure at least one argument follows
	  if (argc < 2)
	    {
	      printUsage();
	      return 0;
	    }
	  surfaceFileName = argv[1];

	  for (int i = 0; i < 2; ++i) 
            {
              argv++; argc--;
            }
	}
      else if (std::string(argv[0]) == std::string("-a"))
	{
	  // anastruct: make sure at least one argument follows
	  if (argc < 2)
	    {
	      printUsage();
	      return 0;
	    }
	  anastructFileName = argv[1];

	  for (int i = 0; i < 2; ++i) 
            {
              argv++; argc--;
            }
	}
      else if (std::string(argv[0]) == std::string("-h"))
	{
	  // hField: make sure at least one argument follows
	  if (argc < 2)
	    {
	      printUsage();
	      return 0;
	    }
	  hFieldFileName = argv[1];
	  argv++; argc--;
	  argv++; argc--;
	}
      else if (std::string(argv[0]) == std::string("-t"))
	{
	  // translation: make sure at least one argument follows
	  if (argc < 4)
	    {
	      printUsage();
	      return 0;
	    }
	  translation.x = atof(argv[1]);
	  translation.y = atof(argv[2]);
	  translation.z = atof(argv[3]);
          doTranslation = true;
	  for (int i = 0; i < 4; ++i) 
            {
              argv++; argc--;
            }
	}
      else if (std::string(argv[0]) == std::string("-n"))
	{
          // nearest neighbor
	  nearestNeighbor = true;
	  argv++; argc--;	  
	}
      else if (std::string(argv[0]) == std::string("-f"))
	{
          // forward
	  forward = true;
          backward = false;
	  argv++; argc--;	  
	}
      else if (std::string(argv[0]) == std::string("-b"))
	{
          // backward
	  forward = false;
	  backward = true;
	  argv++; argc--;	  
	}
      else if (std::string(argv[0]) == std::string("-p"))
	{
          // deformation field interpolation
	  if (argc < 2)
	    {
	      printUsage();
	      return 0;
	    }
          doDeformationInterpolation = true;
          alpha = atof(argv[1]);
	  argv++; argc--;	  
	  argv++; argc--;	  
	}
      else
	{
	  std::cerr << "Invalid argument: " << argv[0] << std::endl;
	  printUsage();
	  return 0;
	}
    }
  
  if ((surfaceFileName == "" && 
       imageFileName == "" && 
       anastructFileName == "") || 
      (!doTranslation &&
       hFieldFileName == "") || 
      outputFilename == "")
    {
      printUsage();
      return 0;
    }

  if ((forward && backward) || !(forward || backward)) {
    std::cerr << "Must specify excatly one of -f and -b." << std::endl;
    printUsage();
    return 0;
  }

  //
  // print parameters
  //
  if (doTranslation)
    {
      std::cerr << "Translation : " << translation << std::endl;
    }
  else
    {
      std::cerr << "hField Filename  : " << hFieldFileName << std::endl;
    }
  if (imageFileName != "")
    {
      std::cerr << "Image Filename   : " << imageFileName << std::endl;
    }
  if (surfaceFileName != "")
    {
      std::cerr << "Surface Filename   : " << surfaceFileName << std::endl;
    }
  if (anastructFileName != "")
    {
      std::cerr << "Anastruct Filename   : " << anastructFileName << std::endl;
    }
  std::cerr << "Output Filename  : " << outputFilename << std::endl;

  //
  // load hField
  //
  Array3D<Vector3D<float> > hField;
  Vector3D<unsigned int> hSize;
  Vector3D<double> hOrigin, hSpacing;
  if (!doTranslation)
    {
      std::cerr << "Loading hField: " << hFieldFileName << "...";
      HField3DIO::readMETAHeader(hSize,hOrigin,hSpacing,
                                 hFieldFileName.c_str());
      HField3DIO::readMETA(hField, hFieldFileName.c_str());
      std::cerr << "DONE" << std::endl;
      std::cerr << "  Dimensions: " << hSize << std::endl;
      std::cerr << "  Origin    : " << hOrigin << std::endl;
      std::cerr << "  Spacing   : " << hSpacing << std::endl;

      if (doDeformationInterpolation) {
        std::cerr << "Interpolating hField: alpha == " << alpha << "...";    
        Array3D<Vector3D<float> > eye(hField.getSize());
        HField3DUtils::setToIdentity(eye);
        eye.scale(1.0-alpha);
        hField.scale(alpha);
        Array3DUtils::sum(hField,eye);
        std::cerr << "DONE" << std::endl;
      }
    }

  if (surfaceFileName != "")
    {
      /*
       //
      // load surface
      //
      Surface surface;
      try 
        {
          surface.readBYU(surfaceFileName.c_str());
        }
      catch (...)
        {
          std::cerr << "Can't load surface: " << surfaceFileName << std::endl;
          exit(0);
        }

      //
      // deform surface
      //
      if (forward)
        {
          if (doTranslation)
            {
              surface.translate(translation);
            }
          else
            {
              // subtract origin and divide by spacing
              SurfaceUtils::worldToImageIndex(surface, srcOrigin, srcSpacing);
              // assumption is that vertices are in image index coords
              HField3DUtils::apply(surface, hField);
              // multiply by spacing and add origin
              SurfaceUtils::imageIndexToWorld(surface, dstOrigin, dstSpacing);
            }
        }
      else
        {
          if (doTranslation)
            {
              surface.translate(-translation);
            }
          else
            {
              SurfaceUtils::worldToImageIndex(surface, srcOrigin, srcSpacing);
              // assumption is that vertices are in image index coords
              HField3DUtils::inverseApply(surface, hField);
              SurfaceUtils::imageIndexToWorld(surface, dstOrigin, dstSpacing);
            }
        }

      //
      // write out surface
      //
      try
        {
          surface.writeBYU(outputFilename.c_str());      
        }
      catch (...)
        {
          std::cerr << "Can't write surface: " << outputFilename << std::endl;
          exit(0);
        }      
    */}
  else if (anastructFileName != "")
    {/*
      //
      // load anastruct
      //
      Surface surface;
      Anastruct ana;
      try 
        {
          AnastructUtils::readPLUNCAnastruct(ana, anastructFileName.c_str());
        }
      catch (...)
        {
          std::cerr << "Can't load anastruct: " << anastructFileName 
                    << std::endl;
          exit(0);
        }

      //
      // create surface from anastruct
      //
      AnastructUtils::anastructToSurface(ana, surface);

      //
      // deform surface
      //
      if (forward)
        {
          if (doTranslation)
            {
              surface.translate(translation);
            }
          else
            {
              SurfaceUtils::worldToImageIndex(surface, srcOrigin, srcSpacing);
              HField3DUtils::apply(surface, hField);
              SurfaceUtils::imageIndexToWorld(surface, dstOrigin, dstSpacing);
            }
        }
      else
        {
          if (doTranslation)
            {
              surface.translate(-translation);
            }
          else
            {
              SurfaceUtils::worldToImageIndex(surface, srcOrigin, srcSpacing);
              HField3DUtils::inverseApply(surface, hField);
              SurfaceUtils::imageIndexToWorld(surface, dstOrigin, dstSpacing);
            }
        }

      //
      // write out surface
      //
      try
        {
          surface.writeBYU(outputFilename.c_str());      
        }
      catch (...)
        {
          std::cerr << "Can't write surface: " << outputFilename << std::endl;
          exit(0);
        }      
   */ }
  else
    {
      //
      // load input image
      //
      Image<float> image;
      std::cerr << "Loading image: " << imageFileName << "...";
      ApplicationUtils::LoadImageITK(imageFileName.c_str(), image);
      std::cerr << "DONE" << std::endl;
      std::cerr << "  Dimensions: " << image.getSize() << std::endl;
      std::cerr << "  Origin    : " << image.getOrigin() << std::endl;
      std::cerr << "  Spacing   : " << image.getSpacing() << std::endl;
      
      //
      // create deformed image and hfield
      //
      Image<VoxelType> def(image);
      
      if (doTranslation)
        {
          def = image;
          if (forward)
            {
              ImageUtils::translate(def, translation);
            }
          else
            {
              ImageUtils::translate(def, -translation);
            }
        }
      else
        {
          //
          // compute starting voxel of roi
          // NB: assume that spacing hfield == spacing image
          //
          if (hSpacing != image.getSpacing()) {
            std::cerr << "WARNING: hField and image spacing do not agree!"
                      << std::endl
                      << "       : this operation is not currently supported"
                      << std::endl;
          }
          Vector3D<double> iOrigin  = image.getOrigin();
          Vector3D<double> iSpacing = image.getSpacing();
          Vector3D<double> roiStart = (hOrigin-iOrigin) / iSpacing;
          Vector3D<int> roi((int)floor(roiStart.x+0.5),
                            (int)floor(roiStart.y+0.5),
                            (int)floor(roiStart.z+0.5));
          std::cerr << "ROI: hField (0,0,0) corresponds to image " << roi 
                    << std::endl 
                    << "     [rounded from " << roiStart << "]" << std::endl;
          std::cerr << "Deforming image...";
          if (forward)
            {
              HField3DUtils::forwardApply(image, hField, def, 
                                          roi.x, roi.y, roi.z, 0.0F);
            }
          else
            {
              HField3DUtils::apply(image, hField, def, roi.x, roi.y, roi.z,
                                   0.0F, nearestNeighbor);
            }
          std::cerr << "DONE" << std::endl;
        }

      //
      // write deformed image
      //
      if (outputFilename != "")
        {
          std::cerr << "Writing Deformed Image...";
          ImageUtils::writeMETA(def, outputFilename.c_str());
          std::cerr << "DONE" << std::endl;
        }
    }      
  return 0;
}
  
  
  
