#include <vector>
#include <iostream>
#include <string>
#include <fstream>
//#include <cstdlib>

#include <itkImage.h>
#include <itkImageFileReader.h> 
#include <itkImageRegionIterator.h>
#include <itkImageFileWriter.h>
#include <itkMetaImageIO.h>
#include <itkDefaultDynamicMeshTraits.h>
#include <itkMinimumMaximumImageCalculator.h> 
#include <itkThresholdImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkDiscreteGaussianImageFilter.h> 
#include <itkCastImageFilter.h>
#include <itkConnectedComponentImageFilter.h>

#include <itkMesh.h>
#include <itkBinaryMask3DMeshSource.h>
#include <itkMeshSpatialObject.h>
#include <itkSceneSpatialObject.h>
#include <itkMetaMeshConverter.h>
#include <itkRelabelComponentImageFilter.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>

#include "vtkActor.h"

#include <vtkCellArray.h>
#include <vtkImageData.h>
#include <vtkImageImport.h>
#include <vtkImageGaussianSmooth.h>
#include <vtkImageThreshold.h>
#include <vtkImageToStructuredPoints.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkStripper.h>
#include <vtkDecimatePro.h>
#include <vtkQuadricDecimation.h>

#include "vtkPolyDataMapper.h"
#include "vtkPolyDataReader.h"
#include "vtkTriangleFilter.h"
#include "vtkActor.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkProperty.h"
#include "vtkCamera.h"
#include "vtkLight.h"
#include "vtkImageMarchingCubes.h"
#include "itkImageToVTKImageFilter.h"
#include "vtkPolyDataWriter.h"
#include "vtkSTLWriter.h"
#include "vtkMassProperties.h"

#include "argio.hh"
#include "FL/Fl_Box.H"
#include "FL/Fl_File_Chooser.H"


//using namespace itk ;
//using namespace std ;

#define Dimension 3 
#define NO_LABEL -1
#define PI 3.141592653589793238462643

class vol2surf : public Fl_Box
{
public:
  void hasThreshold ( bool state );
  void setThreshold (const char *value);
  void setUpperThreshold (const char *value);
  void setThreshold (double value);
  void setUpperThreshold (double value);
  
  vol2surf(int x, int y, int w, int h, const char * = 0) ;
  virtual ~vol2surf();
    
  void setDecimate (bool state);
  void setFilename ( const char *name ) ;
  void setOutputFormat ( bool so, bool vtk, bool byu, bool iv, bool stl ) ;
  void hasGaussian ( bool state ) ;
  void setVariance ( const char *value ) ;
  void setVariance ( double value ) ;
  void hasLabel ( bool state ) ;
  void setLabel ( const char *value ) ;
  void setLabel ( int value ) ;
  void hasConnected (bool state);
  void setConnectedThreshold (const char *value);
  double getVolume();
  double getSurfaceArea();
  double getSphereRadiiRatio();
  char * getVolumeStr();
  char * getSurfaceAreaStr();
  char * getSphereRadiiRatioStr();
  void Run();

  int CommandLine ( int argc, const char *argv[] ) ;
private:
  void readMeshBack();
  void writeSO();
  void VTKSetup();
  void readSOBack();
  void writeIV();
  void writeBYU();
  void writePLY();
  void writeSTL();
  void writeMesh();
  void generateMesh();
  bool gaussian, label, manThreshold, connected ;
  bool outputSO, outputMesh, outputBYU, outputIV, outputSTL ;

  double manThresholdValue ;
  double manUpperThresholdValue ;
  double variance ;
  int connectedThresholdValue;
  int labelValue ;

  bool decimationOn;
  bool quadricDecimationOn;
  int  targetReduction;

  double imageVolume;
  double surfaceArea;
  double sphereRadiiRatio;

  std::string filename;

  // added by joohwi
  // This is to provide specific output filename.
  // The output file's type is determined by a filename from '-outfile'
  // PLY is currently not supported
  enum FileType { UNKNOWN, SO, VTK, BYU, IV, PLY, STL};

  FileType outFileType; 
  std::string outFileName;

  std::vector < double > verts ;
  std::vector < int > tris ;
  double boundingBox[6] ;

  // types used in this routine
  typedef unsigned short ImagePixelType ;
  typedef float SmoothImagePixelType;
  typedef itk::Image < ImagePixelType, Dimension >  ImageType ;
  typedef itk::Image < SmoothImagePixelType, Dimension > SmoothImageType;
  typedef ImageType::SizeType ImageSizeType ;
  typedef itk::ImageFileReader < ImageType > VolumeReaderType ;
  
  typedef itk::MinimumMaximumImageCalculator< ImageType > minMaxCalcType;

  typedef itk::ThresholdImageFilter < ImageType > threshFilterType;
  typedef itk::BinaryThresholdImageFilter < ImageType, ImageType > binThreshFilterType;
  typedef itk::DiscreteGaussianImageFilter < SmoothImageType, SmoothImageType > gaussFilterType;
  typedef itk::CastImageFilter < ImageType, SmoothImageType > castInputFilterType;
  typedef itk::CastImageFilter < SmoothImageType, ImageType > castOutputFilterType;
  typedef itk::ConnectedComponentImageFilter < ImageType, ImageType > conCompFilterType;
  typedef itk::Mesh < float ,3, itk::DefaultDynamicMeshTraits < float,3,3 > > MeshType ;
  typedef itk::BinaryMask3DMeshSource < ImageType, MeshType >  MeshSourceType ;
  typedef itk::MeshSpatialObject < MeshType > MeshObjectType ;
  typedef itk::MetaMeshConverter < Dimension, float, itk::DefaultDynamicMeshTraits < float,3,3 > > MeshConverterType ;
  typedef itk::RelabelComponentImageFilter < ImageType, ImageType > relabelFilterType;

  typedef itk::ImageRegionConstIterator< ImageType > ConstIteratorType;
  typedef itk::ImageRegionIterator< ImageType> IteratorType;
  typedef itk::ImageBase<Dimension>::SpacingType SpacingType;

  typedef itk::SceneSpatialObject<> SceneType ;


  SceneType::Pointer mScene ;
  MeshObjectType::Pointer mMeshObject ;  

  bool commandLine ;  
};

