#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#pragma warning ( disable : 4284 )
#pragma warning ( disable : 4018 )
#endif







#include "vol2surf.h"



vol2surf::vol2surf (int x, int y, int w, int h, const char *l) : Fl_Box (x, y, w, h, l)
{
  gaussian = false ;
  label = false ;
  manThreshold = false ;
  outputSO = true ;
  outputMesh = outputBYU = outputIV = outputSTL = false ;
  connected = false;
  decimationOn = false;
  quadricDecimationOn = false;
  targetReduction = 95;
  imageVolume = 0.0;
  surfaceArea = 0.0;
  sphereRadiiRatio = 0.0;
  commandLine = false ;
}

vol2surf::~vol2surf()
{
}


void vol2surf::setVariance ( double value) 
{
  this->variance = value ;
}

void vol2surf::setThreshold(double value)
{
  manThresholdValue = value ;
}
void vol2surf::setUpperThreshold(double value)
{
  manUpperThresholdValue = value ;
}

void vol2surf::setFilename (const char *name) 
{
  filename = name ;
}

void vol2surf::setOutputFormat ( bool so, bool vtk, bool byu, bool iv, bool stl ) 
{
  outputSO = so ;
  outputMesh = vtk ;
  outputBYU = byu ;
  outputIV = iv ;
  outputSTL = stl;
}


void vol2surf::hasGaussian ( bool state ) 
{
  gaussian = state ;
}


void vol2surf::setVariance (const char *value) 
{
  variance = atof ( value ) ;
}

void vol2surf::setDecimate (bool state) 
{
  decimationOn = state;
}

void vol2surf::hasLabel ( bool state ) 
{
  label = state ;
}


void vol2surf::setLabel (const char *value) 
{
  labelValue = atoi ( value ) ;
}

void vol2surf::setLabel (int value) 
{
  labelValue = value ;
}

void vol2surf::hasConnected (bool state)
{
  connected = state;
}

void vol2surf::setConnectedThreshold (const char* value)
{
  connectedThresholdValue = atoi ( value );
}

double vol2surf::getVolume()
{
  return imageVolume;
}

char * vol2surf::getVolumeStr()
{
  char * string = new char [30];
  sprintf(string, "%f", imageVolume);
  return string;
}

double vol2surf::getSurfaceArea()
{
  return surfaceArea;
}

char * vol2surf::getSurfaceAreaStr()
{
  char * string = new char [30];
  sprintf(string, "%f", surfaceArea);
  return string;
}

double vol2surf::getSphereRadiiRatio()
{
  return sphereRadiiRatio;
}

char * vol2surf::getSphereRadiiRatioStr()
{
  char * string = new char [30];
  sprintf(string, "%f", sphereRadiiRatio);
  return string;
}

void vol2surf::Run()
{
  //Convert the volume into the surface. Then write the vtkpolydata file
  generateMesh () ;
 std::cout<<"done gerenating mesh"<<std::endl;
  //Read the vtk file written in generateMesh, and fill up the triangle structures
  readMeshBack () ;
  std::cout<<"ready to write"<<std::endl;
  //Write the Meta file
  writeSO () ;
  
  MeshConverterType * meshConverter = new MeshConverterType () ;
  
  int i ;
  std::string metaFileName(filename) ;
  i = metaFileName.find_last_of ('.') ;
  metaFileName.replace ( i+1, 4, "meta" ) ;
  mMeshObject = meshConverter->ReadMeta ( metaFileName.c_str() ) ;
  
  if ( !this->commandLine )
    this->VTKSetup() ;

 	if (outFileType == UNKNOWN) { 
		if ( !outputBYU && !outputIV && !outputSTL)
			return ;

		if ( outputBYU ) 
		{
		  writeBYU () ;
		}

		if ( outputIV ) 
		{
		  writeIV () ;
		}
		if ( outputSTL ) 
		{
		  writeSTL () ;
		}
	} else {
		switch (outFileType) {
		case BYU:
		        writeBYU();
			break;
		case IV:
			writeIV();
			break;
		case PLY:
			writePLY();
			break;
		case STL:
		        writeSTL();
			break;
		default:
			break;
		}
	}
}

void vol2surf::generateMesh()
{
  // load the image
  ImageType::Pointer image ;
  VolumeReaderType::Pointer imageReader = VolumeReaderType::New() ;
  try
  {
    imageReader->SetFileName(filename.c_str()) ;
    imageReader->Update() ;
    image = imageReader->GetOutput() ;
  }
  catch (itk::ExceptionObject e)
  {
    e.Print(cout) ;
    exit(0) ;
  }
  
  // Make sure background is 0 and object is 255
  const int LABEL_VAL = 255;
  ImageType::Pointer procImage ;
  binThreshFilterType::Pointer labelFilter;
  binThreshFilterType::Pointer threshFilter;
  binThreshFilterType::Pointer threshFilter2;
  binThreshFilterType::Pointer binTreshFilter;
  conCompFilterType::Pointer conCompFilter;
  relabelFilterType::Pointer relabelFilter;

  SmoothImageType::Pointer smoothImage ;
  
  if (label) 
  {
    // Thresholding at label -> set to 0 otherwise
    labelFilter = binThreshFilterType::New();
    labelFilter->SetInput(image);
    labelFilter->SetLowerThreshold(labelValue);
    labelFilter->SetUpperThreshold(labelValue);
    //labelFilter->ThresholdOutside(labelValue, labelValue);
    labelFilter->SetOutsideValue (0);
    labelFilter->SetInsideValue (LABEL_VAL);
    labelFilter->ReleaseDataFlagOn();
    labelFilter->Update();
    image = labelFilter->GetOutput();
    image->DisconnectPipeline();
  }
  
  // Set object to 255 -> Tresholding all but label
  threshFilter = binThreshFilterType::New();
  threshFilter->SetInput(image);
    
  minMaxCalcType::Pointer minMaxCalc = minMaxCalcType::New();
  minMaxCalc->SetImage(image);
  minMaxCalc->Compute();
  cout << "Min: " << minMaxCalc->GetMinimum() << ";   Max: " << minMaxCalc->GetMaximum() << endl ;
  double max_val = minMaxCalc->GetMaximum();
  double threshold = minMaxCalc->GetMaximum() / 2.0 ; 
  if ( this->manThresholdValue == -1 ) 
  {
    this->manThresholdValue = 0 ;
  }
  if ( this->manUpperThresholdValue == -1 ) 
  {
    this->manUpperThresholdValue = max_val ;
  }

  if ( !manThreshold )
    {
      //threshFilter->ThresholdAbove(threshold);
      threshFilter->SetLowerThreshold((ImagePixelType) threshold);
      threshFilter->SetUpperThreshold((ImagePixelType) max_val);
    }
  else
    {
      cout << manThresholdValue << " ... " << manUpperThresholdValue << endl ;
      //threshFilter->ThresholdAbove(manThresholdValue);
      threshFilter->SetLowerThreshold((ImagePixelType) manThresholdValue);
      threshFilter->SetUpperThreshold((ImagePixelType) manUpperThresholdValue);
    }
  threshFilter->SetOutsideValue (0);
  threshFilter->SetInsideValue(LABEL_VAL);
  threshFilter->Update();
  procImage = threshFilter->GetOutput();
  
  if (connected)
  {
    conCompFilter = conCompFilterType::New();
    conCompFilter->SetInput(procImage); //Connected components computed, in no particular order

    relabelFilter = relabelFilterType::New(); //relabelFilter sorts connected components
    relabelFilter->SetInput(conCompFilter->GetOutput());
    relabelFilter->SetMinimumObjectSize(connectedThresholdValue); //components smaller than our threshold size are discarded

    threshFilter2 = binThreshFilterType::New();
    threshFilter2->SetInput(relabelFilter->GetOutput());
    threshFilter2->SetLowerThreshold(0); //flatten the image, restoring all remaining objects to value 255
    threshFilter2->SetUpperThreshold(0);
    threshFilter2->SetOutsideValue (LABEL_VAL);
    threshFilter2->SetInsideValue(0);
    threshFilter2->Update();
    procImage = threshFilter2->GetOutput();
  }

  //obtain volume
  ImageType::RegionType inputRegion;
  ImageType::RegionType::IndexType inputStart;
  inputStart[0] = 0;
  inputStart[1] = 0;
  inputStart[2] = 0;
  inputRegion.SetSize( procImage->GetLargestPossibleRegion().GetSize() );
  inputRegion.SetIndex( inputStart );
  ConstIteratorType inputIt( procImage, inputRegion );
  SpacingType spacing = imageReader->GetOutput()->GetSpacing();
  for ( inputIt.GoToBegin(); !inputIt.IsAtEnd(); ++inputIt)
  {
    if(inputIt.Get() != 0) imageVolume += spacing[0]*spacing[1]*spacing[2];      
  }
    
  castInputFilterType::Pointer convInFilter = castInputFilterType::New();
  convInFilter->SetInput(procImage);
  convInFilter->Update () ;
  smoothImage = convInFilter->GetOutput () ;
  
  // gaussian smoothing, if requested
  if (gaussian) 
  {
    
    gaussFilterType::Pointer smoothFilter = gaussFilterType::New();  
    smoothFilter->SetInput(smoothImage);
    smoothFilter->SetVariance(variance);
    smoothFilter->Update () ;
    smoothImage = smoothFilter->GetOutput () ;
  }
  
  
  // go to vtk 
  typedef itk::ImageToVTKImageFilter<SmoothImageType> Itk2VtkType;
  Itk2VtkType::Pointer  m_Itk2Vtk = Itk2VtkType::New();
  
  m_Itk2Vtk->SetInput(smoothImage);  
  m_Itk2Vtk->Update();
  
  // create mesh using marching cubes
  vtkImageMarchingCubes *marcher = vtkImageMarchingCubes::New();
  marcher->SetInput(m_Itk2Vtk->GetOutput());
  marcher->SetValue(0, 128);
  marcher->ComputeScalarsOff();
  marcher->ComputeGradientsOff();
  marcher->Update() ;  
  vtkPolyData *pipePolyTail = marcher->GetOutput();

  //TriangleDecimation
  if (decimationOn) {  
		if (!quadricDecimationOn) {
			vtkDecimatePro * DecimateFilter = vtkDecimatePro::New();
			DecimateFilter->SetInput(pipePolyTail);
			DecimateFilter->SetTargetReduction(targetReduction);  
			DecimateFilter->PreserveTopologyOn();
			DecimateFilter->Update();
			pipePolyTail = DecimateFilter->GetOutput();
		} else {
			vtkQuadricDecimation* decimateFilter = vtkQuadricDecimation::New();
			decimateFilter->SetInput(pipePolyTail);
			decimateFilter->SetTargetReduction(targetReduction / 100.0);	
			decimateFilter->Update();
			pipePolyTail = decimateFilter->GetOutput();
		}
  }

  // Get surface area of the mesh
  vtkMassProperties *massprops = vtkMassProperties::New();
  massprops->SetInput( pipePolyTail );
  surfaceArea = massprops->GetSurfaceArea();

  //Calculate Spherical Radii Ratio
  double radSurface = sqrt(surfaceArea / (4.0 * PI) );
  double radVolume = pow(imageVolume / (4.0 * PI / 3.0), 1.0/3.0);
  sphereRadiiRatio = radSurface / radVolume;
  


  // write out the vtk mesh
  std::string outfileName(filename) ;
  int i ;
  i = outfileName.find_last_of ('.') ;
  outfileName.replace ( i+1, 4, "vtk" ) ;
  
  vtkPolyDataWriter *writer = vtkPolyDataWriter::New () ;
  writer->SetInput ( pipePolyTail) ;
  writer->SetFileName ( outfileName.c_str() ) ;
  writer->Update () ;
}

void vol2surf::writeMesh()
{
  ofstream outfile ;
  int nPts, nTris ;
  double normX, normY, normZ ;
  std::string outfileName(filename) ;
  int i ;
  
  i = outfileName.find_last_of ('.') ;
  outfileName.replace ( i+1, 4, "vtk" ) ;
  outfile.open ( outfileName.c_str() ) ;

  nPts = (verts.size () + 1) / 3 ;
  nTris = (tris.size () + 1) / 3 ;

  normX = boundingBox[1] - boundingBox[0] ;
  normY = boundingBox[3] - boundingBox[2] ;
  normZ = boundingBox[5] - boundingBox[4] ;

  // output header
  outfile << "# vtk DataFile Version 1.0" << endl << "vtk output" << endl << "ASCII" << endl << "DATASET POLYDATA" ;
  outfile << endl << "POINTS " << nPts << " float" << endl ;

  for ( i = 0 ; i < nPts ; i++ )
  {
    outfile << (verts[3*i]-boundingBox[0])/normX << " " << (verts[3*i+1]-boundingBox[2])/normY << " " << (verts[3*i+2]-boundingBox[4])/normZ << endl ;
  }
  
  outfile<< endl << "POLYGONS " << nTris << " " << 4*nTris << endl ;
  for ( i = 0 ; i < nTris ; i++)
  {
    outfile << "3 " << tris[3*i] << " " << tris[3*i+1] << " " << tris[3*i+2] << endl ;
  }

  outfile.close () ;
  
}

void vol2surf::writeBYU()
{
  ofstream outfile ;
  int nPts, nTris ;
  std::string outfileName(filename) ;
  int i ;
  
  if (outFileType == UNKNOWN) {
		i = outfileName.find_last_of ('.') ;
		outfileName.replace ( i+1, 4, "byu" ) ;
	} else {
		outfileName = outFileName;
  }
  outfile.open ( outfileName.c_str() ) ;

  nPts = (verts.size () + 1) / 3 ;
  nTris = (tris.size () + 1) / 3 ;

  // output header
  outfile << "1 " << nPts << " " << nTris << " " << 3*nTris << " 0" << endl ;
  outfile << "1 " << nTris << endl ;

  for ( i = 0 ; i < nPts/2 ; i++ )
  {
    outfile << verts[6*i] << " " << verts[6*i+1] << " " << verts[6*i+2] << " " ;
    outfile << verts[6*i+3] << " " << verts[6*i+4] << " " << verts[6*i+5] << endl ;
  }
  
  for ( i = 0 ; i < nTris ; i++)
  {
    outfile << tris[3*i] << " " << tris[3*i+1] << " -" << tris[3*i+2] << endl ;
  }

  outfile.close () ;
}

void vol2surf::writeSTL()
{
  
  std::string outfileName(filename) ;
  int i ;
  
  //Get the vtk file created earlier
  i = outfileName.find_last_of ('.') ;
  outfileName.replace ( i+1, 4, "vtk" ) ;
  
  // read in the vtk polydata file
  vtkPolyDataReader *meshReader = vtkPolyDataReader::New() ;
  meshReader->SetFileName(outfileName.c_str());
  vtkPolyData *mesh = meshReader->GetOutput();

  //Create the stl file name
  i = outfileName.find_last_of ('.') ;
  outfileName.replace ( i+1, 4, "stl" ) ;

  //Create the triangle mesh
  vtkTriangleFilter *tri = vtkTriangleFilter::New();
  vtkSTLWriter *writer = vtkSTLWriter::New();
  tri->SetInput(mesh);

  //Write the stl file
  writer->SetInput(tri->GetOutput());
  writer->SetFileName(outfileName.c_str());
  writer->Update();
  writer->Delete();
  tri->Delete();
}

void vol2surf::writeIV()
{
  ofstream outfile ;
  int nPts, nTris ;
  std::string outfileName(filename) ;
  int i ;

  if (outFileType == UNKNOWN) {
		i = outfileName.find_last_of ('.') ;
    outfileName.replace ( i+1, 4, "iv" ) ;
	} else {
		outfileName = outFileName;
  }
  outfile.open ( outfileName.c_str() ) ;

  nPts = (verts.size () + 1) / 3 ;
  nTris = (tris.size () + 1) / 3 ;

  // output header
  outfile << "#Inventor V2.1 ascii" << endl ;
  outfile << "Separator { " << endl << "Coordinate3 { " << endl << "point [ " << endl ;

  for ( i = 0 ; i < nPts ; i++ )
  {
    outfile << verts[3*i] << " " << verts[3*i+1] << " " << verts[3*i+2] << " ," << endl ;
  }
  
  outfile << "]" << endl << "}" << endl << "IndexedTriangleStripSet { " << endl << "coordIndex[" << endl ;

  for ( i = 0 ; i < nTris ; i++)
  {
    outfile << tris[3*i] << "," << tris[3*i+1] << "," << tris[3*i+2] << ",-1," << endl ;
  }

  outfile << "]" << endl << "}       }  " << endl ;

  outfile.close () ;
  
}

void vol2surf::writePLY()
{
	std::cout << "PLY filetype is currently not supported." << endl;
}

void vol2surf::readSOBack()
{
  // reads in the ITK mesh file we just generated
  ifstream infile ;
  std::string metaFileName (filename);
  char line[200] ;
  double x, y, z ;
  int index, v1, v2, v3 ;
  int nPts, nTris ;
  int i ;
  double minX, maxX, minY, maxY, minZ, maxZ ;

  i = metaFileName.find_last_of ('.') ;
  metaFileName.replace ( i+1, 4, "meta" ) ;

  minX = minY = minZ = 99999 ;
  maxX = maxY = maxZ = -99999 ;


  infile.open ( metaFileName.c_str() ) ;
  
  verts.clear () ;
  tris.clear () ;
  // skip header - all we need here is the number of vertices
  for (int i = 0 ; i < 13 ; i++)
    infile.getline (line, 200 ) ;
  nPts = atoi(strrchr(line,' ')) ;
  infile.getline ( line,200 )  ;

  for ( i = 0 ; i < nPts ; i++ )
  {
    infile >> index ;
    infile >> x >> y >> z ;
    assert ( index == i ) ;
    verts.push_back ( x ) ;
    verts.push_back ( y ) ;
    verts.push_back ( z ) ;

    if ( x > maxX )
      maxX = x ;
    else if ( x < minX )
      minX = x ;
    if ( y > maxY )
      maxY = y ;
    else if ( y < minY )
      minY = y ;
    if ( z > maxZ )
      maxZ = z ;
    else if ( z < minZ )
      minZ = z ;
  }

  boundingBox[0] = minX ;
  boundingBox[1] = maxX ;
  boundingBox[2] = minY ;
  boundingBox[3] = maxY ;
  boundingBox[4] = minZ ;
  boundingBox[5] = maxZ ;

  // skip the next couple of lines - all we need is the number of triangles
  infile.getline ( line, 200 ) ;
  infile.getline ( line, 200 ) ;
  infile.getline ( line, 200 ) ;
  nTris = atoi(strrchr(line,' ')) ;
  infile.getline ( line, 200 ) ;
  
  for ( i = 0 ; i < nTris ; i++)
  {
    infile >> index >> v1 >> v2 >> v3 ;
    assert ( index == i ) ;
    tris.push_back ( v1 ) ;
    tris.push_back ( v2 ) ;
    tris.push_back ( v3 ) ;
  }

  infile.close () ;
}

/*void vol2surf::setupSOV(sov::FlVTKDisplay *sovDisplay)
{
  
  mScene = SceneType::New() ;
  mScene->AddSpatialObject ( mMeshObject ) ;

  // Create a 3D VTK Renderer 
  typedef sov::VTKRenderer3D RendererType;
  RendererType::Pointer myRenderer3D = RendererType::New();

  // Create the Render Method and plug it to the renderer 
  // Surface Render Method 3D 
  typedef sov::SurfaceMeshVTKRenderMethod3D<MeshType> SurfaceRenderMethodType;
  SurfaceRenderMethodType::Pointer surfaceRenderMethod3D = SurfaceRenderMethodType::New();
  
  // Add the object to the renderer
  myRenderer3D->AddRenderMethod(surfaceRenderMethod3D);
  
  myRenderer3D->SetScene(mScene);
  sovDisplay->AddRenderer(myRenderer3D);
  sovDisplay->Update();
  
}*/

void vol2surf::VTKSetup()
{
  std::string outfileName(filename) ;
  int i ;
  
  i = outfileName.find_last_of ('.') ;
  outfileName.replace ( i+1, 4, "vtk" ) ;
  
    // read in the vtk polydata file
  vtkPolyDataReader *meshReader = vtkPolyDataReader::New() ;
  meshReader->SetFileName(outfileName.c_str());

  // The mapper is responsible for pushing the geometry into the graphics
  // library. It may also do color mapping, if scalars or other attributes
  // are defined.
  //
  vtkPolyDataMapper *meshMapper = vtkPolyDataMapper::New();
  meshMapper->SetInput(meshReader->GetOutput());
  meshMapper->Update() ;

  // The actor is a grouping mechanism: besides the geometry (mapper), it
  // also has a property, transformation matrix, and/or texture map.
  // In this example we create eight different spheres (two rows of four
  // spheres) and set the specular lighting coefficients. A little ambient
  // is turned on so the sphere is not completely black on the back side.
  //
  vtkActor *meshActor = vtkActor::New();
  meshActor->SetMapper(meshMapper);
  
  // Create the graphics structure. The renderer renders into the 
  // render window. The render window interactor captures mouse events
  // and will perform appropriate camera or actor manipulation
  // depending on the nature of the events.
  //
  vtkRenderer *ren1 = vtkRenderer::New();
  vtkRenderWindow *renWin = vtkRenderWindow::New();
  renWin->AddRenderer(ren1);
  vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
  iren->SetRenderWindow(renWin);

  // Add the actors to the renderer, set the background and size.
  //
  ren1->AddActor(meshActor);
  ren1->SetBackground(0.1, 0.2, 0.4);
  renWin->SetSize(400, 200);

  // Set up the lighting.
  //
  vtkLight *light = vtkLight::New();
  light->SetFocalPoint(1.875,0.6125,0);
  light->SetPosition(0.875,1.6125,1);
  ren1->AddLight(light);

  // We want to eliminate perspective effects on the apparent lighting.
  // Parallel camera projection will be used. To zoom in parallel projection
  // mode, the ParallelScale is set.
  //
  //ren1->GetActiveCamera()->SetFocalPoint(0,0,0);
  //ren1->GetActiveCamera()->SetPosition(0,0,1);
  //ren1->GetActiveCamera()->SetViewUp(0,1,0);
  //ren1->GetActiveCamera()->ParallelProjectionOn();
  ren1->ResetCamera();
  ren1->GetActiveCamera()->SetParallelScale(1.5);
  
  // This starts the event loop and invokes an initial render.
  //
  iren->Initialize();
  iren->Start();

  // Exiting from here, we have to delete all the instances that
  // have been created.
  //
  meshActor->Delete();
  meshMapper->Delete();
  meshReader->Delete();
  ren1->Delete();
  renWin->Delete();
  iren->Delete();
}

void vol2surf::setThreshold(const char *value)
{
  manThresholdValue = atof ( value ) ;
}
void vol2surf::setUpperThreshold(const char *value)
{
  manUpperThresholdValue = atof ( value ) ;
}

void vol2surf::hasThreshold(bool state)
{
  manThreshold = state ;
}

void vol2surf::writeSO()
{
  ofstream outfile ;
  int nPts, nTris ;
  std::string outfileName(filename) ;
  int i ;
  
  i = outfileName.find_last_of ('.') ;
  outfileName.replace ( i+1, 4, "meta" ) ;
  outfile.open ( outfileName.c_str() ) ;
  nPts = (verts.size () + 1) / 3 ;
  nTris = (tris.size () + 1) / 3 ;

  // output header
  outfile << "ObjectType = Mesh" << endl << "NDims = 3" << endl << "ID = 0" << endl ;
  outfile << "TransformMatrix = 1 0 0 0 1 0 0 0 1" << endl << "Offset = 0 0 0" << endl << "CenterOfRotation = 0 0 0" << endl ;
  outfile << "ElementSpacing = 1 1 1" << endl << "PointType = MET_FLOAT" << endl << "PointDataType = MET_FLOAT" << endl ;
  outfile << "CellDataType = MET_FLOAT" << endl << "NCellTypes = 1" << endl << "PointDim = ID x y ..." << endl ;
  outfile << "NPoints = " << nPts << endl ;

  outfile << "Points = " << endl ;
  
  for ( i = 0 ; i < nPts ; i++ )
  {
    outfile << i << " " << verts[3*i] << " " << verts[3*i+1] << " " << verts[3*i+2] << endl ;
  }
  
  outfile << endl << "CellType = TRI" << endl << "NCells = " << nTris << endl << "Cells = " << endl ;

  for ( i = 0 ; i < nTris ; i++)
  {
    outfile << i << " " << tris[3*i] << " " << tris[3*i+1] << " " << tris[3*i+2] << endl ;
  }

  outfile.close () ;
}

void vol2surf::readMeshBack()
{
  // reads in the VTK mesh file we just generated

  ifstream infile ;
  std::string vtkFileName (filename);
  char line[200] ;
  double x, y, z ;
  int index, v1, v2, v3 ;
  int nPts, nTris ;
  int i ;

  i = vtkFileName.find_last_of ('.') ;
  vtkFileName.replace ( i+1, 4, "vtk" ) ;
  
  infile.open ( vtkFileName.c_str() ) ;
  
  verts.clear () ;
  tris.clear () ;
  // skip header - all we need here is the number of vertices
  for (int i = 0 ; i < 5 ; i++)
    infile.getline (line, 200 ) ;
  nPts = atoi(strchr(line,' ')) ;

  for ( i = 0 ; i < nPts ; i++ )
  {
    infile >> x >> y >> z ;
    verts.push_back ( x ) ;
    verts.push_back ( y ) ;
    verts.push_back ( z ) ;
  }

  // skip the next couple of lines - all we need is the number of triangles
  infile.getline ( line, 200 ) ;
  while ( strncmp (line, "POLYGONS", 8) )
    infile.getline ( line, 200 ) ;
  nTris = atoi(strchr(line,' ')) ;
  
  for ( i = 0 ; i < nTris ; i++)
  {
    infile >> index >> v1 >> v2 >> v3 ;
    assert ( index == 3 ) ;
    tris.push_back ( v1 ) ;
    tris.push_back ( v2 ) ;
    tris.push_back ( v3 ) ;
  }

  infile.close () ;
}
int vol2surf::CommandLine ( int argc, const char *argv[] )
{
  this->commandLine = true ;

  // get the arguments
  char *inputFileName ;
  char *outputFileName;
  bool so, vtk, byu, iv, stl, gaussian, extractLabel, threshold ;
  double var, label, min, max ;

  inputFileName  = ipGetStringArgument(argv, "-input", NULL);  
  outputFileName = ipGetStringArgument(argv, "-outfile", NULL);

  so = ipExistsArgument (argv, "-so") ;
  vtk = ipExistsArgument (argv, "-vtk") ;
  byu = ipExistsArgument (argv, "-byu") ;
  iv = ipExistsArgument (argv, "-iv") ;
  stl = ipExistsArgument (argv, "-stl") ;

  decimationOn = ipExistsArgument(argv, "-decimate");
      quadricDecimationOn = decimationOn && ipExistsArgument(argv, "-quadric");
  targetReduction = ipGetIntArgument(argv, "-targetReduction", targetReduction);

  gaussian = ipExistsArgument (argv, "-gaussian") ;
  extractLabel = ipExistsArgument (argv, "-label") ;
  threshold = ( ipExistsArgument (argv, "-min") || ipExistsArgument (argv, "-max") );

  var = ipGetDoubleArgument (argv, "-gaussian", -1) ;
  label = ipGetDoubleArgument (argv, "-label", -1) ;
  min = ipGetDoubleArgument (argv, "-min", -1) ;
  max = ipGetDoubleArgument (argv, "-max", -1) ;
 
  // make sure the arguments are valid
  if ( !inputFileName || ipExistsArgument(argv, "-help") )
  {
	  cout << "Vol2surf (" << __DATE__ << ")" << endl;
    cout << "Usage:" << endl ;
    cout << "Vol2surf -input fileName [-outfile fileName] [options...]" << endl;
		cout << "		[-so] [-vtk] [-byu] [-iv] [-stl] " << endl;
		cout << "		[-gaussian var] [-label l] [-min a -max b]" << endl ;
		cout << "		[-decimate [-quadric] [-targetReduction ##]]" << endl;
    return 0 ;
  }

  this->setFilename (inputFileName) ;

  if (outputFileName == NULL) {
		outFileType = UNKNOWN;
	} else {
		std::string outfilename = std::string(outputFileName);
	  std::string::size_type npos = outfilename.find_last_of(".");
	  if (npos > 0) {
			std::string extension = outfilename.substr(npos + 1);
			if ("so" == extension) {
				outFileType = SO;
				outFileName = outfilename.substr(0, npos) + ".so";
			} else if ("vtk" == extension) {
				outFileType = VTK;
				outFileName = outfilename.substr(0, npos) + ".vtk";
			} else if ("byu" == extension) {
				outFileType = BYU;
				outFileName = outfilename.substr(0, npos) + ".byu";
			} else if ("iv" == extension) {
				outFileType = IV;
				outFileName = outfilename.substr(0, npos) + ".iv";
			} else if ("stl" == extension) {
				outFileType = STL;
				outFileName = outfilename.substr(0, npos) + ".stl";
			}else if ("ply" == extension) {
				outFileType = PLY;
				outFileName = outfilename.substr(0, npos) + ".ply";
			}
			cout << "Output File: " << outFileName << endl;
	  } else {
			outFileType = UNKNOWN;
		}
  }

  this->setOutputFormat (so, vtk, byu, iv, stl) ;
  this->hasGaussian (gaussian);
  this->hasLabel (extractLabel) ;
  this->hasThreshold (threshold) ;
  this->setThreshold (min) ;
  this->setUpperThreshold (max) ;
  
  if (gaussian) 
    this->setVariance (var) ;
  if ( extractLabel )
    this->setLabel ((int) label) ;

  
  this->Run () ;

  int nPts = (verts.size () + 1) / 3 ;
  int nTris = (tris.size () + 1) / 3 ;

  cout << "Number of Points: " << nPts << endl;
	cout << "Number of Triangles: " << nTris << endl;
  cout << "Surface Area: " << this->getSurfaceAreaStr () << endl ;
  cout << "Volume: " << this->getVolumeStr () << endl ;
  cout << "Sphere Radii Ratio: " << this->getSphereRadiiRatioStr() << endl ;

  return 0 ;
}
