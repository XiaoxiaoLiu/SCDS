#include "PostprocessingTPSCLP.h"

#include <itksys/SystemTools.hxx>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "vtkPolyData.h"
#include "vtkPolyDataReader.h"
#include "vtkPolyDataWriter.h"
#include "vtkFloatArray.h"
#include "vtkGeneralTransform.h"
#include "vtkThinPlateSplineTransform.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkTriangle.h"
#include "vtkPlane.h"
#include "vtkKdTreePointLocator.h"
#include "vtkTriangleFilter.h"

#include "MeshIO.h"


static bool _saveTPSFiles = false;

void Print3 ( float *a, std::string label )
{
  std::cout << label << ": " << a[0] << " " << a[1] << " " << a[2] << std::endl ;
}

int parseParameterFile ( std::string fileName, std::vector < std::string > &particleFileNames, std::vector < std::string > &meshFileNames, std::vector<std::string> &meshOutputNames, std::vector<std::string> &tpsOutputNames )
{
  std::ifstream input ;
  input.open ( fileName.c_str () ) ;

  bool found;
  char * valuePtr;
  char line[1001];
  int nShapes ;

  input.seekg(0,std::ios::beg);
  found = false ;
  while ( !found && !input.eof())
    { input.getline ( line, 1000 ) ;
      if (line[0] != '#' && strstr ( line, "NUMBER_OF_SHAPES" )) found = true;
    }
  valuePtr=strchr(line, '=');
  if (!valuePtr) return false ;
  valuePtr++;
  sscanf(valuePtr, " %d ", &nShapes);
  std::cout << "nShapes: " << nShapes << std::endl ;

  input.seekg(0,std::ios::beg);
  const int numEntries = 1;
  int counter = 0;
  while ( counter < numEntries && !input.eof())
    { input.getline ( line, 1000 ) ;
      if ((line[0] != '#')) counter++;
    }
  
  particleFileNames.resize ( nShapes ) ;
  meshFileNames.resize ( nShapes ) ;
  meshOutputNames.resize(nShapes);
  tpsOutputNames.resize(nShapes);
  for (int i = 0 ; i < nShapes ; i++ )
    {
      if (_saveTPSFiles) {
        input >> particleFileNames[i] >> meshFileNames[i] >> meshOutputNames[i] >> tpsOutputNames[i];
        std::cout << particleFileNames[i] << " : " << meshFileNames[i] << " >> " << meshOutputNames[i] << " & " << tpsOutputNames[i] << std::endl ;
      } else {
        input >> particleFileNames[i] >> meshFileNames[i] >> meshOutputNames[i];
        std::cout << particleFileNames[i] << " : " << meshFileNames[i] << " >> " << meshOutputNames[i] << std::endl ;
      }
    }

  input.close () ;
  return nShapes ;
}

vtkPoints *parseParticleFile ( std::string fileName ) 
{
  int counter = 0;

  std::ifstream in( fileName.c_str() );
  if ( !in )
    {
      std::cout << "Could not open point file for input: " << fileName << std::endl ;
      return NULL ;
    }

  vtkPoints *particles = vtkPoints::New () ;

  float pt[3] ;
  while ( in )
    {
      in >> pt[0] >> pt[1] >> pt[2] ;
      if ( in ) 
	{  
	  particles->InsertNextPoint ( pt ) ;
	  counter++ ;
	}
    }
 
  std::cout << particles->GetNumberOfPoints () << " particles read from " << fileName << std::endl ;

  in.close();

  return particles ;
}

void CheckMeshQuality ( vtkPolyData *mesh )
{
  std::cout << "Checking mesh quality" << std::endl ;
  long int nTriangles = mesh->GetNumberOfCells () ;
  vtkIdList *verts = vtkIdList::New () ;
  double area, a[3], b[3], c[3] ;
  for ( long int i = 0 ; i < nTriangles ; i++ )
    {
      mesh->GetCellPoints ( i, verts ) ;
      mesh->GetPoint ( verts->GetId(0), a ) ;
      mesh->GetPoint ( verts->GetId(1), b ) ;
      mesh->GetPoint ( verts->GetId(2), c ) ;
      area = vtkTriangle::TriangleArea ( a, b, c ) ;
      if ( area < 0.0001 ) std::cout << area << std::endl ;
    }
}

static void ProjectMeshWithLocator(vtkPolyData* m1, vtkPolyData* m2) {
  if (m2->GetCell(0)->GetCellType() == VTK_TRIANGLE_STRIP) {
    vtkTriangleFilter* filter = vtkTriangleFilter::New();
    filter->SetInput(m2);
    filter->Update();
    m2 = filter->GetOutput();
  }

  vtkKdTreePointLocator* pointLocator = vtkKdTreePointLocator::New();
  pointLocator->SetDataSet(m2);
  pointLocator->BuildLocator();

  vtkPoints* pl1 = m1->GetPoints();
  m2->BuildLinks();

  double x[3], cp[3];
  for (int i = 0; i < pl1->GetNumberOfPoints(); i++) {
    pl1->GetPoint(i, x);
    int id = pointLocator->FindClosestPoint(x);
    m2->GetPoint(id, cp);

    vtkIdList* cellIds = vtkIdList::New();
    m2->GetPointCells(id, cellIds);
    
    /**
     * projection a point onto cell surface
     * 1) find cell
     * 2) compute normal 
     * 3) construct implicit plane function
     * 4) projection onto the plane
     * 5) determine the projected point is inside the cell
     * 6) find the point with minimum distance
     */
    double txp[3], xp[3], dist2 = vtkMath::Distance2BetweenPoints(x, cp);
    memcpy(txp, cp, sizeof(cp));

    for (int j = 0; j < cellIds->GetNumberOfIds(); j++) {
      vtkCell* cell = m2->GetCell(j);

      double t1[3], t2[3], t3[3], normal[3];
      cell->GetPoints()->GetPoint(0, t1);
      cell->GetPoints()->GetPoint(1, t2);
      cell->GetPoints()->GetPoint(2, t3);

      vtkTriangle::ComputeNormal(t1, t2, t3, normal);
      vtkPlane::GeneralizedProjectPoint(x, t1, normal, xp);

      int inside = vtkTriangle::PointInTriangle(xp, t1, t2, t3, 0);
      double ndist2 = vtkMath::Distance2BetweenPoints(x, xp);
      if (inside == 1) {
        if (dist2 > ndist2) {
          memcpy(txp, xp, sizeof(xp));
          cout << i << " vertex found closer projection onto a polygon (" << ndist2 << " < " << dist2 << ")" << endl;
          dist2 = ndist2;
        }
      } else {
        if (dist2 > ndist2) {
        }
      }
    }

   m1->GetPoints()->SetPoint(i, txp);
  }
}


void ComputeMeshFromParticles ( vtkPoints *subjectParticles, vtkPoints *templateParticles, vtkPolyData *templateMesh, std::string origMeshFileName, const char *fileName, const char *tpsFileName, bool projectToSurface, const char *debug, float factor )
{
  vtkThinPlateSplineTransform *thin = vtkThinPlateSplineTransform::New () ;
  thin->SetSourceLandmarks ( templateParticles ) ;
  thin->SetTargetLandmarks ( subjectParticles ) ;
  thin->SetBasisToR () ;

  vtkGeneralTransform *t1 = vtkGeneralTransform::New () ;
  t1->SetInput ( thin ) ;

  vtkTransformPolyDataFilter *filter = vtkTransformPolyDataFilter::New () ;
  filter->SetInput ( templateMesh ) ;
  filter->SetTransform ( t1 ) ;
  filter->Update () ;

  vtkPolyData *tpsMesh = filter->GetOutput () ;
  std::cout << "TPS computation finished. " << std::endl ;

  if (_saveTPSFiles) {
    vtkPolyDataWriter *tpsWriter = vtkPolyDataWriter::New () ;
    tpsWriter->SetFileName ( tpsFileName ) ;
    tpsWriter->SetInput ( tpsMesh ) ;
    tpsWriter->Update () ;
  }

  std::string distLogFileName = std::string(fileName) + std::string(".distlog");
  std::cout << "Distance log file: " << distLogFileName << endl;

  ofstream fdistlog;
  fdistlog.open(distLogFileName.c_str());

  if ( projectToSurface ) {
    std::vector < short int > used ;
    std::vector < vtkIdType > closestVerts ;

    std::cout << "Starting projection. " << std::endl ;

    vtkPolyDataReader *meshReader = vtkPolyDataReader::New () ;
    meshReader->SetFileName(origMeshFileName.c_str()) ;
    meshReader->Update () ;
    vtkPolyData *origMesh = meshReader->GetOutput() ;

    ProjectMeshWithLocator(tpsMesh, origMesh);

    /*
    double cp[3], x[3] ;
    for ( int i = 0 ; i < tpsMesh->GetNumberOfPoints () ; i++ ) {
      if ( ! ( i % 1000 ) ) {
        std::cout << i << " vertices out of " << tpsMesh->GetNumberOfPoints () << " completed." << std::endl ;
      }
      tpsMesh->GetPoint ( i, x ) ;

      int subId ;
      double pCoords[3], weights[3], dist ;
      double tol = 0.5 ;
      while ( tol < 2500 ) {
        vtkCell *closestCell = origMesh->FindAndGetCell ( x, NULL, -1, tol, subId, pCoords, (double *) weights ) ;
        if (closestCell) {		    
          int result = closestCell->EvaluatePosition ( x, cp, subId, pCoords, dist, weights ) ;
          if ( result == -1 )  {
            std::cout << "argh" << std::endl ;
          } else {
            tol = -1 ;
            break ;
          }
        }
        tol *= 2 ;
      }


      bool badOrigX = false;
      if ( tol > 0 ) {
        std::cout << i << " vertex : " << tol << std::endl;
        vtkIdType closestVert ;
        double origX[3], temp[3] ;

        closestVert = origMesh->FindPoint ( x ) ;
        origMesh->GetPoint ( closestVert, origX ) ;

        double dist2 = vtkMath::Distance2BetweenPoints(x, origX);
        if (dist2 > 1) {
          badOrigX = true;
          std::cout << i << " vertex has bad origX " << std::endl;
        }
        for ( int j = 0 ; j < 3 ; j++ ) {
          temp[j] = origX[j] - x[j] ; 
          cp[j] = x[j] + temp[j] * factor ;
        }

        vtkCell *closestCell = origMesh->FindAndGetCell ( cp, NULL, -1, tol, subId, pCoords, (double *) weights ) ;
        if ( closestCell ) {		    
          double cp2[3] ;
          int result = closestCell->EvaluatePosition ( cp, cp2, subId, pCoords, dist, weights ) ;
          if ( result == -1 ) {
            std::cout << "argh" << std::endl ;
          } else {
            cp[0] = cp2[0] ;
            cp[1] = cp2[1] ;
            cp[2] = cp2[2] ;
            //std::cout << "yay" << std::endl ;
          }
        }
      //std::cout << "Distance is: " << temp[0] * temp[0] + temp[1] * temp[1] + temp[2] * temp[2] << std::endl ;
      } else if (badOrigX) {
        std::cout << i << " bad origX has no closest cell" << std::endl;
      }

      double dist2 = vtkMath::Distance2BetweenPoints(x, cp);
      fdistlog << dist2 << std::endl;

      if (dist2 > 1) {
        std::cout << i << " distance = " << dist2 << std::endl;
      }

      tpsMesh->GetPoints()->SetPoint(i, cp) ;
    }
    */
  }

  fdistlog.close();

  vtkPolyDataWriter *meshWriter = vtkPolyDataWriter::New () ;
  meshWriter->SetFileName ( fileName ) ;
  meshWriter->SetInput ( tpsMesh ) ;

  meshWriter->Update () ;

  CheckMeshQuality ( tpsMesh ) ;
  return ;
}

int main (int argc, char **argv)
{
  PARSE_ARGS ;

  if (workingDirectory != "") {
    itksys::SystemTools::ChangeDirectory(workingDirectory.c_str());
  }
  _saveTPSFiles = saveTPSFiles;

  std::cout << parameterFileName << std::endl ;
  std::cout << templateID << std::endl ;
  if (outputDirectory.size() != 0 && outputDirectory[outputDirectory.size()-1] != '/')
    outputDirectory.append("/");
  
  // read the parameter file 
  std::vector < std::string > particleFileNames ;
  std::vector < std::string  >  meshFileNames ;
  std::vector < std::string > meshOutputNames ;
  std::vector < std::string > tpsOutputNames ;
  int nShapes = parseParameterFile ( parameterFileName, particleFileNames, meshFileNames, meshOutputNames, tpsOutputNames ) ;

  // STEP 1
  // read template mesh
  vtkPolyData *templateMesh ;
  MeshIO *meshIOtool = new ( MeshIO ) ;
  if ( meshIOtool->Read ( meshFileNames[templateID] ) ) 
    {
      templateMesh = meshIOtool->GetMesh () ;
      std::cout << templateMesh->GetNumberOfPoints () << std::endl ;
    }
  else
    {
      std::cout << "Can't read template mesh: " << meshFileNames[templateID] << std::endl ;
      exit ( 1 ) ;
    }

  // read the particles for the template 
  vtkPoints *templateParticles = parseParticleFile ( particleFileNames[templateID] ) ;
  
  for ( int currentShape = 0 ; currentShape < nShapes ; currentShape++ )
    {
      // STEP 2
      vtkPoints *currentParticles = parseParticleFile ( particleFileNames[currentShape] ) ;

      std::string meshFileName (particleFileNames[currentShape]);
      meshFileName.append(".vtk") ;

      if (outputDirectory.size() != 0)
	{
           int found = meshFileName.find_last_of("/\\");
           meshFileName = meshFileName.substr ( found+1 ) ;
  	   meshFileName.insert(0,outputDirectory) ;
	}

      std::string attrFileName (particleFileNames[currentShape]);
      attrFileName.append(".txt") ;

      if (outputDirectory.size() != 0)
	{
           int found = attrFileName.find_last_of("/\\");
           attrFileName = attrFileName.substr ( found+1 ) ;
  	   attrFileName.insert(0,outputDirectory) ;
	}

      ComputeMeshFromParticles ( currentParticles, templateParticles, templateMesh, meshFileNames[currentShape], meshOutputNames[currentShape].c_str(), tpsOutputNames[currentShape].c_str(), projectToSurface, attrFileName.c_str(), projectionFactor ) ;

      //currentParticles->Delete () ;
    }
  
  // release memory
  //templateParticles->Delete () ;
  //templateMesh->Delete () ;
  //delete ( meshIOtool ) ;

  return 0 ;
}

