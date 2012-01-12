#include "vtkPolyData.h"
#include "vtkPoints.h"
#include "vtkPointData.h"
#include "vtkKdTreePointLocator.h"
#include "vtkPointLocator.h"
#include "vtkIdList.h"
#include "vtkMath.h"
#include "vtkTriangleFilter.h"
#include "vtkThinPlateSplineTransform.h"
#include "vtkGeneralTransform.h"
#include "vtkPolyDataWriter.h"
#include "vtkPolyDataReader.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkCell.h"
#include "vtkTriangle.h"
#include "vtkPlane.h"
#include "vtkXMLPolyDataReader.h"
#include "vtkXMLPolyDataWriter.h"

#include "fstream"
#include "string"

using namespace std;

static bool EndsWith(const string& a, const string& b) {
    if (b.size() > a.size()) return false;
    return std::equal(a.begin() + a.size() - b.size(), a.end(), b.begin());
}

static vtkPolyData* Read(string &name) {
	if (EndsWith(name, ".vtp")) {
		vtkXMLPolyDataReader* r = vtkXMLPolyDataReader::New();
		r->SetFileName(name.c_str());
		r->Update();
		return r->GetOutput();
	} else if (EndsWith(name, ".vtk")) {
		vtkPolyDataReader* r = vtkPolyDataReader::New();
		r->SetFileName(name.c_str());
		r->Update();
		return r->GetOutput();
	}
	return NULL;
}

static void Write(string &name, vtkPolyData* mesh) {
	if (EndsWith(name, ".vtp")) {
		vtkXMLPolyDataWriter* w = vtkXMLPolyDataWriter::New();
		w->SetFileName(name.c_str());
		w->SetInput(mesh);
		w->SetDataModeToAppended();
		w->SetCompressorTypeToZLib();
		w->EncodeAppendedDataOff();
		w->Update();
	} else if (EndsWith(name, ".vtk")) {
		vtkPolyDataWriter* w = vtkPolyDataWriter::New();
		w->SetFileName(name.c_str());
		w->SetInput(mesh);
		w->Write();
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


void ComputeMeshFromParticles ( vtkPoints *subjectParticles, vtkPoints *templateParticles, 
			vtkPolyData *templateMesh, std::string origMeshFileName, 
			const char *fileName, bool projectToSurface)
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

	/*
  if (false) {
    vtkPolyDataWriter *tpsWriter = vtkPolyDataWriter::New () ;
    tpsWriter->SetFileName ( tpsFileName ) ;
    tpsWriter->SetInput ( tpsMesh ) ;
    tpsWriter->Update () ;
  }
	*/

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
  }

  fdistlog.close();

  vtkPolyDataWriter *meshWriter = vtkPolyDataWriter::New () ;
  meshWriter->SetFileName ( fileName ) ;
  meshWriter->SetInput ( tpsMesh ) ;

  meshWriter->Update () ;

  return ;
}

int main(int argc, char* argv[]) {
	if (argc < 6) {
		cout << argv[0] << " " << " source-vtk target-vtk mesh-to-deform mesh-to-project output-vtk" << endl;
		return 0;
	}

	string srcName(argv[1]);
	string dstName(argv[2]);
	string meshToDeformName(argv[3]);

	vtkPolyData* srcPolyData = Read(srcName);
	vtkPolyData* dstPolyData = Read(dstName);
	vtkPolyData* meshToDeform = Read(meshToDeformName);

	string fileMeshToProject(argv[4]);

	vtkPoints* srcPoints = srcPolyData->GetPoints();
	vtkPoints* dstPoints = dstPolyData->GetPoints();

	ComputeMeshFromParticles(srcPoints, dstPoints, meshToDeform, fileMeshToProject, argv[5], true);
}
