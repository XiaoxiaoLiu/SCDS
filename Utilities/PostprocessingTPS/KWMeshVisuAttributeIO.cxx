#include "KWMeshVisuAttributeIO.h"

KWMeshVisuAttributeIO::KWMeshVisuAttributeIO () 
{
  this->m_attributes = 0 ;
  this->m_fileName = "" ;
}

KWMeshVisuAttributeIO::~KWMeshVisuAttributeIO ()
{
}

void KWMeshVisuAttributeIO::ReadAttributes()
{
  if ( this->m_attributes ) 
    {
      this->m_attributes->Delete () ;
    }

  this->m_attributes = vtkFloatArray::New() ;
  
  // read in the file
  std::ifstream attrFile ;
  attrFile.open ( this->m_fileName.c_str() ) ;
  int nPts, nDim ;

  // This is a quite sloppy implementation of reading keyword/value pairs!
  bool found;
  char * valuePtr;
  char typeString[1001], line[1001];

  attrFile.seekg(0,std::ios::beg);
  found = false ;
  while ( !found && !attrFile.eof())
  { attrFile.getline ( line, 1000 ) ;
    if (line[0] != '#' && strstr ( line, "NUMBER_OF_POINTS" )) found = true;
  }
  valuePtr=strchr(line, '=');
  if (!valuePtr) return;
  valuePtr++;
  sscanf(valuePtr, " %d ", &nPts);

  attrFile.seekg(0,std::ios::beg);
  found = false ;
  while ( !found && !attrFile.eof())
  { attrFile.getline ( line, 1000 ) ;
    if (line[0] != '#' && strstr ( line, "DIMENSION" )) found = true;
  }
  valuePtr=strchr(line, '=');
  if (!valuePtr) return;
  valuePtr++;
  sscanf(valuePtr, " %d ", &nDim);

  attrFile.seekg(0,std::ios::beg);
  found = false ;
  while ( !found && !attrFile.eof())
  { attrFile.getline ( line, 1000 ) ;
    if (line[0] != '#' && strstr ( line, "TYPE" )) found = true;
  }
  valuePtr=strchr(line, '=');
  if (!valuePtr) return;
  valuePtr++;
  sscanf(valuePtr, " %s ", typeString);

  attrFile.seekg(0,std::ios::beg);
  const int numEntries = 3;
  int counter = 0;
  while ( counter < numEntries && !attrFile.eof())
  { attrFile.getline ( line, 1000 ) ;
    if ((line[0] != '#')) counter++;
  }

  this->m_attributes->SetNumberOfComponents ( nDim ) ;
  float data[9] ; 
  for (int i = 0 ; i < nPts ; i++ )
  {
    for ( int j = 0 ; j < nDim ; j++ )
      {
	attrFile >> data[j] ;
      }

    this->m_attributes->InsertTupleValue (i, data) ;
  }
  
  attrFile.close () ;
}

void KWMeshVisuAttributeIO::WriteAttributes ()
{
  std::ofstream out ;
  out.open ( this->m_fileName.c_str() ) ;

  // print out header
  unsigned long int nPoints = this->m_attributes->GetNumberOfTuples () ;
  out << "NUMBER_OF_POINTS=" << nPoints << std::endl ;
  short int dim = this->m_attributes->GetNumberOfComponents () ;
  out << "DIMENSION=" << dim << std::endl ;

  if ( dim == 1 ) 
    out << "TYPE=Scalar" << std::endl ;
  else if ( dim == 3 ) 
    out << "TYPE=Vector" << std::endl ;
  else if ( dim == 9 ) 
    out << "TYPE=Tensor" << std::endl ;
  else 
    out << "TYPE=Unknown" << std::endl ;

  float *tuple = new float[dim] ;
  for ( unsigned long int i = 0 ; i < nPoints ; i++ )
    {
      this->m_attributes->GetTupleValue ( i, tuple ) ;
      for ( short int j = 0 ; j < dim ; j++ )
	{
	  out << tuple[j] << " " ;
	}
      out << std::endl ;
    }

  out.close () ;
}
