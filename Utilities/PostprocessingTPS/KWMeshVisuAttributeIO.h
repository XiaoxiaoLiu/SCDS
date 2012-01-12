#ifndef KWMeshVisuAttributeIO_h
#define KWMeshVisuAttributeIO_h

#include <string>
#include <fstream>
#include "vtkFloatArray.h"

class KWMeshVisuAttributeIO
{
 public: 
  KWMeshVisuAttributeIO () ;
  ~KWMeshVisuAttributeIO() ;

  void SetFileName ( std::string fileName ) 
  {
    this->m_fileName = fileName ;
  }

  void SetAttributes ( vtkFloatArray *attributes )
  {
    this->m_attributes = attributes ;
  }

  vtkFloatArray *GetAttributes ()
    {
      return this->m_attributes ;
    } 

  void ReadAttributes () ;
  void WriteAttributes () ;

 private:
  std::string m_fileName ;
  vtkFloatArray *m_attributes ;
};
#endif
