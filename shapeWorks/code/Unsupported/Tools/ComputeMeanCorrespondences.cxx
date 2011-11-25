/*=========================================================================
  Program:   ShapeWorks: Particle-based Shape Correspondence & Visualization
  Module:    $RCSfile: ComputeMeanCorrespondences.cxx,v $
  Date:      $Date: 2009/05/06 21:49:16 $
  Version:   $Revision: 1.1.1.1 $
  Author:    $Author: cates $

  Copyright (c) 2009 Scientific Computing and Imaging Institute.
  See ShapeWorksLicense.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.
=========================================================================*/
#include <iostream>
#include <vector>
#include <string>
#include "itkParticlePositionReader.h"
#include "itkParticlePositionWriter.h"
#include "param.h"
#include <iostream>
#include <fstream>

int main(int argc, char *argv[])
{ 
  if (argc != 2)
    {
    std::cerr << "Use: " << argv[0] << "paramfile" << std::endl;
    return 1;
    }
  typedef itk::ParticlePositionReader<3>::PointType PointType;
  
  param::parameterFile pf(argv[1]);
  std::string outputfile;  
  std::vector< std::string > inputfiles;
  //  std::cout << pf << std::endl;
  std::string tmpa, outfile;
  int ii = 0;
  bool ok = true;
  PARAMSET(pf, outfile, "output", 0, ok, "outputfile.txt");

  
  while (ok == true)
  {
  // Record the point file names.
  PARAMSET(pf, tmpa, "point_files", ii, ok, "");
 
  if (ii==0 && ok != true)
    {
    std::cerr << "No input/output files have been specified" << std::endl;
    throw 1;
    }
  if (ok == true)
    {
    inputfiles.push_back(tmpa);
    //    std::cout << tmpa << std::endl;
    } // if ok == true
  ii++;
  } // while ok == true
  
  itk::ParticlePositionReader<3>::Pointer readerA
    = itk::ParticlePositionReader<3>::New();
  readerA->SetFileName(inputfiles[0].c_str());
  readerA->Update();
  unsigned int N =  readerA->GetOutput().size();


  PointType point;
  std::vector<PointType> meanpoints;
  point[0] = point[1] = point[2] =0.0;
  for (unsigned int i = 0; i < N; i++)
    {
    meanpoints.push_back(point);
    }

  // For all files
  for (unsigned int i = 0; i < inputfiles.size() ; i++)
    {
    itk::ParticlePositionReader<3>::Pointer reader
      = itk::ParticlePositionReader<3>::New();
    reader->SetFileName(inputfiles[i].c_str());
    reader->Update();

    // For all points
    for(int j = 0; j < N; j++)
      {
      meanpoints[j][0] += (reader->GetOutput()[j])[0];
      meanpoints[j][1] += (reader->GetOutput()[j])[1];
      meanpoints[j][2] += (reader->GetOutput()[j])[2];
      }    
    }

  double ninv = 1.0 / static_cast<double>(inputfiles.size());
  for(int j = 0; j < N; j++)
    {
    meanpoints[j][0] *= ninv;
    meanpoints[j][1] *= ninv;
    meanpoints[j][2] *= ninv;
    }

  itk::ParticlePositionWriter<3>::Pointer writer =
    itk::ParticlePositionWriter<3>::New();
  writer->SetInput(meanpoints);
  writer->SetFileName(outfile.c_str());
  writer->Update();
  
  return 0;
}
