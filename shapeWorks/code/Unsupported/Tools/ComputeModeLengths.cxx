/*=========================================================================
  Program:   ShapeWorks: Particle-based Shape Correspondence & Visualization
  Module:    $RCSfile: ComputeModeLengths.cxx,v $
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
#include "CS6350.h"
#include <vector>
#include "param.h"
#include <iostream>
#include <fstream>

#include "itkParticleShapeStatistics.h"

int main(int argc, char *argv[])
{
  try
    {
    if (argc < 2)
      {
      std::cerr << "Use: " << argv[0] << " paramfile [output file]"
                << std::endl;
      return 1;
      }
    param::parameterFile pf(argv[1]);
    ParticleShapeStatistics<3> stats;
    int nummodes;
    bool ok=true;
    PARAMSET(pf, nummodes, "lda_modes", 0, ok, 1);

    if (ok==false)
      {
      std::cerr << "Need to specify the number of modes for LDA" << std::endl;
      return 1;
      }
    
    
    stats.ReadPointFiles(pf);
    stats.ComputeModes();
    stats.PrincipalComponentProjections();
    stats.FisherLinearDiscriminant(nummodes);

		if (argc == 3)
		{
			stats.WriteCSVFile(argv[2]);
		}
		else
		{
			stats.WriteCSVFile("modelengths.txt");
		}

    }
  catch (itk::ExceptionObject &e)
    {
    std::cerr << e << std::endl;
    return 1;
    }
  catch( CS6350::exception & eg )
    {
    std::cerr << eg << std::endl;
    }
  catch( std::exception & ex )
    {
    std::cerr << ex.what() << std::endl;
    }
  catch( ... )
    {
    std::cerr << "Unknown exception" << std::endl;
    return 2;
    }  
  
  return 0;
}

