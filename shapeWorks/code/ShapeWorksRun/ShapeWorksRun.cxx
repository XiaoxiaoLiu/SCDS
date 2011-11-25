/*=========================================================================
  Program:   ShapeWorks: Particle-based Shape Correspondence & Visualization
  Module:    $RCSfile: ShapeWorksRun.cxx,v $
  Date:      $Date: 2009/05/06 21:49:15 $
  Version:   $Revision: 1.1.1.1 $
  Author:    $Author: cates $

  Copyright (c) 2009 Scientific Computing and Imaging Institute.
  See ShapeWorksLicense.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.
=========================================================================*/
#include <iostream>
#include "ShapeWorksRunApp.h"
#include "CS6350.h"

int main(int argc, char *argv[])
{
  try
    {
    if (argc < 2)
      {
      std::cerr << "Use: " << argv[0] << " parameterfile" << std::endl;
      return 1;
      }

    ShapeWorksRunApp<> app(argv[1]);

    app.Run();
    }
  catch (itk::ExceptionObject &e)
    {
    std::cerr << e << std::endl;
    return 1;
    }
  catch( CS6350::exception & eg )
    {
    std::cerr << eg << std::endl;
    return 2;
    }
  catch( std::exception & ex )
    {
    std::cerr << ex.what() << std::endl;
    return 3;
    }
  catch( ... )
    {
    std::cerr << "Unknown exception" << std::endl;
    return 4;
    }  
  
  return 0;
}
