#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <cstring>

int main (int argc, char **argv)
{  
   std::cout << "Tool for generating KWMeshVisu style attribute files containing the spatial location of points from a polydata file. Specify channel 0 for x, 1 for y, 2 for z. " << std::endl ;
   std::cout << "Usage: " << argv[0] << " vtkPolyData channel[0,1,2] outputFile" << std::endl ;
   if ( argc < 4 ) 
      return ( -1 ) ;

   std::ifstream polydata ( argv[1] ) ;
   int channel = atoi ( argv[2] ) ;
   std::ofstream output ( argv[3] ) ;

   char line[500] ;
   double point[3] ;

   while ( strncmp ( line, "POINTS", 6 ) ) 
   { 
      polydata.getline ( line, 500 )  ;
   }
  
   int nPoints = atoi ( line + 6 ) ;
  
   output << "NUMBER_OF_POINTS=" << nPoints << std::endl ;
   output << "DIMENSION=1" << std::endl << "TYPE=Scalar" << std::endl ;
 
   for ( int i = 0 ; i < nPoints ; i++ )
   {
      polydata >> point[0] >> point[1] >> point[2] ;
      output << point[channel] << std::endl ; 
   }   
   polydata.close () ;
   output.close () ;
 
   return 0 ;
}
