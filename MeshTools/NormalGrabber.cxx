#include <fstream>
#include <iostream>
#include <cstring>
#include <stdlib.h>

int main (int argc, char **argv)
{  
   std::cout << "Tool for generating KWMeshVisu style attribute files containing the spatial location of points from a polydata file. " << std::endl ;
   std::cout << "Usage: " << argv[0] << " vtkPolyData  outputFilePrefix" << std::endl ;
   if ( argc < 3 ) 
      return ( -1 ) ;

   std::ifstream polydata ( argv[1] ) ;

    char outputFN1[100];
    strcpy(outputFN1,argv[2]);
    strcat(outputFN1,"1.txt");
    
    char outputFN2[100];
    strcpy(outputFN2,argv[2]);
    strcat(outputFN2,"2.txt");


    char outputFN3[100];
    strcpy(outputFN3,argv[2]);
    strcat(outputFN3,"3.txt");
   

   std::ofstream output1 (outputFN1 ) ;
   std::ofstream output2 (outputFN2);
   std::ofstream output3 (outputFN3);



   char line[500] ;
   double point[3] ;

   while ( strncmp ( line, "POINT_DATA", 10 ) !=0 ) 
   { 
      polydata.getline ( line, 500 )  ;
   }
  
   int nPoints = atoi ( line + 10 ) ;
  
    std::cout << "POINT_DATA=" << nPoints << std::endl ;
    polydata.getline ( line, 500 )  ; //reading: Normals normalsfloat

   output1 << "NUMBER_OF_POINTS=" << nPoints << std::endl ;
   output1 << "DIMENSION=1" << std::endl << "TYPE=Scalar" << std::endl ;


   output2 << "NUMBER_OF_POINTS=" << nPoints << std::endl ;
   output2 << "DIMENSION=1" << std::endl << "TYPE=Scalar" << std::endl ;

   output3 << "NUMBER_OF_POINTS=" << nPoints << std::endl ;
   output3 << "DIMENSION=1" << std::endl << "TYPE=Scalar" << std::endl ;




   for ( int i = 0 ; i < nPoints ; i++ )
   {
      polydata >> point[0] >> point[1] >> point[2] ;
      output1 << point[0] << std::endl ; 
      output2 << point[1] << std::endl ;
      output3 << point[2] << std::endl ;
   }   
   polydata.close () ;
   output1.close () ;
   output2.close () ;
   output3.close () ;
   return 0 ;
}
