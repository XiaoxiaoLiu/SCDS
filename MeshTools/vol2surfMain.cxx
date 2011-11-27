#include "vol2surfGUI.h"

int main(int argc, const char*argv[])
{
  if ( argc == 1 ) 
  {
    vol2surfGUI gui ;
    gui.show(argc, (char **)argv);
    return Fl::run();
  }
  
  else 
  {
    vol2surf cmdObject(0, 0, 800, 600, "Vol2surf") ;
    return cmdObject.CommandLine ( argc, argv ) ;
  }
  
}
