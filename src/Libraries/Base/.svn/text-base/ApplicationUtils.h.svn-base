#ifndef ApplicationUtils_h
#define ApplicationUtils_h

#include <string>
#include "Image.h"

class ApplicationUtils
{
 public:
  static void        ParseOptionVoid  (int& argc, char**& argv);
  static int         ParseOptionInt   (int& argc, char**& argv);
  static double      ParseOptionDouble(int& argc, char**& argv);
  static std::string ParseOptionString(int& argc, char**& argv);
  static bool        ParseOptionBool  (int& argc, char**& argv);

  static bool ITKHasImageFileReader(const char* fileName);

  template <class VoxelType>
    static void        
    LoadImageITK(const char* fileName, Image<VoxelType>& image);

  template <class VoxelType>
    static void        
    SaveImageITK(const char* fileName, const Image<VoxelType>& image);

};

#include "ApplicationUtils.txx"

#endif // ApplicationUtils_h
