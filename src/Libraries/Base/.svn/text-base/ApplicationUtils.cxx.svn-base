#include "ApplicationUtils.h"
#include "itkImageIOBase.h"
#include "itkImageIOFactory.h"
#include <string>

void 
ApplicationUtils::
ParseOptionVoid(int& argc, char**& argv)
{
  argv++; argc--;
}

double 
ApplicationUtils::
ParseOptionDouble(int& argc, char**& argv)
{
  std::string s = ApplicationUtils::ParseOptionString(argc, argv);
  return atof(s.c_str());
}

int
ApplicationUtils::
ParseOptionInt(int& argc, char**& argv)
{
  std::string s = ApplicationUtils::ParseOptionString(argc, argv);
  return atoi(s.c_str());
}

bool
ApplicationUtils::
ParseOptionBool(int& argc, char**& argv)
{
  std::string s = ApplicationUtils::ParseOptionString(argc, argv);
  if (s == "true" || s == "True" || s == "TRUE" || s == "1")
  {
    return true;
  }
  return false;
}

std::string 
ApplicationUtils::
ParseOptionString(int& argc, char**& argv)
{
  if (argc == 0)
  {
    std::stringstream ss;
    ss << "Error parsing argument: no arguments";
    throw std::runtime_error(ss.str().c_str());
  }

  std::string arg(argv[0]);
  std::string::size_type equalPos = arg.find('=');
  std::string val = "";
  if (equalPos != std::string::npos)
  {
    // parse key=value and return value
    val = arg.substr(equalPos+1);
    argv++; argc--;
  }
  else if (argc > 1)
  {
    val = argv[1];
    argv++; argc--;
    argv++; argc--;
  }
  else
  {
    std::stringstream ss;
    ss << "Error parsing argument: " << arg;
    throw std::runtime_error(ss.str().c_str());
  }

  // uncomment this to debug argument parsing
  //std::cerr << "arg=" << arg << ", val=" << val << std::endl;

  return val;
}

bool
ApplicationUtils::
ITKHasImageFileReader(const char* fileName)
{
  itk::ImageIOBase::Pointer imageIO = 
    itk::ImageIOFactory::CreateImageIO(fileName, 
                                       itk::ImageIOFactory::ReadMode);
  return (!imageIO.IsNull());
}
