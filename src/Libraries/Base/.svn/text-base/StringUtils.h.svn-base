#ifndef StringUtils_h
#define StringUtils_h

#include <string>
#include <list>
#include <cstdlib>

#ifdef WIN32
#define STRINGUTILS_PATH_SEPERATOR "\/"  //used to be "\\"
#else
#define STRINGUTILS_PATH_SEPERATOR "/"
#endif

class StringUtils
{
public:
  // make all alpha chars upper case
  static std::string toUpper(const std::string& s);

  // make all alpha chars lower case
  static std::string toLower(const std::string& s);

  // remove whitespace from either end of string
  static std::string trimWhitespace(const std::string& s);

  // fill a list with the tokens in s.  delimiters are thrown out.
  static void tokenize(const std::string& s, 
                       const std::string& delimiters, 
                       std::list<std::string>& tokens);

  // convert to bool
  static bool toBool(const std::string& s)
  {
    return toUpper(trimWhitespace(s)) == "TRUE";
  }

  // convert to integer
  static
  int
  toInt(const std::string& s)
  {
    return atoi(s.c_str());
  }

  // convert to double
  static
  double
  toDouble(const std::string& s)
  {
    return atof(s.c_str());
  }

  // return the path-to-directory part of a path
  // e.g., /afs/radonc/pkg/test.txt yields /afs/radonc/pkg/
  static std::string getPathDirectory(const std::string& s);

  // return the filename part of a path
  // e.g.,  /afs/radonc/pkg/test.txt yields test.txt
  static std::string getPathFile(const std::string& s);

  // return the filename part of a path
  // e.g.,  /afs/radonc/pkg/test.txt yields txt
  static std::string getPathExtension(const std::string& s);

  // remove extension from path
  // e.g.,  /afs/radonc/pkg/test.txt yields /afs/radonc/pkg/test
  static std::string stripPathExtension(const std::string& s);

  // require a particular extension 
  // forceExtension("foo.png", "png") yields "foo.png".
  // forceExtension("foo.eps", "png") yields "foo.eps.png".
  static std::string forcePathExtension(const std::string& filename, 
                                        const std::string& extension);
};

#endif
