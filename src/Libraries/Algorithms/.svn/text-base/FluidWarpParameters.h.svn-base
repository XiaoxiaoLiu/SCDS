#ifndef FluidWarpParameters_h
#define FluidWarpParameters_h

#include <iosfwd>

class FluidWarpParameters
{
public:
  unsigned int    numIterations;
  double          alpha;
  double          beta;
  double          gamma;
  double          maxPerturbation;
  unsigned int    numBasis;  
  bool            jacobianScale;
  bool            divergenceFree;

  FluidWarpParameters();
  void writeASCII(std::ostream& output) const;
  void readASCII(std::istream& input);
};

std::ostream& operator<<(std::ostream& output, 
			 const FluidWarpParameters& parameters);
std::istream& operator>>(std::istream& input, 
			 FluidWarpParameters& parameters);
#endif
