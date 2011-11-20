#ifndef LINE3D_H
#define LINE3D_H

#include "Vector3D.h"
#include <iosfwd>

class Line3D
{
public:
  Vector3D<double> p1, p2;
  void set(const Vector3D<double>& p1In,
           const Vector3D<double>& p2In)
  {
    p1 = p1In;
    p2 = p2In;
  }
};

bool operator<(const Line3D& lhs, const Line3D& rhs);
bool operator==(const Line3D& lhs, const Line3D& rhs);
std::ostream& operator<<(std::ostream& output, const Line3D& line);

#endif
