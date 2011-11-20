#ifndef SURFACE_UTILS_VTK_H
#define SURFACE_UTILS_VTK_H

#include <vtkPolyData.h>

#include <Surface.h>

void SurfaceToVtkPolyData(const Surface& surface, vtkPolyData* polyData);
void VtkPolyDataToSurface(vtkPolyData* polyData, Surface& surface);
void FixSurfaceOrientation(Surface& surface);
void FixSurfaceOrientationByu(const char* in, const char* out);

#endif // ndef SURFACE_UTILS_VTK_H
