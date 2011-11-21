How to Cmake the program in windows:

How to Cmake the program in windows:
1. copy   ./CMake/FindFFTW.cmake to the cmake module folder (C:\Program Files\cmake-2.8.2-win32-x86\share\cmake-2.8\Modules)
2. set the FFTW libs according to "MySetup_win32.cmake" 
3. install  pthreads libarary accordingly:
    . copy the  FindPthreads.cmake to the cmake module folder.
    . right now, the lib contains three .h files and precompiled libraries. 
