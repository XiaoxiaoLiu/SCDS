# sets the following variables
# FFTW3_FOUND         - Set to true when fftw is found.
# FFTW3_INCLUDE_DIR   - Include directory for fftw headers.
# FFTW3_LIBRARY_DIR   - Link directory for fftw libraries.

FIND_PATH(FFTW3_BASE_DIR fftw3.h
   ${RADONC_LIBS_PATH}/fftw
 
  )

IF(FFTW3_BASE_DIR)
  SET(FFTW3_FOUND 1)
  SET(FFTW3_INCLUDE_DIR ${FFTW3_BASE_DIR}/include)
  SET(FFTW3_LIBRARY_DIR ${FFTW3_BASE_DIR}/lib)
ELSE(FFTW3_BASE_DIR)
  SET(FFTW3_FOUND 0)
ENDIF(FFTW3_BASE_DIR)

INCLUDE_DIRECTORIES(${FFTW3_INCLUDE_DIR})
LINK_DIRECTORIES(${FFTW3_LIBRARY_DIR})
