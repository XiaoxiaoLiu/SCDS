# sets the following variables
# FFTW_FOUND         - Set to true when fftw is found.
# FFTW_INCLUDE_DIR   - Include directory for fftw headers.
# FFTW_LIBRARY_DIR   - Link directory for fftw libraries.

FIND_PATH(FFTW_BASE_DIR include/sfftw.h
  ${RADONC_LIBS_PATH}/fftw
  )

IF(FFTW_BASE_DIR)
  SET(FFTW_FOUND 1)
  SET(FFTW_INCLUDE_DIR ${FFTW_BASE_DIR}/include)
  SET(FFTW_LIBRARY_DIR ${FFTW_BASE_DIR}/lib)
ELSE(FFTW_BASE_DIR)
  SET(FFTW_FOUND 0)
ENDIF(FFTW_BASE_DIR)

INCLUDE_DIRECTORIES(${FFTW_INCLUDE_DIR})
LINK_DIRECTORIES(${FFTW_LIBRARY_DIR})
