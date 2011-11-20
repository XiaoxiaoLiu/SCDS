# sets the following variables
# LAPACK_FOUND         - Set to true when fftw is found.
# LAPACK_INCLUDE_DIR   - Include directory for fftw headers.
# LAPACK_LIBRARY_DIR   - Link directory for fftw libraries.
# LAPACK_LIBRARIES     - list of lapack libraries

FIND_PATH(LAPACK_BASE_DIR lib
  ${RADONC_LIBS_PATH}/LAPACK
  )

IF(LAPACK_BASE_DIR)
  SET(LAPACK_FOUND 1)
  SET(LAPACK_INCLUDE_DIR ${LAPACK_BASE_DIR})
  SET(LAPACK_LIBRARY_DIR ${LAPACK_BASE_DIR})
  
  IF(WIN32)
    SET(LAPACK_LIBRARIES
      lapack
      )
  ENDIF(WIN32)
  IF(CMAKE_SYSTEM MATCHES "Darwin*")
    SET(LAPACK_LIBRARIES
      lapack
      blas
      I77
      F77
      )
  ENDIF(CMAKE_SYSTEM MATCHES "Darwin*")
  IF(CMAKE_SYSTEM MATCHES "Sun*")
    SET(LAPACK_LIBRARIES
      lapack
      blas
      g2c
      )
  ENDIF(CMAKE_SYSTEM MATCHES "Sun*")
  IF(CMAKE_SYSTEM MATCHES "Linux*")
    SET(LAPACK_LIBRARIES
      lapack
      blas
      f2c
      )
  ENDIF(CMAKE_SYSTEM MATCHES "Linux*")
  
ELSE(LAPACK_BASE_DIR)
  SET(LAPACK_FOUND 0)
ENDIF(LAPACK_BASE_DIR)

INCLUDE_DIRECTORIES(${LAPACK_INCLUDE_DIR})
LINK_DIRECTORIES(${LAPACK_LIBRARY_DIR})
