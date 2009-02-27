FIND_PATH(FFTW3_INCLUDE fftw3.h /usr/include)

FIND_LIBRARY(FFTW3_LIB fftw3 /usr/lib)

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(FFTW3 DEFAULT_MSG FFTW3_INCLUDE FFTW3_LIB)

