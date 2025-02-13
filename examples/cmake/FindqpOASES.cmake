find_path(qpOASES_INCLUDE_DIR NAMES qpOASES.hpp PATHS /usr/include /usr/local/include)

find_library(qpOASES_LIBRARY NAMES qpOASES PATHS /usr/lib /usr/local/lib)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(qpOASES REQUIRED_VARS qpOASES_INCLUDE_DIR qpOASES_LIBRARY)

if(qpOASES_FOUND)
  set(qpOASES_INCLUDE_DIRS ${qpOASES_INCLUDE_DIR})
  set(qpOASES_LIBRARIES ${qpOASES_LIBRARY})
endif()