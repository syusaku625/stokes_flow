# https://qiita.com/shohirose/items/d9bda00a39a113965c5c

find_path(GLOG_INCLUDE_DIR glog/logging.h
  PATHS
    /usr
    /usr/local
    ${GLOG_DIR}
    $ENV{GLOG_DIR}
  PATH_SUFFIXES
    include
)

find_library(GLOG_LIBRARY 
NAMES glog
PATHS
    /usr
    /usr/local
    ${GLOG_DIR}
    $ENV{GLOG_DIR}
  PATH_SUFFIXES
    lib
) 

mark_as_advanced(
  GLOG_INCLUDE_DIR
  GLOG_LIBRARY     # ヘッダーのみのライブラリの場合は不要
  )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GLOG
  REQUIRED_VARS
    GLOG_INCLUDE_DIR
    GLOG_LIBRARY      # ヘッダーのみのライブラリの場合は不要
  )

if(GLOG_FOUND AND NOT TARGET glog)
  add_library(glog UNKNOWN IMPORTED)
  set_target_properties(glog PROPERTIES
    IMPORTED_LINK_INTERFACE_LANGUAGES ["C"|"CXX"]  # ヘッダーのみのライブラリの場合は不要
    IMPORTED_LOCATION "${GLOG_LIBRARY}"       # ヘッダーのみのライブラリの場合は不要
    INTERFACE_INCLUDE_DIRECTORIES "${GLOG_INCLUDE_DIR}"
    )
endif()