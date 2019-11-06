project (mcd-recurse VERSION "0.1" DESCRIPTION "mcd-recurse")

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set_source_files_properties(easyloggingpp/src/easylogging++.cc PROPERTIES COMPILE_FLAGS "-Wno-float-equal -Wno-missing-noreturn -Wno-unknown-pragmas -Wno-shadow")

add_executable(mcd-recurse
                       easyloggingpp/src/easylogging++.cc
                       io.cpp
                       geom.cpp
                       main.cpp
                       )

configure_file (
  "${PROJECT_SOURCE_DIR}/config.h.in"
  "${PROJECT_BINARY_DIR}/config.h"
  @ONLY
)
