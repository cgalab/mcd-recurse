project (mcd-recurse VERSION "0.1" LANGUAGES "C")

set(_target "mcd-recurse")

set(CMAKE_INCLUDE_CURRENT_DIR ON)

add_definitions(-DELPP_CUSTOM_COUT=std::cerr)
set_source_files_properties(easyloggingpp/src/easylogging++.cc PROPERTIES COMPILE_FLAGS "-Wno-float-equal -Wno-missing-noreturn -Wno-unknown-pragmas -Wno-shadow")
add_definitions(-DTRILIBRARY)
add_definitions(-DANSI_DECLARATORS)

add_executable("${_target}"
                       easyloggingpp/src/easylogging++.cc
                       io.cpp
                       geom.cpp
                       main.cpp
                       triangle/triangle.c
                       "${PROJECT_BINARY_DIR}/gitversion.h"
                       )

add_custom_command(
  OUTPUT "${PROJECT_BINARY_DIR}/gitversion.h"
  DEPENDS "${CMAKE_SOURCE_DIR}/.git/logs/HEAD"
  COMMAND ${CMAKE_COMMAND}
    -Dlocal_dir="${CMAKE_SOURCE_DIR}"
    -Doutput_dir="${PROJECT_BINARY_DIR}"
    -P "${CMAKE_SOURCE_DIR}/cmake/gitversion.cmake"
)
