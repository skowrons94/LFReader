cmake_minimum_required(VERSION 3.0)
project(LFReader VERSION 1.0 LANGUAGES CXX)

set(CMAKE_INSTALL_DIR /usr/bin/LFReader)

set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(HEADER_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)

find_package(ROOT REQUIRED COMPONENTS RIO Net)

include(${ROOT_USE_FILE})

include_directories(
  ${HEADER_DIR}
  ${ROOT_INCLUDE_DIR}
  ${CAENSLC_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
)

add_executable(LFReader
  ${SRC_DIR}/main.cc
  ${SRC_DIR}/LFReader.cc
)

target_compile_features(LFReader PUBLIC cxx_std_14)

target_link_libraries(LFReader PRIVATE ${ROOT_LIBRARIES})

install(TARGETS LFReader
  DESTINATION ${CMAKE_INSTALL_DIR})
