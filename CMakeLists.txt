cmake_minimum_required(VERSION 3.21)

# Only set the cxx_standard if it is not set by someone else
if (NOT DEFINED CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 17)
endif()

# Set the project name and language
project(
  QuESo
  VERSION 0.0.1
  DESCRIPTION "Qadrature for Embedded SOlids"
  LANGUAGES CXX C
)

add_subdirectory(queso)