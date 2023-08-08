#!/bin/bash
# Please do not modify this script

# Change directory to queso/
cd "queso/"

export APP_SOURCE="${APP_SOURCE:-"$( cd "$(dirname "$0")" ; pwd -P )"}"
export APP_BUILD="${APP_SOURCE}/../build/"

# Set basic configuration
export CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE:-"Release"}
export PYTHON_EXECUTABLE=${PYTHON_EXECUTABLE:-"/usr/bin/python3"}

# Clean
clear
rm -rf "${APP_BUILD}/${CMAKE_BUILD_TYPE}/cmake_install.cmake"
rm -rf "${APP_BUILD}/${CMAKE_BUILD_TYPE}/CMakeCache.txt"
rm -rf "${APP_BUILD}/${CMAKE_BUILD_TYPE}/CMakeFiles"

# Configure
cmake -H"${APP_SOURCE}" -B"${APP_BUILD}/${CMAKE_BUILD_TYPE}" \
-DQUESO_BUILD_TESTING=ON

# Buid
cmake --build "${APP_BUILD}/${CMAKE_BUILD_TYPE}" --target install -- -j2
