# Install script for directory: /Users/daniel/Dropbox/projects/camb/camb/src/build_scripts/indigo-renderer

# Set the install prefix
IF(NOT DEFINED CMAKE_INSTALL_PREFIX)
  SET(CMAKE_INSTALL_PREFIX "/usr/local")
ENDIF(NOT DEFINED CMAKE_INSTALL_PREFIX)
STRING(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
IF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  IF(BUILD_TYPE)
    STRING(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  ELSE(BUILD_TYPE)
    SET(CMAKE_INSTALL_CONFIG_NAME "Release")
  ENDIF(BUILD_TYPE)
  MESSAGE(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
ENDIF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)

# Set the component getting installed.
IF(NOT CMAKE_INSTALL_COMPONENT)
  IF(COMPONENT)
    MESSAGE(STATUS "Install component: \"${COMPONENT}\"")
    SET(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  ELSE(COMPONENT)
    SET(CMAKE_INSTALL_COMPONENT)
  ENDIF(COMPONENT)
ENDIF(NOT CMAKE_INSTALL_COMPONENT)

IF(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  INCLUDE("/Users/daniel/Dropbox/projects/camb/camb/src/build/indigo-renderer/indigo/cmake_install.cmake")
  INCLUDE("/Users/daniel/Dropbox/projects/camb/camb/src/build/indigo-renderer/png/cmake_install.cmake")
  INCLUDE("/Users/daniel/Dropbox/projects/camb/camb/src/build/indigo-renderer/pixman/cmake_install.cmake")
  INCLUDE("/Users/daniel/Dropbox/projects/camb/camb/src/build/indigo-renderer/cairo/cmake_install.cmake")
  INCLUDE("/Users/daniel/Dropbox/projects/camb/camb/src/build/indigo-renderer/render2d/cmake_install.cmake")
  INCLUDE("/Users/daniel/Dropbox/projects/camb/camb/src/build/indigo-renderer/renderer/cmake_install.cmake")

ENDIF(NOT CMAKE_INSTALL_LOCAL_ONLY)

