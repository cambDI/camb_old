# CMake generated Testfile for 
# Source directory: /Users/daniel/Dropbox/projects/camb/camb/src/build_scripts/indigo-all
# Build directory: /Users/daniel/Dropbox/projects/camb/camb/src/build
# 
# This file includes the relevent testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
ADD_TEST(dlopen-test "/Users/daniel/Dropbox/projects/camb/camb/src/build/dist/Mac/10.7/shared/dlopen-test" "/Users/daniel/Dropbox/projects/camb/camb/src/api/libs/shared/Mac/10.7/libindigo.dylib" "/Users/daniel/Dropbox/projects/camb/camb/src/api/libs/shared/Mac/10.7/libindigo-inchi.dylib" "/Users/daniel/Dropbox/projects/camb/camb/src/api/libs/shared/Mac/10.7/libindigo-renderer.dylib")
SUBDIRS(indigo)
SUBDIRS(indigo-inchi)
SUBDIRS(indigo-renderer)
