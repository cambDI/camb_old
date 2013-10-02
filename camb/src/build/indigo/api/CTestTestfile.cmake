# CMake generated Testfile for 
# Source directory: /Users/daniel/Dropbox/projects/camb/camb/src/api
# Build directory: /Users/daniel/Dropbox/projects/camb/camb/src/build/indigo/api
# 
# This file includes the relevent testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
ADD_TEST(indigo-c-test-static "/Users/daniel/Dropbox/projects/camb/camb/src/build/dist/Mac/10.7/shared/indigo-c-test-static")
ADD_TEST(indigo-c-test-shared "/Users/daniel/Dropbox/projects/camb/camb/src/build/dist/Mac/10.7/shared/indigo-c-test-shared")
ADD_TEST(dlopen-indigo-test "/Users/daniel/Dropbox/projects/camb/camb/src/build/dist/Mac/10.7/shared/dlopen-test" "/Users/daniel/Dropbox/projects/camb/camb/src/api/libs/shared/Mac/10.7/libindigo.dylib")
