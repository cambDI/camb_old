# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/Cellar/cmake/2.8.10.1/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/2.8.10.1/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/local/Cellar/cmake/2.8.10.1/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/daniel/Dropbox/projects/camb/camb/src/build_scripts/indigo-all

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/daniel/Dropbox/projects/camb/camb/src/build

# Include any dependencies generated for this target.
include indigo/api/CMakeFiles/indigo-c-test-shared.dir/depend.make

# Include the progress variables for this target.
include indigo/api/CMakeFiles/indigo-c-test-shared.dir/progress.make

# Include the compile flags for this target's objects.
include indigo/api/CMakeFiles/indigo-c-test-shared.dir/flags.make

indigo/api/CMakeFiles/indigo-c-test-shared.dir/tests/c/indigo-test.c.o: indigo/api/CMakeFiles/indigo-c-test-shared.dir/flags.make
indigo/api/CMakeFiles/indigo-c-test-shared.dir/tests/c/indigo-test.c.o: /Users/daniel/Dropbox/projects/camb/camb/src/api/tests/c/indigo-test.c
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/daniel/Dropbox/projects/camb/camb/src/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object indigo/api/CMakeFiles/indigo-c-test-shared.dir/tests/c/indigo-test.c.o"
	cd /Users/daniel/Dropbox/projects/camb/camb/src/build/indigo/api && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/indigo-c-test-shared.dir/tests/c/indigo-test.c.o   -c /Users/daniel/Dropbox/projects/camb/camb/src/api/tests/c/indigo-test.c

indigo/api/CMakeFiles/indigo-c-test-shared.dir/tests/c/indigo-test.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/indigo-c-test-shared.dir/tests/c/indigo-test.c.i"
	cd /Users/daniel/Dropbox/projects/camb/camb/src/build/indigo/api && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /Users/daniel/Dropbox/projects/camb/camb/src/api/tests/c/indigo-test.c > CMakeFiles/indigo-c-test-shared.dir/tests/c/indigo-test.c.i

indigo/api/CMakeFiles/indigo-c-test-shared.dir/tests/c/indigo-test.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/indigo-c-test-shared.dir/tests/c/indigo-test.c.s"
	cd /Users/daniel/Dropbox/projects/camb/camb/src/build/indigo/api && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /Users/daniel/Dropbox/projects/camb/camb/src/api/tests/c/indigo-test.c -o CMakeFiles/indigo-c-test-shared.dir/tests/c/indigo-test.c.s

indigo/api/CMakeFiles/indigo-c-test-shared.dir/tests/c/indigo-test.c.o.requires:
.PHONY : indigo/api/CMakeFiles/indigo-c-test-shared.dir/tests/c/indigo-test.c.o.requires

indigo/api/CMakeFiles/indigo-c-test-shared.dir/tests/c/indigo-test.c.o.provides: indigo/api/CMakeFiles/indigo-c-test-shared.dir/tests/c/indigo-test.c.o.requires
	$(MAKE) -f indigo/api/CMakeFiles/indigo-c-test-shared.dir/build.make indigo/api/CMakeFiles/indigo-c-test-shared.dir/tests/c/indigo-test.c.o.provides.build
.PHONY : indigo/api/CMakeFiles/indigo-c-test-shared.dir/tests/c/indigo-test.c.o.provides

indigo/api/CMakeFiles/indigo-c-test-shared.dir/tests/c/indigo-test.c.o.provides.build: indigo/api/CMakeFiles/indigo-c-test-shared.dir/tests/c/indigo-test.c.o

# Object files for target indigo-c-test-shared
indigo__c__test__shared_OBJECTS = \
"CMakeFiles/indigo-c-test-shared.dir/tests/c/indigo-test.c.o"

# External object files for target indigo-c-test-shared
indigo__c__test__shared_EXTERNAL_OBJECTS =

dist/Mac/10.7/shared/indigo-c-test-shared: indigo/api/CMakeFiles/indigo-c-test-shared.dir/tests/c/indigo-test.c.o
dist/Mac/10.7/shared/indigo-c-test-shared: indigo/api/CMakeFiles/indigo-c-test-shared.dir/build.make
dist/Mac/10.7/shared/indigo-c-test-shared: dist/Mac/10.7/lib/libindigo.dylib
dist/Mac/10.7/shared/indigo-c-test-shared: indigo/api/CMakeFiles/indigo-c-test-shared.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking C executable ../../dist/Mac/10.7/shared/indigo-c-test-shared"
	cd /Users/daniel/Dropbox/projects/camb/camb/src/build/indigo/api && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/indigo-c-test-shared.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
indigo/api/CMakeFiles/indigo-c-test-shared.dir/build: dist/Mac/10.7/shared/indigo-c-test-shared
.PHONY : indigo/api/CMakeFiles/indigo-c-test-shared.dir/build

indigo/api/CMakeFiles/indigo-c-test-shared.dir/requires: indigo/api/CMakeFiles/indigo-c-test-shared.dir/tests/c/indigo-test.c.o.requires
.PHONY : indigo/api/CMakeFiles/indigo-c-test-shared.dir/requires

indigo/api/CMakeFiles/indigo-c-test-shared.dir/clean:
	cd /Users/daniel/Dropbox/projects/camb/camb/src/build/indigo/api && $(CMAKE_COMMAND) -P CMakeFiles/indigo-c-test-shared.dir/cmake_clean.cmake
.PHONY : indigo/api/CMakeFiles/indigo-c-test-shared.dir/clean

indigo/api/CMakeFiles/indigo-c-test-shared.dir/depend:
	cd /Users/daniel/Dropbox/projects/camb/camb/src/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/daniel/Dropbox/projects/camb/camb/src/build_scripts/indigo-all /Users/daniel/Dropbox/projects/camb/camb/src/api /Users/daniel/Dropbox/projects/camb/camb/src/build /Users/daniel/Dropbox/projects/camb/camb/src/build/indigo/api /Users/daniel/Dropbox/projects/camb/camb/src/build/indigo/api/CMakeFiles/indigo-c-test-shared.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : indigo/api/CMakeFiles/indigo-c-test-shared.dir/depend

