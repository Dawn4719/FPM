# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.27

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /opt/clion-2023.3.2/bin/cmake/linux/x64/bin/cmake

# The command to remove a file.
RM = /opt/clion-2023.3.2/bin/cmake/linux/x64/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/dbia/qsl/FINISHED/PSM/nm

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/dbia/qsl/FINISHED/PSM/nm/cmake-build-debug

# Include any dependencies generated for this target.
include matching/CMakeFiles/GraphConverter.out.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include matching/CMakeFiles/GraphConverter.out.dir/compiler_depend.make

# Include the progress variables for this target.
include matching/CMakeFiles/GraphConverter.out.dir/progress.make

# Include the compile flags for this target's objects.
include matching/CMakeFiles/GraphConverter.out.dir/flags.make

matching/CMakeFiles/GraphConverter.out.dir/GraphConverter.cpp.o: matching/CMakeFiles/GraphConverter.out.dir/flags.make
matching/CMakeFiles/GraphConverter.out.dir/GraphConverter.cpp.o: /home/dbia/qsl/FINISHED/PSM/nm/matching/GraphConverter.cpp
matching/CMakeFiles/GraphConverter.out.dir/GraphConverter.cpp.o: matching/CMakeFiles/GraphConverter.out.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/dbia/qsl/FINISHED/PSM/nm/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object matching/CMakeFiles/GraphConverter.out.dir/GraphConverter.cpp.o"
	cd /home/dbia/qsl/FINISHED/PSM/nm/cmake-build-debug/matching && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT matching/CMakeFiles/GraphConverter.out.dir/GraphConverter.cpp.o -MF CMakeFiles/GraphConverter.out.dir/GraphConverter.cpp.o.d -o CMakeFiles/GraphConverter.out.dir/GraphConverter.cpp.o -c /home/dbia/qsl/FINISHED/PSM/nm/matching/GraphConverter.cpp

matching/CMakeFiles/GraphConverter.out.dir/GraphConverter.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/GraphConverter.out.dir/GraphConverter.cpp.i"
	cd /home/dbia/qsl/FINISHED/PSM/nm/cmake-build-debug/matching && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/dbia/qsl/FINISHED/PSM/nm/matching/GraphConverter.cpp > CMakeFiles/GraphConverter.out.dir/GraphConverter.cpp.i

matching/CMakeFiles/GraphConverter.out.dir/GraphConverter.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/GraphConverter.out.dir/GraphConverter.cpp.s"
	cd /home/dbia/qsl/FINISHED/PSM/nm/cmake-build-debug/matching && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/dbia/qsl/FINISHED/PSM/nm/matching/GraphConverter.cpp -o CMakeFiles/GraphConverter.out.dir/GraphConverter.cpp.s

# Object files for target GraphConverter.out
GraphConverter_out_OBJECTS = \
"CMakeFiles/GraphConverter.out.dir/GraphConverter.cpp.o"

# External object files for target GraphConverter.out
GraphConverter_out_EXTERNAL_OBJECTS =

matching/GraphConverter.out: matching/CMakeFiles/GraphConverter.out.dir/GraphConverter.cpp.o
matching/GraphConverter.out: matching/CMakeFiles/GraphConverter.out.dir/build.make
matching/GraphConverter.out: graph/libgraph.so
matching/GraphConverter.out: utility/libutility.so
matching/GraphConverter.out: matching/CMakeFiles/GraphConverter.out.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/home/dbia/qsl/FINISHED/PSM/nm/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable GraphConverter.out"
	cd /home/dbia/qsl/FINISHED/PSM/nm/cmake-build-debug/matching && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/GraphConverter.out.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
matching/CMakeFiles/GraphConverter.out.dir/build: matching/GraphConverter.out
.PHONY : matching/CMakeFiles/GraphConverter.out.dir/build

matching/CMakeFiles/GraphConverter.out.dir/clean:
	cd /home/dbia/qsl/FINISHED/PSM/nm/cmake-build-debug/matching && $(CMAKE_COMMAND) -P CMakeFiles/GraphConverter.out.dir/cmake_clean.cmake
.PHONY : matching/CMakeFiles/GraphConverter.out.dir/clean

matching/CMakeFiles/GraphConverter.out.dir/depend:
	cd /home/dbia/qsl/FINISHED/PSM/nm/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/dbia/qsl/FINISHED/PSM/nm /home/dbia/qsl/FINISHED/PSM/nm/matching /home/dbia/qsl/FINISHED/PSM/nm/cmake-build-debug /home/dbia/qsl/FINISHED/PSM/nm/cmake-build-debug/matching /home/dbia/qsl/FINISHED/PSM/nm/cmake-build-debug/matching/CMakeFiles/GraphConverter.out.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : matching/CMakeFiles/GraphConverter.out.dir/depend

