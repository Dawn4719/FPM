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
include matching/CMakeFiles/SubgraphMatching.out.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include matching/CMakeFiles/SubgraphMatching.out.dir/compiler_depend.make

# Include the progress variables for this target.
include matching/CMakeFiles/SubgraphMatching.out.dir/progress.make

# Include the compile flags for this target's objects.
include matching/CMakeFiles/SubgraphMatching.out.dir/flags.make

matching/CMakeFiles/SubgraphMatching.out.dir/matchingcommand.cpp.o: matching/CMakeFiles/SubgraphMatching.out.dir/flags.make
matching/CMakeFiles/SubgraphMatching.out.dir/matchingcommand.cpp.o: /home/dbia/qsl/FINISHED/PSM/nm/matching/matchingcommand.cpp
matching/CMakeFiles/SubgraphMatching.out.dir/matchingcommand.cpp.o: matching/CMakeFiles/SubgraphMatching.out.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/dbia/qsl/FINISHED/PSM/nm/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object matching/CMakeFiles/SubgraphMatching.out.dir/matchingcommand.cpp.o"
	cd /home/dbia/qsl/FINISHED/PSM/nm/cmake-build-debug/matching && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT matching/CMakeFiles/SubgraphMatching.out.dir/matchingcommand.cpp.o -MF CMakeFiles/SubgraphMatching.out.dir/matchingcommand.cpp.o.d -o CMakeFiles/SubgraphMatching.out.dir/matchingcommand.cpp.o -c /home/dbia/qsl/FINISHED/PSM/nm/matching/matchingcommand.cpp

matching/CMakeFiles/SubgraphMatching.out.dir/matchingcommand.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/SubgraphMatching.out.dir/matchingcommand.cpp.i"
	cd /home/dbia/qsl/FINISHED/PSM/nm/cmake-build-debug/matching && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/dbia/qsl/FINISHED/PSM/nm/matching/matchingcommand.cpp > CMakeFiles/SubgraphMatching.out.dir/matchingcommand.cpp.i

matching/CMakeFiles/SubgraphMatching.out.dir/matchingcommand.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/SubgraphMatching.out.dir/matchingcommand.cpp.s"
	cd /home/dbia/qsl/FINISHED/PSM/nm/cmake-build-debug/matching && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/dbia/qsl/FINISHED/PSM/nm/matching/matchingcommand.cpp -o CMakeFiles/SubgraphMatching.out.dir/matchingcommand.cpp.s

matching/CMakeFiles/SubgraphMatching.out.dir/FilterVertices.cpp.o: matching/CMakeFiles/SubgraphMatching.out.dir/flags.make
matching/CMakeFiles/SubgraphMatching.out.dir/FilterVertices.cpp.o: /home/dbia/qsl/FINISHED/PSM/nm/matching/FilterVertices.cpp
matching/CMakeFiles/SubgraphMatching.out.dir/FilterVertices.cpp.o: matching/CMakeFiles/SubgraphMatching.out.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/dbia/qsl/FINISHED/PSM/nm/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object matching/CMakeFiles/SubgraphMatching.out.dir/FilterVertices.cpp.o"
	cd /home/dbia/qsl/FINISHED/PSM/nm/cmake-build-debug/matching && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT matching/CMakeFiles/SubgraphMatching.out.dir/FilterVertices.cpp.o -MF CMakeFiles/SubgraphMatching.out.dir/FilterVertices.cpp.o.d -o CMakeFiles/SubgraphMatching.out.dir/FilterVertices.cpp.o -c /home/dbia/qsl/FINISHED/PSM/nm/matching/FilterVertices.cpp

matching/CMakeFiles/SubgraphMatching.out.dir/FilterVertices.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/SubgraphMatching.out.dir/FilterVertices.cpp.i"
	cd /home/dbia/qsl/FINISHED/PSM/nm/cmake-build-debug/matching && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/dbia/qsl/FINISHED/PSM/nm/matching/FilterVertices.cpp > CMakeFiles/SubgraphMatching.out.dir/FilterVertices.cpp.i

matching/CMakeFiles/SubgraphMatching.out.dir/FilterVertices.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/SubgraphMatching.out.dir/FilterVertices.cpp.s"
	cd /home/dbia/qsl/FINISHED/PSM/nm/cmake-build-debug/matching && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/dbia/qsl/FINISHED/PSM/nm/matching/FilterVertices.cpp -o CMakeFiles/SubgraphMatching.out.dir/FilterVertices.cpp.s

matching/CMakeFiles/SubgraphMatching.out.dir/BuildTable.cpp.o: matching/CMakeFiles/SubgraphMatching.out.dir/flags.make
matching/CMakeFiles/SubgraphMatching.out.dir/BuildTable.cpp.o: /home/dbia/qsl/FINISHED/PSM/nm/matching/BuildTable.cpp
matching/CMakeFiles/SubgraphMatching.out.dir/BuildTable.cpp.o: matching/CMakeFiles/SubgraphMatching.out.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/dbia/qsl/FINISHED/PSM/nm/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object matching/CMakeFiles/SubgraphMatching.out.dir/BuildTable.cpp.o"
	cd /home/dbia/qsl/FINISHED/PSM/nm/cmake-build-debug/matching && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT matching/CMakeFiles/SubgraphMatching.out.dir/BuildTable.cpp.o -MF CMakeFiles/SubgraphMatching.out.dir/BuildTable.cpp.o.d -o CMakeFiles/SubgraphMatching.out.dir/BuildTable.cpp.o -c /home/dbia/qsl/FINISHED/PSM/nm/matching/BuildTable.cpp

matching/CMakeFiles/SubgraphMatching.out.dir/BuildTable.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/SubgraphMatching.out.dir/BuildTable.cpp.i"
	cd /home/dbia/qsl/FINISHED/PSM/nm/cmake-build-debug/matching && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/dbia/qsl/FINISHED/PSM/nm/matching/BuildTable.cpp > CMakeFiles/SubgraphMatching.out.dir/BuildTable.cpp.i

matching/CMakeFiles/SubgraphMatching.out.dir/BuildTable.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/SubgraphMatching.out.dir/BuildTable.cpp.s"
	cd /home/dbia/qsl/FINISHED/PSM/nm/cmake-build-debug/matching && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/dbia/qsl/FINISHED/PSM/nm/matching/BuildTable.cpp -o CMakeFiles/SubgraphMatching.out.dir/BuildTable.cpp.s

matching/CMakeFiles/SubgraphMatching.out.dir/GenerateQueryPlan.cpp.o: matching/CMakeFiles/SubgraphMatching.out.dir/flags.make
matching/CMakeFiles/SubgraphMatching.out.dir/GenerateQueryPlan.cpp.o: /home/dbia/qsl/FINISHED/PSM/nm/matching/GenerateQueryPlan.cpp
matching/CMakeFiles/SubgraphMatching.out.dir/GenerateQueryPlan.cpp.o: matching/CMakeFiles/SubgraphMatching.out.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/dbia/qsl/FINISHED/PSM/nm/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object matching/CMakeFiles/SubgraphMatching.out.dir/GenerateQueryPlan.cpp.o"
	cd /home/dbia/qsl/FINISHED/PSM/nm/cmake-build-debug/matching && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT matching/CMakeFiles/SubgraphMatching.out.dir/GenerateQueryPlan.cpp.o -MF CMakeFiles/SubgraphMatching.out.dir/GenerateQueryPlan.cpp.o.d -o CMakeFiles/SubgraphMatching.out.dir/GenerateQueryPlan.cpp.o -c /home/dbia/qsl/FINISHED/PSM/nm/matching/GenerateQueryPlan.cpp

matching/CMakeFiles/SubgraphMatching.out.dir/GenerateQueryPlan.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/SubgraphMatching.out.dir/GenerateQueryPlan.cpp.i"
	cd /home/dbia/qsl/FINISHED/PSM/nm/cmake-build-debug/matching && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/dbia/qsl/FINISHED/PSM/nm/matching/GenerateQueryPlan.cpp > CMakeFiles/SubgraphMatching.out.dir/GenerateQueryPlan.cpp.i

matching/CMakeFiles/SubgraphMatching.out.dir/GenerateQueryPlan.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/SubgraphMatching.out.dir/GenerateQueryPlan.cpp.s"
	cd /home/dbia/qsl/FINISHED/PSM/nm/cmake-build-debug/matching && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/dbia/qsl/FINISHED/PSM/nm/matching/GenerateQueryPlan.cpp -o CMakeFiles/SubgraphMatching.out.dir/GenerateQueryPlan.cpp.s

matching/CMakeFiles/SubgraphMatching.out.dir/EvaluateQuery.cpp.o: matching/CMakeFiles/SubgraphMatching.out.dir/flags.make
matching/CMakeFiles/SubgraphMatching.out.dir/EvaluateQuery.cpp.o: /home/dbia/qsl/FINISHED/PSM/nm/matching/EvaluateQuery.cpp
matching/CMakeFiles/SubgraphMatching.out.dir/EvaluateQuery.cpp.o: matching/CMakeFiles/SubgraphMatching.out.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/dbia/qsl/FINISHED/PSM/nm/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object matching/CMakeFiles/SubgraphMatching.out.dir/EvaluateQuery.cpp.o"
	cd /home/dbia/qsl/FINISHED/PSM/nm/cmake-build-debug/matching && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT matching/CMakeFiles/SubgraphMatching.out.dir/EvaluateQuery.cpp.o -MF CMakeFiles/SubgraphMatching.out.dir/EvaluateQuery.cpp.o.d -o CMakeFiles/SubgraphMatching.out.dir/EvaluateQuery.cpp.o -c /home/dbia/qsl/FINISHED/PSM/nm/matching/EvaluateQuery.cpp

matching/CMakeFiles/SubgraphMatching.out.dir/EvaluateQuery.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/SubgraphMatching.out.dir/EvaluateQuery.cpp.i"
	cd /home/dbia/qsl/FINISHED/PSM/nm/cmake-build-debug/matching && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/dbia/qsl/FINISHED/PSM/nm/matching/EvaluateQuery.cpp > CMakeFiles/SubgraphMatching.out.dir/EvaluateQuery.cpp.i

matching/CMakeFiles/SubgraphMatching.out.dir/EvaluateQuery.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/SubgraphMatching.out.dir/EvaluateQuery.cpp.s"
	cd /home/dbia/qsl/FINISHED/PSM/nm/cmake-build-debug/matching && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/dbia/qsl/FINISHED/PSM/nm/matching/EvaluateQuery.cpp -o CMakeFiles/SubgraphMatching.out.dir/EvaluateQuery.cpp.s

matching/CMakeFiles/SubgraphMatching.out.dir/GenerateFilteringPlan.cpp.o: matching/CMakeFiles/SubgraphMatching.out.dir/flags.make
matching/CMakeFiles/SubgraphMatching.out.dir/GenerateFilteringPlan.cpp.o: /home/dbia/qsl/FINISHED/PSM/nm/matching/GenerateFilteringPlan.cpp
matching/CMakeFiles/SubgraphMatching.out.dir/GenerateFilteringPlan.cpp.o: matching/CMakeFiles/SubgraphMatching.out.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/dbia/qsl/FINISHED/PSM/nm/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object matching/CMakeFiles/SubgraphMatching.out.dir/GenerateFilteringPlan.cpp.o"
	cd /home/dbia/qsl/FINISHED/PSM/nm/cmake-build-debug/matching && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT matching/CMakeFiles/SubgraphMatching.out.dir/GenerateFilteringPlan.cpp.o -MF CMakeFiles/SubgraphMatching.out.dir/GenerateFilteringPlan.cpp.o.d -o CMakeFiles/SubgraphMatching.out.dir/GenerateFilteringPlan.cpp.o -c /home/dbia/qsl/FINISHED/PSM/nm/matching/GenerateFilteringPlan.cpp

matching/CMakeFiles/SubgraphMatching.out.dir/GenerateFilteringPlan.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/SubgraphMatching.out.dir/GenerateFilteringPlan.cpp.i"
	cd /home/dbia/qsl/FINISHED/PSM/nm/cmake-build-debug/matching && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/dbia/qsl/FINISHED/PSM/nm/matching/GenerateFilteringPlan.cpp > CMakeFiles/SubgraphMatching.out.dir/GenerateFilteringPlan.cpp.i

matching/CMakeFiles/SubgraphMatching.out.dir/GenerateFilteringPlan.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/SubgraphMatching.out.dir/GenerateFilteringPlan.cpp.s"
	cd /home/dbia/qsl/FINISHED/PSM/nm/cmake-build-debug/matching && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/dbia/qsl/FINISHED/PSM/nm/matching/GenerateFilteringPlan.cpp -o CMakeFiles/SubgraphMatching.out.dir/GenerateFilteringPlan.cpp.s

matching/CMakeFiles/SubgraphMatching.out.dir/StudyPerformance.cpp.o: matching/CMakeFiles/SubgraphMatching.out.dir/flags.make
matching/CMakeFiles/SubgraphMatching.out.dir/StudyPerformance.cpp.o: /home/dbia/qsl/FINISHED/PSM/nm/matching/StudyPerformance.cpp
matching/CMakeFiles/SubgraphMatching.out.dir/StudyPerformance.cpp.o: matching/CMakeFiles/SubgraphMatching.out.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/dbia/qsl/FINISHED/PSM/nm/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object matching/CMakeFiles/SubgraphMatching.out.dir/StudyPerformance.cpp.o"
	cd /home/dbia/qsl/FINISHED/PSM/nm/cmake-build-debug/matching && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT matching/CMakeFiles/SubgraphMatching.out.dir/StudyPerformance.cpp.o -MF CMakeFiles/SubgraphMatching.out.dir/StudyPerformance.cpp.o.d -o CMakeFiles/SubgraphMatching.out.dir/StudyPerformance.cpp.o -c /home/dbia/qsl/FINISHED/PSM/nm/matching/StudyPerformance.cpp

matching/CMakeFiles/SubgraphMatching.out.dir/StudyPerformance.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/SubgraphMatching.out.dir/StudyPerformance.cpp.i"
	cd /home/dbia/qsl/FINISHED/PSM/nm/cmake-build-debug/matching && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/dbia/qsl/FINISHED/PSM/nm/matching/StudyPerformance.cpp > CMakeFiles/SubgraphMatching.out.dir/StudyPerformance.cpp.i

matching/CMakeFiles/SubgraphMatching.out.dir/StudyPerformance.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/SubgraphMatching.out.dir/StudyPerformance.cpp.s"
	cd /home/dbia/qsl/FINISHED/PSM/nm/cmake-build-debug/matching && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/dbia/qsl/FINISHED/PSM/nm/matching/StudyPerformance.cpp -o CMakeFiles/SubgraphMatching.out.dir/StudyPerformance.cpp.s

# Object files for target SubgraphMatching.out
SubgraphMatching_out_OBJECTS = \
"CMakeFiles/SubgraphMatching.out.dir/matchingcommand.cpp.o" \
"CMakeFiles/SubgraphMatching.out.dir/FilterVertices.cpp.o" \
"CMakeFiles/SubgraphMatching.out.dir/BuildTable.cpp.o" \
"CMakeFiles/SubgraphMatching.out.dir/GenerateQueryPlan.cpp.o" \
"CMakeFiles/SubgraphMatching.out.dir/EvaluateQuery.cpp.o" \
"CMakeFiles/SubgraphMatching.out.dir/GenerateFilteringPlan.cpp.o" \
"CMakeFiles/SubgraphMatching.out.dir/StudyPerformance.cpp.o"

# External object files for target SubgraphMatching.out
SubgraphMatching_out_EXTERNAL_OBJECTS =

matching/SubgraphMatching.out: matching/CMakeFiles/SubgraphMatching.out.dir/matchingcommand.cpp.o
matching/SubgraphMatching.out: matching/CMakeFiles/SubgraphMatching.out.dir/FilterVertices.cpp.o
matching/SubgraphMatching.out: matching/CMakeFiles/SubgraphMatching.out.dir/BuildTable.cpp.o
matching/SubgraphMatching.out: matching/CMakeFiles/SubgraphMatching.out.dir/GenerateQueryPlan.cpp.o
matching/SubgraphMatching.out: matching/CMakeFiles/SubgraphMatching.out.dir/EvaluateQuery.cpp.o
matching/SubgraphMatching.out: matching/CMakeFiles/SubgraphMatching.out.dir/GenerateFilteringPlan.cpp.o
matching/SubgraphMatching.out: matching/CMakeFiles/SubgraphMatching.out.dir/StudyPerformance.cpp.o
matching/SubgraphMatching.out: matching/CMakeFiles/SubgraphMatching.out.dir/build.make
matching/SubgraphMatching.out: /usr/lib/x86_64-linux-gnu/libtbb.so.2
matching/SubgraphMatching.out: graph/libgraph.so
matching/SubgraphMatching.out: utility/libutility.so
matching/SubgraphMatching.out: matching/CMakeFiles/SubgraphMatching.out.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/home/dbia/qsl/FINISHED/PSM/nm/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Linking CXX executable SubgraphMatching.out"
	cd /home/dbia/qsl/FINISHED/PSM/nm/cmake-build-debug/matching && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/SubgraphMatching.out.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
matching/CMakeFiles/SubgraphMatching.out.dir/build: matching/SubgraphMatching.out
.PHONY : matching/CMakeFiles/SubgraphMatching.out.dir/build

matching/CMakeFiles/SubgraphMatching.out.dir/clean:
	cd /home/dbia/qsl/FINISHED/PSM/nm/cmake-build-debug/matching && $(CMAKE_COMMAND) -P CMakeFiles/SubgraphMatching.out.dir/cmake_clean.cmake
.PHONY : matching/CMakeFiles/SubgraphMatching.out.dir/clean

matching/CMakeFiles/SubgraphMatching.out.dir/depend:
	cd /home/dbia/qsl/FINISHED/PSM/nm/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/dbia/qsl/FINISHED/PSM/nm /home/dbia/qsl/FINISHED/PSM/nm/matching /home/dbia/qsl/FINISHED/PSM/nm/cmake-build-debug /home/dbia/qsl/FINISHED/PSM/nm/cmake-build-debug/matching /home/dbia/qsl/FINISHED/PSM/nm/cmake-build-debug/matching/CMakeFiles/SubgraphMatching.out.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : matching/CMakeFiles/SubgraphMatching.out.dir/depend
