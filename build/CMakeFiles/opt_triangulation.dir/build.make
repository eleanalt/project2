# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/eleanalt/project

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/eleanalt/project/build

# Include any dependencies generated for this target.
include CMakeFiles/opt_triangulation.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/opt_triangulation.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/opt_triangulation.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/opt_triangulation.dir/flags.make

CMakeFiles/opt_triangulation.dir/src/main.cpp.o: CMakeFiles/opt_triangulation.dir/flags.make
CMakeFiles/opt_triangulation.dir/src/main.cpp.o: ../src/main.cpp
CMakeFiles/opt_triangulation.dir/src/main.cpp.o: CMakeFiles/opt_triangulation.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/eleanalt/project/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/opt_triangulation.dir/src/main.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/opt_triangulation.dir/src/main.cpp.o -MF CMakeFiles/opt_triangulation.dir/src/main.cpp.o.d -o CMakeFiles/opt_triangulation.dir/src/main.cpp.o -c /home/eleanalt/project/src/main.cpp

CMakeFiles/opt_triangulation.dir/src/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/opt_triangulation.dir/src/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/eleanalt/project/src/main.cpp > CMakeFiles/opt_triangulation.dir/src/main.cpp.i

CMakeFiles/opt_triangulation.dir/src/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/opt_triangulation.dir/src/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/eleanalt/project/src/main.cpp -o CMakeFiles/opt_triangulation.dir/src/main.cpp.s

CMakeFiles/opt_triangulation.dir/src/TriangulationOptimizer.cpp.o: CMakeFiles/opt_triangulation.dir/flags.make
CMakeFiles/opt_triangulation.dir/src/TriangulationOptimizer.cpp.o: ../src/TriangulationOptimizer.cpp
CMakeFiles/opt_triangulation.dir/src/TriangulationOptimizer.cpp.o: CMakeFiles/opt_triangulation.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/eleanalt/project/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/opt_triangulation.dir/src/TriangulationOptimizer.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/opt_triangulation.dir/src/TriangulationOptimizer.cpp.o -MF CMakeFiles/opt_triangulation.dir/src/TriangulationOptimizer.cpp.o.d -o CMakeFiles/opt_triangulation.dir/src/TriangulationOptimizer.cpp.o -c /home/eleanalt/project/src/TriangulationOptimizer.cpp

CMakeFiles/opt_triangulation.dir/src/TriangulationOptimizer.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/opt_triangulation.dir/src/TriangulationOptimizer.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/eleanalt/project/src/TriangulationOptimizer.cpp > CMakeFiles/opt_triangulation.dir/src/TriangulationOptimizer.cpp.i

CMakeFiles/opt_triangulation.dir/src/TriangulationOptimizer.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/opt_triangulation.dir/src/TriangulationOptimizer.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/eleanalt/project/src/TriangulationOptimizer.cpp -o CMakeFiles/opt_triangulation.dir/src/TriangulationOptimizer.cpp.s

CMakeFiles/opt_triangulation.dir/src/InputParser.cpp.o: CMakeFiles/opt_triangulation.dir/flags.make
CMakeFiles/opt_triangulation.dir/src/InputParser.cpp.o: ../src/InputParser.cpp
CMakeFiles/opt_triangulation.dir/src/InputParser.cpp.o: CMakeFiles/opt_triangulation.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/eleanalt/project/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/opt_triangulation.dir/src/InputParser.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/opt_triangulation.dir/src/InputParser.cpp.o -MF CMakeFiles/opt_triangulation.dir/src/InputParser.cpp.o.d -o CMakeFiles/opt_triangulation.dir/src/InputParser.cpp.o -c /home/eleanalt/project/src/InputParser.cpp

CMakeFiles/opt_triangulation.dir/src/InputParser.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/opt_triangulation.dir/src/InputParser.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/eleanalt/project/src/InputParser.cpp > CMakeFiles/opt_triangulation.dir/src/InputParser.cpp.i

CMakeFiles/opt_triangulation.dir/src/InputParser.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/opt_triangulation.dir/src/InputParser.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/eleanalt/project/src/InputParser.cpp -o CMakeFiles/opt_triangulation.dir/src/InputParser.cpp.s

CMakeFiles/opt_triangulation.dir/src/OutputGenerator.cpp.o: CMakeFiles/opt_triangulation.dir/flags.make
CMakeFiles/opt_triangulation.dir/src/OutputGenerator.cpp.o: ../src/OutputGenerator.cpp
CMakeFiles/opt_triangulation.dir/src/OutputGenerator.cpp.o: CMakeFiles/opt_triangulation.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/eleanalt/project/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/opt_triangulation.dir/src/OutputGenerator.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/opt_triangulation.dir/src/OutputGenerator.cpp.o -MF CMakeFiles/opt_triangulation.dir/src/OutputGenerator.cpp.o.d -o CMakeFiles/opt_triangulation.dir/src/OutputGenerator.cpp.o -c /home/eleanalt/project/src/OutputGenerator.cpp

CMakeFiles/opt_triangulation.dir/src/OutputGenerator.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/opt_triangulation.dir/src/OutputGenerator.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/eleanalt/project/src/OutputGenerator.cpp > CMakeFiles/opt_triangulation.dir/src/OutputGenerator.cpp.i

CMakeFiles/opt_triangulation.dir/src/OutputGenerator.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/opt_triangulation.dir/src/OutputGenerator.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/eleanalt/project/src/OutputGenerator.cpp -o CMakeFiles/opt_triangulation.dir/src/OutputGenerator.cpp.s

# Object files for target opt_triangulation
opt_triangulation_OBJECTS = \
"CMakeFiles/opt_triangulation.dir/src/main.cpp.o" \
"CMakeFiles/opt_triangulation.dir/src/TriangulationOptimizer.cpp.o" \
"CMakeFiles/opt_triangulation.dir/src/InputParser.cpp.o" \
"CMakeFiles/opt_triangulation.dir/src/OutputGenerator.cpp.o"

# External object files for target opt_triangulation
opt_triangulation_EXTERNAL_OBJECTS =

opt_triangulation: CMakeFiles/opt_triangulation.dir/src/main.cpp.o
opt_triangulation: CMakeFiles/opt_triangulation.dir/src/TriangulationOptimizer.cpp.o
opt_triangulation: CMakeFiles/opt_triangulation.dir/src/InputParser.cpp.o
opt_triangulation: CMakeFiles/opt_triangulation.dir/src/OutputGenerator.cpp.o
opt_triangulation: CMakeFiles/opt_triangulation.dir/build.make
opt_triangulation: /usr/lib/x86_64-linux-gnu/libgmpxx.so
opt_triangulation: /usr/lib/x86_64-linux-gnu/libmpfr.so
opt_triangulation: /usr/lib/x86_64-linux-gnu/libgmp.so
opt_triangulation: CMakeFiles/opt_triangulation.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/eleanalt/project/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Linking CXX executable opt_triangulation"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/opt_triangulation.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/opt_triangulation.dir/build: opt_triangulation
.PHONY : CMakeFiles/opt_triangulation.dir/build

CMakeFiles/opt_triangulation.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/opt_triangulation.dir/cmake_clean.cmake
.PHONY : CMakeFiles/opt_triangulation.dir/clean

CMakeFiles/opt_triangulation.dir/depend:
	cd /home/eleanalt/project/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/eleanalt/project /home/eleanalt/project /home/eleanalt/project/build /home/eleanalt/project/build /home/eleanalt/project/build/CMakeFiles/opt_triangulation.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/opt_triangulation.dir/depend

