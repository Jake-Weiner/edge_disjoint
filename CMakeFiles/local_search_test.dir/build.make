# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.12

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


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
CMAKE_COMMAND = /usr/local/bin/cmake

# The command to remove a file.
RM = /usr/local/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/jake/PhD/Edge_Disjoint/c++/Lagrangian_Relax/git

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/jake/PhD/Edge_Disjoint/c++/Lagrangian_Relax/git

# Include any dependencies generated for this target.
include CMakeFiles/local_search_test.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/local_search_test.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/local_search_test.dir/flags.make

CMakeFiles/local_search_test.dir/main_test.cpp.o: CMakeFiles/local_search_test.dir/flags.make
CMakeFiles/local_search_test.dir/main_test.cpp.o: main_test.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jake/PhD/Edge_Disjoint/c++/Lagrangian_Relax/git/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/local_search_test.dir/main_test.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/local_search_test.dir/main_test.cpp.o -c /home/jake/PhD/Edge_Disjoint/c++/Lagrangian_Relax/git/main_test.cpp

CMakeFiles/local_search_test.dir/main_test.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/local_search_test.dir/main_test.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jake/PhD/Edge_Disjoint/c++/Lagrangian_Relax/git/main_test.cpp > CMakeFiles/local_search_test.dir/main_test.cpp.i

CMakeFiles/local_search_test.dir/main_test.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/local_search_test.dir/main_test.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jake/PhD/Edge_Disjoint/c++/Lagrangian_Relax/git/main_test.cpp -o CMakeFiles/local_search_test.dir/main_test.cpp.s

CMakeFiles/local_search_test.dir/ED.cpp.o: CMakeFiles/local_search_test.dir/flags.make
CMakeFiles/local_search_test.dir/ED.cpp.o: ED.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jake/PhD/Edge_Disjoint/c++/Lagrangian_Relax/git/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/local_search_test.dir/ED.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/local_search_test.dir/ED.cpp.o -c /home/jake/PhD/Edge_Disjoint/c++/Lagrangian_Relax/git/ED.cpp

CMakeFiles/local_search_test.dir/ED.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/local_search_test.dir/ED.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jake/PhD/Edge_Disjoint/c++/Lagrangian_Relax/git/ED.cpp > CMakeFiles/local_search_test.dir/ED.cpp.i

CMakeFiles/local_search_test.dir/ED.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/local_search_test.dir/ED.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jake/PhD/Edge_Disjoint/c++/Lagrangian_Relax/git/ED.cpp -o CMakeFiles/local_search_test.dir/ED.cpp.s

CMakeFiles/local_search_test.dir/LaPSO.cpp.o: CMakeFiles/local_search_test.dir/flags.make
CMakeFiles/local_search_test.dir/LaPSO.cpp.o: LaPSO.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jake/PhD/Edge_Disjoint/c++/Lagrangian_Relax/git/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/local_search_test.dir/LaPSO.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/local_search_test.dir/LaPSO.cpp.o -c /home/jake/PhD/Edge_Disjoint/c++/Lagrangian_Relax/git/LaPSO.cpp

CMakeFiles/local_search_test.dir/LaPSO.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/local_search_test.dir/LaPSO.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jake/PhD/Edge_Disjoint/c++/Lagrangian_Relax/git/LaPSO.cpp > CMakeFiles/local_search_test.dir/LaPSO.cpp.i

CMakeFiles/local_search_test.dir/LaPSO.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/local_search_test.dir/LaPSO.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jake/PhD/Edge_Disjoint/c++/Lagrangian_Relax/git/LaPSO.cpp -o CMakeFiles/local_search_test.dir/LaPSO.cpp.s

CMakeFiles/local_search_test.dir/CpuTimer.cpp.o: CMakeFiles/local_search_test.dir/flags.make
CMakeFiles/local_search_test.dir/CpuTimer.cpp.o: CpuTimer.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jake/PhD/Edge_Disjoint/c++/Lagrangian_Relax/git/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/local_search_test.dir/CpuTimer.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/local_search_test.dir/CpuTimer.cpp.o -c /home/jake/PhD/Edge_Disjoint/c++/Lagrangian_Relax/git/CpuTimer.cpp

CMakeFiles/local_search_test.dir/CpuTimer.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/local_search_test.dir/CpuTimer.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jake/PhD/Edge_Disjoint/c++/Lagrangian_Relax/git/CpuTimer.cpp > CMakeFiles/local_search_test.dir/CpuTimer.cpp.i

CMakeFiles/local_search_test.dir/CpuTimer.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/local_search_test.dir/CpuTimer.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jake/PhD/Edge_Disjoint/c++/Lagrangian_Relax/git/CpuTimer.cpp -o CMakeFiles/local_search_test.dir/CpuTimer.cpp.s

CMakeFiles/local_search_test.dir/anyoption.cpp.o: CMakeFiles/local_search_test.dir/flags.make
CMakeFiles/local_search_test.dir/anyoption.cpp.o: anyoption.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jake/PhD/Edge_Disjoint/c++/Lagrangian_Relax/git/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/local_search_test.dir/anyoption.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/local_search_test.dir/anyoption.cpp.o -c /home/jake/PhD/Edge_Disjoint/c++/Lagrangian_Relax/git/anyoption.cpp

CMakeFiles/local_search_test.dir/anyoption.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/local_search_test.dir/anyoption.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jake/PhD/Edge_Disjoint/c++/Lagrangian_Relax/git/anyoption.cpp > CMakeFiles/local_search_test.dir/anyoption.cpp.i

CMakeFiles/local_search_test.dir/anyoption.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/local_search_test.dir/anyoption.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jake/PhD/Edge_Disjoint/c++/Lagrangian_Relax/git/anyoption.cpp -o CMakeFiles/local_search_test.dir/anyoption.cpp.s

CMakeFiles/local_search_test.dir/VolVolume.cpp.o: CMakeFiles/local_search_test.dir/flags.make
CMakeFiles/local_search_test.dir/VolVolume.cpp.o: VolVolume.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jake/PhD/Edge_Disjoint/c++/Lagrangian_Relax/git/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/local_search_test.dir/VolVolume.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/local_search_test.dir/VolVolume.cpp.o -c /home/jake/PhD/Edge_Disjoint/c++/Lagrangian_Relax/git/VolVolume.cpp

CMakeFiles/local_search_test.dir/VolVolume.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/local_search_test.dir/VolVolume.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jake/PhD/Edge_Disjoint/c++/Lagrangian_Relax/git/VolVolume.cpp > CMakeFiles/local_search_test.dir/VolVolume.cpp.i

CMakeFiles/local_search_test.dir/VolVolume.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/local_search_test.dir/VolVolume.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jake/PhD/Edge_Disjoint/c++/Lagrangian_Relax/git/VolVolume.cpp -o CMakeFiles/local_search_test.dir/VolVolume.cpp.s

CMakeFiles/local_search_test.dir/prep_mip.cpp.o: CMakeFiles/local_search_test.dir/flags.make
CMakeFiles/local_search_test.dir/prep_mip.cpp.o: prep_mip.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jake/PhD/Edge_Disjoint/c++/Lagrangian_Relax/git/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object CMakeFiles/local_search_test.dir/prep_mip.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/local_search_test.dir/prep_mip.cpp.o -c /home/jake/PhD/Edge_Disjoint/c++/Lagrangian_Relax/git/prep_mip.cpp

CMakeFiles/local_search_test.dir/prep_mip.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/local_search_test.dir/prep_mip.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jake/PhD/Edge_Disjoint/c++/Lagrangian_Relax/git/prep_mip.cpp > CMakeFiles/local_search_test.dir/prep_mip.cpp.i

CMakeFiles/local_search_test.dir/prep_mip.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/local_search_test.dir/prep_mip.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jake/PhD/Edge_Disjoint/c++/Lagrangian_Relax/git/prep_mip.cpp -o CMakeFiles/local_search_test.dir/prep_mip.cpp.s

CMakeFiles/local_search_test.dir/djikstra.cpp.o: CMakeFiles/local_search_test.dir/flags.make
CMakeFiles/local_search_test.dir/djikstra.cpp.o: djikstra.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jake/PhD/Edge_Disjoint/c++/Lagrangian_Relax/git/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object CMakeFiles/local_search_test.dir/djikstra.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/local_search_test.dir/djikstra.cpp.o -c /home/jake/PhD/Edge_Disjoint/c++/Lagrangian_Relax/git/djikstra.cpp

CMakeFiles/local_search_test.dir/djikstra.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/local_search_test.dir/djikstra.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jake/PhD/Edge_Disjoint/c++/Lagrangian_Relax/git/djikstra.cpp > CMakeFiles/local_search_test.dir/djikstra.cpp.i

CMakeFiles/local_search_test.dir/djikstra.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/local_search_test.dir/djikstra.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jake/PhD/Edge_Disjoint/c++/Lagrangian_Relax/git/djikstra.cpp -o CMakeFiles/local_search_test.dir/djikstra.cpp.s

# Object files for target local_search_test
local_search_test_OBJECTS = \
"CMakeFiles/local_search_test.dir/main_test.cpp.o" \
"CMakeFiles/local_search_test.dir/ED.cpp.o" \
"CMakeFiles/local_search_test.dir/LaPSO.cpp.o" \
"CMakeFiles/local_search_test.dir/CpuTimer.cpp.o" \
"CMakeFiles/local_search_test.dir/anyoption.cpp.o" \
"CMakeFiles/local_search_test.dir/VolVolume.cpp.o" \
"CMakeFiles/local_search_test.dir/prep_mip.cpp.o" \
"CMakeFiles/local_search_test.dir/djikstra.cpp.o"

# External object files for target local_search_test
local_search_test_EXTERNAL_OBJECTS =

local_search_test: CMakeFiles/local_search_test.dir/main_test.cpp.o
local_search_test: CMakeFiles/local_search_test.dir/ED.cpp.o
local_search_test: CMakeFiles/local_search_test.dir/LaPSO.cpp.o
local_search_test: CMakeFiles/local_search_test.dir/CpuTimer.cpp.o
local_search_test: CMakeFiles/local_search_test.dir/anyoption.cpp.o
local_search_test: CMakeFiles/local_search_test.dir/VolVolume.cpp.o
local_search_test: CMakeFiles/local_search_test.dir/prep_mip.cpp.o
local_search_test: CMakeFiles/local_search_test.dir/djikstra.cpp.o
local_search_test: CMakeFiles/local_search_test.dir/build.make
local_search_test: CMakeFiles/local_search_test.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/jake/PhD/Edge_Disjoint/c++/Lagrangian_Relax/git/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Linking CXX executable local_search_test"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/local_search_test.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/local_search_test.dir/build: local_search_test

.PHONY : CMakeFiles/local_search_test.dir/build

CMakeFiles/local_search_test.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/local_search_test.dir/cmake_clean.cmake
.PHONY : CMakeFiles/local_search_test.dir/clean

CMakeFiles/local_search_test.dir/depend:
	cd /home/jake/PhD/Edge_Disjoint/c++/Lagrangian_Relax/git && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/jake/PhD/Edge_Disjoint/c++/Lagrangian_Relax/git /home/jake/PhD/Edge_Disjoint/c++/Lagrangian_Relax/git /home/jake/PhD/Edge_Disjoint/c++/Lagrangian_Relax/git /home/jake/PhD/Edge_Disjoint/c++/Lagrangian_Relax/git /home/jake/PhD/Edge_Disjoint/c++/Lagrangian_Relax/git/CMakeFiles/local_search_test.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/local_search_test.dir/depend

