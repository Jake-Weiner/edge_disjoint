# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.12

# Default target executed when no arguments are given to make.
default_target: all

.PHONY : default_target

# Allow only one "make -f Makefile2" at a time, but pass parallelism.
.NOTPARALLEL:


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

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/usr/local/bin/cmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache

.PHONY : rebuild_cache/fast

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "No interactive CMake dialog available..."
	/usr/local/bin/cmake -E echo No\ interactive\ CMake\ dialog\ available.
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache

.PHONY : edit_cache/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /home/jake/PhD/Edge_Disjoint/c++/Lagrangian_Relax/git/CMakeFiles /home/jake/PhD/Edge_Disjoint/c++/Lagrangian_Relax/git/CMakeFiles/progress.marks
	$(MAKE) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /home/jake/PhD/Edge_Disjoint/c++/Lagrangian_Relax/git/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean

.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named main_test

# Build rule for target.
main_test: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 main_test
.PHONY : main_test

# fast build rule for target.
main_test/fast:
	$(MAKE) -f CMakeFiles/main_test.dir/build.make CMakeFiles/main_test.dir/build
.PHONY : main_test/fast

CpuTimer.o: CpuTimer.cpp.o

.PHONY : CpuTimer.o

# target to build an object file
CpuTimer.cpp.o:
	$(MAKE) -f CMakeFiles/main_test.dir/build.make CMakeFiles/main_test.dir/CpuTimer.cpp.o
.PHONY : CpuTimer.cpp.o

CpuTimer.i: CpuTimer.cpp.i

.PHONY : CpuTimer.i

# target to preprocess a source file
CpuTimer.cpp.i:
	$(MAKE) -f CMakeFiles/main_test.dir/build.make CMakeFiles/main_test.dir/CpuTimer.cpp.i
.PHONY : CpuTimer.cpp.i

CpuTimer.s: CpuTimer.cpp.s

.PHONY : CpuTimer.s

# target to generate assembly for a file
CpuTimer.cpp.s:
	$(MAKE) -f CMakeFiles/main_test.dir/build.make CMakeFiles/main_test.dir/CpuTimer.cpp.s
.PHONY : CpuTimer.cpp.s

ED.o: ED.cpp.o

.PHONY : ED.o

# target to build an object file
ED.cpp.o:
	$(MAKE) -f CMakeFiles/main_test.dir/build.make CMakeFiles/main_test.dir/ED.cpp.o
.PHONY : ED.cpp.o

ED.i: ED.cpp.i

.PHONY : ED.i

# target to preprocess a source file
ED.cpp.i:
	$(MAKE) -f CMakeFiles/main_test.dir/build.make CMakeFiles/main_test.dir/ED.cpp.i
.PHONY : ED.cpp.i

ED.s: ED.cpp.s

.PHONY : ED.s

# target to generate assembly for a file
ED.cpp.s:
	$(MAKE) -f CMakeFiles/main_test.dir/build.make CMakeFiles/main_test.dir/ED.cpp.s
.PHONY : ED.cpp.s

LaPSO.o: LaPSO.cpp.o

.PHONY : LaPSO.o

# target to build an object file
LaPSO.cpp.o:
	$(MAKE) -f CMakeFiles/main_test.dir/build.make CMakeFiles/main_test.dir/LaPSO.cpp.o
.PHONY : LaPSO.cpp.o

LaPSO.i: LaPSO.cpp.i

.PHONY : LaPSO.i

# target to preprocess a source file
LaPSO.cpp.i:
	$(MAKE) -f CMakeFiles/main_test.dir/build.make CMakeFiles/main_test.dir/LaPSO.cpp.i
.PHONY : LaPSO.cpp.i

LaPSO.s: LaPSO.cpp.s

.PHONY : LaPSO.s

# target to generate assembly for a file
LaPSO.cpp.s:
	$(MAKE) -f CMakeFiles/main_test.dir/build.make CMakeFiles/main_test.dir/LaPSO.cpp.s
.PHONY : LaPSO.cpp.s

VolVolume.o: VolVolume.cpp.o

.PHONY : VolVolume.o

# target to build an object file
VolVolume.cpp.o:
	$(MAKE) -f CMakeFiles/main_test.dir/build.make CMakeFiles/main_test.dir/VolVolume.cpp.o
.PHONY : VolVolume.cpp.o

VolVolume.i: VolVolume.cpp.i

.PHONY : VolVolume.i

# target to preprocess a source file
VolVolume.cpp.i:
	$(MAKE) -f CMakeFiles/main_test.dir/build.make CMakeFiles/main_test.dir/VolVolume.cpp.i
.PHONY : VolVolume.cpp.i

VolVolume.s: VolVolume.cpp.s

.PHONY : VolVolume.s

# target to generate assembly for a file
VolVolume.cpp.s:
	$(MAKE) -f CMakeFiles/main_test.dir/build.make CMakeFiles/main_test.dir/VolVolume.cpp.s
.PHONY : VolVolume.cpp.s

anyoption.o: anyoption.cpp.o

.PHONY : anyoption.o

# target to build an object file
anyoption.cpp.o:
	$(MAKE) -f CMakeFiles/main_test.dir/build.make CMakeFiles/main_test.dir/anyoption.cpp.o
.PHONY : anyoption.cpp.o

anyoption.i: anyoption.cpp.i

.PHONY : anyoption.i

# target to preprocess a source file
anyoption.cpp.i:
	$(MAKE) -f CMakeFiles/main_test.dir/build.make CMakeFiles/main_test.dir/anyoption.cpp.i
.PHONY : anyoption.cpp.i

anyoption.s: anyoption.cpp.s

.PHONY : anyoption.s

# target to generate assembly for a file
anyoption.cpp.s:
	$(MAKE) -f CMakeFiles/main_test.dir/build.make CMakeFiles/main_test.dir/anyoption.cpp.s
.PHONY : anyoption.cpp.s

main_test.o: main_test.cpp.o

.PHONY : main_test.o

# target to build an object file
main_test.cpp.o:
	$(MAKE) -f CMakeFiles/main_test.dir/build.make CMakeFiles/main_test.dir/main_test.cpp.o
.PHONY : main_test.cpp.o

main_test.i: main_test.cpp.i

.PHONY : main_test.i

# target to preprocess a source file
main_test.cpp.i:
	$(MAKE) -f CMakeFiles/main_test.dir/build.make CMakeFiles/main_test.dir/main_test.cpp.i
.PHONY : main_test.cpp.i

main_test.s: main_test.cpp.s

.PHONY : main_test.s

# target to generate assembly for a file
main_test.cpp.s:
	$(MAKE) -f CMakeFiles/main_test.dir/build.make CMakeFiles/main_test.dir/main_test.cpp.s
.PHONY : main_test.cpp.s

prep_mip.o: prep_mip.cpp.o

.PHONY : prep_mip.o

# target to build an object file
prep_mip.cpp.o:
	$(MAKE) -f CMakeFiles/main_test.dir/build.make CMakeFiles/main_test.dir/prep_mip.cpp.o
.PHONY : prep_mip.cpp.o

prep_mip.i: prep_mip.cpp.i

.PHONY : prep_mip.i

# target to preprocess a source file
prep_mip.cpp.i:
	$(MAKE) -f CMakeFiles/main_test.dir/build.make CMakeFiles/main_test.dir/prep_mip.cpp.i
.PHONY : prep_mip.cpp.i

prep_mip.s: prep_mip.cpp.s

.PHONY : prep_mip.s

# target to generate assembly for a file
prep_mip.cpp.s:
	$(MAKE) -f CMakeFiles/main_test.dir/build.make CMakeFiles/main_test.dir/prep_mip.cpp.s
.PHONY : prep_mip.cpp.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... rebuild_cache"
	@echo "... edit_cache"
	@echo "... main_test"
	@echo "... CpuTimer.o"
	@echo "... CpuTimer.i"
	@echo "... CpuTimer.s"
	@echo "... ED.o"
	@echo "... ED.i"
	@echo "... ED.s"
	@echo "... LaPSO.o"
	@echo "... LaPSO.i"
	@echo "... LaPSO.s"
	@echo "... VolVolume.o"
	@echo "... VolVolume.i"
	@echo "... VolVolume.s"
	@echo "... anyoption.o"
	@echo "... anyoption.i"
	@echo "... anyoption.s"
	@echo "... main_test.o"
	@echo "... main_test.i"
	@echo "... main_test.s"
	@echo "... prep_mip.o"
	@echo "... prep_mip.i"
	@echo "... prep_mip.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system
