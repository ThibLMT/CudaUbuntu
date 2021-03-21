# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.17

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
CMAKE_COMMAND = /snap/clion/145/bin/cmake/linux/bin/cmake

# The command to remove a file.
RM = /snap/clion/145/bin/cmake/linux/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/thiblmt/CLionProjects/CudaUbuntu

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/thiblmt/CLionProjects/CudaUbuntu/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/CudaUbuntu.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/CudaUbuntu.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/CudaUbuntu.dir/flags.make

CMakeFiles/CudaUbuntu.dir/main.cu.o: CMakeFiles/CudaUbuntu.dir/flags.make
CMakeFiles/CudaUbuntu.dir/main.cu.o: ../main.cu
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/thiblmt/CLionProjects/CudaUbuntu/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CUDA object CMakeFiles/CudaUbuntu.dir/main.cu.o"
	/usr/local/cuda-11.2/bin/nvcc -forward-unknown-to-host-compiler $(CUDA_DEFINES) $(CUDA_INCLUDES) $(CUDA_FLAGS) -x cu -dc /home/thiblmt/CLionProjects/CudaUbuntu/main.cu -o CMakeFiles/CudaUbuntu.dir/main.cu.o

CMakeFiles/CudaUbuntu.dir/main.cu.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CUDA source to CMakeFiles/CudaUbuntu.dir/main.cu.i"
	$(CMAKE_COMMAND) -E cmake_unimplemented_variable CMAKE_CUDA_CREATE_PREPROCESSED_SOURCE

CMakeFiles/CudaUbuntu.dir/main.cu.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CUDA source to assembly CMakeFiles/CudaUbuntu.dir/main.cu.s"
	$(CMAKE_COMMAND) -E cmake_unimplemented_variable CMAKE_CUDA_CREATE_ASSEMBLY_SOURCE

# Object files for target CudaUbuntu
CudaUbuntu_OBJECTS = \
"CMakeFiles/CudaUbuntu.dir/main.cu.o"

# External object files for target CudaUbuntu
CudaUbuntu_EXTERNAL_OBJECTS =

CMakeFiles/CudaUbuntu.dir/cmake_device_link.o: CMakeFiles/CudaUbuntu.dir/main.cu.o
CMakeFiles/CudaUbuntu.dir/cmake_device_link.o: CMakeFiles/CudaUbuntu.dir/build.make
CMakeFiles/CudaUbuntu.dir/cmake_device_link.o: CMakeFiles/CudaUbuntu.dir/dlink.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/thiblmt/CLionProjects/CudaUbuntu/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CUDA device code CMakeFiles/CudaUbuntu.dir/cmake_device_link.o"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/CudaUbuntu.dir/dlink.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/CudaUbuntu.dir/build: CMakeFiles/CudaUbuntu.dir/cmake_device_link.o

.PHONY : CMakeFiles/CudaUbuntu.dir/build

# Object files for target CudaUbuntu
CudaUbuntu_OBJECTS = \
"CMakeFiles/CudaUbuntu.dir/main.cu.o"

# External object files for target CudaUbuntu
CudaUbuntu_EXTERNAL_OBJECTS =

CudaUbuntu: CMakeFiles/CudaUbuntu.dir/main.cu.o
CudaUbuntu: CMakeFiles/CudaUbuntu.dir/build.make
CudaUbuntu: CMakeFiles/CudaUbuntu.dir/cmake_device_link.o
CudaUbuntu: CMakeFiles/CudaUbuntu.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/thiblmt/CLionProjects/CudaUbuntu/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CUDA executable CudaUbuntu"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/CudaUbuntu.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/CudaUbuntu.dir/build: CudaUbuntu

.PHONY : CMakeFiles/CudaUbuntu.dir/build

CMakeFiles/CudaUbuntu.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/CudaUbuntu.dir/cmake_clean.cmake
.PHONY : CMakeFiles/CudaUbuntu.dir/clean

CMakeFiles/CudaUbuntu.dir/depend:
	cd /home/thiblmt/CLionProjects/CudaUbuntu/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/thiblmt/CLionProjects/CudaUbuntu /home/thiblmt/CLionProjects/CudaUbuntu /home/thiblmt/CLionProjects/CudaUbuntu/cmake-build-debug /home/thiblmt/CLionProjects/CudaUbuntu/cmake-build-debug /home/thiblmt/CLionProjects/CudaUbuntu/cmake-build-debug/CMakeFiles/CudaUbuntu.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/CudaUbuntu.dir/depend

