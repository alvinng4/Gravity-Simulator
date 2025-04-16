# Installation in C
Below are the dependencies for the package:

- C compiler, preferably GCC or Clang
- CMake

Optional dependencies:

- HDF5 (For data storage. Optional but recommended)
- FFTW3 (Only required for cosmological simulations)
- OpenMP

Below is an example for compiling a C program with grav_sim.

1. Clone the repository to your local machine:
```
git clone https://github.com/alvinng4/grav_sim
```

2. Create a project directory. Lets call it `test/`
```
mkdir test
cd test
```

3. Create `test.c` and `CMakeLists.txt` in `test/`:
```C title="test.c"
#include <grav_sim.h>

int main(void)
{
    print_compilation_info();
    return 0;
}
```
```CMake title="CMakeLists.txt"
cmake_minimum_required(VERSION 3.10)
project(test C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -Wall -Wextra -Wpedantic")

# Define source and library directories
set(SRC_DIR ../grav_sim)    # <-- Make sure the path is correct
include_directories(${SRC_DIR})

# Define the executable
add_executable(test test.c)

# Link against the grav_sim library
add_subdirectory(${SRC_DIR} build)
target_link_libraries(test PRIVATE grav_sim)
```

4. Create a build directory and compile the code:
```
mkdir build
cd build
cmake [-DUSE_OPENMP=ON -DUSE_FFTW3=ON -DUSE_HDF5=ON] [-DCMAKE_C_COMPILER=gcc] ..
cmake --build .
```

    * `-DUSE_OPENMP=ON -DUSE_FFTW3=ON -DUSE_HDF5=ON`: optional flags for the dependencies.
    * `-DCMAKE_C_COMPILER=gcc`: optional flag to specify the C compiler.

5. Run the code:
```
./test
```
You should see something like this:
```
-----------------------------------------------------------------
                                              __                   
    __   _ __    __     __  __           ____/\_\    ___ ___       
  /'_ `\/\`'__\/'__`\  /\ \/\ \         /',__\/\ \ /' __` __`\     
 /\ \L\ \ \ \//\ \L\.\_\ \ \_/ |       /\__, `\ \ \/\ \/\ \/\ \    
 \ \____ \ \_\\ \__/.\_\\ \___/        \/\____/\ \_\ \_\ \_\ \_\   
  \/___L\ \/_/ \/__/\/_/ \/__/   _______\/___/  \/_/\/_/\/_/\/_/   
    /\____/                     /\______\                          
    \_/__/                      \/______/                          


grav_sim version 0.0.4

Operating System: MacOS
Compilation Info:
  Compiled with OpenMP: true
  Compiled with HDF5: true
    Version: 1.14.6
  Compiled with FFTW3: true
    Version: fftw-3.3.10

Build time: Apr 16 2025 18:31:10
Compiler: GCC (version: 14)
-----------------------------------------------------------------
```

Congrats! :partying_face: Now you have successfully built a project
with the grav_sim package in C.
