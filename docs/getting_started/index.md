`grav_sim` is a package written in C provided with a Python wrapper.
Therefore, you can choose to use `grav_sim` in either C or Python.
The performance difference between the two is negligible.

Here we provide a quick guide for building `grav_sim` on Linux / macOS.

/// tab | C

Below are the dependencies for the package:

- C compiler, preferably GCC or Clang
- CMake

Optional dependencies:

- OpenMP (For parallel computing. Optional but recommended)
- HDF5 (For data storage. Optional but recommended)
- FFTW3 (Only required for cosmological simulations)

Below is an example for compiling a C program with grav_sim.

1. Clone the repository to your local machine:
```
git clone https://github.com/alvinng4/grav_sim
cd grav_sim
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
cmake_minimum_required(VERSION 3.14)

project(test
    VERSION 0.0.1
    LANGUAGES C
    DESCRIPTION "test example using grav_sim"
)

### grav_sim ###
# Directory containing grav_sim's top-level CMakeLists.txt.
set(GRAV_SIM_SOURCE_DIR
    "${CMAKE_CURRENT_SOURCE_DIR}/../"
    CACHE PATH
    "Path to the grav_sim source directory"
)

if(NOT EXISTS "${GRAV_SIM_SOURCE_DIR}/CMakeLists.txt")
    message(FATAL_ERROR
        "grav_sim was not found at:\n"
        "  ${GRAV_SIM_SOURCE_DIR}\n"
        "Set GRAV_SIM_SOURCE_DIR to the grav_sim source directory."
    )
endif()

add_subdirectory(
    "${GRAV_SIM_SOURCE_DIR}"
    "${CMAKE_CURRENT_BINARY_DIR}/grav_sim"
)

### Executable ###
add_executable(test
    test.c
)

set_target_properties(test PROPERTIES
    C_STANDARD 99
    C_STANDARD_REQUIRED ON
    C_EXTENSIONS OFF
)

### Dependencies ###
target_link_libraries(test PRIVATE
    grav_sim::grav_sim
)

### Warnings ###
if(MSVC)
    target_compile_options(test PRIVATE
        /W4
    )
elseif(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
    target_compile_options(test PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Wshadow
        -Wpointer-arith
        -Wcast-qual
        -Wstrict-prototypes
        -Wmissing-prototypes
    )
endif()
```

4. Create a build directory and compile the code:
```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release --fresh
cmake --build build
```

5. Run the code:
```
./build/test
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


grav_sim version 1.0.1

Operating System: Linux
Compilation Info:
  Compiled with OpenMP: false
  Compiled with HDF5: false
  Compiled with FFTW3: false

Build time: Oct 19 2026 12:04:43
Compiler: GCC (version: 13)
-----------------------------------------------------------------
```

6. Some options can be enabled if needed
```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DUSE_OPENMP=ON -DUSE_FFTW3=ON -DUSE_HDF5=ON -DCMAKE_C_COMPILER=gcc --fresh
cmake --build build
```

    * `-DUSE_OPENMP=ON -DUSE_FFTW3=ON -DUSE_HDF5=ON`: optional flags for the dependencies.
    * `-DCMAKE_C_COMPILER=gcc`: optional flag to specify the C compiler.
Now rerun the code
```
./build/test
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


grav_sim version 1.0.1

Operating System: Linux
Compilation Info:
  Compiled with OpenMP: true
  Compiled with HDF5: true
    Version: 1.10.10
  Compiled with FFTW3: true
    Version: fftw-3.3.10-sse2-avx

Build time: Oct 19 2026 12:05:38
Compiler: GCC (version: 13)
-----------------------------------------------------------------
```

Congrats! :partying_face: You have successfully built a project
with the grav_sim package in C.

Check out our example projects for more usage with the C API.
///

/// tab | Python
To use the Python wrapper, the C library needs to be compiled first.

Below are the dependencies for the package:

- C compiler, preferably GCC or Clang
- CMake

Optional dependencies:

- OpenMP (For parallel computing. Optional but recommended)
- HDF5 (For data storage. Optional but recommended)
- FFTW3 (Only required for cosmological simulations)

Below is an example for compiling a C program with grav_sim.

1. Clone the repository to your local machine:
```
git clone https://github.com/alvinng4/grav_sim
cd grav_sim
```

2. Create a build directory and compile the code:
```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release --fresh
cmake --build build
```
3. Some options can be enabled if needed
```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DUSE_OPENMP=ON -DUSE_FFTW3=ON -DUSE_HDF5=ON -DCMAKE_C_COMPILER=gcc --fresh
cmake --build build
```

    * `-DUSE_OPENMP=ON -DUSE_FFTW3=ON -DUSE_HDF5=ON`: optional flags for the dependencies.
    * `-DCMAKE_C_COMPILER=gcc`: optional flag to specify the C compiler.

The C library should be compiled. Next, we install the Python dependencies.
/// tab | pip
First, let us create a virtual environment
```
python -m venv test_grav
source test_grav/bin/activate
```
Then, install the package via
```
pip install .
```
Now, you may test it by
```
python -m grav_sim
```
You should see something like
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


grav_sim version 1.0.1

Operating System: Linux
Compilation Info:
  Compiled with OpenMP: true
  Compiled with HDF5: true
    Version: 1.10.10
  Compiled with FFTW3: true
    Version: fftw-3.3.10-sse2-avx

Build time: Oct 19 2026 12:05:38
Compiler: GCC (version: 13)
-----------------------------------------------------------------
```
///

/// tab | uv
Simply run
```
uv sync
```
Now, you may test it by
```
uv run -m grav_sim
```
You should see something like
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


grav_sim version 1.0.1

Operating System: Linux
Compilation Info:
  Compiled with OpenMP: true
  Compiled with HDF5: true
    Version: 1.10.10
  Compiled with FFTW3: true
    Version: fftw-3.3.10-sse2-avx

Build time: Oct 19 2026 12:05:38
Compiler: GCC (version: 13)
-----------------------------------------------------------------
```
///
Done! :partying_face: Check out our [Python tutorial](http://alvinng4.com/grav_sim/examples/tutorial/) for more usage.
///
