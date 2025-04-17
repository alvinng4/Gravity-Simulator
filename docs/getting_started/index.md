grav_sim is a package written in C provided with a Python wrapper.
Therefore, you can choose to install in either C or Python.
The performance difference between the two is negligible.

Here we provide a quick installation guide for Python. For a more
detailed guide, please refer to [Installation in C](installation_in_c.md) and
[Installation in Python](installation_in_python.md).


/// tab | MacOS / Linux
The package is available on PyPI, and could be installed directly using pip:
```
pip install grav_sim
```
To check whether it is successfully installed, run
```
python -m grav_sim
```
You should see the compilation information and the path to the
compiled library. If not, you may need to refer
to [Installation in Python](installation_in_python.md).
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
  Compiled with FFTW3: false

Build time: Apr 16 2025 13:12:32
Compiler: GCC (version: 14)
-----------------------------------------------------------------
C library location: /Library/Frameworks/Python.framework/Versions/3.11/lib/python3.11/site-packages/libgrav_sim.dylib
```
///

/// tab | Windows
The package is available on PyPI. However, I am not sure if it would work
on Windows. Try following the guides for MacOS / Linux.
If it does not work, please refer to [Installation in Python](installation_in_python.md),
or install WSL (Windows Subsystem for Linux).
///