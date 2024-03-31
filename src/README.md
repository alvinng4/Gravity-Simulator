## Compilation
Below is the command I used to compile the source files.

On macbook air M1:
```
gcc-13 -shared -fPIC -O3 c_lib.c -o c_lib.so
```

On windows:
```
gcc -shared -fPIC -O3 c_lib.c -o c_lib.dll
```

After compilation, put the files into the `gravity_plot` or `gravity_sim` folder.
