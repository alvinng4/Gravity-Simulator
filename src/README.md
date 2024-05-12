## Compilation
Below is the command I used to compile the source files.

On mac:
```
gcc -shared -fPIC -O3 -std=c99 -Wall -Wextra -Wpedantic -Wmisleading-indentation c_lib.c -o c_lib.dylib
```

On windows:
```
gcc -shared -fPIC -O3 -std=c99 -Wall -Wextra -Wpedantic -Wmisleading-indentation c_lib.c -o c_lib.dll
```

On linux:
```
gcc -shared -fPIC -O3 -std=c99 -Wall -Wextra -Wpedantic -Wmisleading-indentation c_lib.c -o c_lib.so
```

After compilation, put the files into the `gravity_plot` or `gravity_sim` folder.
