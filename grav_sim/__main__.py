"""
Print compilation information for the grav_sim C library.

Usage:
    python -m grav_sim

This script searches for the compiled C library in the parent directory of the current file,
loads it using ctypes, and calls the `print_compilation_info` function defined in the C library.
User can run this script to verify the installation of the grav_sim library.
"""

import ctypes
from pathlib import Path


def main() -> None:
    search_path = Path(__file__).parent.parent
    c_lib_files = [str(p) for p in search_path.rglob("*libgrav_sim*")]
    if len(c_lib_files) == 0:
        raise FileNotFoundError(f"C library not found from path: {search_path}")

    c_lib_path = c_lib_files[0]
    c_lib: ctypes.CDLL = ctypes.cdll.LoadLibrary(c_lib_path)
    c_lib.print_compilation_info()
    print(f"C library location: {c_lib_path}")


if __name__ == "__main__":
    main()
