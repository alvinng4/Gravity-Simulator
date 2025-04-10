import ctypes
import glob
from pathlib import Path


def main() -> None:
    search_path = Path(__file__).parent.parent
    c_lib_files = [str(p) for p in search_path.rglob('*libgrav_sim*')]
    if len(c_lib_files) == 0:
        raise FileNotFoundError(f"C library not found from path: {search_path}")

    c_lib_path = c_lib_files[0]
    print(c_lib_path)
    c_lib: ctypes.CDLL = ctypes.cdll.LoadLibrary(c_lib_path)
    c_lib.print_compilation_info()

if __name__ == "__main__":
    main()
