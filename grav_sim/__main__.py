import argparse
import ctypes
import glob
from pathlib import Path


def parse_args() -> str:
    parser = argparse.ArgumentParser(description="Gravity-Simulator")
    parser.add_argument(
        "path",
        type=str,
        help="Path to the C library",
        default="../src",
    )
    args = parser.parse_args()
    return args.path

def main() -> None:
    c_lib_path: Path = Path(parse_args())
    if not c_lib_path.exists():
        raise FileNotFoundError(f"Path {c_lib_path} does not exist.")
    













if __name__ == "__main__":
    main()