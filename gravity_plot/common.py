import re

import numpy as np


def get_bool(msg: str) -> bool:
    """
    Get boolean input from user

    Args: msg (str): Message to display to user
    Display: "msg (y/n):"
    Returns: bool
    """
    while True:
        if matches := re.search(
            r"^\s*(yes|no|y|n)\s*$", input(f"{msg} (y/n): "), re.IGNORECASE
        ):
            if matches.group(1).lower() in ["y", "yes"]:
                return True
            elif matches.group(1).lower() in ["n", "no"]:
                return False

        print("Invalid input. Please try again.\n")


def acceleration(objects_count, x, m, G):
    """
    Calculate acceleration by a = - GM/r^3 vec{r}
    """
    # Allocate memory
    temp_a = np.zeros((objects_count * objects_count, 3))

    # Calculations
    for j in range(objects_count):
        for k in range(j + 1, objects_count):
            R = x[j] - x[k]
            temp_value = G * R / np.linalg.norm(R) ** 3
            temp_a[j * objects_count + k] = -temp_value * m[k]
            temp_a[k * objects_count + j] = temp_value * m[j]

    temp_a = temp_a.reshape((objects_count, objects_count, 3))
    a = np.sum(temp_a, axis=1)

    return a
