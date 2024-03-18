import numpy as np

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