import numpy as np
import itertools
from typing import Tuple


def lower_backwards(A: np.ndarray, b: np.ndarray) -> np.ndarray:
    x = np.zeros((len(b), 1))
    for i in range(0, len(b)):
        x[i] = (1 / A[i, i] * (b[i] - sum([A[i, j] * x[j] for j in range(0, i)])))
    return np.array(x)


def upper_backwards(A: np.ndarray, b: np.ndarray) -> np.ndarray:
    x = np.zeros((len(b), 1))
    for i in reversed(range(0, len(b))):
        x[i] = (1 / A[i, i] * (b[i] - sum([A[i, j] * x[j]
                                           for j in range(1 + i, len(b))])))
    return np.array(x)

if __name__ == "__main__":
    L = np.array([[1, 0, 0],
                  [1, 3, 0],
                  [1, 2, 1]])

    U = np.array([[1, 2, 1],
                  [0, 3, 1],
                  [0, 0, 1]])

    b = np.array([[1], [5], [3]])

    print(lower_backwards(L, b))
    
