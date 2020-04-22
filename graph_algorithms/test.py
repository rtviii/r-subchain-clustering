import functools
import numpy as np


def test():

    ar1 = np.array([
        [1, 0, 1],
        [1, 0, 1],
        [1, 0, 1],
    ])
    ar2 = np.array([
        [1, 0, 1],
        [1, 0, 1],
        [1, 0, 1],
    ])
    ar3 = np.array([
        [1, 0, 1],
        [1, 7, 1],
        [1, 0, 1],
    ]);
    tot = [ar1, ar2, ar3];

    w =functools.reduce(lambda x,y: np.add(x,y), tot);
    s = np.divide(w, len(w))
    print(s)

test()
