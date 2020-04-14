import numpy as np


def main():
    xs = [
        np.array([[3, 3, 0],
                  [3, 4, 0],
                  [0, 0, 0]]),

        np.array([[3, 0, 3],
                  [0, 0, 0],
                  [3, 0, 4]]),

        np.array([[4, 0, 3],
                  [0, 0, 0],
                  [3, 0, 3]]),

        np.array([[4, 3, 0],
                  [3, 3, 0],
                  [0, 0, 0]]),

        np.array([[0, 0, 0],
                  [0, 3, 3],
                  [0, 3, 4]]),

        np.array([[0, 0, 0],
                  [0, 4, 3],
                  [0, 3, 3]]),
    ]
    Y = np.array([[5, 1, 4], [1, 1, 0], [4, 0, 3]])
    print(Y)

    def norm(x, y):
        return np.linalg.norm(x-y, ord=2)

    norms = [norm(x, Y) for x in xs]
    print(norms)


main()
