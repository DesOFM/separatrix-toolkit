import argparse
import numpy as np
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("file", help="saddles.dat")
    args = parser.parse_args()

    data = np.loadtxt(args.file, comments="#")

    if data.ndim == 1:
        data = data.reshape(1, -1)

    x = data[:, 2]
    y = data[:, 3]
    ids = data[:, 1].astype(int)

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(x, y)

    for sid, xs, ys in zip(ids, x, y):
        ax.annotate(f"S{sid}", (xs, ys))

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("Detected saddle points")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True)

    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
