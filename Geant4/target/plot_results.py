import sys
from glob import glob
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import mplhep as hep

hep.style.use("CMS")

from wcsv import read_wcsv


def main():
    files = sys.argv[1:]
    print(f"Found {len(files)} files")

    df = []
    for f in files:
        #Evt  Edep   EProton  ENeutron         X         Y             Z
        data = read_wcsv(f)
        df.append(data)
    df = pd.concat(df, ignore_index=True)
    print(df)


    plt.figure(figsize=(8, 6))
    plt.hist(
        df["ENeutron"],
        bins=np.linspace(0, 700, 700//5),
        histtype="step",
        density=True,
    )
    plt.xlabel("Neutron energy (keV)")
    plt.ylabel("Normalized counts")
    plt.legend()
    plt.tight_layout()


    plt.figure(figsize=(8, 6))
    plt.hist(
        df["EProton"],
        bins=np.linspace(0, 700, 700//5),
        histtype="step",
        density=True,
    )
    plt.xlabel("Proton energy (keV)")
    plt.ylabel("Normalized counts")
    plt.legend()
    plt.tight_layout()


    # plot xy 2d heatmap
    plt.figure(figsize=(8, 6))
    plt.hist2d(
        df["X"],
        df["Y"],
        bins=100,
        range=[[-1, 1], [-1, 1]],
        density=True,
        cmap="viridis",
    )
    plt.colorbar(label="Normalized counts")
    plt.xlabel("X (mm)")
    plt.ylabel("Y (mm)")
    plt.axis("equal")
    plt.tight_layout()


    plt.figure(figsize=(8, 6))
    plt.hist(
        df["Z"],
        bins=np.linspace(-1, 1, 100),
        histtype="step",
        density=True,
    )
    plt.xlabel("Z (mm)")
    plt.ylabel("Normalized counts")
    plt.tight_layout()

    plt.show()


if __name__ == "__main__":
    main()
