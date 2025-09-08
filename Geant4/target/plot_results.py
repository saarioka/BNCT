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
    df['X'] *= 1e3  # mm to um
    df['Y'] *= 1e3  # mm to um
    df['Z'] *= 1e3  # mm to um

    df['pX'] /= 1e3
    df['pY'] /= 1e3
    df['pZ'] /= 1e3
    df['pT'] = np.sqrt(df['pX']**2 + df['pY']**2)

    print(df)


    # plot xy 2d heatmap on log z scale
    plt.figure(figsize=(8, 6))
    plt.hist2d(
        df["X"],
        df["Y"],
        bins=100,
        range=1e-1*np.array(((-1, 1), (-1, 1))),
        cmap="viridis",
        norm=plt.matplotlib.colors.LogNorm(),
    )
    plt.colorbar()
    plt.xlabel("X (um)")
    plt.ylabel("Y (um)")
    plt.axis("equal")
    plt.tight_layout()


    plt.figure(figsize=(8, 6))
    plt.hist(
        df["Z"],
        bins=100,
        histtype="step",
    )
    plt.xlabel("Z (um)")
    #plt.ylabel("Normalized counts")
    plt.tight_layout()


    # radial distribution
    r = np.sqrt(df["X"]**2 + df["Y"]**2)
    plt.figure(figsize=(8, 6))
    plt.hist(
        r,
        bins=50,
        histtype="step",
    )
    plt.xlabel("Radial distance r (um)")
    #plt.ylabel("Normalized counts")
    plt.tight_layout()


    ## cylindrical binning

    # Define binning
    nbins_z = 50
    nbins_r = 50
    #z_edges = np.linspace(0, np.max(df["Z"]), nbins_z + 1)
    z_edges = np.linspace(0, 1.5, nbins_z + 1)
    r_edges = np.linspace(0, 0.1, nbins_r + 1)

    # Histogram counts (no weights yet)
    H, z_edges, r_edges = np.histogram2d(df["Z"], r, bins=[z_edges, r_edges])

    # Bin centers
    z_centers = 0.5 * (z_edges[:-1] + z_edges[1:])
    r_centers = 0.5 * (r_edges[:-1] + r_edges[1:])

    # Bin widths
    dz = np.diff(z_edges)
    dr = np.diff(r_edges)

    # 2D meshgrid of bin widths and centers
    Zc, Rc = np.meshgrid(z_centers, r_centers, indexing="ij")
    DZ, DR = np.meshgrid(dz, dr, indexing="ij")

    # Cylindrical shell volume for each bin
    dV = 2 * np.pi * Rc * DR * DZ

    # Density = counts / volume
    rho = H / dV

    # Normalize if desired (probability density)
    rho /= rho.sum() * dV.mean()  # ensures ∫ρ dV = 1

    # Plot
    plt.figure(figsize=(8, 6))
    plt.pcolormesh(z_edges, r_edges, rho.T, cmap="viridis", norm=plt.matplotlib.colors.LogNorm())
    plt.colorbar(label="Density (1/mm³)")
    plt.xlabel("Z (um)")
    plt.ylabel("Radial distance r (um)")
    plt.tight_layout()



    # Energy spectrum with cut on momentum direction: onyl forward neutrons
    pz = df["pZ"]
    p = np.sqrt(df["pX"]**2 + df["pY"]**2 + df["pZ"]**2)
    cos_theta = pz / p
    #print(f"Cos(theta) min: {cos_theta.min()}, max: {cos_theta.max()}")

    # different cuts on cos(theta)
    plt.figure(figsize=(8, 6))
    thetacuts = (-1, 0, 0.6, 0.7, 0.8, 0.9, 1 - round(1/13.8, 3))
    for c in thetacuts:
        plt.hist(
            df["ENeutron"][cos_theta > c],
            bins=np.linspace(0, 700, 700//5),
            histtype="step",
            label=c,
        )
    plt.xlabel("Neutron energy (keV), forward")
    plt.legend(title=r"$\cos(\theta)$ lower cut")
    #plt.ylabel("Normalized counts")
    plt.tight_layout()


    # forward neutrons only, acceptance cut
    plt.figure(figsize=(8, 6))
    plt.hist(
        df["ENeutron"][cos_theta > thetacuts[-1]],
        bins=np.linspace(0, 600, 600//5),
        histtype="step",
    )
    plt.xlabel("Neutron energy (keV)")
    #plt.ylabel("Normalized counts")
    plt.title(r"$\cos(\theta) >$" + f" {thetacuts[-1]} (device acceptance)")
    plt.xlim(0, 600)
    plt.tight_layout()

    # momentum direction distribution
    plt.figure(figsize=(8, 6))
    plt.hist(
        df["pX"],
        bins=50,
        histtype="step",
    )
    plt.hist(
        df["pY"],
        bins=50,
        histtype="step",
    )
    plt.hist(
        df["pZ"],
        bins=50,
        histtype="step",
    )
    plt.xlabel("Momentum components")
    plt.legend(["pX", "pY", "pZ"])
    plt.tight_layout()


    # momenta co-variance
    plt.figure(figsize=(8, 6))
    plt.hist2d(
        np.sqrt(df["pX"]**2 + df["pY"]**2),
        df["pZ"],
        bins=100,
        cmap="viridis",
        norm=plt.matplotlib.colors.LogNorm(),
    )
    plt.colorbar()
    plt.xlabel("pT")
    plt.ylabel("pZ")
    plt.axis("equal")
    plt.tight_layout()

    # pz vs energy
    plt.figure(figsize=(8, 6))
    d = df[cos_theta > thetacuts[-1]]
    plt.hist2d(
        d["ENeutron"],
        d["pZ"],
        bins=100,
        cmap="viridis",
        norm=plt.matplotlib.colors.LogNorm(),
    )
    plt.colorbar()
    plt.xlabel("Neutron energy (keV)")
    plt.ylabel("pZ")
    plt.title(r"$\cos(\theta) >$" + f" {thetacuts[-1]} (device acceptance)")
    plt.tight_layout()

    
    # Distributions of X, Y, Z, pX, pY, pZ in same figure
    fig, axs = plt.subplots(2, 3, figsize=(15, 8))
    axs = axs.flatten()
    variables = ["X", "Y", "Z", "pX", "pY", "pZ"]
    for i, var in enumerate(variables):
        axs[i].hist(df[var], bins=50, histtype="step")
        axs[i].set_xlabel(var)
        axs[i].set_ylabel("Counts")
        #axs[i].set_yscale("log")
    plt.tight_layout()


    # Covariances of X, Y, Z, pX, pY, pZ in same figure
    fig, axs = plt.subplots(2, 3, figsize=(15, 8))
    axs = axs.flatten()
    cov_pairs = [("X", "Y"), ("X", "Z"), ("Y", "Z"), ("pX", "pY"), ("pX", "pZ"), ("pY", "pZ")]
    for i, (var1, var2) in enumerate(cov_pairs):
        axs[i].hist2d(
            df[var1],
            df[var2],
            bins=100,
            cmap="viridis",
            norm=plt.matplotlib.colors.LogNorm(),
        )
        axs[i].set_xlabel(var1)
        axs[i].set_ylabel(var2)
        axs[i].axis("equal")
    plt.tight_layout()


    plt.show()


if __name__ == "__main__":
    main()
