import os
import sys
import argparse
from glob import glob
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import mplhep as hep

hep.style.use("CMS")

from wcsv import read_wcsv

THETA_CUT = 0.95
MOMEMTUM_SPREAD = 3.0  # % RMS
ENERGY_SPREAD = 1.0  # % RMS
MULTIPLY_ENTRIES = 5  # how many times to multiply the entries by smearing

def main():
    parser = argparse.ArgumentParser(description="Skim and plot neutron data")
    parser.add_argument("folder", type=str, help="Input folder containing CSV files")
    parser.add_argument("pattern", type=str, help="Filename pattern to match")
    parser.add_argument("-t", "--threads", type=int, default=22, help="Number of threads for Geant4 run. The files are split into this many parts.")
    args = parser.parse_args()

    foldername = args.folder
    filename_pattern = args.pattern

    files = glob(f"{foldername}/{filename_pattern}*.csv")
    print(f"Found {len(files)} files")

    output_dir = foldername.replace("results", "results_skimmed")
    os.makedirs(output_dir, exist_ok=True)

    df = []
    header = ""  # assume same header for all files
    for f in files:
        data, header = read_wcsv(f, header=True)
        df.append(data)
    df = pd.concat(df, ignore_index=True)
    
    df['pT'] = np.sqrt(df['pX']**2 + df['pY']**2)
    df['costheta'] = df['pZ'] / np.sqrt(df['pX']**2 + df['pY']**2 + df['pZ']**2)

    df_out = df[df['costheta'] > THETA_CUT]

    for t in range(args.threads):
        df_out_thread = df_out[df_out.index % args.threads == t]
        outfile = f"{output_dir}/{args.pattern}_skimmed_nt_Kinematics_t{t}.csv"
        with open(outfile, 'w') as file:
            file.write(header + "\n")
            df_out_thread.to_csv(file, index=False, header=None)
    
    print(f"Skimmed data written to {output_dir}")
    print(f'Amount of data per file (for {args.threads} threads): {len(df_out)/args.threads:.0f} events/file')
    print(f"Total neutrons: {len(df)}, forward neutrons (cos(theta)>{THETA_CUT}): {len(df_out)}")
    print(f"Corresponding angle: {np.arccos(THETA_CUT)*180/np.pi:.1f} deg")
    print(f"Forward neutrons fraction: {len(df_out)/len(df)*100:.1f} %")
    print(f"Forward neutrons mean energy: {df_out['ENeutron'].mean():.1f} keV")

    # Create more entries by introducing some spread to momentum and energy
    np.random.seed(42)
    df_out = pd.concat([df_out] * MULTIPLY_ENTRIES, ignore_index=True)
    df_out['pX'] *= np.random.normal(1, MOMEMTUM_SPREAD/100, size=len(df_out))
    df_out['pY'] *= np.random.normal(1, MOMEMTUM_SPREAD/100, size=len(df_out))
    df_out['pZ'] *= np.random.normal(1, MOMEMTUM_SPREAD/100, size=len(df_out))
    df_out['ENeutron'] *= np.random.normal(1, ENERGY_SPREAD/100, size=len(df_out))

    print(f"After smearing, amount of forward neutrons: {len(df_out)}")

    df['X'] *= 1e3  # mm to um
    df['Y'] *= 1e3  # mm to um
    df['Z'] *= 1e3  # mm to um

    df['pX'] /= 1e3
    df['pY'] /= 1e3
    df['pZ'] /= 1e3

    #print(df)

    # Energy spectrum with cut on momentum direction: onyl forward neutrons
    # different cuts on cos(theta)
    plt.figure(figsize=(8, 6))
    thetacuts = (-1, 0, 0.6, 0.7, 0.8, 0.9, 1 - round(1/13.8, 3))
    for c in thetacuts:
        plt.hist(
            df["ENeutron"][df["costheta"] > c],
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
        df_out["ENeutron"],
        bins=np.linspace(0, 600, 600//5),
        histtype="step",
    )
    plt.xlabel("Neutron energy (keV)")
    #plt.ylabel("Normalized counts")
    plt.title(r"$\cos(\theta) >$" + f" {THETA_CUT}")
    plt.xlim(0, 600)
    plt.tight_layout()

    # momentum direction distribution
    plt.figure(figsize=(8, 6))
    plt.hist(
        df_out["pX"],
        bins=50,
        histtype="step",
    )
    plt.hist(
        df_out["pY"],
        bins=50,
        histtype="step",
    )
    plt.hist(
        df_out["pZ"],
        bins=50,
        histtype="step",
    )
    plt.xlabel("Momentum components")
    plt.legend(["pX", "pY", "pZ"])
    plt.tight_layout()


    plt.show()


if __name__ == "__main__":
    main()
