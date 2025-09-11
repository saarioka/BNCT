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

THETA_CUT = -1
MOMEMTUM_SPREAD = 0  # % RMS
ENERGY_SPREAD = 0  # % RMS
MULTIPLY_ENTRIES = 1  # how many times to multiply the entries by smearing

def main():
    parser = argparse.ArgumentParser(description="Skim and plot neutron data")
    parser.add_argument('-f', "--folders", type=str, help="Input folders containing CSV files")
    parser.add_argument('-p', "--pattern", type=str, help="Filename pattern to match")
    parser.add_argument("-t", "--threads", type=int, required=False, help="Number of threads for Geant4 run. The files are split into this many parts.")
    args = parser.parse_args()

    foldername = args.folders
    filename_pattern = args.pattern

    print(foldername)

    files = glob(f"{foldername}/{filename_pattern}*.csv")

    if len(files) == 0:
        print("No files found matching the pattern.")
        sys.exit(1)

    print(f"Found {len(files)} files")

    threads = args.threads or len(files)

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

    print(f"Total neutrons: {len(df)}, forward neutrons (cos(theta)>{THETA_CUT}): {len(df_out)}")
    print(f"Corresponding angle: {np.arccos(THETA_CUT)*180/np.pi:.1f} deg")
    print(f"Forward neutrons fraction: {len(df_out)/len(df)*100:.1f} %")
    print(f"Forward neutrons mean energy: {df_out['ENeutron'].mean():.1f} keV")

    if MULTIPLY_ENTRIES > 1:
        print(f"Multiplying entries by {MULTIPLY_ENTRIES} with smearing: momentum {MOMEMTUM_SPREAD} % RMS, energy {ENERGY_SPREAD} % RMS")
        # Create more entries by introducing some spread to momentum and energy
        np.random.seed(42)
        df_out = pd.concat([df_out] * MULTIPLY_ENTRIES, ignore_index=True)
        df_out['pX'] *= np.random.normal(1, MOMEMTUM_SPREAD/100, size=len(df_out))
        df_out['pY'] *= np.random.normal(1, MOMEMTUM_SPREAD/100, size=len(df_out))
        df_out['pZ'] *= np.random.normal(1, MOMEMTUM_SPREAD/100, size=len(df_out))
        if ENERGY_SPREAD > 0:
            df_out['ENeutron'] *= np.random.normal(1, ENERGY_SPREAD/100, size=len(df_out))

    for t in range(threads):
        df_out_thread = df_out[df_out.index % threads == t]
        outfile = f"{output_dir}/{args.pattern}_skimmed_{round(THETA_CUT*100)}_{MULTIPLY_ENTRIES}_nt_Kinematics_t{t}.csv"
        with open(outfile, 'w') as file:
            file.write(header + "\n")
            df_out_thread.to_csv(file, index=False, header=None)
    
    print(f"Skimmed data written to {output_dir}/{args.pattern}_skimmed_{round(THETA_CUT*100)}_{MULTIPLY_ENTRIES}_nt_Kinematics_t?.csv")
    print(f'Data split into {threads} files for {threads} threads')
    print(f"Total forward neutrons after smearing: {len(df_out)}")
    print(f'Amount of data per file (for {threads} threads): {len(df_out)/threads:.0f} events/file')

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
    plt.xlabel("Neutron energy (keV)")
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
