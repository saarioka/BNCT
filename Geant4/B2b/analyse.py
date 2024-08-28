import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def read_csv(fname):
    with open(fname, 'r') as f:
        f.readline()
        f.readline()
        line3 = f.readline()
    #cols = line3.replace('# ', '').replace(' [', '-[').split()
    cols = line3.split(',')
    cols = [l.replace('# ', '') for l in cols]
    cols = [l.replace(' ', '') for l in cols]
    cols = [l.replace('\n', '') for l in cols]
    print(cols)
    return pd.read_csv(fname, skiprows=3, names=cols, sep=',')


def main():
    current = read_csv('current.csv')
    print(current)

    nofstep = read_csv('nofstep.csv')
    print(nofstep)
    print(nofstep['entry'].sum())

    noftrack = read_csv('noftrack.csv')
    print(noftrack)
    print(noftrack['entry'].sum())

    flux = read_csv('flux.csv')
    print(flux)
    print(flux['total(value)[percm2]'].sum())

    ax1 = plt.contourf(


if __name__ == '__main__':
    main()
