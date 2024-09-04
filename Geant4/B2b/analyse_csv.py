import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import yaml
import glob

def read_histo(file_name):
    print('Reading file: ', file_name)
    with open(file_name, 'r') as f:
        l = 0
        for line in f.readlines():
            ll = line.replace('\n', '')
            if ll.startswith('#title'):
                title = ll.replace('#title', '')
            if ll.startswith('#axis'):
                _, _, bins, lower, upper = ll.split()
            if not line.startswith('#'):
                break
            l += 1
    df = pd.read_csv(file_name, skiprows=l)

    info = {
        'Bins:': bins,
        'Lower lim:': lower,
        'Upper lim:': upper,
        'Underflow:': df['entries'].iloc[0],
        'Overflow:': df['entries'].iloc[-1],
        'Total entries:': df['entries'].sum()
    }
    info_str = '\n'.join([f'{k} {v}' for k, v in info.items()])

    x = np.linspace(float(lower), float(upper), int(bins))

    plt.figure(figsize=(10, 6))
    plt.title(title)
    plt.step(x, df['entries'][1:-1], label='Entries per bin')
    plt.legend(title=info_str)
    plt.savefig(file_name.replace('.csv', '.pdf'))


def read_ntuples(file_names):
    data_start = 0
    dfs = []
    for i, fn in enumerate(file_names):
        if i == 0:
            print('Reading file: ', fn)

            with open(fn, 'r') as f:
                cols = []
                for line in f.readlines():
                    ll = line.replace('\n', '')
                    if ll.startswith('#title'):
                        title = ll.replace('#title', '')
                    elif ll.startswith('#column'):
                        cols.append(ll.split()[-1])
                    if not line.startswith('#'):
                        break
                    data_start += 1

        dfs.append(pd.read_csv(fn, skiprows=data_start, names=cols))
    df = pd.concat(dfs)
    print(df.describe())
    plt.figure()
    Es = df[df.Edep > 0].groupby('Evt').sum()
    print(Es.describe())

    plt.figure(figsize=(10, 6))
    plt.hist(Es.Edep, bins=100, histtype='step')


def main():
    read_histo('Run0_h1_E.csv')
    read_histo('Run0_h1_Edep.csv')
    read_histo('Run0_h1_X.csv')
    read_histo('Run0_h1_Y.csv')
    read_histo('Run0_h1_Z.csv')

    read_ntuples(glob.glob('Run0_nt_Ntuple_t*.csv'))


    plt.show()


if __name__ == '__main__':
    main()
