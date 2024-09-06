import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import yaml
import glob
import mplhep as hep

hep.style.use('CMS')

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

    return pd.concat(dfs)

    bin_edges = np.logspace(0, 4, 100)
    print('Bin edges', bin_edges)

    bin_middles = (bin_edges[1:] + bin_edges[:-1]) / 2
    print('Bin middles', bin_middles)

    #hist, bin_edges = np.histogram(Es.E, bins=bin_edges)
    #print('Hist', hist)
    #print('Bin edges', bin_edges)

    #plt.plot(bin_middles, hist * bin_middles)
    #plt.hist(Es.E, bins=np.logspace(0, 4, 100), **kwargs)
    #plt.hist(Es.E, bins=bin_edges, **kwargs)

    #_, bin_edges = np.histogram(np.log10(Es.E), bins=100)
    #plt.hist(Es.E, bins=10**bin_edges, **kwargs)

    #hist, bin_edges = np.histogram(Es.E, bins=bin_edges)
    #print('Hist', hist)
    #print('Bin edges', bin_edges)

    #a,b,_ = plt.hist(hist*bin_middles, bins=bin_edges, **kwargs)
    #print(a)
    #print(b)

    # https://indico.cern.ch/event/677259/contributions/2815150/attachments/1585765/2507230/Vassilis_Zermatt2018.pdf
    xmin = np.log10(1e-5)
    xmax = np.log10(10000)
    N = 100
    s = (xmax - xmin) / N
    lowE = 10**np.linspace(xmin, xmax, N)
    f = np.sqrt(10**s) / (10**s - 1)


def main():
    #read_histo('Run0_h1_E.csv')
    #read_histo('Run0_h1_Edep.csv')
    #read_histo('Run0_h1_X.csv')
    #read_histo('Run0_h1_Y.csv')
    #read_histo('Run0_h1_Z.csv')

    #read_ntuples(glob.glob('Run0_nt_Ntuple_t*.csv'))

    files = glob.glob('Run0_*mm_nt_Ntuple_t*.csv')
    ds = [int(f.split('_')[1].replace('mm', '')) for f in files]
    ds = np.flip(np.unique(ds))

    N_primaries = 1e7
    uA = N_primaries / 6.25e12

    plt.figure(figsize=(16, 8))
    Fs = []
    for i,d in enumerate(ds):
        df = read_ntuples(glob.glob(f'Run0_{d}mm_nt_Ntuple_t*.csv'))
        Es = df[df.Detector == 4].groupby('Evt').sum()

        _, bin_edges = np.histogram(np.log10(Es.E), bins=100)
        plt.hist(Es.E, bins=10**bin_edges, label=f'{d} mm' if d in (2, 80) else None, histtype='step', color=[0.9*i/len(ds)]*3, linewidth=2)

        Fs.append(len(Es))
    plt.semilogx()

    x = np.array(ds)
    y = np.array(Fs) / uA
    plt.xlabel('Energy [keV]')
    plt.ylabel('Neutrons exiting the moderator')
    plt.legend(title='Moderator thickness')
    plt.tight_layout()

    plt.figure(figsize=(16, 8))
    plt.scatter(x, y)

    #fit a line
    m, b = np.polyfit(x, y, 1)
    plt.plot(x, m*x+b, label=f'Fit: {m:.4f}x+{b:.4f}', color='black')

    plt.xlabel('Moderator thickness [mm]')
    plt.ylabel('Neutrons exiting the moderator')
    plt.tight_layout()

    plt.show()


if __name__ == '__main__':
    main()
