import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import mplhep as hep

hep.style.use(hep.style.CMS)

# It seems this = primaries x cycles
N_PRIMARIES_FLUKA = 20e6
N_PRIMARIES_MCNP = 1e9

colors = ['deeppink', 'deepskyblue', 'orangered', 'limegreen']

def read_MCNP():
    df = pd.read_csv('8cm_600kev_out.csv')
    df['Counts'] /= N_PRIMARIES_MCNP
    df['bin_centers'] = (df['Energy'] + np.roll(df['Energy'], 1)) / 2
    return df

def read_FLUKA():
    fnames = [f'bnct_{i}_tab.lis' for i in [34, 42, 43, 51]]
    data = {}
    for f in fnames:
        with open(f) as ff:
            lines = ff.readlines()
            indata = []
            for l,line in enumerate(lines):
                if l == 0:
                    name = line.split()[4]
                elif l == 1:
                    bins = int(line.split()[-1])
                elif line.strip() == '':
                    break
                else:
                    b1, b2, c, u = line.split()
                    b1, b2, c, u = float(b1), float(b2), float(c), float(u)
                    b1 *= 1e6
                    b2 *= 1e6
                    u *= 1e6
                    mb = (b1 + b2) / 2
                    c /= N_PRIMARIES_FLUKA
                    u /= N_PRIMARIES_FLUKA
                    indata.append([b1, b2, c, u, mb])
            print('Name:', name, ', bins:', bins)
            df = pd.DataFrame(indata, columns=['bin_left', 'bin_right', 'Counts', 'Uncertainty', 'mean_bin'])
            df.drop(columns=['bin_left'], inplace=True)
            df.rename(columns={'bin_right': 'Energy'}, inplace=True)
            df.loc[-1] = [0, 0, 0, 0]
            df.index = df.index + 1
            df.sort_index(inplace=True)
            data[name] = df
            print(len(df))
    return data


def main():
    dfm = read_MCNP()
    dff = read_FLUKA()

    kwargs_fluka = {'label': 'FLUKA <X>Y', 'linestyle': '-', 'linewidth': 2, 'color': colors[0], 'where': 'pre'}
    kwargs_mcnp = {'label': 'MCNP', 'linestyle': '-', 'linewidth': 2, 'color': colors[1], 'where': 'pre'}

    plt.figure(figsize=(13, 7))
    plt.title('Neutron energy spectra IN\ntwo-way current')

    df = dff['I2_1']
    plt.step(df['Energy'] * 1e3 , df['Counts'] * df['mean_bin'], **kwargs_fluka)
    plt.errorbar(df['mean_bin'] * 1e3, df['Counts'] * df['mean_bin'], yerr=df['Counts'] * df['mean_bin'] * df['Uncertainty'], fmt='none', color=colors[0])

    df = dfm[dfm['Surface'] == 7]
    plt.step(df['Energy'] * 1e3, df['Counts'], **kwargs_mcnp)
    plt.errorbar(df['bin_centers'] * 1e3, df['Counts'], yerr=df['Counts'] * df['Uncertainty'], linewidth=2, fmt='none', color=colors[1])

    plt.legend()
    plt.loglog()
    plt.xlabel('Energy (eV)')
    plt.xlim([1e-3, 1e6])
    plt.ylabel('$I_n$ (1/primary/400cm$^2$)')
    plt.tight_layout()


    plt.figure(figsize=(13, 7))
    plt.title('Neutron energy spectra OUT\ntwo-way current')

    df = dff['I2_9']
    plt.step(df['Energy'] * 1e3 , df['Counts'] * df['mean_bin'], **kwargs_fluka)
    plt.errorbar(df['mean_bin'] * 1e3, df['Counts'] * df['mean_bin'], yerr=df['Counts'] * df['mean_bin'] * df['Uncertainty'], fmt='none', color=colors[0])

    df = dfm[dfm['Surface'] == 15]
    plt.step(df['Energy'] * 1e3, df['Counts'], **kwargs_mcnp)
    plt.errorbar(df['bin_centers'] * 1e3, df['Counts'], yerr=df['Counts'] * df['Uncertainty'], linewidth=2, fmt='none', color=colors[1])

    plt.legend()
    plt.loglog()
    plt.xlabel('Energy (eV)')
    plt.xlim([1e-3, 1e6])
    plt.ylabel('$I_n$ (1/primary/400cm$^2$)')
    plt.tight_layout()


    plt.figure(figsize=(13, 7))
    plt.title('Neutron energy spectra\none-way current')

    df = dff['I1_1']
    kwargs_fluka['label'] = 'FLUKA IN'
    kwargs_fluka['color'] = colors[2]
    plt.step(df['Energy'] * 1e3 , df['Counts'] * df['mean_bin'], **kwargs_fluka)
    plt.errorbar(df['mean_bin'] * 1e3, df['Counts'] * df['mean_bin'], yerr=df['Counts'] * df['mean_bin'] * df['Uncertainty'], fmt='none', color=colors[2])

    df = dff['I1_9']
    kwargs_fluka['label'] = 'FLUKA OUT'
    kwargs_fluka['color'] = colors[3]
    plt.step(df['Energy'] * 1e3 , df['Counts'] * df['mean_bin'], **kwargs_fluka)
    plt.errorbar(df['mean_bin'] * 1e3, df['Counts'] * df['mean_bin'], yerr=df['Counts'] * df['mean_bin'] * df['Uncertainty'], fmt='none', color=colors[3])

    plt.legend()
    plt.loglog()
    plt.xlabel('Energy (eV)')
    plt.xlim([1e-3, 1e6])
    plt.ylabel('$I_n$ (1/primary/400cm$^2$)')
    plt.tight_layout()

    plt.show()


if __name__ == '__main__':
    main()
