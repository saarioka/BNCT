import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import mplhep as hep

hep.style.use(hep.style.CMS)

# It seems this = primaries x cycles
N_PRIMARIES_FLUKA = 20e6
N_PRIMARIES_MCNP = 1e9
N_PRIMARIES_GEANT4 = 1e8

colors = ['deeppink', 'deepskyblue', 'orangered', 'limegreen', 'black', 'blue', 'red', 'green']

def read_MCNP():
    df = pd.read_csv('8cm_600kev_out.csv')
    df['Counts'] /= N_PRIMARIES_MCNP
    df['bin_centers'] = (df['Energy'] + np.roll(df['Energy'], 1)) / 2
    return df

def read_Geant4():
    df = pd.read_csv('600keV_8cm_n_out_timo.csv', header=None, names=['Energy', 'Counts'])
    df['Counts'] /= N_PRIMARIES_GEANT4
    df['bin_centers'] = (df['Energy'] + np.roll(df['Energy'], 1)) / 2
    df['bin_widths'] = df['Energy'] - np.roll(df['Energy'], 1)
    return df

def read_FLUKA(file_index):
    fnames = [f'data/bnct-{file_index}_{i}_tab.lis' for i in [34, 42, 43, 51]]
    if file_index == 0:
        print('Special FLUKA case')
        fnames = [f'data/bnct_{i}_tab.lis' for i in [34, 42, 43, 51]]
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
            #df.index = df.index + 1
            #df.sort_index(inplace=True)
            data[name] = df
            #print(len(df))
    return data


def main():
    kwargs_fluka = {'label': 'FLUKA <X>Y', 'linestyle': '-', 'linewidth': 2, 'color': colors[0], 'where': 'pre'}
    kwargs_mcnp = {'label': 'MCNP', 'linestyle': '-', 'linewidth': 2, 'color': colors[1], 'where': 'pre'}
    kwargs_geant4 = {'label': 'Geant4', 'linestyle': '-', 'linewidth': 2, 'color': colors[2], 'where': 'post'}

    plt.figure(figsize=(13, 7))
    plt.title('Neutron energy spectra IN\ntwo-way current')

    dff = read_FLUKA(0)

    df = dff['I2_1']
    plt.step(df['Energy'] * 1e3 , df['Counts'] * df['mean_bin'], **kwargs_fluka)
    plt.errorbar(df['mean_bin'] * 1e3, df['Counts'] * df['mean_bin'], yerr=df['Counts'] * df['mean_bin'] * df['Uncertainty'], fmt='none', color=colors[0])

    dfm = read_MCNP()

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
    #plt.step(df['Energy'] * 1e3 , df['Counts'] * df['mean_bin'], **kwargs_fluka)
    plt.step(df['Energy'] * 1e3 , df['Counts'], **kwargs_fluka)
    #plt.errorbar(df['mean_bin'] * 1e3, df['Counts'] * df['mean_bin'], yerr=df['Counts'] * df['mean_bin'] * df['Uncertainty'], fmt='none', color=colors[0])

    df = dfm[dfm['Surface'] == 15]
    #plt.step(df['Energy'] * 1e3, df['Counts'], **kwargs_mcnp)
    plt.step(df['Energy'] * 1e3, df['Counts'] / df['bin_centers'], **kwargs_mcnp)
    #plt.errorbar(df['bin_centers'] * 1e3, df['Counts'], yerr=df['Counts'] * df['Uncertainty'], linewidth=2, fmt='none', color=colors[1])

    dfg = read_Geant4()

    df = dfg
    #plt.step(df['Energy'] * 1e3, df['Counts'] / df['bin_widths'] * df['bin_centers'], **kwargs_geant4)
    #plt.step(df['Energy'] * 1e3, df['Counts'] * df['bin_centers'], **kwargs_geant4)
    #plt.step(df['Energy'] * 1e3, df['Counts'] / df['bin_widths'] * df['bin_centers'], **kwargs_geant4)
    plt.step(df['Energy'] * 1e3, df['Counts'] / df['bin_centers'], **kwargs_geant4)

    plt.gca().axvline(12.75e-3, color='black', linestyle='--', label='12.75 meV')

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


    plt.figure(figsize=(13, 7))
    plt.title('Neutron energy spectra IN\ntwo-way current')

    for i,file_index in enumerate((1,2,3,4,5,6)):
        dff = read_FLUKA(file_index)
        df = dff['I2_1']
        plt.step(df['Energy'] * 1e3 , df['Counts'] * df['mean_bin'], color=colors[i], label=f'FLUKA <X>Y {file_index*100} keV')
        plt.errorbar(df['mean_bin'] * 1e3, df['Counts'] * df['mean_bin'], yerr=df['Counts'] * df['mean_bin'] * df['Uncertainty'], fmt='none', color=colors[i])

    dfm = read_MCNP()
    df = dfm[dfm['Surface'] == 7]
    #plt.step(df['Energy'] * 1e3, df['Counts'], label='MCNP', color=colors[1], linewidth=2, linestyle='--', where='pre')
    #plt.errorbar(df['bin_centers'] * 1e3, df['Counts'], yerr=df['Counts'] * df['Uncertainty'], linewidth=2, fmt='none', color=colors[1])

    plt.legend()
    plt.loglog()
    plt.xlabel('Energy (eV)')
    plt.xlim([1e-3, 1e6])
    plt.ylabel('$I_n$ (1/primary/400cm$^2$)')
    plt.tight_layout()


    plt.figure(figsize=(13, 7))
    plt.title('Neutron energy spectra OUT\ntwo-way current')

    for i,file_index in enumerate((1,2,3,4,5,6)):
        dff = read_FLUKA(file_index)
        df = dff['I2_9']
        plt.step(df['Energy'] * 1e3 , df['Counts'] * df['mean_bin'], color=colors[i], label=f'FLUKA <X>Y {file_index*100} keV')
        plt.errorbar(df['mean_bin'] * 1e3, df['Counts'] * df['mean_bin'], yerr=df['Counts'] * df['mean_bin'] * df['Uncertainty'], fmt='none', color=colors[i])

    dfm = read_MCNP()
    df = dfm[dfm['Surface'] == 15]
    #plt.step(df['Energy'] * 1e3, df['Counts'], label='MCNP', color=colors[1], linewidth=2, linestyle='--', where='pre')
    #plt.errorbar(df['bin_centers'] * 1e3, df['Counts'], yerr=df['Counts'] * df['Uncertainty'], linewidth=2, fmt='none', color=colors[1])

    plt.legend()
    plt.loglog()
    plt.xlabel('Energy (eV)')
    plt.xlim([1e-3, 1e6])
    plt.ylabel('$I_n$ (1/primary/400cm$^2$)')
    plt.tight_layout()


    plt.show()


if __name__ == '__main__':
    main()
