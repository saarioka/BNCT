import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import mplhep as hep

hep.style.use(hep.style.CMS)

TIME_FMT = r'%H:%M:%S.%f'

def main():
    df_current = pd.read_csv('HIP neutronitesti virrat.csv', skiprows=6, header=None, names=['Row', 'Time', 'Current'])
    df_current['Time'] = pd.to_datetime(df_current['Time'], format=TIME_FMT)
    df_current = df_current[['Time', 'Current']]
    df_current['ts'] = pd.to_datetime(df_current['Time']).astype(np.int64)
    df_current['Current'] -= df_current['Current'].min()
    df_current.dropna(inplace=True)
    df_current['Current'] *= 1e9
    
    df_doserate = pd.read_csv('HIP neutronitesti neutroniannosnopeudet.csv', skiprows=8, header=None, names=['Row', 'Time', 'DoseRate'])
    df_doserate['Time'] = pd.to_datetime(df_doserate['Time'], format=TIME_FMT)
    df_doserate = df_doserate[['Time', 'DoseRate']]
    df_doserate['ts'] = pd.to_datetime(df_doserate['Time']).astype(np.int64)
    df_doserate.dropna(inplace=True)

    ts_min = df_current['ts'].min()
    df_current['ts'] = (df_current['ts'] - ts_min) / 1e9
    df_doserate['ts'] = (df_doserate['ts'] - ts_min) / 1e9

    # find intervals where the beam was on, group them, calculate the mean current and dose rate and plot them as lines
    df_current['beam_on'] = df_current['Current'] > 4

    df_out = pd.DataFrame(columns=['ts', 'Current', 'DoseRate', 'DosePerCurrent'])

    plt.figure(figsize=(13, 7))
    ax1 = plt.gca()
    ax2 = ax1.twinx()
    ax1.scatter(df_current['ts'] / 60, df_current['Current'], 4, label='Current', color='dodgerblue')
    ax2.scatter(df_doserate['ts'] / 60, df_doserate['DoseRate'], 4, label='Dose rate', color='coral')
    # show periods with beam on
    for i, group in df_current.groupby((df_current['beam_on'] != df_current['beam_on'].shift()).cumsum()):
        if group['beam_on'].iloc[0]:
            t1 = group['ts'].iloc[0]
            t2 = group['ts'].iloc[-1]

            t1 += 15  # setting time for dose measurement

            # calculate mean current and dose rate
            mean_current = group['Current'].mean()
            mean_doserate = df_doserate[(df_doserate['ts'] >= t1) & (df_doserate['ts'] <= t2)]['DoseRate'].mean()

            # calculate dose per current
            dose_per_current = mean_doserate / mean_current

            df_out.loc[len(df_out)] = [t1, mean_current, mean_doserate, dose_per_current]

            ax1.axvspan(t1 / 60, t2 / 60, color='green', alpha=0.2, zorder=-1)
            ax1.plot([t1 / 60, t2 / 60], [mean_current, mean_current], color='blue', lw=3)
            ax2.plot([t1 / 60, t2 / 60], [mean_doserate, mean_doserate], color='red', lw=3)

    df_out.to_csv('acclab_output.csv', index=False)

    ax1.set_xlabel('Time [min]')
    ax1.set_ylabel('Current [nA]', color='blue')
    ax2.set_ylabel('Dose rate [µSv/h]', color='red')
    plt.tight_layout()

    plt.figure(figsize=(13, 7))
    plt.scatter(df_out['ts'], df_out['DosePerCurrent'])
    plt.xlabel('Current [nA]')
    plt.ylabel('Dose rate per current [µSv/h/nA]')
    plt.tight_layout()

    print(df_out)
    in_beam = df_out[3:7]
    print(in_beam)

    behind_moderator = df_out[7:10]
    print(behind_moderator)

    on_side = df_out[10:]
    print(on_side)

    fraction = behind_moderator['DosePerCurrent'].mean() / in_beam['DosePerCurrent'].mean()
    fraction_err = fraction * np.sqrt((behind_moderator['DosePerCurrent'].std() / behind_moderator['DosePerCurrent'].mean())**2 + (in_beam['DosePerCurrent'].std() / in_beam['DosePerCurrent'].mean())**2)
    print(f'Adding moderator plates: {fraction:.4f} ± {fraction_err:.4f}')

    fraction2 = on_side['DosePerCurrent'].mean() / in_beam['DosePerCurrent'].mean()
    fraction2_err = fraction2 * np.sqrt((on_side['DosePerCurrent'].std() / on_side['DosePerCurrent'].mean())**2 + (in_beam['DosePerCurrent'].std() / in_beam['DosePerCurrent'].mean())**2)
    print(f'Moving the possu to side WRT in beam without moderator: {fraction2:.4f} ± {fraction2_err:.4f}')

    ibpe_side = on_side['DoseRate'] / fraction2
    ibpe_mod = behind_moderator['DoseRate'] / fraction

    #possu_response_90deg = 1/1.33
    
    plt.figure(figsize=(13, 7))
    ax1 = plt.gca()
    ax2 = ax1.twinx()
    ax1.scatter(on_side['ts'], ibpe_side, label='Measured', color='turquoise')
    ax1.scatter(behind_moderator['ts'], ibpe_mod, label='Measured', color='turquoise')
    ax2.scatter(on_side['ts'], ibpe_side * 0.75, label='Fluence', color='magenta')
    ax2.scatter(behind_moderator['ts'], ibpe_mod * 0.75, label='Measured', color='magenta')

    ax1.set_ylim(0, None)
    ax2.set_ylim(0, None)
    ax1.set_xlabel('Time [s]')
    ax1.set_ylabel('Dose rate [µSv/h]')
    ax2.set_ylabel('Fluence [n/cm²]')
    plt.tight_layout()

    plt.show()

if __name__ == '__main__':
    main()
