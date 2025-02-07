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
    
    df_doserate = pd.read_csv('HIP neutronitesti neutroniannosnopeudet.csv', skiprows=8, header=None, names=['Row', 'Time', 'DoseRate'])
    df_doserate['Time'] = pd.to_datetime(df_doserate['Time'], format=TIME_FMT)
    df_doserate = df_doserate[['Time', 'DoseRate']]
    df_doserate['ts'] = pd.to_datetime(df_doserate['Time']).astype(np.int64)
    df_doserate.dropna(inplace=True)

    ts_min = df_current['ts'].min()
    df_current['ts'] = (df_current['ts'] - ts_min) / 1e9
    df_doserate['ts'] = (df_doserate['ts'] - ts_min) / 1e9

    df_current['Current_smooth'] = df_current['Current'].rolling(window=20).mean()
    df_doserate['DoseRate_smooth'] = df_doserate['DoseRate'].rolling(window=20).mean()

    # find intervals where the beam was on, group them, calculate the mean current and dose rate and plot them as lines
    df_current['beam_on'] = df_current['Current'] > 4e-9

    plt.figure(figsize=(13, 7))
    ax1 = plt.gca()
    ax2 = ax1.twinx()
    ax1.scatter(df_current['ts'] / 60, df_current['Current'] * 1e9, 3, label='Current', color='blue')
    # show periods with beam on
    for i, group in df_current.groupby((df_current['beam_on'] != df_current['beam_on'].shift()).cumsum()):
        if group['beam_on'].iloc[0]:
            ax1.axvspan(group['ts'].iloc[0] / 60, group['ts'].iloc[-1] / 60, color='green', alpha=0.3)

    ax2.scatter(df_doserate['ts'] / 60, df_doserate['DoseRate'], 3, label='Dose rate', color='red')
    ax1.set_xlabel('Time [min]')
    ax1.set_ylabel('Current [nA]', color='blue')
    ax2.set_ylabel('Dose rate [µSv/h]', color='red')

    # shared legend for all axes
    #lines, labels = ax1.get_legend_handles_labels()
    #lines2, labels2 = ax2.get_legend_handles_labels()
    #ax2.legend(lines + lines2, labels + labels2, loc='upper center')

    plt.show()

if __name__ == '__main__':
    main()
