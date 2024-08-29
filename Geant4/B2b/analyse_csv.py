import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import yaml

def read_csv(file_name):
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

def main():
    read_csv('Run0_h1_E.csv')
    read_csv('Run0_h1_Edep.csv')
    read_csv('Run0_h1_X.csv')
    read_csv('Run0_h1_Y.csv')
    read_csv('Run0_h1_Z.csv')
    plt.show()

if __name__ == '__main__':
    main()
