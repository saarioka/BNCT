import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

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

    print('Title:', title)
    print('Bins:', bins)
    print('Lower:', lower)
    print('Upper:', upper)
    print('Underflow:', df['entries'].iloc[0])
    print('Overflow:', df['entries'].iloc[-1])
    print('Entries:', df['entries'].sum())

    plt.figure(figsize=(10, 6))
    plt.title(title)
    x = np.linspace(float(lower), float(upper), int(bins))
    plt.step(x, df['entries'][1:-1])
    plt.savefig(file_name.replace('.csv', '.pdf'))

def main():
    read_csv('Run0_h1_E.csv')
    read_csv('Run0_h1_X.csv')
    plt.show()

if __name__ == '__main__':
    main()
