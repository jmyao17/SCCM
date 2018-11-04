import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from sys import argv

def main():
    script,filename,figname=argv
    print(filename,figname)
    df=pd.read_csv(filename,delimiter=r'\s+')
    plt.figure(figsize=(6,6))
    plt.plot(df['bet2t'],df['Etot'])
    plt.savefig(figname)

if __name__ == '__main__':
      main()




