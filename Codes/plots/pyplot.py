import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


df=pd.read_csv('../../data/Ti/Ti48_EHFB.dat',delimiter=r'\s+')

plt.figure(figsize=(6,6))

plt.plot(df['bet2t'],df['Etot'])

plt.savefig('Ti48_EHFB.eps')

