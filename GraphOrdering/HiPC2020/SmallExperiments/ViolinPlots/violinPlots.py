# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 00:38:37 2020

@author: reetb
"""

import os
os.chdir('C:/Local/Coursework/Research-HPCBIO-Lab/GraphOrdering/HiPC2020/Small Experiments/Violin Plots/')

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
plt.style.use('ggplot')


filename = 'arenas-email'
v = 1132
e = 5451

title = 'Input: ' + filename + ', n: ' + str(v) + ', m: ' + str(e)

df = pd.read_csv(filename + '.csv')

df = df.melt(var_name='Reordering Schemes', value_name='Gaps')

ax = sns.violinplot(x="Reordering Schemes", y="Gaps", data=df)

ax.set_title(title)

ax.get_figure().savefig(filename + '.png', dpi = 100)