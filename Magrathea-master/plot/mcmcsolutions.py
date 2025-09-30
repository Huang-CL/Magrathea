# -*- coding: utf-8 -*-
"""
Created on Mon Jul  7 11:45:21 2025

@author: drice
"""

import pandas as pd
import corner
import matplotlib.pyplot as plt

# Load the MCMC output file
df = pd.read_csv("mcmcoutput.txt", delim_whitespace=True)

# OPTIONAL: Drop duplicates if MCMC got stuck
df = df.drop_duplicates()

# Choose the columns to plot
params = ['Mass', 'fCore', 'fMantle', 'fWater', 'fAtm', 'RPlanet']
#params = ['Mass', 'fCore', 'fMantle', 'fWater', 'RPlanet']
data = df[params]

# Make the corner plot
figure = corner.corner(
    data,
    labels=params,
    show_titles=True,
    title_fmt=".3f",         # number format
    title_kwargs={"fontsize": 12},
    label_kwargs={"fontsize": 12}
)

plt.savefig('cornerplot.pdf',bbox_inches='tight')
plt.show()
plt.close()

# Define chain length+1
chain_length = 1001

# Create a chain index and step number for plotting
df['chain'] = df.index // chain_length
df['step'] = df.groupby('chain').cumcount()


# Plot tracer plots
fig, axes = plt.subplots(len(params), 1, figsize=(10, 2.5 * len(params)), sharex=True)

for i, param in enumerate(params):
    ax = axes[i]
    for chain_id, chain_df in df.groupby('chain'):
        ax.plot(chain_df['step'], chain_df[param], label=f'Chain {chain_id}', alpha=0.6)
    ax.set_ylabel(param)
    ax.legend(loc='upper right', fontsize='small')

axes[-1].set_xlabel("Step")
plt.tight_layout(rect=[0, 0, 1, 0.97])
plt.savefig('traceplot.pdf',bbox_inches='tight')
plt.show()
plt.close()