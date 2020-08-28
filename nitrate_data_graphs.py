import os
import numpy as np
import matplotlib.pyplot as plt

data_cow = np.genfromtxt("nl_cows.txt", delimiter = ',', skip_header = 1)
x_cow = data_cow[:,0]
y_cow = data_cow[:,1]

data_nconc = np.genfromtxt("nl_n.txt", delimiter = ',', skip_header = 1)
x_nconc = data_nconc[:,0]
y_nconc = data_nconc[:,1]

plt.figure()
plt.plot(x_cow, y_cow)
plt.xlabel('Years')
plt.ylabel('Number of Cows')
plt.title('Number of cows in Southland over the years')
plt.savefig('cows.png')

plt.figure()
plt.plot(x_nconc, y_nconc)
plt.xlabel('Years')
plt.ylabel('Nitrate Concentration(mg/L)')
plt.title('Nitrate Concentration in Southland over the years')
plt.savefig('nitrate_concentration.png')

fig = plt.figure(figsize = (20,10))
ax = fig.add_subplot(111)
cows = ax.plot(x_cow, y_cow, color = 'red', label = 'Cows')
ax.set_xlabel('Years', fontsize = 20)
ax.set_ylabel('Number of Cows', fontsize = 20)
ax2 = ax.twinx()
n_conc = ax2.plot(x_nconc, y_nconc, color = 'blue', label = 'Nitrate Concentration')
ax2.set_ylabel('Nitrate Concentration(mg/L)', fontsize = 20)
plt.title('Annual Cattle Numbers in Southland and Nitrate Concentrations of Edendale Farm', fontsize = 20)
lines = cows + n_conc
labels = [l.get_label() for l in lines]
ax.legend(lines, labels,loc = 0)
fig.savefig('nitrate_conc_vs_cows.png', dpi = 200)