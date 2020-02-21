import numpy as np
import matplotlib.colors as col
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sea
import matplotlib.cm as cm
import plotSimTracks
import h5py as h5

def plot():
  f = h5.File("../data/fields/e2-000066.h5","r")
  datasetNames = [n for n in f.keys()] #Two Datasets: AXIS and e2
  E = f['e2'][:]
  E = E.astype(float)
  a1_bounds = f['AXIS']['AXIS1']#.astype(float)
  a2_bounds = f['AXIS']['AXIS2']#.astype(float)

  t0 = 858.95 
  xi = np.linspace(a1_bounds[0] - t0,a1_bounds[1] - t0,len(E[0]))
  r = np.linspace(a2_bounds[0],a2_bounds[1],len(E))
  plt.style.use('seaborn-poster')

  fig, ax = plt.subplots()
  #Make color axis of electric field
  colors = ax.pcolormesh(xi ,r,E,norm=col.SymLogNorm(linthresh=0.03,linscale=0.03,vmin=-E.max(),vmax=E.max()),cmap="RdBu_r")
    
  cbar = fig.colorbar(colors,ax=ax)
  cbar.set_label('Transverse Electric Field ($m_e c\omega_p / e$)')
  ax.set_xlabel("$\\xi$ ($c/\omega_p$)")
  ax.set_ylabel('r ($c/\omega_p$)')
  ax.set_title('Transverse Electric Field from OSIRIS')

  plt.xlim(xi[0], xi[-1])
  plt.savefig("../plots/trans_efield.png",transparent=True)
  plt.show()
plot()
