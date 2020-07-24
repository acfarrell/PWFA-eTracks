import numpy as np
import matplotlib.colors as col
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
import plotSimTracks
import h5py as h5
import pandas

def plot(fname,t0):
  f = h5.File(fname,"r")
  print(list(f.keys()))
  datasetNames = [n for n in f.keys()] #Two Datasets: AXIS and e2
  E = f['exslicexz'][:]
  E = E.astype(float)
  a1_bounds = f['AXIS']['AXIS1']#.astype(float)
  a2_bounds = f['AXIS']['AXIS2']#.astype(float)

  xi = np.linspace(a2_bounds[0],a2_bounds[1],len(E))
  r = np.linspace(a1_bounds[0],a1_bounds[1],len(E[0]))

  fig, ax = plt.subplots()
  #Make color axis of electric field
  colors = ax.pcolormesh(xi,r,E.transpose(),norm=col.SymLogNorm(linthresh=0.03,linscale=0.03,vmin=-E.max(),vmax=E.max()),cmap="RdBu_r")
    
  cbar = fig.colorbar(colors,ax=ax)
  cbar.set_label('Electric Field ($m_e c\omega_p / e$)')
  ax.set_xlabel("$\\xi$ ($c/\omega_p$)")
  ax.set_ylabel('r ($c/\omega_p$)')
  ax.set_title('Longitudinal Electric Field from OSIRIS')

  #plt.xlim(xi[0], xi[-1])
  plt.savefig("transverse_efield.png",transparent=True)
  plt.show()
plot("data/fields/exslicexz_00000100.h5",100*14.7)
