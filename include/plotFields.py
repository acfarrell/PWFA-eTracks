import numpy as np
import matplotlib.colors as col
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib as mpl
import matplotlib.cm as cm
import plotSimTracks
import h5py as h5

def plot(field):
  fig, ax = plt.subplots()
  if field =="er":
    f = h5.File("../data/fields/e2-000066.h5","r")
    cbarname = ('$E_r$, Transverse Electric Field ($m_e c\omega_p / e$)')
    ax.set_title('Transverse Electric Field from OSIRIS')
    fname="../plots/trans_efield.png"
    dat = 'e2'
  if field =="ez":
    f = h5.File("../data/fields/e1-000066.h5","r")
    cbarname = ('$E_z$, Longitudinal Electric Field ($m_e c\omega_p / e$)')
    ax.set_title('Longitudinal Electric Field from OSIRIS')
    fname="../plots/long_efield.png"
    dat='e1'
  if field =="b":
    f = h5.File("../data/fields/b3-000066.h5","r")
    cbarname = ('$B_\phi$, Azimuthal Magnetic Field ($m_e c\omega_p / e$)')
    ax.set_title('Azimuthal Magnetic Field from OSIRIS')
    fname="../plots/phi_bfield.png"
    dat='b3'
  datasetNames = [n for n in f.keys()] #Two Datasets: AXIS and e2
  E = f[dat][:]
  E = E.astype(float)
  a1_bounds = f['AXIS']['AXIS1']#.astype(float)
  a2_bounds = f['AXIS']['AXIS2']#.astype(float)

  t0 = 858.95 
  xi = np.linspace(a1_bounds[0] - t0,a1_bounds[1] - t0,len(E[0]))
  r = np.linspace(a2_bounds[0],a2_bounds[1],len(E))

  #Make color axis of electric field
  colors = ax.pcolormesh(xi ,r,E,norm=col.SymLogNorm(linthresh=0.03,linscale=0.03,vmin=-E.max(),vmax=E.max()),cmap="RdBu_r")
  
  tick_locations=[x*0.01 for x in range(2,10)]+ [x*0.01 for x in range(-10,-1)] + [x*0.1 for x in range(-10,10)] +[ x for x in range(-10,10)]
  cbar = fig.colorbar(colors,ax=ax,ticks=tick_locations, format=ticker.LogFormatterMathtext())
  cbar.set_label(cbarname)
  ax.set_xlabel("$\\xi$ ($c/\omega_p$)")
  ax.set_ylabel('r ($c/\omega_p$)')

  plt.xlim(xi[0], xi[-1])
  plt.savefig(fname,transparent=True)
  plt.show()
plot("er")
plot("ez")
plot("b")
