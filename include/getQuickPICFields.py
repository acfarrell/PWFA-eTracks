# This file retrieves the fields from OSIRIS data files stored within the data/ folder.

import h5py as h5
import pandas
import numpy as np
import matplotlib.colors as col
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm

def getField(fname):
  f = h5.File('data/'+fname,"r")
  datasetNames = [n for n in f.keys()] #Two Datasets: AXIS and e2
  field = datasetNames[-1]
  Field_dat = f[field][:].astype(float)
  return Field_dat.transpose()

def axes(fname): 
  f = h5.File('data/'+fname,"r")
  datasetNames = [n for n in f.keys()] #Two Datasets: AXIS and e2
  field = datasetNames[-1]
  Field_dat = f[field][:].astype(float)
  a1_bounds = f['AXIS']['AXIS1']
  a2_bounds = f['AXIS']['AXIS2']

  xi_dat = np.linspace(a2_bounds[0] ,a2_bounds[1] ,len(Field_dat))
  r_dat = np.linspace(a1_bounds[0],a1_bounds[1],len(Field_dat[0]))
  

  return r_dat, xi_dat 

def transE(fname):
  return np.flip(getField(fname),1)

def longE(fname):
  return np.flip(getField(fname),1)

def phiB(fname):
  return np.flip(getField(fname),1)

def plotFields():
  r, xi = axes("fields/exslicexz_00000100.h5")
  Er = transE("fields/exslicexz_00000100.h5")
  Ez = longE("fields/ezslicexz_00000100.h5")
  Bphi = phiB("fields/byslicexz_00000100.h5")


  fig, ax = plt.subplots(3,sharex=True,figsize=(10, 15))
  #Make color axis of electric field
  colors = ax[0].pcolormesh(xi,r,Er,norm=col.SymLogNorm(linthresh=0.03,linscale=0.03,vmin=-Er.max(),vmax=Er.max()),cmap="RdBu_r")
    
  cbar = fig.colorbar(colors,ax=ax[0])
  cbar.set_label('Transverse Electric Field ($m_e c\omega_p / e$)')
  #ax[0].set_xlabel("$\\xi$ ($c/\omega_p$)")
  ax[0].set_ylabel('r ($c/\omega_p$)')
  ax[0].set_title('Fields from QuickPIC')
  ax[0].set_ylim(0, r[-1])

  colors = ax[1].pcolormesh(xi,r,Ez,norm=col.SymLogNorm(linthresh=0.03,linscale=0.03,vmin=-Er.max(),vmax=Er.max()),cmap="RdBu_r")
    
  cbar = fig.colorbar(colors,ax=ax[1])
  cbar.set_label('Longitudinal Electric Field ($m_e c\omega_p / e$)')
  #ax[1].set_xlabel("$\\xi$ ($c/\omega_p$)")
  ax[1].set_ylabel('r ($c/\omega_p$)')
  ax[1].set_ylim(0, r[-1])
  
  colors = ax[2].pcolormesh(xi,r,Bphi,norm=col.SymLogNorm(linthresh=0.03,linscale=0.03,vmin=-Er.max(),vmax=Er.max()),cmap="RdBu_r")
    
  cbar = fig.colorbar(colors,ax=ax[2])
  cbar.set_label('Azimuthal Magnetic Field ($m_e c\omega_p / e$)')
  ax[2].set_xlabel("$\\xi$ ($c/\omega_p$)")
  ax[2].set_ylabel('r ($c/\omega_p$)')
  ax[2].set_ylim(0, r[-1])
  plt.savefig("quickPIC_Fields.png",transparent=True)
  plt.show()
#plotFields()
