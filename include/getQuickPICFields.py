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

def getDepth(): 
  f = h5.File('data/eyslicexz_00000000.h5',"r")
  datasetNames = [n for n in f.keys()] #Two Datasets: AXIS and e2
  field = datasetNames[-1]
  Field_dat = f[field][:].astype(float)
  a1_bounds = f['AXIS']['AXIS1']
  a2_bounds = f['AXIS']['AXIS2']

  xi_dat = np.linspace(a2_bounds[0] ,a2_bounds[1] ,len(Field_dat))
  y_dat = np.linspace(a1_bounds[0],a1_bounds[1],len(Field_dat[0]))
  return (y_dat[1] - y_dat[0])
def spliceLowRes(fname):
  r, xi = axes(fname)
  field = getField(fname)
  #print("r-axis has ", len(r)," bins")
  #print("xi-axis has ", len(xi)," bins")

  lowRes_fname = fname[:8] + '_lowRes' + fname[8:]
  lowRes_r, lowRes_xi = axes(lowRes_fname)
  lowRes_field = getField(lowRes_fname)
  #print("Low res r-axis has ", len(lowRes_r)," bins")
  #print("Low res xi-axis has ", len(lowRes_xi)," bins")
  
  xi_IDX = (np.abs(lowRes_xi - 9)).argmin()
  for j in range(xi_IDX, len(lowRes_xi)):
    xi = np.append(xi, lowRes_xi[j])
    xi_slice = field.sum(1)[...,None]
    xi_slice.shape
    for i in range(len(field)):
      lowRes_IDX = int(i//4)
      xi_slice[i] = lowRes_field[lowRes_IDX,j]
    #print(len(field)," rows, " ,len(xi_slice))
    field = np.append(field, xi_slice,1)
    field.shape
  #print(field[:][512:])

  #xi = np.append(xi, lowRes_xi[xi_IDX:])
  #field = np.append(field, lowRes_field[:,xi_IDX:])
  
  return r, xi, np.flip(field,1)

def plotFields():
#r, xi = axes("fields/exslicexz_00000139.h5")
  r,xi,Er = spliceLowRes("quickPIC/exslicexz_00000139.h5")
  r,xi,Ez = spliceLowRes("quickPIC/ezslicexz_00000139.h5")
  r,xi,Bphi = spliceLowRes("quickPIC/byslicexz_00000139.h5")


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
