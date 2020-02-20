# This file retrieves the fields from OSIRIS data files stored within the data/ folder.

import h5py as h5
import numpy as np

def getField(fname):
  f = h5.File(fname,"r")
  datasetNames = [n for n in f.keys()] #Two Datasets: AXIS and e2
  field = datasetNames[-1]
  Field_dat = f[field][:].astype(float)
  return Field_dat

def axes(): 
  f = h5.File('data/EField_r.h5',"r")
  datasetNames = [n for n in f.keys()] #Two Datasets: AXIS and e2
  field = datasetNames[-1]
  Field_dat = f[field][:].astype(float)
  a1_bounds = f['AXIS']['AXIS1']
  a2_bounds = f['AXIS']['AXIS2']
  
  z_dat = np.linspace(a1_bounds[0],a1_bounds[1],len(Field_dat[0]))
  r_dat = np.linspace(a2_bounds[0],a2_bounds[1],len(Field_dat))
  
  t0 = 858.95 #time at which field data was simulated, constant for all fields

  return r_dat, z_dat, t0 

def transE():
  return getField('data/EField_r.h5')

def longE():
  return getField('data/EField_z.h5')

def phiB():
  return getField('data/BField_phi.h5')
