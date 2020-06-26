# This file retrieves the fields from OSIRIS data files stored within the data/ folder.

import h5py as h5
import numpy as np

def getField(fname):
  f = h5.File('data/'+fname,"r")
  datasetNames = [n for n in f.keys()] #Two Datasets: AXIS and e2
  field = datasetNames[-1]
  Field_dat = f[field][:].astype(float)
  return Field_dat

def axes(fname,t0): 
  f = h5.File('data/'+fname,"r")
  datasetNames = [n for n in f.keys()] #Two Datasets: AXIS and e2
  field = datasetNames[-1]
  Field_dat = f[field][:].astype(float)
  a1_bounds = f['AXIS']['AXIS1']
  a2_bounds = f['AXIS']['AXIS2']

  xi_dat = np.linspace(a1_bounds[0] - t0,a1_bounds[1] - t0,len(Field_dat[0]))
  r_dat = np.linspace(a2_bounds[0],a2_bounds[1],len(Field_dat))
  

  return r_dat, xi_dat 

def transE(fname):
  return getField(fname)

def longE(fname):
  return getField(fname)

def phiB(fname):
  return getField(fname)
