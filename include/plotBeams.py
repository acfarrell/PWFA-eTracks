import sys
import csv
import numpy as np
import random
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as col
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
from tempfile import TemporaryFile as tmp

import getOsirisFields as osiris

Er = osiris.transE()
r,xi,t0 = osiris.axes()

def plot():
  dat = np.load('data.npz')
  xi = dat['xi']
  escaped = dat['esc']
  trail = dat['beam']
  fig, ax  = plt.subplots()
  E = osiris.transE()

  #define binary color map
  ternary_cmaplist = [(1.0,0.,0.,1.0),(0.,0.,0.,0.0),(0.0,1.0,0.0,1.0)]
  ternary_cmap = mpl.colors.LinearSegmentedColormap.from_list('Custom cmap', ternary_cmaplist, 3)
  
  colors = ax.pcolormesh(xi ,r,E,norm=col.SymLogNorm(linthresh=0.03,linscale=0.03,vmin=-E.max(),vmax=E.max()),cmap="RdBu_r")
  
  colors2 = ax.pcolormesh(xi ,r,escaped,cmap=ternary_cmap)
  tick_locations=[x*0.01 for x in range(2,10)]+ [x*0.01 for x in range(-10,-1)] + [x*0.1 for x in range(-10,10)] +[ x for x in range(-10,10)]
  cbar = fig.colorbar(colors,ax=ax,ticks=tick_locations, format=ticker.LogFormatterMathtext())
  cbar.set_label('$E_r$, Transverse Electric Field ($m_e c\omega_p / e$)')

  ax.set_ylabel('r ($c/\omega_p$)')
  ax.set_title('Captured Electron Ionization Positions')
  
  plt.xlim(xi[0], xi[-1])
  fn = "capturedRegion.png"
  plt.savefig(fn,dpi=300,transparent=True)
  plt.show()
  return

def plotBeams():
  dat = np.load('data.npz')
  xi = dat['xi']
  escaped = dat['esc']
  trail = dat['beam']
  drive = []
  n = 0
  for i in range(len(r)):
    for j in range(len(xi)):
      if escaped[i,j] == 1:
        drive.append(xi[j])
        n += 1
  print("Captured Driving Electrons = ", n)
  plt.style.use('seaborn-poster')
  fig, axs  = plt.subplots(2,sharex=True)
  E = osiris.transE()

  #define binary color map
  ternary_cmaplist = [(1.0,0.,0.,1.0),(0.,0.,0.,0.0),(0.0,1.0,0.0,1.0)]
  ternary_cmap = mpl.colors.LinearSegmentedColormap.from_list('Custom cmap', ternary_cmaplist, 3)
  
  colors = axs[0].pcolormesh(xi ,r,E,norm=col.SymLogNorm(linthresh=0.03,linscale=0.03,vmin=-E.max(),vmax=E.max()),cmap="RdBu_r")
  
  colors2 = axs[0].pcolormesh(xi ,r,escaped,cmap=ternary_cmap)
  
  tick_locations=[x*0.01 for x in range(2,10)]+ [x*0.01 for x in range(-10,-1)] + [x*0.1 for x in range(-10,10)] +[ x for x in range(-10,10)]
  cbar = fig.colorbar(colors,ax=axs.ravel().tolist(),ticks=tick_locations, format=ticker.LogFormatterMathtext())
  cbar.set_label('$E_r$, Transverse Electric Field ($m_e c\omega_p / e$)')

  axs[0].set_ylabel('r ($c/\omega_p$)')
  axs[0].set_title('Captured Electron Beam Test')
  
  nbins = int(len(E[0])/20)
  axs[1].hist(drive, bins = nbins,color='red',label="Driving Beam")
  axs[1].hist(trail, bins = nbins,color='blue',label="Trailing Beam")
  counts, bin_edges = np.histogram(trail,bins=nbins)
  axs[1].set_xlabel("$\\xi$ ($c/\omega_p$)")
  axs[1].legend()
  axs[1].set_ylim(None, counts.max()+20 )
  plt.xlim(xi[0], xi[-1])
  fn = "capturedBeam.png"
  plt.savefig(fn,dpi=300,transparent=True)
  plt.show()
  return


def plotBeamOverlaps():
  #Plot a closer look at the beam profiles with a square driving beam
  dat = np.load('data.npz')
  xi = dat['xi']
  escaped = dat['esc']
  trail = dat['beam']
  drive = []
  square_drive = [[],[],[],[],[],[],[],[]]
  square_trail = [[],[],[],[],[],[],[],[]]
  counts = [0,0,0,0,0,0,0,0]
  n = 0
  for i in range(len(r)):
    for j in range(len(xi)):
      if escaped[i,j] == 1:
        drive.append(xi[j])
        n += 1
  while True:
    if counts[0] >= 200 and counts[5] >= 200 and counts[1] >= 200 and counts[2]>=200 and counts[3]>=200 and counts[4]>=200:
      break
    i = random.randint(0,len(drive)-1)
    start = 6.8
    end = 7.1
    if start <= drive[i] <= end:
      binw = (end - start)/6.0
      for nbin in range(6):
        if start + nbin*binw < drive[i] < start + (nbin+1)*binw and counts[nbin] <= 200:
          square_drive[nbin].append(drive[i])
          square_trail[nbin].append(trail[i])
          counts[nbin] += 1
  print("Captured Driving Electrons = ", n)
  plt.style.use('seaborn-poster')
  fig, ax  = plt.subplots()#2,sharex=True)
  E = osiris.transE()

  #Find number of bins for trailing beam hist
  maxes = []
  mins = []
  for nbin in range(6):
    maxes.append(max(square_trail[nbin]))
    mins.append(min(square_trail[nbin]))

  trail_range = max(maxes) - min(mins)
  trailBins = int(trail_range / (binw))
  #define binary color map
  ternary_cmaplist = [(1.0,0.,0.,1.0),(0.,0.,0.,0.0),(0.0,1.0,0.0,1.0)]
  ternary_cmap = mpl.colors.LinearSegmentedColormap.from_list('Custom cmap', ternary_cmaplist, 3)
  
  nbins = int(len(E[0])/20)
  counts, bin_edges = np.histogram(drive,bins=nbins)
  
  bin_colors = ['b','g','r','c','m','y','k','0.75']
  
  ax.hist(square_drive,6,color=bin_colors, stacked = True)
  ax.hist(square_trail,trailBins,color=bin_colors, stacked=True)#histtype='step', linewidth=1.5)
  ax.set_xlabel("$\\xi$ ($c/\omega_p$)")
  plt.title("Trailing Beam Longitudinal Distribution")
  plt.xlim(xi[0], xi[-1])
  plt.ylim(0, 400)
  fn = "plots/overlap.png"
  plt.savefig(fn,dpi=300,transparent=True)
  plt.show()

  return
plot()
plotBeams()
plotBeamOverlaps()
