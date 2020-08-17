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

import getQuickPICFields as quickPIC
import getOsirisFields as osiris

#Definition of Constants
M_E = 9.109e-31                   #electron rest mass in kg
EC = 1.60217662e-19               #electron charge in C
EP_0 = 8.854187817e-12                #vacuum permittivity in C/(V m)
C = 299892458                     #speed of light in vacuum in m/s
N = 5e22                          #electron number density in 1/m^3
WP = np.sqrt(N*EC**2/(M_E*EP_0))   #plasma frequency in 1/s

r,xi, Er = quickPIC.spliceLowRes("quickPIC/exslicexz_00000139.h5")
r,xi, Ez = quickPIC.spliceLowRes("quickPIC/ezslicexz_00000139.h5")
#Er = osiris.transE('fields/exslicexz_00000100.h5')
#r,xi = osiris.axes('fields/exslicexz_00000100.h5')

def plot():
  dat = np.load('data.npz')
  xi = dat['xi']
  escaped = dat['esc']
  trail = dat['beam']
  fig, ax  = plt.subplots()
  E = Er

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
  ax.set_xlabel("$\\xi$ ($c/\omega_p$)")
  
  plt.xlim(xi[0], xi[-1])
  fn = "capturedRegion.png"
  plt.savefig(fn,dpi=300,transparent=True)
  plt.show()
  return


def plotBeams():
  dat = np.load('quickPIC/quickPIC_139.npz')
  #xi = dat['xi']
  escaped = dat['esc']
  trail = dat['trail']
  drive = dat['drive']
  n = 0
  print("Captured Electrons = ", len(trail))
  fig, axs  = plt.subplots(2,sharex=True)
  E = Er 
  
  xiWidth = (xi[1] - xi[0]) / WP
  spread = (max(trail) - min(trail))
  spread_SI = spread / WP

  print("Spread of trailing beam = ",round(spread,3),"c/wp, ",spread_SI," seconds")
  print("Bin width in xi = ", xiWidth)
  #define binary color map
  ternary_cmaplist = [(1.0,0.,0.,1.0),(0.,0.,0.,0.0),(1.0,1.0,0.0,1.0),(0.0,1.0,0.0,1.0)]
  ternary_cmap = mpl.colors.LinearSegmentedColormap.from_list('Custom cmap', ternary_cmaplist, 4)
  
  colors = axs[0].pcolormesh(xi,r,E,norm=col.SymLogNorm(linthresh=0.03,linscale=0.03,vmin=-E.max(),vmax=E.max()),cmap="RdBu_r")
  
  colors2 = axs[0].pcolormesh(xi ,r,escaped,cmap=ternary_cmap)
  
  tick_locations=[x*0.01 for x in range(2,10)]+ [x*0.01 for x in range(-10,-1)] + [x*0.1 for x in range(-10,10)] +[ x for x in range(-10,10)]
  cbar = fig.colorbar(colors,ax=axs.ravel().tolist(),ticks=tick_locations, format=ticker.LogFormatterMathtext())
  cbar.set_label('$E_r$, Transverse Electric Field ($m_e c\omega_p / e$)')

  axs[0].set_ylabel('r ($c/\omega_p$)')
  axs[0].set_title('Captured Electron Beam Test')
  axs[0].set_xlabel("$\\xi$ ($c/\omega_p$)")
  axs[0].set_ylim(0,3)
  #axs[0].set_xlim(0,8)
  #axs[1].set_xlim(0,8)
  plt.xlim(0,8)
  
  binwidth = (max(drive) - min(drive))/15.0
  nbins = int(round( (max(trail) - min(trail))/binwidth))
  axs[1].hist(drive, bins = 15,color='red',label="Driving Beam")
  axs[1].hist(trail, bins = nbins,color='blue',label="Trailing Beam")
  counts, bin_edges = np.histogram(trail,bins=nbins)
  axs[1].set_xlabel("$\\xi$ ($c/\omega_p$)")
  axs[1].legend()
  axs[1].set_ylim(None, counts.max()+20 )
  #plt.xlim(xi[0], xi[-1])
  fn = "capturedIonizedBeam_quickPIC_139.png"
  plt.savefig(fn,dpi=300,transparent=True)
  plt.show()
  return
plotBeams()

def plotBeamOverlaps():
  #Plot a closer look at the beam profiles with a square driving beam
  dat = np.load('data.npz')
  xi = dat['xi']
  escaped = dat['esc']
  trail = []#dat['beam']
  drive = []
  square_drive = [[],[],[],[],[],[],[],[]]
  square_trail = [[],[],[],[],[],[],[],[]]
  counts = [0,0,0,0,0,0,0,0]
  n = 0
  for i in range(len(r)):
    for j in range(len(xi)):
      if escaped[i,j] == 1 and r[i] < 0.2:
        drive.append(xi[j])
        trail.append(dat['beam'][j])
        n += 1
  while True:
    if counts[0] >= 200 and counts[5] >= 200 and counts[1] >= 200 and counts[2]>=200 and counts[3]>=200 and counts[4]>=200:
      break
    if counts[0]>=10:
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

def getSquareBeam():
  #Get initial and final electron location mean and error for 6 xi-slices of a square driving beam
  dat = np.load('data.npz')
  xi = dat['xi']
  escaped = dat['esc']
  trail = dat['beam']
  drive = []
  square_drive = [[],[],[],[],[],[],[],[]]
  square_trail = [[],[],[],[],[],[],[],[]]
  drive_tot=[]
  trail_tot=[]
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
          drive_tot.append(drive[i])
          trail_tot.append(trail[i])
          counts[nbin] += 1
  #Get Stats
  drive_means = []
  drive_err = []
  trail_means = []
  trail_err = []
  for nbin in range(6):
    trail_means.append(np.mean(square_trail[nbin]))
    drive_means.append(np.mean(square_drive[nbin]))
    trail_err.append(np.std(square_trail[nbin]))
    drive_err.append(np.std(square_drive[nbin]))
  return drive_means, drive_err, trail_means, trail_err

def plotBeamError():
  #Plot initial v. final xi position for square driving beam over a number of random trials
  fig, ax  = plt.subplots()#2,sharex=True)

  #Get Stats
  drive = []
  drive_err = []
  trail = []
  trail_err = []
  for nbin in range(100):
    d,derr,t,terr = getSquareBeam()
    trail.append(t)
    drive.append(d)
    trail_err.append(terr)
    drive_err.append(derr)
  #define binary color map
  ax.errorbar(np.mean(drive,axis=0),np.mean(trail,axis=0),np.mean(trail_err,axis=0),np.mean(drive_err,axis=0), fmt='ro',ecolor='b', ls='none',capsize=5,markersize=5,elinewidth=2, markeredgewidth=2)
  ax.set_xlabel("$\\xi_i$ ($c/\omega_p$)")
  ax.set_ylabel("$\\xi_f$ ($c/\omega_p$)")
  plt.title("Initial v. Final Electron Position")
  #plt.xlim(xi[0], xi[-1])
  #plt.ylim(xi[0],xi[-1])
  fn = "plots/error.png"
  plt.savefig(fn,dpi=300,transparent=True)
  plt.show()

  return

def plotBeamOverlapsLowR():
  #Plot a closer look at the beam profiles with a square driving beam
  dat = np.load('data.npz')
  xi = dat['xi']
  escaped = dat['esc']
  trail = []#dat['beam']
  drive = []
  square_drive = [[],[],[],[],[],[],[],[]]
  square_trail = [[],[],[],[],[],[],[],[]]
  counts = [0,0,0,0,0,0,0,0]
  n = 0
  for i in range(len(r)):
    for j in range(len(xi)):
      if escaped[i,j] == 1 and r[i] < 0.5:
        drive.append(xi[j])
        trail.append(dat['beam'][j])
        n += 1
  for i in range(len(drive)):
    start = min(drive)
    end = max(drive)
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
  trailBins = int(trail_range / (binw))*2
  #define binary color map
  ternary_cmaplist = [(1.0,0.,0.,1.0),(0.,0.,0.,0.0),(0.0,1.0,0.0,1.0)]
  ternary_cmap = mpl.colors.LinearSegmentedColormap.from_list('Custom cmap', ternary_cmaplist, 3)
  
  nbins = int(len(E[0])/20)
  counts, bin_edges = np.histogram(drive,bins=nbins)
  
  bin_colors = ['b','g','r','c','m','y','k','0.75']
  
  ax.hist(square_drive,12,color=bin_colors, stacked = True)
  ax.hist(square_trail,trailBins,color=bin_colors, stacked=True)#False,histtype='step', linewidth=1.5)
  ax.set_xlabel("$\\xi$ ($c/\omega_p$)")
  plt.title("Trailing Beam Longitudinal Distribution")
  plt.xlim(xi[0], xi[-1])
  plt.ylim(0, 500)
  fn = "plots/overlapLowR_DriveTrail.png"
  plt.savefig(fn,dpi=300,transparent=True)
  plt.show()

  return

def plotCorrespondence():
  dat = np.load('data.npz')
  xi = dat['xi']
  escaped = dat['esc']
  trail = []#dat['beam']
  drive = []
  correspondence = np.zeros((len(r),len(xi)))
  n = 0
  for i in range(len(r)):
    for j in range(len(xi)):
      if escaped[i,j] == 1:
        ri = r[i]
        xii = xi[j]
        xif = dat['beam'][j]
        drive.append(xi[j])
        trail.append(dat['beam'][j])
        cor = xif**2 - (xii**2 + ri**2)
        correspondence[i][j] = cor
        n += 1
  
  fig, ax  = plt.subplots()
  E = osiris.transE()
  correspondence = np.ma.masked_where(correspondence == 0, correspondence)
  cmap = plt.cm.OrRd
  cmap.set_bad(color = (1,1,1,0))

  colors = ax.pcolormesh(xi ,r,E,norm=col.SymLogNorm(linthresh=0.03,linscale=0.03,vmin=-E.max(),vmax=E.max()),cmap="RdBu_r")
  
  colors2 = ax.pcolormesh(xi ,r,correspondence, cmap=cmap)
  tick_locations=[x*0.01 for x in range(2,10)]+ [x*0.01 for x in range(-10,-1)] + [x*0.1 for x in range(-10,10)] +[ x for x in range(-10,10)]
  cbar = fig.colorbar(colors,ax=ax,ticks=tick_locations, format=ticker.LogFormatterMathtext())
  cbar.set_label('$E_r$, Transverse Electric Field ($m_e c\omega_p / e$)')
  cbar2 = fig.colorbar(colors2,ax=ax)
  cbar2.set_label('$\\xi_f^2 + r_f^2 - (\\xi_i^2 + r_i^2)$, Xinlu\'s Correspondence ($c^2/\omega_p^2$)')

  ax.set_ylabel('r ($c/\omega_p$)')
  ax.set_title('Captured Electron Ionization Positions')
  ax.set_xlabel("$\\xi$ ($c/\omega_p$)")
  
  plt.xlim(xi[0], xi[-1])
  fn = "PositionCorrespondence.png"
  plt.savefig(fn,dpi=300,transparent=True)
  plt.show()
  return
#plot()
#plotBeams()
#plotBeamOverlapsLowR()
#plotCorrespondence()
#plotBeamError()
