import sys
import csv
import numpy as np
import random
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as col
from mpl_toolkits.axes_grid1 import make_axes_locatable
from tempfile import TemporaryFile as tmp

import include.eTracks as eTracks 
import include.getOsirisFields as osiris

Er = osiris.transE()
r,xi,t0 = osiris.axes()

def plot():

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
    if counts[0] >= 500 and counts[5] >= 500 and counts[1] >= 500 and counts[2]>=500 and counts[3]>=500 and counts[4]>=500:
      break
    i = random.randint(0,len(drive)-1)
    start = 6.6
    end = 7.2
    if start <= drive[i] <= end:
      binw = (end - start)/6.0
      for nbin in range(6):
        if start + nbin*binw < drive[i] < start + (nbin+1)*binw and counts[nbin] <= 500:
          square_drive[nbin].append(drive[i])
          square_trail[nbin].append(trail[i])
          counts[nbin] += 1
  print("Captured Driving Electrons = ", n)
  plt.style.use('seaborn-poster')
  fig, ax  = plt.subplots()#2,sharex=True)
  E = osiris.transE()

  #define binary color map
  ternary_cmaplist = [(1.0,0.,0.,1.0),(0.,0.,0.,0.0),(0.0,1.0,0.0,1.0)]
  ternary_cmap = mpl.colors.LinearSegmentedColormap.from_list('Custom cmap', ternary_cmaplist, 3)
  
  #colors = axs[0].pcolormesh(z - t0 ,r,E,norm=col.SymLogNorm(linthresh=0.03,linscale=0.03,vmin=-E.max(),vmax=E.max()),cmap="RdBu_r")
  
  #colors2 = axs[0].pcolormesh(z - t0 ,r,escaped,cmap=ternary_cmap)
  
  #cbar = fig.colorbar(colors,ax=axs.ravel().tolist())
  #cbar.set_label('Transverse Electric Field ($m_e c\omega_p / e$)')

  #axs[0].set_ylabel('r ($c/\omega_p$)')
  #axs[0].set_title('Captured Electron Beam Test')
  
  nbins = int(len(E[0])/20)
  #axs[1].hist(drive, bins = nbins,color='red',label="Driving Beam")
  #axs[1].hist(trail, bins = nbins,color='blue',label="Trailing Beam")
  counts, bin_edges = np.histogram(drive,bins=nbins)
  
  bin_colors = ['b','g','r','c','m','y','k','0.75']
  
  ax.hist(square_drive,12,color=bin_colors, stacked = True)
  ax.hist(square_trail,60,color=bin_colors, stacked=True)#histtype='step', linewidth=1.5)
  ax.set_xlabel("$\\xi$ ($c/\omega_p$)")
  plt.title("Trailing Beam Longitudinal Distribution")
  plt.xlim(z[0]- t0, z[-1]-t0)
  plt.ylim(0, 600)
  fn = "escaped.png"
  plt.savefig(fn,transparent=True)
  plt.show()

  return

  

def main():
  nrows = len(Er)
  ncols = len(Er[0])
  escaped = np.zeros((nrows,ncols))

  driveBeamProf = []
  trailBeamProf = []
  
  fname = 'data.npz'
  eCount = 0
  num = 0
  for i in range(nrows):
    for j in range(int(ncols/2), ncols):
      if Er[i,j] < -0.5:
        eCount += 1
        plotTrack = False
        if eCount % 50 == 0:
          plotTrack = True
          num +=1
        
        escaped[i,j], xiPos = eTracks.GetTrajectory(r[i],0,0,xi[j],0,0,plotTrack, num)
        if escaped[i,j] == 1:
          trailBeamProf.append(xiPos)
      print('Row ',i, "/",nrows,", Column ", int(j - ncols/2),"/",int(ncols/2)," : ",eCount," electrons", end="\r", flush=True)
  np.savez(fname, r=r,xi=xi, esc=escaped, beam=trailBeamProf)
  plot()
main()
#plot()
