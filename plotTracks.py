#Script for generating plots of electron trajectories

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sea
import matplotlib.cm as cm
def plot(r,t,E):
        
        plt.style.use('seaborn-poster')
        #ax1 = sea.heatmap(E, cmap = "YlGnBu")
        
        #ax2 = ax1.twinx()

        plt.plot(t,r)
        plt.xlabel("time (1/w_p")
        plt.ylabel("radius (c/w_p)")
        plt.title("Electron Radial Trajectory")
        plt.show()
