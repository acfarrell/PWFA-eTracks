#Script for generating plots of electron trajectories

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sea
import matplotlib.cm as cm

def plot(r,t,E):

        plt.style.use('seaborn-poster')
        cmap = mpl.cm.get_cmap("YlGnBu")
        normalize = mpl.colors.Normalize(vmin=min(E),vmax=max(E))
        colors = [cmap(normalize(value)) for value in E]

        fig, ax = plt.subplots(figsize=(10,10))
        ax.plot(t,r)
        #cax, _ = mpl.colorbar.make_axes(ax)
        #cax.set_title("E Field (m_e c \omega_p / c)")
        #cbar = mpl.colorbar.ColorbarBase(cax,cmap=cmap, norm=normalize)
        ax.set_xlabel("time ($\omega_p^{-1}$)")
        ax.set_ylabel("radius ($c/\omega_p$)")
        ax.set_title("Electron Radial Trajectory")
        plt.show()
