#Script for generating plots of electron trajectories

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sea
import matplotlib.cm as cm
def plot(E):
        
        plt.style.use('seaborn-poster')
        sea.heatmap(E, cmap = "YlGnBu")
        plt.show()
