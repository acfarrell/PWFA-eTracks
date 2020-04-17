import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import csv
import matplotlib.animation as animation

t = []
xi = []
z = []
r = []
p_z = []
p_r = []
E_z = []
E_r = []
B = []

def get_data(fname):
  with open(fname, 'r') as csvfile:
    plots = csv.reader(csvfile,delimiter='\t')
    next(plots)
    next(plots)
    for row in plots:
      t.append(float(row[0]))
      xi.append(float(row[1]))
      z.append(float(row[2]))
      r.append(float(row[3]))
      p_z.append(float(row[4]))
      p_r.append(float(row[4]))
      E_z.append(float(row[5]))
      E_r.append(float(row[6]))
      B.append(float(row[7]))

def get_xir(track):
  if track == 'min':
    get_data('data/TRACKS_data_Min_X2.txt')
  elif track == 'max':
    get_data('data/TRACKS_data_Max_X2.txt')
  else:
    get_data('data/TRACKS_data_Med_X2.txt')
    
  return xi,r

def animate(i, x, y, line, text):
  x_dat = x[:i*10] #select data range to display in frame
  y_dat = y[:i*10]
  #update line
  line.set_data(x_dat,y_dat)
  #update time stamp
  text.set_text("t = "+str(t[i*10])+" $1/\omega_p$")
def plot_r_v_xi():
  fig, ax = plt.subplots()
  ax.plot(xi,r,'k')
  plt.xlabel("$\\xi$ ($c/\omega_p$)")
  plt.ylabel("$r$ ($c/\omega_p$)")
  plt.title("Radial Trajectory of Simulated Tracks")
  plt.show()

def plot_E_v_r():
  fig, ax = plt.subplots()
  plt.ylim(-5,1)
  plt.xlim(-.5,2)
  plt.xlabel("$r$ ($c/\omega_p$)")
  plt.ylabel("$E_r$ ($m_e c \omega_p/e$)")
  plt.title("Radial Electric Field of Simulated Tracks")
  line, = ax.plot([0],[0],label="Simulated Data")
  text = ax.text(.75,.95,"",transform=ax.transAxes, ha="left",va="center")
  
  #plot theoretical curve
  x = np.linspace(.001,2,100)
  ax.plot(x,-1.0/4.0* x, 'k--',label="Theory")
  ax.legend(loc="lower right")
  def animateEvr(i):
    animate(i,r,E_r,line,text)

#  ani = animation.FuncAnimation(fig,animateEvr,frames=int(len(r)/10.0),interval=1,repeat=True)
#  plt.show()

