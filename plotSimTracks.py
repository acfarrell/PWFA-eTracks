import numpy as np
import matplotlib.pyplot as plt
import csv

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
    plots = csv.reader(csvfile,delimiter-'  ')
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
                    
