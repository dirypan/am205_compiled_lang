import matplotlib.pyplot as plt
import numpy as np
import math
#load data and plot them using imshow

data = np.loadtxt('data_file.txt')
xs = data[:,0]
ys = data[:,1]
roots = data[:,2]
L = max(xs)
N = int(math.sqrt(len(xs)))
#resize roots to be a 2D array to use plt.imshow()
roots = np.reshape(roots,(-1,N))
# print(roots)
plt.imshow(roots,extent = [-L,L,-L,L])
plt.show()