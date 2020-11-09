# create a function that displays the output when run this way:
#	python layerview.py ouput.dat

import sys
import os
from time import time
import subprocess
import struct
import numpy as np
import matplotlib
import math
import matplotlib.pyplot as plt

from mpl_toolkits.axes_grid1 import ImageGrid

def intensity(E):
    Econj = np.conj(E)
    I = np.sum(E*Econj, axis=-1)
    return np.real(I)

#evaluate a solved homogeneous substrate
# Returns a complex NxMx3 array representing the cross section of the field at Y=0
def evaluate(Depths, k, d, n0, sz, Pt, Pr, X, Y, Z):
    Depths = np.array(Depths)
    sz = np.array(sz)
    Pt = np.array(Pt)
    Pr = np.array(Pr)
    s = np.array(d) * n0
    #allocate space for layer indices
    LI = np.zeros(Z.shape, dtype=np.int)
            
    #find the layer index for each sample point
    L = len(Depths)
    LI[Z < Depths[0]] = 0
    for l in range(L-1):
        idx = np.logical_and(Z > Depths[l], Z <= Depths[l+1])
        LI[idx] = l
        LI[Z > Depths[-1]] = L - 1
    
    #calculate the appropriate phase shift for the wave transmitted through the layer
    Ph_t = np.exp(1j * k * sz[LI] * (Z - Depths[LI]))
    
    #calculate the appropriate phase shift for the wave reflected off of the layer boundary
    LIp = LI + 1
    LIp[LIp >= L] = 0
    Ph_r = np.exp(-1j * k * sz[LI] * (Z - Depths[LIp]))
    Ph_r[LI >= L-1] = 0
    
    #calculate the phase shift based on the X and Y positions
    Ph_xy = np.exp(1j * k * (s[0] * X + s[1] * Y))
    
    #apply the phase shifts
    Et = Pt[:, LI] * Ph_t[:, :]
    Er = Pr[:, LI] * Ph_r[:, :]
    
    #add everything together coherently
    E = (Et + Er) * Ph_xy[:, :]
    
    #return the electric field
    return np.moveaxis(E, 0, -1)
        
class planewave:
    def __int__(self):
        self.LAYERS = 0                          #Number of layers.                          int
        self.depths = []                         #z positions of layers. [1, 5, ..., 10]     double
        self.k0 = 0.0                              #wavenumber at free space.                  double
        self.d = []                              #direction of propogation. [0.5, 0]         double
        self.n0 = 0.0+0.0j                           #the refractive index of the first layer.   complex<double>
        self.sz = []                             #z-component of propagation for each layer. complex<double>
        self.Pt = [[] for i in range(3)]         #transmission                               complex<double>
        self.Pr = [[],[],[]]                     #refraction                                 complex<double>
        
# display a binary file produced using the coupled wave C code
def layer(strc):    
    f = open(strc, "rb")

    # create an empty plane wave structure
    L = planewave()
    L.depths = []
    L.d = []
    L.sz = []
    L.Pt = [[],[],[]]
    L.Pr = [[],[],[]]

    # open the input file for reading
    file_bytes = os.path.getsize(strc)

    # calculate the number of layers in the sample
    L.LAYERS = int((file_bytes/8-5)/15)

    # load the raw layer data into the plane wave structure
    data_raw = struct.unpack('d' * (15*L.LAYERS+5), f.read((15*L.LAYERS+5)* 8))
    data = np.asarray(data_raw)    
    L.k0 = data[0]
    L.d.append(data[1])
    L.d.append(data[2])
    L.n0 = complex(data[3], data[4])

    # load each layer's plane waves from the binary file
    for i in range(L.LAYERS):
        L.depths.append(data[5+15*i])
        L.sz.append(complex(data[6+15*i], data[7+15*i]))
        L.Pt[0].append(complex(data[8+15*i], data[9+15*i]))
        L.Pt[1].append(complex(data[15*i+10], data[15*i+11]))
        L.Pt[2].append(complex(data[15*i+12], data[15*i+13]))
        L.Pr[0].append(complex(data[15*i+14], data[15*i+15]))
        L.Pr[1].append(complex(data[15*i+16], data[15*i+17]))
        L.Pr[2].append(complex(data[15*i+18], data[15*i+19]))

    N = 512															# simulation resolution NxM
    M = 1024
    #DAVID: Don't hard-code the dimensions - you'll have to calculate them based on the sample information in the file
    D = [-110, 110, 0, 60]											# dimensions of the simulation
    x = np.linspace(D[2], D[3], N)									# set the sample points for the simulation
    z = np.linspace(D[0], D[1], M)
    [X, Z] = np.meshgrid(x, z)										# create a mesh grid to evaluate layers
    Y = np.zeros(X.shape)

    # evaluate the field across all layers
    E = evaluate(L.depths, L.k0, L.d, L.n0, L.sz, L.Pt, L.Pr, X, Y, Z)
    Er = np.real(E)
    I = intensity(E)

    plt.set_cmap("afmhot")											# set the color map
    plt.subplot(1, 4, 1)
    plt.imshow(Er[:, :, 0], extent=(D[3], D[2], D[1], D[0]))
    #plt.colorbar()
    plt.title("Ex")

    plt.subplot(1, 4, 2)
    plt.imshow(Er[:, :, 1], extent=(D[3], D[2], D[1], D[0]))
    #plt.colorbar()
    plt.title("Ey")

    plt.subplot(1, 4, 3)
    plt.imshow(Er[:, :, 2], extent=(D[3], D[2], D[1], D[0]))
    #plt.colorbar()
    plt.title("Ez")

    plt.subplot(1, 4, 4)
    plt.imshow(I, extent=(D[3], D[2], D[1], D[0]))
    plt.colorbar()
    plt.title("I")
    
    #fig = plt.figure(1, (5, 10))
    #plt.set_cmap("afmhot")
    #matplotlib.rcParams.update({'font.size': 10})
    #grid = ImageGrid(fig, rect = 211, nrows_ncols = (1, 3), axes_pad = 0.2, label_mode = "1", cbar_mode = "single", cbar_size = "18%")
    #Title = ["Ex", "Ey", "Ez"]
    #for i in range(3):
        # grid[i].axis('off')
    #    im = grid[i].imshow(Er[..., i], extent=(D[3], D[2], D[1], D[0]), interpolation="nearest")
    #    grid[i].set_title(Title[i])
    #grid.cbar_axes[0].colorbar(im)
    #plt.title("E")
    #plt.subplot(212)
    #plt.imshow(I, extent=(D[3], D[2], D[1], D[0]))
    #plt.title("I")
    #plt.colorbar()
    plt.show()

# function displays usage text to the console
def usage():
	print("Usage:")
	print("     layerview input.dat")

if __name__ == '__main__':
    start = time()
    if len(sys.argv) < 2:				# if there are no command line arguments
    	usage()							# display the usage text
    	exit()							# exit
    else:
    	layer(sys.argv[1])				# otherwise display the given data file

    end = time()
    print("The elapsed time is " + str(end - start) + " s. ")