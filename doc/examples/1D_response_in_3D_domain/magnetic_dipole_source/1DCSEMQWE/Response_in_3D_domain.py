import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from mpl_toolkits.mplot3d import axes3d

#================================================================================================
#load the (transmitter-receiver)TX-RX configurations
configs = pd.read_csv("TX_RX_configurations.csv")
#check the dipole type of the transmitter
print("Transmitter is: ")
print(configs.source_type[0])
#now get the TX-RX configurations 
cols = [x for x in configs.columns if x != 'source_type']
TX_RX = configs[cols].astype('float64')

#test the TX-RX data
#print(TX_RX)
#print(TX_RX.RX_y)

#================================================================================================
#load the EM fields at the first frequency
EMfield = pd.read_csv("CSEM_1D_RESPONSES_at_frequency_index_1.csv",dtype='float64')



#================================================================================================
mm = 1/25.4  # milimeters in inches
fig = plt.figure(figsize=(190*mm, 100*mm), facecolor='w')
#fig.subplots_adjust(wspace=.4, hspace=.4)

#In this example, the x-, y-, and z-positions are arranged in the same range
#from -3000 to 3000 meter
pos = np.arange(-3000., 3100., 500.)/1000.

x, y, z = np.meshgrid(pos,
                      pos,
                      pos, indexing='ij')
nx = len(pos)
ny = len(pos)
nz = len(pos)

#pre-allocate the direction data for the real vector (inphase)
Bx_re = np.zeros((nx,ny,nz))
By_re = np.zeros((nx,ny,nz))
Bz_re = np.zeros((nx,ny,nz))

Ex_re = np.zeros((nx,ny,nz))
Ey_re = np.zeros((nx,ny,nz))
Ez_re = np.zeros((nx,ny,nz))

#pre-allocate the direction data for the imaginary vector (quadrature phase)
Bx_im = np.zeros((nx,ny,nz))
By_im = np.zeros((nx,ny,nz))
Bz_im = np.zeros((nx,ny,nz))

Ex_im = np.zeros((nx,ny,nz))
Ey_im = np.zeros((nx,ny,nz))
Ez_im = np.zeros((nx,ny,nz))

#Get the data
for i in range(nx):
    for j in range(ny):
        for k in range(nz):
          index = i * ny * nz + j * nz + k
          #data[i,j][k] = E[index]
          Bx_re[i,j][k] = EMfield.Bx_re[index]
          By_re[i,j][k] = EMfield.By_re[index]
          Bz_re[i,j][k] = EMfield.Bz_re[index]
          
          Bx_im[i,j][k] = EMfield.Bx_im[index]
          By_im[i,j][k] = EMfield.By_im[index]
          Bz_im[i,j][k] = EMfield.Bz_im[index]
          
          Ex_re[i,j][k] = EMfield.Ex_re[index]
          Ey_re[i,j][k] = EMfield.Ey_re[index]
          Ez_re[i,j][k] = EMfield.Ez_re[index]

# set up the axes for the first plot
ax1 = fig.add_subplot(1, 2, 1, projection='3d')          
ax1.quiver(x, y, z, Bx_re, By_re, Bz_re, length=0.3, normalize=True, label='Real(B) (T)')
ax1.invert_zaxis()
ax1.legend(loc='best', fontsize=9)
ax1.set_xlabel('$X (km)$', fontsize=9)
ax1.set_ylabel('$Y (km)$', fontsize=9)
ax1.set_zlabel('$Z (km)$', fontsize=9)
plt.xticks(np.arange(-3,4,step=2))
plt.yticks(np.arange(-3,4,step=2))
ax1.set_zticks(np.arange(-3,4,step=2))

# set up the axes for the second plot
ax2 = fig.add_subplot(1, 2, 2, projection='3d')
ax2.quiver(x, y, z, Ex_re, Ey_re, Ez_re, length=0.3, normalize=True, label='Real(E) (V/m)')
ax2.invert_zaxis()
ax2.legend(loc='best', fontsize=9)
ax2.set_xlabel('$X (km)$', fontsize=9)
ax2.set_ylabel('$Y (km)$', fontsize=9)
ax2.set_zlabel('$Z (km)$', fontsize=9)
plt.xticks(np.arange(-3,4,step=2))
plt.yticks(np.arange(-3,4,step=2))
ax2.set_zticks(np.arange(-3,4,step=2))

plt.show()
# Save files in pdf and eps format
fig.savefig("Fig 6.eps",bbox_inches = 'tight', dpi=1000)
fig.savefig("Fig 6.pdf",bbox_inches = 'tight', dpi=1000)