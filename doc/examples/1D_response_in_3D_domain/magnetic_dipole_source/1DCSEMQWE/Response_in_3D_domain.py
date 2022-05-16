from turtle import width
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
fig = plt.figure(figsize=(190*mm, 190*mm), facecolor='w')
fig.subplots_adjust(wspace=.4, hspace=.4)

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

#======================================================================================================================
# set up the axes for the first plot (upper left)
ax1 = fig.add_subplot(2, 2, 1, projection='3d')          

# 3D plot has problem of the width of the vectors.
# I change to plot each vector instead of the whole data
x2 = x.reshape(np.product(x.shape))
y2 = y.reshape(np.product(y.shape))
z2 = z.reshape(np.product(z.shape))
Bx = Bx_re.reshape(np.product(x.shape))
By = By_re.reshape(np.product(x.shape))
Bz = Bz_re.reshape(np.product(x.shape))

for x1,y1,z1,Bx1,By1,Bz1 in zip(x2,y2,z2,Bx,By,Bz):
    ax1.quiver(x1, y1, z1, Bx1, By1, Bz1, pivot = 'middle', normalize=True, length=0.36, arrow_length_ratio = 0.45)

#ax1.quiver(x, y, z, Bx_re, By_re, Bz_re, arrow_length_ratio = 0.8, length=0.3, normalize=True, label='Real(B) (T)')
ax1.invert_zaxis()
#ax1.legend(loc='best', fontsize=9)
ax1.set_xlabel('$X (km)$', fontsize=9)
ax1.set_ylabel('$Y (km)$', fontsize=9)
ax1.set_zlabel('$Z (km)$', fontsize=9)
plt.xticks(np.arange(-3,4,step=2))
plt.yticks(np.arange(-3,4,step=2))
ax1.set_zticks(np.arange(-3,4,step=2))
plt.title('a. Inphase magnetic flux density (B) vectors',fontsize=9, fontweight="bold")

# set up the axes for the second plot (upper right)
ax2 = fig.add_subplot(2, 2, 2, projection='3d')

Ex = Ex_re.reshape(np.product(x.shape))
Ey = Ey_re.reshape(np.product(x.shape))
Ez = Ez_re.reshape(np.product(x.shape))

for x1,y1,z1,Ex1,Ey1,Ez1 in zip(x2,y2,z2,Ex,Ey,Ez):
    ax2.quiver(x1, y1, z1, Ex1, Ey1, Ez1, pivot = 'middle', normalize=True, length=0.36, arrow_length_ratio = 0.45)

#ax2.quiver(x, y, z, Ex_re, Ey_re, Ez_re, arrow_length_ratio = 0.45, length=0.36, normalize=True, label='Real(E) (V/m)',)
ax2.invert_zaxis()
#ax2.legend(loc='best', fontsize=9)
ax2.set_xlabel('$X (km)$', fontsize=9)
ax2.set_ylabel('$Y (km)$', fontsize=9)
ax2.set_zlabel('$Z (km)$', fontsize=9)
plt.xticks(np.arange(-3,4,step=2))
plt.yticks(np.arange(-3,4,step=2))
ax2.set_zticks(np.arange(-3,4,step=2))
plt.title('b. Inphase electric field (E) vectors',fontsize=9,fontweight="bold")

# cross section of a at y = 0 (B only has Bx and Bz components)
ax3 = fig.add_subplot(2, 2, 3)    

x3 = x[y==0]
z3 = z[y==0]
Bx = Bx_re[y==0]
Bz = Bz_re[y==0]
x3 = x3.reshape(np.product(x3.shape))
z3 = z3.reshape(np.product(z3.shape))
Bx = Bx.reshape(np.product(x3.shape))
# Since we invert y axis, multiply Bz with -1
Bz = Bz.reshape(np.product(x3.shape)) * (-1.)

lengths = np.sqrt(Bx**2 + Bz**2)
for x1,z1,Bx1,Bz1,l in zip(x3,z3,Bx,Bz,lengths):
    ax3.quiver(x1, z1, Bx1, Bz1, pivot = 'middle')
    
ax3.invert_yaxis()
ax3.set_xlabel('$X (km)$', fontsize=9)
ax3.set_ylabel('$Z (km)$', fontsize=9)
plt.xticks(np.arange(-3,4,step=1))
plt.yticks(np.arange(-3,4,step=1))
plt.title('c. Inphase B vectors at y = 0 km',fontsize=9, fontweight="bold")
plt.plot(0, 2.05, 'v',color='red', markersize=5, label='VMD source')
ax3.legend(fontsize=9, loc='upper right')

for ly in (0, 1, 2, 2.1):
    ax3.axhline(y=ly, color='black', linestyle='-', lw=0.5)

#ax3.quiver(x[y==0], y[y==0], z[y==0], Bx_re[y==0], By_re[y==0], Bz_re[y==0], arrow_length_ratio = 0.76, length=0.3, normalize=True, label='Real(B) (T)')
#ax3.quiver(x3, z3, Bx, Bz, label='Real(B) (T)')
#ax3.invert_zaxis()
#ax3.legend(loc='best', fontsize=9)
#ax3.set_xlabel('$X (km)$', fontsize=9)
#ax3.set_ylabel('$Y (km)$', fontsize=9)
#ax3.set_zlabel('$Z (km)$', fontsize=9)
#ax3.set_ylabel('$Z (km)$', fontsize=9)
#plt.xticks(np.arange(-3,4,step=2))
#plt.yticks(np.arange(-3,4,step=2))
#ax3.set_zticks(np.arange(-3,4,step=2))
#plt.title('c. B vectors at y = 0 m',fontsize=9, fontweight="bold")


# cross section of E vectors at z = 1 km, E(Ex,Ey,0)
ax4 = fig.add_subplot(2, 2, 4)   
x4 = x[z==1]
y4 = y[z==1]
Ex = Ex_re[z==1]
Ey = Ey_re[z==1]
x4 = x4.reshape(np.product(x4.shape))
y4 = y4.reshape(np.product(y4.shape))
Ex = Ex.reshape(np.product(x4.shape))
Ey = Ey.reshape(np.product(x4.shape))
for x1,y1,Ex1,Ey1 in zip(x4,y4,Ex,Ey):
    ax4.quiver(x1, y1, Ex1, Ey1, pivot = 'middle')

ax4.set_xlabel('$X (km)$', fontsize=9)
ax4.set_ylabel('$Y (km)$', fontsize=9)
plt.xticks(np.arange(-3,4,step=1))
plt.yticks(np.arange(-3,4,step=1))
plt.title('d. Inphase E vectors at z = 1 km',fontsize=9, fontweight="bold")

'''
ax4 = fig.add_subplot(2, 2, 4, projection='3d')          
ax4.quiver(x[z==0], y[z==0], z[z==0], Ex_re[z==0], Ey_re[z==0], Ez_re[z==0],arrow_length_ratio = 0.5, length=0.5, normalize=True,  label='Real(B) (T)')
ax4.invert_zaxis()
ax4.legend(loc='best', fontsize=9)
ax4.set_xlabel('$X (km)$', fontsize=9)
ax4.set_ylabel('$Y (km)$', fontsize=9)
ax4.set_zlabel('$Z (km)$', fontsize=9)
plt.xticks(np.arange(-3,4,step=2))
plt.yticks(np.arange(-3,4,step=2))
ax4.set_zticks(np.arange(-3,4,step=2))
plt.title('d. E vectors at z = 0 km',fontsize=9, fontweight="bold")
'''

#print(x[y==0].reshape((len(pos),len(pos))))
#print(z[y==0].reshape((len(pos),len(pos))))
plt.show()

# Save files in pdf and eps format
fig.savefig("Fig 6.eps",bbox_inches = 'tight', dpi=1000)
fig.savefig("Fig 6.pdf",bbox_inches = 'tight', dpi=1000)

'''
x = x.reshape(np.product(x.shape))
y = y.reshape(np.product(y.shape))
z = z.reshape(np.product(z.shape))

scale = 0.02
u = np.sin(np.pi * x) * np.cos(np.pi * y) * np.cos(np.pi * z)
v = -np.cos(np.pi * x) * np.sin(np.pi * y) * np.cos(np.pi * z)
w = np.sqrt(2.0 / 3.0) * np.cos(np.pi * x) * np.cos(np.pi * y) * np.sin(np.pi * z)
lengths = np.sqrt(x**2+y**2+z**2)

for x1,y1,z1,u1,v1,w1,l in zip(x,y,z,u,v,w,lengths):
    ax.quiver(x1, y1, z1, u1, v1, w1, pivot = 'middle', length=l*0.5)
'''