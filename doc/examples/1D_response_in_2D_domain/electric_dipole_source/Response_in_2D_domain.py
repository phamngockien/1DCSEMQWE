from locale import normalize
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

#In this example, y = 0 meters, 
#NOTE here we scaled by 1000 to convert meter to kilometer
#from -3000 to 3000 meter, spacing 100 meter, in the x-direction
posx = np.arange(-3000., 3100., 100.)/1000.
#from -1000 to 5000 meter, spacing 100 meter, in the z-direction
posz = np.arange(-1000., 5100., 100.)/1000.

# A two-dimensional meshgrid
x, z = np.meshgrid(posx,
                   posz, indexing='ij')

nx = len(posx)
nz = len(posz)

#pre-allocate the direction data for the real vector (inphase)
Bx_re = np.zeros((nx,nz))
By_re = np.zeros((nx,nz))
Bz_re = np.zeros((nx,nz))

Ex_re = np.zeros((nx,nz))
Ey_re = np.zeros((nx,nz))
Ez_re = np.zeros((nx,nz))

#pre-allocate the direction data for the imaginary vector (quadrature phase)
Bx_im = np.zeros((nx,nz))
By_im = np.zeros((nx,nz))
Bz_im = np.zeros((nx,nz))

Ex_im = np.zeros((nx,nz))
Ey_im = np.zeros((nx,nz))
Ez_im = np.zeros((nx,nz))

#Get the data
#In this example, Bx, Bz, and Ey are zero (see the response CSV file)
for i in range(nx):
    for k in range(nz):
        index = i * nz + k
        #data[i,k] = E or H at [index]
        Bx_re[i,k] = EMfield.Bx_re[index]
        By_re[i,k] = EMfield.By_re[index]
        Bz_re[i,k] = EMfield.Bz_re[index]
         
        Bx_im[i,k] = EMfield.Bx_im[index]
        By_im[i,k] = EMfield.By_im[index]
        Bz_im[i,k] = EMfield.Bz_im[index]
          
        Ex_re[i,k] = EMfield.Ex_re[index]
        Ey_re[i,k] = EMfield.Ey_re[index]
        Ez_re[i,k] = EMfield.Ez_re[index]
          
        Ex_im[i,k] = EMfield.Ex_im[index]
        Ey_im[i,k] = EMfield.Ey_im[index]
        Ez_im[i,k] = EMfield.Ez_im[index]

#We then plot the field in (Oxz) section

args = {'cmap': 'viridis', 'levels': 500, 'vmin': -18, 'vmax': -6}

# The first plot is for real E(Ex,0,Ez)
ax1 = fig.add_subplot(2, 2, 1)  
ax1.set_title('real E',fontsize=9)
#amplitude of the real E field
cp = plt.contourf(x, z, np.log10(np.sqrt(Ex_re**2 + Ez_re**2)))

#normalize the real E vectors      
u = Ex_re/np.sqrt(Ex_re**2 + Ez_re**2)
v = Ez_re/np.sqrt(Ex_re**2 + Ez_re**2)

#plot the vectors
skip=(slice(None,None,3),slice(None,None,3))
plt.quiver(x[skip], z[skip], u[skip], v[skip])

#inner interface boundaries
#plt.axhline(y=1., color='black', linestyle='-')
#plt.axhline(y=2., color='black', linestyle='-')
#plt.axhline(y=2.1, color='black', linestyle='-')

plt.gca().invert_yaxis()
ax1.set_xlabel('$X (km)$', fontsize=9)
ax1.set_ylabel('$Z (km)$', fontsize=9)
plt.xticks(np.arange(-3,4,step=1))
plt.yticks(np.arange(-1,6,step=1))


# The second plot is for imaginary E(Ex,0,Ez)
ax2 = fig.add_subplot(2, 2, 2)  
ax2.set_title('imaginary E',fontsize=9)
#amplitude of the imaginary E field
cp = plt.contourf(x, z, np.log10(np.sqrt(Ex_im**2 + Ez_im**2)))

#normalize the imaginary E vectors      
u = Ex_im/np.sqrt(Ex_im**2 + Ez_im**2)
v = Ez_im/np.sqrt(Ex_im**2 + Ez_im**2)

#plot the vectors
skip=(slice(None,None,3),slice(None,None,3))
plt.quiver(x[skip], z[skip], u[skip], v[skip])

plt.gca().invert_yaxis()
ax2.set_xlabel('$X (km)$', fontsize=9)
ax2.set_ylabel('$Z (km)$', fontsize=9)
plt.xticks(np.arange(-3,4,step=1))
plt.yticks(np.arange(-1,6,step=1))

clb = fig.colorbar(cp, ax=(ax1, ax2), orientation="horizontal",
ticks=np.arange(args['vmin']+2, args['vmax'], 2),
pad=0.2, aspect=30, fraction=0.1, shrink=0.7)
clb.set_label(label='Amplitude [log10(V/m)]')

# The third plot is for real B(0,By,0)
ax3 = fig.add_subplot(2, 2, 3)  
ax3.set_title('real B',fontsize=9)
#amplitude of the real B field
cp = plt.contourf(x, z, np.log10(np.abs(By_re)))

plt.gca().invert_yaxis()
ax3.set_xlabel('$X (km)$', fontsize=9)
ax3.set_ylabel('$Z (km)$', fontsize=9)
plt.xticks(np.arange(-3,4,step=1))
plt.yticks(np.arange(-1,6,step=1))

# The fourth plot is for imaginary B(0,By,0)
ax4 = fig.add_subplot(2, 2, 4)  
ax4.set_title('imaginary B',fontsize=9)
#amplitude of the imaginary B field
cp = plt.contourf(x, z, np.log10(np.abs(By_im)))

plt.gca().invert_yaxis()
ax4.set_xlabel('$X (km)$', fontsize=9)
ax4.set_ylabel('$Z (km)$', fontsize=9)
plt.xticks(np.arange(-3,4,step=1))
plt.yticks(np.arange(-1,6,step=1))


clb = fig.colorbar(cp, ax=(ax3, ax4), orientation="horizontal",
ticks=np.arange(args['vmin']-4, args['vmax']-4, 2),
pad=0.2, aspect=30, fraction=0.1, shrink=0.7)
clb.set_label(label='Amplitude [log10(T)]')

plt.show()
# Save figure in pdf format
fig.savefig("EM_field_of_a_HED_in_2D_domain.pdf",bbox_inches = 'tight', dpi=1000)