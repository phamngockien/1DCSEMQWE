import numpy as np

file = open("1D_CSEM_INPUT.txt","w")

############################################################################
#type 0 for electric source, 1 for magnetic source
file.write("# DIPOLE TYPE: 0 \n")

file.write('\n')

############################################################################
#Number of quadrature points to evaluate the integral
# 101 is a very accurate value, you can change it for faster computation
# for e.g, 9
file.write("# NQUAD: 101 \n")

file.write('\n')

############################################################################
#Relative and absolute tolerances for convergence criteria of the QWE technique
# Herebelow are the value for the "Truth" EM field (Key, 2012)
# In practice, you may want to change them to a lesser accuracy for faster computation
# For e.g, RTOL: 1e-3, ATOL: 1e-6
file.write("# RTOL: 1e-12 \n")
file.write("# ATOL: 1e-30 \n")

file.write('\n')

############################################################################
# the number of transmitters, here is 1
nTX = 1
file.write("# TRANSMITTERS: ")
file.write(str(nTX)+'\n')

#generate the transmitters' parameters
''' 
    The azimuth ranges from [0, 360], being 0 in x-direction 
    and has positive value w.r.t counter clockwise direction.
    The dip ranges from [-90, 90] degrees with -90 is upward direction and vise versa
'''
#transmitters' positions, in this example [0, 0, 950]
TX_x = np.zeros(nTX).reshape(nTX,1) 
TX_y = np.zeros(nTX).reshape(nTX,1) 
TX_z = (950*np.ones(nTX)).reshape(nTX,1)
#moments, in this example all TX have a moment of 1
TX_m_list = np.ones(nTX).reshape(nTX,1)
#azimuths, in this example the azimuth is 45 degrees
TX_azith_list = (45*np.ones(nTX)).reshape(nTX,1)
#dips, in this example all TX are horizontal (i.e., dip = 45)
TX_dip_list = (45*np.ones(nTX)).reshape(nTX,1)

TX = np.concatenate((TX_x,TX_y,TX_z,TX_m_list,TX_azith_list,TX_dip_list),axis=1)

for row in TX:
    file.write(' '.join([str(a) for a in row]) + '\n')

file.write('\n')
    
############################################################################
# the number of frequencies and their values
nfreqs = 1
file.write("# FREQUENCIES: ")
file.write(str(nfreqs)+'\n')
freqs=[1]
for f in freqs:
    file.write(str(f)+'\n')

file.write('\n')

############################################################################
# the number of layers, their top depths (m), and resistivities (ohm-m).
'''
! This example is a three-layer model
! Air layer (10^12 ohm-m), it top depth should be set to a negative infinity as default (-1e150).  
! Seawater layer (0.3 ohm-m), extending from 0 - 1000 m depth.  
! The seafloor layer (1 ohm-m), extending from 1000-2000 m depth.
'''
nLayers = 3
file.write("# LAYERS: ")
file.write(str(nLayers)+'\n')

#top depths
depths = np.array([-1., 0, 1000]).reshape(nLayers,1)
#resistivities (ohm-m)
res = np.array([1e12, 0.3, 1.]).reshape(nLayers,1)
model = np.concatenate((depths,res),axis=1)
for row in model:
    file.write(' '.join([str(a) for a in row]) + '\n')
    
file.write('\n')
    
############################################################################
# the number of receivers, here is 1
nRX = 11
file.write("# RECEIVERS: ")
file.write(str(nRX)+'\n')
#positions, in this example [range from 0 to 10000, 0, 1000]
RX_x = np.arange(0., 11000.,1000.).reshape(nRX,1)
RX_y = np.zeros(nRX).reshape(nRX,1)
RX_z = (1000*np.ones(nRX)).reshape(nRX,1)
RX = np.concatenate((RX_x,RX_y,RX_z),axis=1)

for row in RX:
    file.write(' '.join([str(a) for a in row]) + '\n')
