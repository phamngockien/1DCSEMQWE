'''
We validate the improvement of our implementation of the QWE method with the previous attemp
in empymod (Werthmüller, 2017).
It is noted that empymod use digital filter method (DLF) for computing the EM field, 
and the use of QWE is to validate the result from the DLF.
The algorithm of QWE in empymod is based on the work of Key (2012)

We use a full space model of 1 Ohm.m for the validation.
The TX-RX configuration is: 
+ A x-oriented electric diplole locates at the origin (0,0,0)
+ The receivers are at x = 300, y = 0, and z = 0 to 3 km.
'''


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
# In practice, you may want to change them to a lesser accuracy
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
#position: a x-oriented electric dipole at the origin
TX_x = np.zeros(nTX).reshape(nTX,1) 
TX_y = np.zeros(nTX).reshape(nTX,1) 
TX_z = np.zeros(nTX).reshape(nTX,1) 
#moments, in this example all TX have a moment of 1
TX_m_list = np.ones(nTX).reshape(nTX,1)
#azimuths, in this example all TX are directed in y direction
TX_azith_list = (0.*np.ones(nTX)).reshape(nTX,1)
#dips, in this example all TX are horizontal (i.e., dip = 0)
TX_dip_list = (0*np.ones(nTX)).reshape(nTX,1)

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
nLayers = 1
file.write("# LAYERS: ")
file.write(str(nLayers)+'\n')

'''
! This example is a full space model of 1 Ohm.m
'''
#top depths
depths = np.array([-1e150]).reshape(nLayers,1)
#resistivities (ohm-m)
res = np.array([ 1.]).reshape(nLayers,1)
model = np.concatenate((depths,res),axis=1)
for row in model:
    file.write(' '.join([str(a) for a in row]) + '\n')
    
file.write('\n')
    
############################################################################
# the number of receivers, here is 31
nRX = 31
file.write("# RECEIVERS: ")
file.write(str(nRX)+'\n')
#positions, in this example [300, 0, range from zero to 3000] meter
RX_x = (300.*np.ones(nRX)).reshape(nRX,1) 
RX_y = np.zeros(nRX).reshape(nRX,1)
RX_z = (np.arange(0,3100,100)).reshape(nRX,1)
RX = np.concatenate((RX_x,RX_y,RX_z),axis=1)

for row in RX:
    file.write(' '.join([str(a) for a in row]) + '\n')