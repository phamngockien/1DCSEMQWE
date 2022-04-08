import empymod
import numpy as np
import pandas as pd 

x=300
y=0
z = np.arange(0., 3100.,100.)

#x-oriented electric dipole source at the origin
# [x, y, z, azimuth, dip] (dipole, infinitesimal small)
src = [0, 0, 0, 0, 0]

# layer resistivities
res = 1.
# layer boundaries
depth = []
# Frequency
freq = 1

#==========================================================================================================================
#THIS PART USING QWE
#msrc = False for ED, true for MD; mrec = False return E (V/m)
#xdirect true only use when src and rec in the same layer, else the program compute the field in wavenumer domain

Resutls=[]

#loop over the depth
for d in z:
    
    # x-directed dipole receiver-array: x, y, z, azimuth, dip
    rec_x = [x, y, d, 0, 0]

    inpdat_x = {'src': src, 'rec': rec_x,
                'depth': depth, 'freqtime': freq,
                'verb': 0, 'msrc':False,
                'ht':'qwe', 'htarg': {'nquad': 101}}


    E_x = empymod.bipole(**inpdat_x, res=res, mrec = False)

    #output results
    Resutls.append([E_x.real, E_x.imag])

pd.DataFrame(Resutls).to_csv('empymod_QWE.csv')

#==========================================================================================================================
#THIS PART USING fast Hankel transform (FHT : the digital filter method)
Resutls=[]

#loop over the depth
for d in z:
    
    # x-directed dipole receiver-array: x, y, z, azimuth, dip
    rec_x = [x, y, d, 0, 0]

    inpdat_x = {'src': src, 'rec': rec_x,
                'depth': depth, 'freqtime': freq,
                'verb': 0, 'msrc':False,
                'ht': 'dlf', 'htarg': {'pts_per_dec': -1}} #DLF settings

    E_x = empymod.bipole(**inpdat_x, res=res, mrec = False)

    #output results
    Resutls.append([E_x.real, E_x.imag])

pd.DataFrame(Resutls).to_csv('empymod_FHT.csv')

#==========================================================================================================================
#This part is for analytical solution
ab = 11 #both TX and RX are x-directed ED
Resutls=[]
#loop over the depth
for d in z:
    inp = {'src': [0, 0, 0], 'rec': [x, y, d], 'depth': [],
            'res': res, 'freqtime': freq, 'verb': 0}
    E_x = empymod.analytical(ab=ab, xdirect=True,**inp)
    #output results
    Resutls.append([E_x.real, E_x.imag])

pd.DataFrame(Resutls).to_csv('empymod_analytical.csv')