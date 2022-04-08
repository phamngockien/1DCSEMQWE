import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

#================================================================================================
#load the (transmitter-receiver)TX-RX configurations
configs = pd.read_csv("1DCSEMQWE/TX_RX_configurations.csv")
#check the dipole type of the transmitter
print("Transmitter is: ")
print(configs.source_type[0])
#now get the TX-RX configurations 
cols = [x for x in configs.columns if x != 'source_type']
TX_RX = configs[cols].astype('float64')

#test the TX-RX data
#print(TX_RX)

# In this example we only consider the amplitude of Ex component.
# The validation for other components can be derived easily by changing the following code
#================================================================================================
#load the EM fields at the first frequency
our_prog_dat = pd.read_csv("1DCSEMQWE/CSEM_1D_RESPONSES_at_frequency_index_1.csv",dtype='float64')
empymod_qwe_dat = pd.read_csv("empymod/empymod_QWE.csv",usecols=["0" , "1"])
analytic_dat = pd.read_csv("analytical_solution_full_space/analytical_solution.csv")

#test the data
#print(empymod_qwe_dat['0'])
#print(analytic_dat.Ex_re)

#================================================================================================
#plot a profile of EM fields w.r.t the depths of receivers
#where offset = sqrt[(TX_x - RX_x)^2 + (TX_y - RX_y)^2]
#Here we divide by 1000. to convert the unit from meter to kilometer
depths = TX_RX.RX_z

#---------------------------------------------------------------------------
#the EM field in frequency domain is often illustrated by its amplitude and phase
# amplitude =  sqrt(real_component^2 + imaginary_component^2)  (V/m for electric field, A/m for magnetic field, T for magnetic flux density)
# phase =  arctan(imaginary_component/real_component)  (rad)
# note in this example, we use a unit moment dipole, thus there is no need to implement the normalization of the EM field w.r.t the dipole moment

# Now compute for the amplitude of Ex component
Ex_amp_our_prog = np.sqrt(np.add(np.square(our_prog_dat.Ex_re),np.square(our_prog_dat.Ex_im)))
Ex_amp_empymod_qwe = np.sqrt(np.add(np.square(empymod_qwe_dat['0']),np.square(empymod_qwe_dat['1'])))
Ex_amp_analytic = np.sqrt(np.add(np.square(analytic_dat.Ex_re),np.square(analytic_dat.Ex_im)))

#---------------------------------------------------------------------------
#Now generate an illustrative figure
mm = 1/25.4  # milimeters in inches
fig1 = plt.figure(figsize=(90*mm, 48*mm), facecolor='w')
#fig1.subplots_adjust(wspace=.4, hspace=.4)

#all left panels are for the amplitude
#plot the amplitude of the Ey at the two frequencies
#ax1 = plt.subplot(1, 1, 1)
plt.semilogy(depths, Ex_amp_analytic,'-',color='black', linewidth=1, label='analytic')
plt.semilogy(depths, Ex_amp_our_prog,'--o',color='black', markersize=5, markerfacecolor='None', linewidth=1, label='our_qwe')
plt.semilogy(depths, Ex_amp_empymod_qwe,':*',color='black', markersize=5, linewidth=1, label='previous_qwe')
plt.legend(loc='best', fontsize=9)
#ax1.set_ylabel('$log_{10}(V/m)$',fontsize=9)
#ax1.set_ylim([-14,-6])
plt.xticks(fontsize=9)
plt.yticks(fontsize=9)
plt.xlabel('Depth (m)',fontsize=9)
plt.ylabel('Amplitude of Ex(V/m)',fontsize=9)
#plt.xticks(np.arange(0,11,step=2),fontsize=9)
#plt.yticks(np.arange(-14,-5,step=2),fontsize=9)

plt.show()

fig1.savefig("Fig 5.eps",bbox_inches = 'tight', dpi=1000)
fig1.savefig("Fig 5.pdf",bbox_inches = 'tight', dpi=1000)