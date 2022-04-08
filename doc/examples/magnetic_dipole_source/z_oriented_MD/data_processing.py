from this import d
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

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
print(TX_RX)
#print(TX_RX.RX_x)

#================================================================================================
#load the EM fields at the first frequency
EMfield_freq1 = pd.read_csv("CSEM_1D_RESPONSES_at_frequency_index_1.csv",dtype='float64')

#load the EM fields at the second frequency
EMfield_freq2 = pd.read_csv("CSEM_1D_RESPONSES_at_frequency_index_2.csv",dtype='float64')

#test the loaded EM fields
print(EMfield_freq1)
#print(EMfield_freq1.Ex_re)
#print(EMfield_freq2)
#print(EMfield_freq2.Ex_re)

#Using the test print above, we can see that the main components of the EM fields are:
#Ey, Bx, and Bz.
#Other components are smaller than the main component for at least 10 order of magnitude.
#Thus, the main signal that observed at the receivers depends mostly on the main components of the EM field.
#Below we only illustrate the main components. Your can also try for the other components similarly.
#NOTE: for a different TX-RX configuration, the main components of the EM field may change accordingly.

#================================================================================================
#plot a profile of EM fields w.r.t the depth of RX
depth = TX_RX.RX_z

#---------------------------------------------------------------------------
#the EM field in frequency domain is often illustrated by its amplitude and phase
# amplitude =  sqrt(real_component^2 + imaginary_component^2)  (V/m for electric field, A/m for magnetic field, T for magnetic flux density)
# phase =  arctan(imaginary_component/real_component)  (rad)
# note in this example, we use a unit moment dipole, thus there is no need to implement the normalization of the EM field w.r.t the dipole moment

# Now compute for the first observation frequency (Ey, Bx, and Bz)
Ey_amplitude_freq_1 = np.sqrt(np.add(np.square(EMfield_freq1.Ey_re),np.square(EMfield_freq1.Ey_im)))
Ey_phase_freq_1 = np.arctan2(EMfield_freq1.Ey_im,EMfield_freq1.Ey_re)


Bx_amplitude_freq_1 = np.sqrt(np.add(np.square(EMfield_freq1.Bx_re),np.square(EMfield_freq1.Bx_im)))
Bx_phase_freq_1 = np.arctan2(EMfield_freq1.Bx_im,EMfield_freq1.Bx_re)

Bz_amplitude_freq_1 = np.sqrt(np.add(np.square(EMfield_freq1.Bz_re),np.square(EMfield_freq1.Bz_im)))
Bz_phase_freq_1 = np.arctan2(EMfield_freq1.Bz_im,EMfield_freq1.Bz_re)

#---------------------------------------------------------------------------
#similarly for the second observation frequency and so on...
#in this example, we only use 2 frequencies. The following code completes the rest part of computation
Ey_amplitude_freq_2 = np.sqrt(np.add(np.square(EMfield_freq2.Ey_re),np.square(EMfield_freq2.Ey_im)))
Ey_phase_freq_2 = np.arctan2(EMfield_freq2.Ey_im,EMfield_freq2.Ey_re)

Bx_amplitude_freq_2 = np.sqrt(np.add(np.square(EMfield_freq2.Bx_re),np.square(EMfield_freq2.Bx_im)))
Bx_phase_freq_2 = np.arctan2(EMfield_freq2.Bx_im,EMfield_freq2.Bx_re)

Bz_amplitude_freq_2 = np.sqrt(np.add(np.square(EMfield_freq2.Bz_re),np.square(EMfield_freq2.Bz_im)))
Bz_phase_freq_2 = np.arctan2(EMfield_freq2.Bz_im,EMfield_freq2.Bz_re)

#---------------------------------------------------------------------------
#Now generate an illustrative figure
mm = 1/25.4  # milimeters in inches
fig1 = plt.figure(figsize=(190*mm, 120*mm), facecolor='w')
fig1.subplots_adjust(wspace=.4, hspace=.4)

#all left panels are for the amplitude
#plot the amplitude of the Ey at the two frequencies
ax1 = plt.subplot(3, 2, 1)
plt.plot(depth, np.log10(Ey_amplitude_freq_1),'-',color='black', linewidth=1, label='Ey @ 0.25 Hz')
plt.plot(depth, np.log10(Ey_amplitude_freq_2),'--',color='black', linewidth=1, label='Ey @ 1. Hz')
plt.legend(loc='best', fontsize=7)
ax1.set_ylabel('$log_{10}(V/m)$',fontsize=10)
#ax1.set_ylim([-14,-6])
plt.xticks(fontsize=9)
plt.yticks(fontsize=9)
#plt.xticks(np.arange(0,11,step=2),fontsize=9)
#plt.yticks(np.arange(-14,-5,step=2),fontsize=9)

#plot the amplitude of the Bx at the two frequencies
ax3 = plt.subplot(3, 2, 3)
plt.plot(depth, np.log10(Bx_amplitude_freq_1),'-',color='black', linewidth=1, label='Bx @ 0.25 Hz')
plt.plot(depth, np.log10(Bx_amplitude_freq_2),'--',color='black', linewidth=1, label='Bx @ 1. Hz')
plt.legend(loc='best', fontsize=7)
ax3.set_ylabel('$log_{10}(T)$',fontsize=10)
plt.xticks(fontsize=9)
plt.yticks(fontsize=9)
#plt.xticks(np.arange(0,11,step=2),fontsize=9)
#plt.yticks(np.arange(-14,-7,step=2),fontsize=9)

#plot the amplitude of the Bz at the two frequencies
ax5 = plt.subplot(3, 2, 5)
plt.plot(depth, np.log10(Bz_amplitude_freq_1),'-',color='black', linewidth=1, label='Bz @ 0.25 Hz')
plt.plot(depth, np.log10(Bz_amplitude_freq_2),'--',color='black', linewidth=1, label='Bz @ 1. Hz')
plt.legend(loc='best', fontsize=7)
ax5.set_ylabel('$log_{10}(T)$',fontsize=10)
ax5.set_xlabel('Depth (m)',fontsize=10)
#plt.xticks(np.arange(0,11,step=2),fontsize=9)
#plt.yticks(np.arange(-18,-9,step=2),fontsize=9)

#---------------------------------------------------------------------------
#all right panels are for the phase
#plot the phase of the Ey at the two frequencies
ax2 = plt.subplot(3, 2, 2)
plt.plot(depth, Ey_phase_freq_1,'-',color='black', linewidth=1, label='Ey @ 0.25 Hz')
plt.plot(depth, Ey_phase_freq_2,'--',color='black', linewidth=1, label='Ey @ 1. Hz')
plt.legend(loc='best', fontsize=7)
plt.ylabel('Phase (rad)',fontsize=10)
plt.xticks(fontsize=9)
plt.yticks(fontsize=9)
#plt.xticks(np.arange(0,11,step=2),fontsize=9)
#plt.yticks(np.arange(-3,4,step=2),fontsize=9)

#plot the phase of the Bx at the two frequencies
ax4 = plt.subplot(3, 2, 4)
plt.plot(depth, Bx_phase_freq_1,'-',color='black', linewidth=1, label='Bx @ 0.25 Hz')
plt.plot(depth, Bx_phase_freq_2,'--',color='black', linewidth=1, label='Bx @ 1. Hz')
plt.legend(loc='best', fontsize=7)
plt.ylabel('Phase (rad)',fontsize=10)
plt.xticks(fontsize=9)
plt.yticks(fontsize=9)
#plt.xticks(np.arange(0,11,step=2),fontsize=9)
#plt.yticks(np.arange(-0.5,0.4,0.2),fontsize=9)

#plot the phase of the Bz at the two frequencies
ax6 = plt.subplot(3, 2, 6)
plt.plot(depth, Bz_phase_freq_1,'-',color='black', linewidth=1, label='Bz @ 0.25 Hz')
plt.plot(depth, Bz_phase_freq_2,'--',color='black', linewidth=1, label='Bz @ 1. Hz')
plt.legend(loc='best', fontsize=7)
plt.ylabel('Phase (rad)',fontsize=10)
ax6.set_xlabel('Depth (m)',fontsize=10)
plt.xticks(fontsize=9)
plt.yticks(fontsize=9)
#plt.xticks(np.arange(0,11,step=2),fontsize=9)
#plt.yticks(np.arange(-3,4,step=2),fontsize=9)

plt.show()

# Save files in pdf format
fig1.savefig("z_oriented_MD.pdf",bbox_inches = 'tight', dpi=1000)