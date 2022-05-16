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
#print(TX_RX)

#================================================================================================
#load the EM fields at the first frequency
EMfield_freq1 = pd.read_csv("CSEM_1D_RESPONSES_at_frequency_index_1.csv",dtype='float64')

#load the EM fields at the second frequency
EMfield_freq2 = pd.read_csv("CSEM_1D_RESPONSES_at_frequency_index_2.csv",dtype='float64')

#test the loaded EM fields
#print(EMfield_freq1)

#================================================================================================
#plot a profile of EM fields w.r.t the depth of RX
depth = TX_RX.RX_z

#---------------------------------------------------------------------------
#the EM field in frequency domain is often illustrated by its amplitude and phase
# amplitude =  sqrt(real_component^2 + imaginary_component^2)  (V/m for electric field, A/m for magnetic field, T for magnetic flux density)
# phase =  arctan(imaginary_component/real_component)  (rad)
# note in this example, we use a unit moment dipole, thus there is no need to implement the normalization of the EM field w.r.t the dipole moment

# Now compute for the first observation frequency 
Ex_amplitude_freq_1 = np.sqrt(np.add(np.square(EMfield_freq1.Ex_re),np.square(EMfield_freq1.Ex_im)))
Ex_phase_freq_1 = np.arctan2(EMfield_freq1.Ex_im,EMfield_freq1.Ex_re)

Ey_amplitude_freq_1 = np.sqrt(np.add(np.square(EMfield_freq1.Ey_re),np.square(EMfield_freq1.Ey_im)))
Ey_phase_freq_1 = np.arctan2(EMfield_freq1.Ey_im,EMfield_freq1.Ey_re)

Ez_amplitude_freq_1 = np.sqrt(np.add(np.square(EMfield_freq1.Ez_re),np.square(EMfield_freq1.Ez_im)))
Ez_phase_freq_1 = np.arctan2(EMfield_freq1.Ez_im,EMfield_freq1.Ez_re)

Bx_amplitude_freq_1 = np.sqrt(np.add(np.square(EMfield_freq1.Bx_re),np.square(EMfield_freq1.Bx_im)))
Bx_phase_freq_1 = np.arctan2(EMfield_freq1.Bx_im,EMfield_freq1.Bx_re)

By_amplitude_freq_1 = np.sqrt(np.add(np.square(EMfield_freq1.By_re),np.square(EMfield_freq1.By_im)))
By_phase_freq_1 = np.arctan2(EMfield_freq1.By_im,EMfield_freq1.By_re)

Bz_amplitude_freq_1 = np.sqrt(np.add(np.square(EMfield_freq1.Bz_re),np.square(EMfield_freq1.Bz_im)))
Bz_phase_freq_1 = np.arctan2(EMfield_freq1.Bz_im,EMfield_freq1.Bz_re)

# Now compute for the second observation frequency 
Ex_amplitude_freq_2 = np.sqrt(np.add(np.square(EMfield_freq2.Ex_re),np.square(EMfield_freq2.Ex_im)))
Ex_phase_freq_2 = np.arctan2(EMfield_freq2.Ex_im,EMfield_freq2.Ex_re)

Ey_amplitude_freq_2 = np.sqrt(np.add(np.square(EMfield_freq2.Ey_re),np.square(EMfield_freq2.Ey_im)))
Ey_phase_freq_2 = np.arctan2(EMfield_freq2.Ey_im,EMfield_freq2.Ey_re)

Ez_amplitude_freq_2 = np.sqrt(np.add(np.square(EMfield_freq2.Ez_re),np.square(EMfield_freq2.Ez_im)))
Ez_phase_freq_2 = np.arctan2(EMfield_freq2.Ez_im,EMfield_freq2.Ez_re)

Bx_amplitude_freq_2 = np.sqrt(np.add(np.square(EMfield_freq2.Bx_re),np.square(EMfield_freq2.Bx_im)))
Bx_phase_freq_2 = np.arctan2(EMfield_freq2.Bx_im,EMfield_freq2.Bx_re)

By_amplitude_freq_2 = np.sqrt(np.add(np.square(EMfield_freq2.By_re),np.square(EMfield_freq2.By_im)))
By_phase_freq_2 = np.arctan2(EMfield_freq2.By_im,EMfield_freq2.By_re)

Bz_amplitude_freq_2 = np.sqrt(np.add(np.square(EMfield_freq2.Bz_re),np.square(EMfield_freq2.Bz_im)))
Bz_phase_freq_2 = np.arctan2(EMfield_freq2.Bz_im,EMfield_freq2.Bz_re)

#---------------------------------------------------------------------------
#Now generate an illustrative figure
mm = 1/25.4  # milimeters in inches
fig1 = plt.figure(figsize=(190*mm, 190*mm), facecolor='w')
fig1.subplots_adjust(wspace=.4, hspace=.4)

#all left panels are for the amplitude
#plot the amplitude of E at the two frequencies (see input file)
ax1 = plt.subplot(2, 2, 1)
plt.plot(depth, np.log10(Ex_amplitude_freq_1),'-',color='black', linewidth=1, label='Ex @ 0.25 Hz')
plt.plot(depth, np.log10(Ey_amplitude_freq_1),':',color='black', linewidth=1, label='Ey @ 0.25 Hz')
plt.plot(depth, np.log10(Ez_amplitude_freq_1),'--',color='black', linewidth=1, label='Ez @ 0.25 Hz')
plt.plot(depth, np.log10(Ex_amplitude_freq_2),'-',color='blue', linewidth=1, label='Ex @ 1 Hz')
plt.plot(depth, np.log10(Ey_amplitude_freq_2),':',color='blue', linewidth=1, label='Ey @ 1 Hz')
plt.plot(depth, np.log10(Ez_amplitude_freq_2),'--',color='blue', linewidth=1, label='Ez @ 1 Hz')
plt.legend(loc='best', fontsize=7)
ax1.set_ylabel('$log_{10}(V/m)$',fontsize=10)
#ax1.set_ylim([-14,-6])
plt.xticks(fontsize=9)
plt.yticks(fontsize=9)
#plt.xticks(np.arange(0,11,step=2),fontsize=9)
#plt.yticks(np.arange(-14,-5,step=2),fontsize=9)

#plot the amplitude of B at frequency of 1 Hz (see input file)
ax3 = plt.subplot(2, 2, 3)
plt.plot(depth, np.log10(Bx_amplitude_freq_1),'-',color='black', linewidth=1, label='Bx @ 0.25 Hz')
plt.plot(depth, np.log10(By_amplitude_freq_1),':',color='black', linewidth=1, label='By @ 0.25 Hz')
plt.plot(depth, np.log10(Bz_amplitude_freq_1),'--',color='black', linewidth=1, label='Bz @ 0.25 Hz')
plt.plot(depth, np.log10(Bx_amplitude_freq_2),'-',color='blue', linewidth=1, label='Bx @ 1 Hz')
plt.plot(depth, np.log10(By_amplitude_freq_2),':',color='blue', linewidth=1, label='By @ 1 Hz')
plt.plot(depth, np.log10(Bz_amplitude_freq_2),'--',color='blue', linewidth=1, label='Bz @ 1 Hz')
plt.legend(loc='best', fontsize=7)
ax3.set_ylabel('$log_{10}(T)$',fontsize=10)
plt.xticks(fontsize=9)
plt.yticks(fontsize=9)
#plt.xticks(np.arange(0,11,step=2),fontsize=9)
#plt.yticks(np.arange(-14,-7,step=2),fontsize=9)
ax3.set_xlabel('Depth (m)',fontsize=10)

#---------------------------------------------------------------------------
#all right panels are for the phase
#plot the phase of  E at frequency of 1 Hz (see input file)
ax2 = plt.subplot(2, 2, 2)
plt.plot(depth, Ex_phase_freq_1,'-',color='black', linewidth=1, label='Ex @ 0.25 Hz')
plt.plot(depth, Ey_phase_freq_1,':',color='black', linewidth=1, label='Ey @ 0.25 Hz')
plt.plot(depth, Ez_phase_freq_1,'--',color='black', linewidth=1, label='Ez @ 0.25 Hz')
plt.plot(depth, Ex_phase_freq_2,'-',color='blue', linewidth=1, label='Ex @ 1 Hz')
plt.plot(depth, Ey_phase_freq_2,':',color='blue', linewidth=1, label='Ey @ 1 Hz')
plt.plot(depth, Ez_phase_freq_2,'--',color='blue', linewidth=1, label='Ez @ 1 Hz')
plt.legend(loc='best', fontsize=7)
plt.ylabel('Phase (rad)',fontsize=10)
plt.xticks(fontsize=9)
plt.yticks(fontsize=9)
#plt.xticks(np.arange(0,11,step=2),fontsize=9)
#plt.yticks(np.arange(-3,4,step=2),fontsize=9)

#plot the phase of B at frequency of 1 Hz (see input file)
ax4 = plt.subplot(2, 2, 4)
plt.plot(depth, Bx_phase_freq_1,'-',color='black', linewidth=1, label='Bx @ 0.25 Hz')
plt.plot(depth, By_phase_freq_1,':',color='black', linewidth=1, label='By @ 0.25 Hz')
plt.plot(depth, Bz_phase_freq_1,'--',color='black', linewidth=1, label='Bz @ 0.25 Hz')
plt.plot(depth, Bx_phase_freq_2,'-',color='blue', linewidth=1, label='Bx @ 1 Hz')
plt.plot(depth, By_phase_freq_2,':',color='blue', linewidth=1, label='By @ 1 Hz')
plt.plot(depth, Bz_phase_freq_2,'--',color='blue', linewidth=1, label='Bz @ 1 Hz')
plt.legend(loc='best', fontsize=7)
plt.ylabel('Phase (rad)',fontsize=10)
plt.xticks(fontsize=9)
plt.yticks(fontsize=9)
#plt.xticks(np.arange(0,11,step=2),fontsize=9)
#plt.yticks(np.arange(-0.5,0.4,0.2),fontsize=9)
ax4.set_xlabel('Depth (m)',fontsize=10)

plt.show()

# Save files in pdf format
fig1.savefig("Arbitrary_oriented_MD.pdf",bbox_inches = 'tight', dpi=1000)