import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm
import matplotlib.ticker as ticker
import sys

INF = sys.argv[1]
df = int(sys.argv[2]) #Downsample in frequency 
toffset = float(sys.argv[3]) #Start time offset in seconds, e.g. if the input npz starts at a non-zero time stamp
INFcands = sys.argv[4] #The singlepulse file from presto


loaded_data = np.load(INF, allow_pickle=True)
array1 = loaded_data['array1']
array1 = array1.item()

starttime = 0  #In seconds
endtime = 20  #In seconds
duration = endtime - starttime
freqmin = 39.746856689453125 #In MHz
freqmax = 75.68435668945312 #In MHz

factor_dt = 4 #the time downsample factor used for generating the npz.
dt = (10485.75994e-6)*factor_dt 
nbinlim = np.int64(duration/dt)
print('nbinlim=', nbinlim)
data = array1[..., :nbinlim]
print('Shape',data.shape)


def downsample_freq(arr, k):
    a, b = arr.shape
    if a % k != 0:
        raise ValueError(f"Number of rows ({a}) must be divisible by downsampling factor ({k}).")
    # Reshape into shape (a/k, k, b), then take the mean along axis=1
    return arr.reshape(a // k, k, b).mean(axis=1)


def mask_chan(data):
    #Only work on data with full frequency resolution
    data_masked = np.ma.masked_array(data)
    masked_chans = np.zeros(24576,dtype=bool)
    maskfromfile = np.loadtxt("mask_zap_chans.txt",dtype=int)
    maskfromfile = 24576 - 1 - maskfromfile  # Flip index positions
    masked_chans[maskfromfile] = 1
    data_masked[masked_chans] = np.ma.masked
    return data_masked

#Mask bad frequency channels with a predetermined zap list
if data.shape[0] == 24576:
    print('Got full resolution data, can apply mask')
    data = mask_chan(data)

#Downsample in frequency
downsampled_arr = downsample_freq(data, k=df)

#Plot
fig = plt.figure(figsize=(10,6))
gs = gridspec.GridSpec(2, 1, height_ratios=[1,2])
plt.subplots_adjust(wspace=0.05, hspace=0.05)

ax = plt.subplot(gs[1])  #FT data
plt.imshow(downsampled_arr[..., :nbinlim], aspect='auto',
           cmap='gist_yarg',
           interpolation='nearest', origin='upper',
           extent=(starttime, endtime,
                   freqmin, freqmax))
ax.set_xlim(starttime, endtime)
ax.set_xlabel("Time in sec")
ax.set_ylabel("Frequency channel")

#Show the presto detections
ax2 = plt.subplot(gs[0]) #Label candidates found from presto

SPcands = np.loadtxt(INFcands,dtype=float) #in seconds
SNRs = SPcands[:, 1]
Times = SPcands[:, 2]

ax2.scatter(Times, SNRs, marker='x',c='k')
for xi, yi in zip(Times, SNRs):
    ax2.vlines(xi, 0, yi, color='black', zorder=2)  # Vertical line from 0 to y-value

#Say we know 17.259561 is good, plot the known pulsar period in blue
"""for i in range(500):
    ax2.axvline(x=17.259561+1.29227*i, c='b', ls=':', zorder=0)
    ax2.axvline(x=17.259561-1.29227*i, c='b', ls=':', zorder=0)    
"""

ax2.set_xlim(starttime, endtime)
ax2.set_ylim(0, 19)
ax2.set_xlabel("")
ax2.set_ylabel("SNR from presto")
ax2.xaxis.set_major_locator(ticker.MaxNLocator(nbins=10)) 
ax2.grid(True)


ax2.xaxis.set_label_position('top') 
ax2.xaxis.tick_top() 
plt.show()

            
