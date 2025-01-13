import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm

loaded_data = np.load('data_compressed.npz', allow_pickle=True)
array1 = loaded_data['array1']
array1 = array1.item()

starttime = 0  #In seconds
endtime = 100 #In seconds
freqmin = 39.746856689453125 #In MHz
freqmax = 75.68435668945312 #In MHz
nbinlim = 298

#downsampled_data = data[::2, ::2]

plt.imshow(array1[..., :nbinlim], aspect='auto',
           cmap='gist_yarg',
           interpolation='nearest', origin='upper',
           extent=(starttime, endtime,
                   freqmin, freqmax))

plt.show()

            
