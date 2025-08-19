from presto import rfifind
import sys
import numpy as np

INF = sys.argv[1] #A .mask file
print('going to read', INF)
a = rfifind.rfifind(INF)

zchan = a.mask_zap_chans
zchan_array = np.array(list(zchan))
np.savetxt('mask_zap_chans.txt', zchan_array, fmt='%d')

zint = a.mask_zap_ints
zint_array = np.array(list(zint))
np.savetxt('mask_zap_ints.txt', zint_array, fmt='%d')


data = a.mask_zap_chans_per_int
with open("mask_zap_chans_per_int.txt", "w") as f:
    for item in data:
        f.write(f"{item}\n")  # Convert each element to a string and write it

print('Done writing zap chan file')
exit
