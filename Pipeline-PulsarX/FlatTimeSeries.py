# -----------------------------------------------------------------------------------------------------------------------

#               FLAT TIME SERIES

#   Script to flat the time series extract from a filterbank file.

#       Author :    Mark BRIONNE
#       Date :      29/08/2023
#       Version :   4b2

#   Comments :  Fit replaced by a normalisation with a running average computed by convolution.
#           Modify the division by a soustraction to remove the running average.
#           Better control of multiprocessing chunksizes to be lower than 1000.

# -----------------------------------------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------
#   IMPORTING PACKAGES

import time
import warnings
import glob
import argparse
import NegSNR as neg
import astropy.stats as st
import multiprocessing as mp
import scipy.signal as sig
import sys
import scipy.optimize as opt
import pylab as plt
import numpy as np
import filterbank as fil
import matplotlib
matplotlib.use('Agg')


# -------------------------------------------------------------------------------
#   PREPROCESSING PARAMETERS

warnings.simplefilter("ignore", opt.OptimizeWarning)

findPeak = False
debugPeak = False
debugFit = False
debugNorm = False
debugTimeRFI = False


# -------------------------------------------------------------------------------
#   FUNCTIONS

def FitChunk(tser0):

    tser0 = (tser0 - tser0.min() + 1) * 100
    conv = sig.convolve(tser0, win, mode='same')

    conv *= win.size
    conv /= normconv

    newtser = tser0 - conv

    if newtser.std() == 0.:
        return newtser - newtser.mean()
    else:
        return (newtser - newtser.mean()) / newtser.std()


def NormalisedArray2d(array2d, axis):

    if axis == 0:
        ch = iter(array2d)
        cs = array2d.shape[0] // args.ncpus + 1
        if cs >= 1000:
            cs //= 2
    else:
        ch = iter(array2d.transpose())
        cs = array2d.shape[1] // args.ncpus + 1
        if cs >= 1000:
            cs //= 2

    pool = mp.Pool(processes=args.ncpus)
    normPool = pool.imap(NegNorm, ch, chunksize=cs)
    pool.close()
    pool.join()

#   return np.transpose( np.array( list( normPool ) ) )

    if axis == 0:
        return np.array(list(normPool))
    else:
        return np.transpose(np.array(list(normPool)))


def NegNorm(data):

    # return neg.EstimSNR( data , verb=False , quiet=True , onlyprof=True ,
    # rfiErasing=True )
    if data.std() != 0.:
        return (data - np.ma.median(data)) / data.std()
    else:
        return data - np.ma.median(data)


def FillMed1d(data):

    if data.count() > 0:
        return data.filled(np.ma.median(data))
    else:
        return data.filled(np.mean(data))


def FillMed2d(array2d, axis=1):

    if axis == 0:
        ch = iter(array2d)
        cs = array2d.shape[0] // args.ncpus + 1
        if cs >= 1000:
            cs //= 2
    else:
        ch = iter(array2d.transpose())
        cs = array2d.shape[1] // args.ncpus + 1
        if cs >= 1000:
            cs //= 2

    pool = mp.Pool(processes=args.ncpus)
    fillPool = pool.imap(FillMed1d, ch, chunksize=cs)
    pool.close()
    pool.join()

    if axis == 0:
        return np.array(list(fillPool))
    else:
        return np.transpose(np.array(list(fillPool)))


def Mean1d(data):

    return data.mean()


def Mean2d(array2d, axis=2):

    if axis == 0:
        ch = iter(array2d)
        cs = array2d.shape[0] // args.ncpus + 1
        if cs >= 1000:
            cs //= 2
    else:
        ch = iter(array2d.transpose())
        cs = array2d.shape[1] // args.ncpus + 1
        if cs >= 1000:
            cs //= 2

    pool = mp.Pool(processes=args.ncpus)
    meanPool = pool.imap(Mean1d, ch, chunksize=cs)
    pool.close()
    pool.join()

    if axis == 0:
        return np.array(list(meanPool))
    elif axis == 1:
        return np.transpose(np.array(list(meanPool)))
    else:
        return Mean1d(np.transpose(np.array(list(meanPool))))


def Median1d(data):

    return np.ma.median(data)


def Median2d(array2d, axis=2):

    if axis == 0:
        ch = iter(array2d)
        cs = array2d.shape[0] / args.ncpus + 1
        if cs >= 1000:
            cs /= 2
    else:
        ch = iter(array2d.transpose())
        cs = array2d.shape[1] / args.ncpus + 1
        if cs >= 1000:
            cs /= 2

    pool = mp.Pool(processes=args.ncpus)
    meanPool = pool.imap(Median1d, ch, chunksize=cs)
    pool.close()
    pool.join()

    if axis == 0:
        return np.array(list(meanPool))
    elif axis == 1:
        return np.transpose(np.array(list(meanPool)))
    else:
        return Mean1d(np.transpose(np.array(list(meanPool))))


def DynSpec():

    spc = fbk.get_spectra(0, 1e4)
    normspc = NormalisedArray2d(spc, 1)

    plt.figure('Dyn spec', (20, 8))
    plt.subplot(211)
    plt.imshow(np.transpose(normspc), aspect='auto', origin='lower')
    plt.subplot(212)
    plt.plot(normspc.mean(1), '+')
    plt.show()


def FindChunk(tser):

    t_chunk = 200
    t_chunk = int(t_chunk)

    tser = np.reshape(tser[: len(tser) // t_chunk * t_chunk],
                      (len(tser) // t_chunk, t_chunk))
    avg = np.median(tser, 1)

    grad = (np.gradient(avg[:-1]) + np.gradient(avg[1:])) / 2.
    grad /= -np.std(grad)

    peaks = sig.argrelmax(grad, order=100)
    maj_peaks = grad[peaks] > np.std(grad)
    peaks_grad = np.argsort(grad[peaks][maj_peaks])
    t = np.arange(len(grad))

    if debugPeak:
        ax1.plot(np.arange(len(avg)), avg)
        ax1.vlines(t[peaks][maj_peaks][peaks_grad], min(avg), max(avg), 'r')

    if findPeak:
        ax4.plot(t * t_chunk + t_chunk // 2, grad)
        ax4.vlines(
            t[peaks][maj_peaks][peaks_grad] *
            t_chunk +
            t_chunk //
            2,
            min(grad),
            max(grad),
            'r')

    ind_jump = []
    for indc in t[peaks][maj_peaks][peaks_grad]:
        ind_jump.append(FindJump(tser, indc) + 1)

    return ind_jump


def FindJump(tser, ind_chunk):

    t_chunk = 20
    t_chunk = int(t_chunk)
    jump = tser[ind_chunk - 1][-t_chunk // 2:]
    jump = np.append(jump, tser[ind_chunk])
    jump = np.append(jump, tser[ind_chunk + 1])
    jump = np.append(jump, tser[ind_chunk + 2][: t_chunk // 2])

    avg = []
    i = t_chunk // 2
    while i < len(jump) - t_chunk / 2:
        avg.append(np.median(jump[i - t_chunk // 2: i + t_chunk // 2 + 1]))
        i += 1

    grad = (np.gradient(avg[:-1]) + np.gradient(avg[1:])) / 2.
    grad /= -np.std(grad)

    if findPeak or debugPeak:
        ax3.plot(np.arange(len(jump)), jump)
        ax3.vlines(np.argmax(grad) + t_chunk / 2, min(jump), max(jump), 'r')

    if debugPeak:
        ax2.plot(np.arange(len(grad)) + t_chunk / 2, grad)
        ax2.vlines(np.argmax(grad) + t_chunk / 2, min(grad), max(grad), 'r')

    return ind_chunk * tser.shape[1] + np.argmax(grad)


def Conv8b(data):

    data8b = np.round(data * scalePtp + deltaMin)
    data8b = data8b.clip(0, 255)

    return data8b.astype(int)


def Dowsampling(ser, ds):

    nbchunk = ser.size / ds
    newser = np.zeros((nbchunk + 1, ))

    meanser = ser[: nbchunk * ds]
    newser[: nbchunk] = np.mean(meanser.reshape(nbchunk, ds), 1)
    newser[-1] = ser[nbchunk * ds:].mean()

    return newser


# -------------------------------------------------------------------------------
#   MAIN PROGRAM

parser = argparse.ArgumentParser(
    'Script to flat the filterbank NenuFAR Blind Survey data.')
parser.add_argument(dest='fbk', type=str, help='Filterbank file to use.')
parser.add_argument('-c', '--ncpus', type=int, default=50,
                    help='Number of CPUs to use (default = 50).')
parser.add_argument(
    '-t',
    '--thr',
    type=float,
    default=3.,
    help='Threshold to consider as an RFI (default = 3 sigmas).')
parser.add_argument(
    '-o',
    '--out',
    type=str,
    help='Output name for the flattened created filterbank file and the results PNG file (default is the fbk filename).')
args = parser.parse_args()


# -------------------------------------------------------------------------------
print("EXTRACTING DATA")
print((time.ctime()))

print("\nExtracting data\n")

fbk = fil.FilterbankFile(args.fbk)
spc0 = fbk.get_spectra()
fbk.close()

if not args.out:
    args.out = args.fbk.split('.')[0]

# spc0 = spc0[:,:10000]  #[32000:38000,:5000]    #[68475:102780,:]
# fbk.times = fbk.times[32000:38000] #[68475:102780]


# -------------------------------------------------------------------------------
# ANALOGICAL JUMPS RESEARCH
print("ANALOGICAL JUMPS RESEARCH")
print((time.ctime()))

print("Jump research")

if findPeak:
    plt.figure('Jump finding', (20, 8))
    ax4 = plt.subplot(212)
    ax3 = plt.subplot(211, sharex=ax4)

INDJ = FindChunk(spc0.sum(1))
INDJ.sort()

print("Jumps found at :")
for i in INDJ:
    print(("\t{:d} s".format(int(i * fbk.dt))))

if findPeak:
    ax3.plot(np.arange(len(spc0.sum(1))), spc0.sum(1), '+')
    ax3.vlines(INDJ, min(spc0.sum(1)), max(spc0.sum(1)), 'r')
    plt.show()


# -------------------------------------------------------------------------------
# MASKING SATURATED VALUES

print("Mask bad values\n")

maskSatur = (spc0 == 0) + (spc0 == 255)
spc = np.ma.masked_where(maskSatur, spc0)
del spc0


# -------------------------------------------------------------------------------
# NORMALISATION ON FREQUENCY
print("NORMALISATION ON FREQUENCY")
print((time.ctime()))

print("Normalisation\n")
# print NegNorm( spc[:,0] )
spc = NormalisedArray2d(spc, 1)
print("Array normalized\n")


# -------------------------------------------------------------------------------
# FREQUENCY RFI RESEARCH
print("FREQUENCY RFI RESEARCH")
print((time.ctime()))

print(("Bandpass RFI erasing to {:.1f} sigmas".format(args.thr)))
print("Computing bandpass")

bandp = Mean2d(spc, 1)

print("Continuum research")

snrBp, offBp, maskBandpass = neg.EstimSNR(
    bandp, verb=False, quiet=True, onlyprof=False, rfiErasing=False)
mask5sig = np.abs((bandp - offBp[2])) > args.thr * offBp[3]
maskBandpass = np.resize(mask5sig, spc.shape)

print("RFI erasing on each channel")

spc = np.ma.masked_where(maskSatur + maskBandpass, spc)
del maskSatur
del maskBandpass
del mask5sig
del bandp
spc = FillMed2d(spc, 1)

print("Bandpass RFI erased\n")


# -------------------------------------------------------------------------------
# TIME RFI RESEARCH
print("TIME RFI RESEARCH")
print((time.ctime()))

print(("Time series RFI erasing to {:.1f} sigmas".format(args.thr)))
print("Computing mean time series")

timeser = Mean2d(spc, 0)

print("Continuum research")

snrTs, offTs, maskTimeser = neg.EstimSNR(
    timeser, verb=False, quiet=True, onlyprof=False, rfiErasing=False)
mask5sig = (timeser - offTs[2]) > args.thr * offTs[3]
maskTimeser = np.repeat(mask5sig, spc.shape[1]).reshape(spc.shape)

print("RFI erasing on each subintegration")

spc = np.ma.masked_where(maskTimeser, spc)
del timeser
del mask5sig
del maskTimeser
spc = FillMed2d(spc, 0)

if debugTimeRFI:
    plt.figure('Time RFI erasing', (20, 8))
    plt.subplot(221)
    plt.hist(timeser, density=True, bins=100)
    plt.vlines(offTs[2] + args.thr * offTs[3], 0, 1, color='Red')
    plt.subplot(222)
    plt.hist(spc1.mean(1), density=True, bins=100)
    plt.vlines(offTs[2] + args.thr * offTs[3], 0, 1, color='Red')
    plt.subplot(223)
    plt.plot(fbk.times, timeser, '+')
    plt.plot(fbk.times[maskTimeser[:, 0]], timeser[maskTimeser[:, 0]], '+')
    plt.subplot(224)
    plt.plot(fbk.times, spc1.mean(1), '+')
    plt.savefig('Diag_rfi_time.png')
    plt.show()
    sys.exit()

print("Time series RFI erased\n")

# np.save( 'Test66_reduced_data' , spc3 )
# sys.exit()


# -------------------------------------------------------------------------------
# NORMALISATION ON TIME
print("NORMALISATION ON TIME")
print((time.ctime()))

print("\nNormalisation of the time series channel by channel")
print("Runing average using a gaussian window of size 15 samples and with a std of 2 samples.")

fitSer = np.zeros((0, spc.shape[1]))
nbJump = 0

win = sig.windows.gaussian(176, 17)

for splitSpc in np.split(spc, INDJ, 0):

    nbJump += 1
    print(("Processing of the analogical jump {:d}".format(nbJump)))

    normconv = sig.convolve(
        np.ones(
            splitSpc.shape[0]) *
        win.size,
        win,
        mode='same')

    pool = mp.Pool(processes=args.ncpus)
    chanser = iter(splitSpc.transpose())
    fitSerChan = pool.imap(
        FitChunk,
        chanser,
        chunksize=spc.shape[1] //
        args.ncpus +
        1)
    print("Waiting result\n")

    pool.close()
    pool.join()

    fitSerSplit = np.array(list(fitSerChan))
    fitSerSplit = fitSerSplit.transpose()

    fitSer = np.vstack((fitSer, fitSerSplit))

del spc
del splitSpc
del pool
del chanser
del fitSerChan
del fitSerSplit

if debugFit:
    print(("Chi-2 = {:.2f}".format(np.sum(np.abs(Mean2d(fitSer, 0))))))
    plt.figure('Time series fit', (20, 8))
    plt.suptitle("Chi-2 = {:.2f}".format(np.sum(np.abs(Mean2d(fitSer, 0)))))
    plt.subplot(311)
    plt.plot(fbk.times, Mean2d(spc3, 0), '+')
    plt.plot(fbk.times, Mean2d(spc3 + fitSer, 0), '+')
    plt.subplot(312)
    plt.plot(fbk.times, Mean2d(fitSer, 0), '+')
    plt.subplot(313)
    plt.plot(fbk.times, abs(np.gradient(Mean2d(fitSer, 0))))
    plt.savefig('Diag_fit_log.png')
    plt.show()
    sys.exit()


# -------------------------------------------------------------------------------
# CONVERSION ON 8 BITS
print("CONVERSION ON 8 BITS")
print((time.ctime()))

print("Conversion to 8 bits")

print("Scale and offset calculating")

scalePtp = 256 / fitSer.ptp()
deltaMin = fitSer.min()

print("Rescaling and correct saturations")

fitSer8b = (fitSer - deltaMin) * scalePtp
fitSer8b = np.round(fitSer8b)
fitSer8b = fitSer8b.clip(0, 255)
fitSer8b = fitSer8b.astype(int)

# -------------------------------------------------------------------------------
# WRITING THE NEW FILTERBANK FILE
print("WRITING THE NEW FILTERBANK FILE")
print((time.ctime()))

print("Reading original fbk header\n")

head_data, head_size = fil.read_header(args.fbk)

print("Writing new flattened fbk\n")

fil.create_filterbank_file(
    '{:s}.flat.fbk'.format(
        args.out),
    header=head_data,
    spectra=fitSer8b,
    verbose=True)

sys.exit()

# -------------------------------------------------------------------------------
# CREATING THE RESULTS PLOT

print("Observation data extraction\n")

_obsFile = open(glob.glob('*.readfile')[0], 'r')
obsData = _obsFile.readlines()
_obsFile.close()

col1 = obsData[2:10]
col2 = obsData[14:21]
col2.append(obsData[26])
cols = np.column_stack((col1, col2))

print("Plotting results\n")

plt.figure('Flattening', (18, 20))
plt.suptitle(
    "Observation file : {:s}".format(
        obsData[1].split("'")[1]),
    fontsize=16,
    fontweight='bold')

ax0 = plt.subplot2grid((4, 2), (0, 0), colspan=2)
tablePlot = ax0.table(cellText=cols, loc=8, cellLoc='right', edges='open')
tablePlot.set_fontsize(16)
tablePlot.scale(1, 2)
ax0.set_axis_off()

ax1 = plt.subplot(423)
ax1.set_title('Raw data')
ax1.imshow(np.transpose(spc0), aspect='auto', origin='lower', extent=[
           fbk.times[0], fbk.times[-1], fbk.frequencies[-1], fbk.frequencies[0]])
ax1.set_xlabel('Time (s)')
ax1.set_ylabel('Frequencies (MHz)')

ax2 = plt.subplot(424, sharex=ax1, sharey=ax1)
ax2.set_title('Flat data')
ax2.imshow(np.transpose(fitSer8b), aspect='auto', origin='lower', extent=[
           fbk.times[0], fbk.times[-1], fbk.frequencies[-1], fbk.frequencies[0]])
ax2.set_xlabel('Time (s)')
ax2.set_ylabel('Frequencies (MHz)')

ax3 = plt.subplot(425)
ax3.plot(fbk.times, Mean2d(spc0, 0), '+')
ax3.set_xlabel('Time (s)')
ax3.set_ylabel('Amplitude (A.U.)')

ax4 = plt.subplot(426, sharex=ax3)
ax4.plot(fbk.times, Mean2d(fitSer8b, 0), '+')
ax4.set_xlabel('Time (s)')
ax4.set_ylabel('Amplitude (A.U.)')

ax5 = plt.subplot(427)
ax5.plot(fbk.frequencies, Mean2d(spc0, 1), '+')
ax5.set_xlabel('Frequencies (MHz)')
ax5.set_ylabel('Amplitude (A.U.)')

ax6 = plt.subplot(428, sharex=ax5)
ax6.plot(fbk.frequencies, Mean2d(fitSer8b, 1), '+')
ax6.set_xlabel('Frequencies (MHz)')
ax6.set_ylabel('Amplitude (A.U.)')

plt.savefig('Flat_{:s}.png'.format(args.out))
plt.show()
