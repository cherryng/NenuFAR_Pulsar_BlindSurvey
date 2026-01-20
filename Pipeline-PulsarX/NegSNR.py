import numpy as np
# import psrchive as psr
import pylab as plt
import sys
import argparse
import filterbank as fil


'''def LoadArch ( arName , compar=True ) :

	ar = psr.Archive_load( arName )

	ar.pscrunch()
	ar.tscrunch()
	if not ar.get_dedispersed() :
		ar.dedisperse()
	ar.fscrunch()

	data = ar.get_data()[0,0,0,:]

	if compar == True :
		ar.remove_baseline()
		data0 = ar.get_data()[0,0,0,:]

	if compar == False :
		return data
	else :
		return data , data0'''


def LoadFbk(fbkName):

    fbk = fil.FilterbankFile(fbkName)
    spc = fbk.get_spectra()

    return spc.mean(0)


def EstimSNR(
        data,
        data0=None,
        verb=False,
        quiet=False,
        onlyprof=False,
        rfiErasing=False,
        nsig=3.0):

    flagMask = False
    if np.ma.is_masked(data):

        if data.count() == 0:
            print("All data are masked !")

            if onlyprof:
                if rfiErasing:
                    return data.filled(0.)
                else:
                    return data

            else:
                snrInfo = [0, 0, 0]
                offInfo = [100, 0, 0, 0, 0]
                return snrInfo, offInfo, np.ones(data.shape).astype(bool)

        else:
            flagMask = True
            dataMasked = data.copy()
            data = data.data[~data.mask]

    medGood = np.median(data)
    avgGood = data.mean()

    centGood = data.min() + data.ptp() / 2.
    p = 1

    if (centGood < medGood) and verb:
        print("\nNeed to remove points extremely low !")
        print("\nPerc , Cent_Good , Med_tot , Cent_Good - Med_tot")

    percGood = data.min()
    while centGood < medGood:

        percGood = np.percentile(data, p)
        maskGood = data >= percGood
        dataGood = data[maskGood]
        centGood = dataGood.min() + dataGood.ptp() / 2.
        p += 1

        if verb:
            print(p - 1, centGood, medGood, centGood - medGood)

    if p == 1:
        dataGood = data

    medOff = np.median(dataGood)
    avgOff = dataGood.mean()

    if verb:
        print("\nPerc , Avg_Off , Med_Off , Avg_Off - Med_Off")

    p = 99
    if avgOff - medOff < 0.:
        avgOff = medOff + 1

    percOff = dataGood.max()
    while avgOff - medOff > 0.:

        percOff = np.percentile(dataGood, p)
        maskOff = dataGood <= percOff
        dataOff = dataGood[maskOff]
        medOff = np.median(dataOff)
        avgOff = dataOff.mean()

        if verb:
            print(p, avgOff, medOff, avgOff - medOff)

        if p == 0:
            print('Do not converge !\n')
            break
        else:
            p -= 1

    if p < 99:
        stdOff = dataOff.std()
        del dataOff
    else:
        stdOff = dataGood.std()
    del dataGood

    if stdOff == 0.:
        stdOff = 1.

    if flagMask:
        prof = (dataMasked - medOff) / stdOff
    else:
        prof = (data - medOff) / stdOff

    if not onlyprof:

        maskOff = data <= percOff
        pulse = prof[~maskOff]
        if p < 99:
            snrInfo = [prof.max(),
                       prof.sum() * np.sqrt(pulse.max() / pulse.sum()),
                       np.sum(prof**2) / (len(prof) - 1)]
        else:
            snrInfo = [0., 0., np.sum(prof**2) / (len(prof) - 1)]
        offInfo = [p + 1, avgOff, medOff, stdOff, avgOff - medOff]

    if not quiet:
        if data0 is None:
            print("\nPeak SNR : {:.3f}".format(snrInfo[0]))
            print("\nSNR = {:.3f}".format(snrInfo[1]))
            print("Chi2 = {:.3f}\n".format(snrInfo[2]))
        else:
            print(
                "\nPeak SNR : \t{:.3f}\t{:.3f}".format(
                    max(data0), snrInfo[0]))
            print(
                "\nSNR = \t{:.3f}\t{:.3f}".format(
                    data0.sum() *
                    np.sqrt(
                        max(pulse) /
                        pulse.sum()),
                    snrInfo[1]))
            print("Chi2 = \t{:.3f}\t{:.3f}\n".format(
                np.sum(data0**2) / (len(prof) - 1), snrInfo[2]))

    if rfiErasing:
        prof = RfiErasing(
            prof,
            pgood=(
                percGood -
                medOff) /
            stdOff,
            poff=(
                percOff -
                medOff) /
            stdOff)
        prof = RfiErasing(prof, pgood=-nsig, poff=nsig)

    if onlyprof:
        return prof
    else:
        return snrInfo, offInfo, maskOff


def PlotProfile(offInfo, maskOff, data, data0=None):

    bins = np.arange(len(data))
    binsOff = bins[maskOff]
    dataOff = data[maskOff]
    prof = (data - offInfo[2]) / offInfo[3]

    plt.figure('SNR comparison', (20, 8))
    plt.subplot(121)
    plt.plot(bins, data, '+-')
    plt.plot(binsOff, dataOff, '+-')
    plt.hlines(np.median(data), 0, max(bins), color='g')
    plt.hlines(data.mean(), 0, max(bins), color='r')
    plt.hlines(offInfo[1], 0, max(bins), color='Purple')
    plt.hlines(
        offInfo[1] +
        offInfo[3],
        0,
        max(bins),
        color='Purple',
        linestyle='dashed')
    plt.hlines(
        offInfo[1] -
        offInfo[3],
        0,
        max(bins),
        color='Purple',
        linestyle='dashed')
    plt.grid()
    plt.subplot(122)
    plt.plot(prof)
    if data0 is not None:
        plt.plot(data0)
    plt.hlines(0, 0, max(bins), color='Purple')
    plt.hlines(-1, 0, max(bins), color='Purple')
    plt.hlines(1, 0, max(bins), color='Purple')
    plt.grid()
    plt.show()


def RfiErasing(data, pgood=None, poff=None):

    med = np.ma.median(data)
    if not pgood:
        pgood = data.min()
    if not poff:
        poff = data.max()

    if np.ma.is_masked(data):
        data.filled(med)

    data = np.ma.masked_outside(data, pgood, poff)

    return data.filled(med)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        'Script to estimates the SNR of a profile with the negative profil method.')
    parser.add_argument(
        dest='file',
        help='Archive file name to calculates the SNR.')
    parser.add_argument(
        '-c',
        '--compar',
        action='store_true',
        help='Option to compare with the SNR of Psrchive.')
    parser.add_argument(
        '-p',
        '--plot',
        action='store_true',
        help='Option to plot profiles with ON-OFF selected points.')
    parser.add_argument(
        '-v',
        '--verb',
        action='store_true',
        help='Option to print information on the OFF research.')
    parser.add_argument(
        '-q',
        '--quiet',
        action='store_true',
        help="Option to don't print SNR result.")
    parser.add_argument('-f', '--fbk', action='store_true',
                        help='Option to work on a filterbank file.')
    args = parser.parse_args()

    if args.compar:
        data, data0 = LoadArch(args.file, args.compar)
        snrInfo, offInfo, maskOff = EstimSNR(
            data, data0, verb=args.verb, quiet=args.quiet)
        if args.plot:
            PlotProfile(offInfo, maskOff, data, data0=data0)

    else:
        if args.fbk:
            data = LoadFbk(args.file)
        else:
            data = LoadArch(args.file, args.compar)
        snrInfo, offInfo, maskOff = EstimSNR(
            data, verb=args.verb, quiet=args.quiet)
        if args.plot:
            PlotProfile(offInfo, maskOff, data, data0=None)
