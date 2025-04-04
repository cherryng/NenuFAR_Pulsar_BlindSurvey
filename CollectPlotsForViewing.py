import numpy as np
import argparse
import sys
import pylab as plt
import os

parser = argparse.ArgumentParser( description='Script to read and to sort the significances of the candidates.' )

parser.add_argument( '-t' , dest='thr' , type=float , nargs=2 , default=[8.,100.] , help='Events with significance b/w the two given values.' )
parser.add_argument( '-p' , dest='periods' , type=float , nargs=2 , default=[0.,300000.] , help='Events with a period b/w the two given values in ms.' )
parser.add_argument( '-d' , dest='dms', type=float , nargs=2 , default=[0.,300.] , help='Events with DM b/w the two given values.' )
parser.add_argument( '-month', dest='date', type=str, nargs=1, help='Specific month to view in YYYYMM') 
parser.add_argument( '-m' , dest='harm' , type=int , default=20 , help='Search for harmonics of known pulsars up to the given value.' )

args = parser.parse_args()

# --------------------------------------------------------------------------------------------------------------------------------------
# EXTRACTING DATA FROM TABLES

psrTable = np.loadtxt( '/home/mbrionne/Tables_summary/ATNF_per_dm.csv' , usecols=(2,3) , dtype=float , delimiter=';' )
psrNames = np.loadtxt( '/home/mbrionne/Tables_summary/ATNF_per_dm.csv' , usecols=(1,) , dtype=str , delimiter=';' )
psrNames = psrNames[ psrTable[:,1] != 0.00 ]
psrTable = psrTable[ psrTable[:,1] != 0.00 ] #[Period in sec, DM]

#Read in the database
candData = np.loadtxt( '/home/mbrionne/Tables_summary/CanTableBlindSurvey.csv' , dtype=str , delimiter=',' , skiprows=1 )
candFil = candData[:,3] #Filename
candDms = candData[:,6].astype( float ) #DM
candPer = candData[:,7].astype( float ) #Period
raw_sig = candData[:,9:].astype( float ) #Sigmas[1:8]
split_data = np.char.split(candFil, '_')
obs_date = [item[3][0:6] for item in split_data] #YYYYMM

del candData

# COMPUTING GLOBAL SIGNIFICANCE
# Optimized coefficients to calculate the global significance
coeffs = [ 0.30728 , 0.18576 , 0.09124 , 0.09124 , 0.09124 , 0.22177 , 0.00573 , 0.00573 ]
sig = np.dot( raw_sig , coeffs )

#----------------------
print("Length of data before any masking:     ", len(candFil))
print("Filtering based on input requested parameters:")
print("    Observation date range (YYYYMM): {}".format(str(args.date[0])))
print("    Threshold range: {:.2f} < Sig < {:.2f}".format( args.thr[0] , args.thr[1] ))
print("    Period range: {:.2f} < P0 < {:.2f}".format( args.periods[0] , args.periods[1] ))
print("    DMs range: {:.2f} < DM < {:.2f}".format( args.dms[0] , args.dms[1] ))

#Masking events that are not the requested month
obs_date = np.array(obs_date)
maskDate = obs_date == args.date[0]
sig = sig[maskDate]
candFil = candFil[ maskDate ]
candDms = candDms[ maskDate ]
candPer = candPer[ maskDate ]
print("Length of data after masking for YYYYMM:", len(sig))

# ------------------------------------------------------------------------------------------------
#Sort arrays by sigma
ord = np.argsort( sig )
ord = np.flip( ord , 0 )

sig = sig[ ord ]
candFil = candFil[ ord ]
candDms = candDms[ ord ]
candPer = candPer[ ord ]

# ------------------------------------------------------------------------------------------

#Masking events that are out of the requested sig range
maskClass = np.ones( sig.shape )
maskClass = maskClass.astype('bool')
if ( args.thr[0] > 0. ) or ( args.thr[1] < 100. ) :
    maskClass = ( sig > args.thr[0] ) * ( sig < args.thr[1] )
sig = sig[ maskClass ]
candFil = candFil[ maskClass ]
candDms = candDms[ maskClass ]
candPer = candPer[ maskClass ]
print("Length of data after masking for sig:    ", len(sig))

#Masking events that are out of the requested period range
maskPer = ( candPer > args.periods[0] ) * ( candPer < args.periods[1] )
sig = sig[ maskPer ]
candFil = candFil[ maskPer ]
candDms = candDms[ maskPer ]
candPer = candPer[ maskPer ]
print("Length of data after masking for period: ", len(sig))

#Masking events out of DM range
maskDms = ( candDms > args.dms[0] ) * ( candDms < args.dms[1] )
sig = sig[ maskDms ]
candFil = candFil[ maskDms ]
candDms = candDms[ maskDms ]
candPer = candPer[ maskDms ]
print("Length of data after masking for DM:     ", len(sig))
print("==============================")
ct_psr = 0 #Counter for known psr
ct_cand = 0 #Counter for candidates
detPSR = list() #List of known psr names

listCAND = list()
listPSR = list()
for cf , cp , cd   in zip( candFil , candPer , candDms ) :
    print("Working on ", cf.split("_A")[0], cp, "DM", cd)
    nbPSR = 0
    maskDM = abs( cd - psrTable[:,1] ) < 0.85 #Find pulsars with <0.1 diff DM units
    #test = psrTable[maskDM]
    #print('test', test)
    if maskDM.sum() > 0. :
	perPSR = psrTable[:,0][ maskDM ]
        #print('perPSR', perPSR)
	potPSR = psrNames[ maskDM ]
	maskPer = abs( cp - perPSR * 1e3 ) < 1. #Find pulsars with period diff < 1e3 
	if args.harm > 1 :
	    for harm in range(1,args.harm+1) :
		maskPer += abs( cp / harm - perPSR * 1e3 ) < 1.
		maskPer += abs( cp * harm - perPSR * 1e3 ) < 1.
	    nbPSR = maskPer.sum()

    if nbPSR > 0. :  #Consider this event a real pulsar 
	knownPSR = potPSR[ maskPer ]
	indClosest = np.argmin( abs( cp - perPSR[ maskPer ] * 1e3 ) )
	knownPSR = knownPSR[ indClosest ]
	detPSR.append( knownPSR )
	ct_psr += 1
	listPSR.append( cf )
        index = np.where(psrNames == knownPSR)
        print("Associating with", knownPSR, 'period=',psrTable[index][0,0], 'DM=',psrTable[index][0,1])

    else: #No known pulsar found
	ct_cand += 1
	listCAND.append( cf )
print("==============================")
print("Number of candidates:", ct_cand)
print("Number of removed events due to potential pulsar associations", ct_psr,"of which", len(np.unique(detPSR)),"are unique")


if ct_psr>0:
    print(detPSR)
    #psrCount = np.unique( detPSR , return_counts=True )
    #print("Number of unique pulsars:", len(psrCount))
#    print("         PSR name      |  Period (sec)| DM ")    
#    for i in detPSR:
#        index = np.where(psrNames == i)
#        print("   ",i,    "|",psrTable[index][0,0], "|", psrTable[index][0,1])
            
# ----------------------------------------------------------------------------------------------------
# Collect plots
import glob
import tarfile

YYYY = str(args.date)[2:6]
MM = str(args.date)[6:8]
path1 = glob.glob('/data?/BS_'+YYYY+'_'+MM+'_PROC')[0]



def make_tarfile( inputlist, label):
    PDFs=[]
    for i in inputlist:
        path2 = i.split('_DM')[0]
        PDFs.append(path1+'/'+path2+'/Candidate_files_'+path2+'/'+i+'.pfd.pdf')

    tarball_name = "Plots_"+label+"_"+YYYY+MM+".tar.gz"
    with tarfile.open(tarball_name, "w:gz") as tar:
        for file in PDFs:
            tar.add(file, arcname=os.path.basename(file))  # Add each file
    print("Tarball ",tarball_name," created successfully with",len(inputlist),"entries")


if len(listCAND)>0:
    make_tarfile(listCAND, 'CAND')
if len(listPSR)>0:
    make_tarfile(listPSR, 'PSR')
