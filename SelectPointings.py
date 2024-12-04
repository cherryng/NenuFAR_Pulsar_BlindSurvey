#-------------
#
# Script to select the 4 closest pointings for a given time range on given dates.
# Based on original code by Mark Brionne
#
#---------------

import numpy as np
import pylab as plt
import sys
from astropy.time import Time
import astropy.coordinates as co
import astropy.units as u
import argparse
from astropy.utils import iers
import datetime
import matplotlib.colors as mcolors
import matplotlib.animation as animation
from meridian import Meridian

LOC_NANCAY = co.EarthLocation( lat=47.376511 , lon=2.1924002 , height=153.834 )         # Nancay coordinates
OBS_DUR = 30                                                                            # Observation duration in minutes
MAX_SEP = 5.                                                                            # Maximum separation angle between the different pointings
MERID_SEP = 2.                                                                          # Maximum separation angle with the meridian
MIN_ELEV = 30.                                                                          # Minimum elevation to consider as a transit

def FindPointing ( flagsGrid , date , slot ) :
    print("this is FindPointing")
    HourStart = int(slot[0])
    HourEnd = int(slot[1])
    print('date is ', date)
    start = datetime.datetime( int( date[0] ) , int( date[1] ) , int( date[2] ) , HourStart , 0 ) #Assuming always start top of hour
    end = datetime.datetime( int( date[0] ) , int( date[1] ) , int( date[2] ) , HourEnd, 0)

    if HourEnd < HourStart: #If ending next day
        end = datetime.datetime( int( date[0] ) , int( date[1] ) , int( date[2] ) + 1 , HourEnd , 0 )
    dur = datetime.timedelta( minutes=OBS_DUR )
    t = start + datetime.timedelta( minutes=OBS_DUR/2. )

    timeList = list()
    obsList = list()
    obsAr1d = np.zeros( (0,) , dtype=int )
    while t < end :
        print('planning for time=',t,'--------')
        #Remove already obs ($6 == 1) and not selected phase ($7)
        maskNotObs = ( flagsGrid == False ) * ( passage == args.numpass ) #An array of boolean
        if maskNotObs.sum() == 0 :
            break 
        gridCoord = co.SkyCoord( radec[ maskNotObs ] , obstime = Time( t ) , location=LOC_NANCAY , unit=( u.hourangle , u.deg ) , frame='icrs' )
        gridAltaz = gridCoord.altaz
        ptNotObs = numpt[ maskNotObs ] #Array of line numbers

        #Remove circumpolar sources, VCR bug lead to ambiguity in meridian transit position
        maskCircum = ( gridCoord.dec.deg > ( 90. - LOC_NANCAY.lat.deg ) ) * ( gridAltaz.alt.deg < LOC_NANCAY.lat.deg ) #Array of boo
        gridCoord = gridCoord[ ~maskCircum ]
        gridAltaz = gridAltaz[ ~maskCircum ]
        ptNotObs = ptNotObs[ ~maskCircum ]

        #Remove the lines below visibility horizon
        maskAlt = ( gridCoord.dec.deg < ( 90. - LOC_NANCAY.lat.deg ) ) * ( gridAltaz.alt.deg < MIN_ELEV )
        gridCoord = gridCoord[ ~maskAlt ]
        gridAltaz = gridAltaz[ ~maskAlt ]
        ptNotObs = ptNotObs[ ~maskAlt ]
        print('len ptNotObs', len(ptNotObs))

        #Sources are only obsered at meridian transit (Azimuth 0 and 180)
        meridSouth = co.SkyCoord( np.zeros( ( len( gridAltaz ), ) ) , gridAltaz.alt.deg , obstime = Time( t ) , location=LOC_NANCAY , unit=u.deg , frame='altaz' )
        meridNorth = co.SkyCoord( 180 * np.ones( ( len( gridAltaz ), ) ) , gridAltaz.alt.deg , obstime = Time( t ) , location=LOC_NANCAY , unit=u.deg , frame='altaz')
        daz = np.fmin( meridSouth.separation( gridAltaz ) , meridNorth.separation( gridAltaz )) #The smaller separation
        maskTransit = daz.value < MERID_SEP
        gridCoord = gridCoord[ maskTransit ]
        ptNotObs = ptNotObs[ maskTransit ]
        print('len ptNotObs after masking Transit', len(ptNotObs), maskTransit[0:5], daz[0:5], 'gridCoord', len(gridCoord))

        if len( gridCoord ) == 0 :
            t += dur #Increment time
            print("Incremented time to t=",t)            
            continue

        #The one pointing with the highest declination
        maxpt = gridCoord[ np.argmax( gridCoord.dec.deg ) ]
        sortsep = np.argsort( maxpt.separation( gridCoord ) )
        ptlist = sortsep[:4]  #the 4 pts with smallest separation
        sep = maxpt.separation( gridCoord )[ ptlist ]
        ptlist = ptlist.astype( int )[ sep.value < MAX_SEP ]

        #Function to find holes to observe if there are not at least 3 possible pointings
        lenptl = len( ptlist )
        while len( ptlist ) < 3 :
            print("---in the ptlist<3 loop")
            ptlist = sortsep[lenptl:lenptl+4]
            try:
                maxpt = gridCoord[ ptlist[0] ]
            except IndexError :
                ptlist = FindHole( sortsep , gridCoord )
                break
            sep = maxpt.separation( gridCoord )[ ptlist ]
            ptlist = ptlist.astype( int )[ sep.value < MAX_SEP ]
            lenptl += len( ptlist )
            
        if ptlist.size == 0 :
            pass
        elif ptlist.size == 1 :
            if ptNotObs[ ptlist ] not in obsAr1d :
                obsAr1d = np.hstack( ( obsAr1d , ptNotObs[ ptlist ] ) )
                obsList.append( ptNotObs[ ptlist ] )
                flagsGrid[ ptNotObs[ ptlist ] ] = True
                if args.verb :
                    print(gridCoord.obstime , gridCoord[ ptlist ].to_string( 'hmsdms' , sep=':' ))
            else:
                pass
        else:
            ptlist = ptlist[ ~ np.in1d( ptNotObs[ ptlist ] , obsAr1d ) ]
            obsAr1d = np.hstack( ( obsAr1d , ptNotObs[ ptlist ] ) )
            obsList.append( ptNotObs[ ptlist ] )
            print('obsList', len(obsList), obsList)
            flagsGrid[ ptNotObs[ ptlist ] ] = True
            if args.verb :
                print(gridCoord.obstime , gridCoord[ ptlist ].to_string( 'hmsdms' , sep=':' ))

        timeList.append( t )
        t += dur
        print('second loop incrementing time to t=', t)
                
    return timeList , obsList


def FindHole ( sortpts , coordPts ) :
    lenptl = 0
    pts = []
    while len( pts ) < 2 :
        pts = sortpts[lenptl:lenptl+2]
        try :
            maxpt = coordPts[ pts[0] ]
        except IndexError :
            break
        sep = maxpt.separation( coordPts )[ pts ]
        pts = pts.astype( int )[ sep.value < MAX_SEP ]
        lenptl += len( pts )

        try :
            pts = np.array( sortpts[0].astype( int ) )
        except IndexError :
            return np.array( [] )

        return pts


def MeanPosition ( Coord ) :
    sc = co.SkyCoord( ra=Coord.ra , dec=Coord.dec , frame='icrs' , unit=(u.hourangle,u.deg) )          # Coordinates object in equatorial HMS$
    raMean = sc[0].ra.deg
    decMean = sc[0].dec.deg

    i = 1
    while i < len(Coord) :
        if abs( raMean - sc[i].ra.deg ) > 180. :
            raMean = ( i * raMean + 360. + sc[i].ra.deg ) / ( 1 + i )
        else :
            raMean = ( i * raMean + sc[i].ra.deg ) / ( 1 + i )
            
        if abs( decMean - sc[i].dec.deg ) > 90. :
            decMean = ( i * decMean + 180. + sc[i].dec.deg ) / ( 1 + i )
        else :
            decMean = ( i * decMean + sc[i].dec.deg ) / ( 1 + i )
        i += 1
    return co.SkyCoord( ra=raMean , dec=decMean , frame='icrs' , unit=u.deg )

    
def ParsetUserWriting ( numobs , dateobs , numpt , coordpt ) :
    if coordpt.size > 1 : #It should be 4 with the 4 beams
        avgpos = MeanPosition( coordpt )
    else :
        avgpos = coordpt
    merid_beam = Meridian( avgpos.to_string( 'hmsdms' , sep=':' ).split()[0] , avgpos.to_string( 'hmsdms' , sep=':' ).split()[1] , dateobs.iso.split()[0] , dateobs.iso.split()[1])
    parsetName = "BS_{:04d}_{:s}.parset_user".format( numobs , dateobs.datetime.strftime('%Y%m%dT%H%M') )
    _file0 = open( parsetName , 'w' )
    print("Writing parset_user file :\t{:s}".format( parsetName ))

    _file0.write( 'Observation.contactName=cng\n' )
    _file0.write( 'Observation.name="BS_{:04d}_{:s}"\n'.format( numobs , avgpos.to_string('hmsdms' , precision=0 , sep="").replace( ' ' , '' ) ) )
    if numpt.size == 4 :
        _file0.write( 'Observation.title="BS_{:04d}_{:04d}_{:04d}_{:04d}"\n'.format( numpt[0] , numpt[1] , numpt[2] , numpt[3] ) )
    elif numpt.size == 3 :
        _file0.write( 'Observation.title="BS_{:04d}_{:04d}_{:04d}"\n'.format( numpt[0] , numpt[1] , numpt[2] ) )
    elif numpt.size == 2 :
        _file0.write( 'Observation.title="BS_{:04d}_{:04d}"\n'.format( numpt[0] , numpt[1] ) )
    else :
        _file0.write( 'Observation.title="BS_{:04d}"\n'.format( numpt ) )

    _file0.write( 'Observation.contactEmail=cherry.ng-guiheneuf@cnrs-orleans.fr\n' )
    _file0.write( 'Observation.nrBeams={:d}\n'.format( numpt.size ) )
    _file0.write( 'Observation.topic=lt03_pulsars\n\n' )

    _file0.write( 'Anabeam[0].target="J2000_{:s}_BS_TRACKING"\n'.format( avgpos.to_string('hmsdms' , precision=0 , sep="").replace( ' ' , '' ) ) )
    _file0.write( 'Anabeam[0].trackingType=tracking\n' )
    _file0.write( 'Anabeam[0].transitDate={:s}Z\n'.format( merid_beam.isot ) )
    _file0.write( "Anabeam[0].ra='{:s}'\n".format( avgpos.to_string( 'hmsdms' , sep=':' ).split()[0] ) )
    _file0.write( "Anabeam[0].dec='{:s}'\n".format( avgpos.to_string( 'hmsdms' , sep=':' ).split()[1] ) )
    _file0.write( 'Anabeam[0].startTime={:s}Z\n'.format( dateobs.isot ) )
    _file0.write( 'Anabeam[0].duration=27m\n' )
    _file0.write( 'Anabeam[0].maList=[0,1,6,8..16,22..26,32,33,41..43,45,47,48]\n\n' )
    i = 0
    
    if coordpt.size == 1 :
        print("Only got one pointing, abort")
        exit
    else:
        for c in coordpt :
            merid_beam = Meridian( c.to_string( 'hmsdms' , sep=':' ).split()[0] , c.to_string( 'hmsdms' , sep=':' ).split()[1] , dateobs.iso.split()[0] , dateobs.iso.split()[1] )

            _file0.write( 'Beam[{:d}].target="PT_{:04d}_BS_TRACKING"\n'.format( i , numpt[i] ) )
            _file0.write( 'Beam[{:d}].trackingType=tracking\n'.format( i ) )
            _file0.write( 'Beam[{:d}].transitDate={:s}Z\n'.format( i , merid_beam.isot ) )
            _file0.write( "Beam[{:d}].ra='{:s}'\n".format( i , c.to_string( 'hmsdms' , sep=':' ).split()[0] ) )
            _file0.write( "Beam[{:d}].dec='{:s}'\n".format( i , c.to_string( 'hmsdms' , sep=':' ).split()[1] ) )
            _file0.write( 'Beam[{:d}].startTime={:s}Z\n'.format( i , dateobs.isot ) )
            _file0.write( 'Beam[{:d}].duration=27m\n'.format( i ) )
            _file0.write( 'Beam[{:d}].subbandList=[200..391]\n'.format( i ) )
            _file0.write( 'Beam[{:d}].toDo=dynamicspectrum\n'.format( i ) )
            _file0.write( 'Beam[{:d}].parameters="tf: df=1.52 dt=10.0 hamm"\n\n'.format( i ) )
            i += 1
    
    _file0.write( 'Output.hd_receivers=[undysputed]' )
    _file0.close()

def StartTimeComputing ( dateobs ) :
    #For time stamps to be written in parset files
    dt = datetime.timedelta( seconds = 830 )
    
    i = 0
    while i < len( dateobs ) :
        dateobs[i] = dateobs[i] - dt
        i += 1    
    return Time( dateobs , precision=0 )


    
def GetNumObs ( obsnum ) :
    #Fint out the obs ID of the last session 
    if len( obsnum[ obsnum != 9999 ] ) == 0 : #These are the not yet observed ones. 
        return 0
    else :
        return max( obsnum[ obsnum != 9999 ] ) + 1


if __name__ == '__main__' :
    parser = argparse.ArgumentParser( description='Script to select the 4 closest pointings for all observations for a given time range on a given date.' )
    parser.add_argument( '-g' , dest='gridPath' , help='grid file with all the pointings')
    parser.add_argument( '-f' , dest='datefile' , help='File with observation dates and time slots.' )
    parser.add_argument( '-plot' , type=str , help="Option to plot the figure for each day or observation [expected value : 'days', 'obs']." )
    parser.add_argument( '-save' , action='store_true' , help='Option to save the figure.' )
    parser.add_argument( '-write' , action='store_true' , help='Option to write the parset_user file(s).' )
    parser.add_argument( '-flag' , action='store_true' , help='Option to write a new grid file with new flags.' )
    parser.add_argument( '-numpass' , type=int , dest='numpass' , default=2 , help='Option to set the number of the pass to search (default = 2).' )
    parser.add_argument( '-v' , dest='verb' , action='store_true' , help='Verbose option to print information of each observation.' )
    args = parser.parse_args()

    dateData = np.loadtxt( args.datefile , dtype=str ,ndmin=2)

    #Load big table    
    gridData = np.loadtxt( args.gridPath , dtype=str ) 
    numpt = gridData[:,0].astype( 'int' )
    radec = gridData[:,1:3]
    sizes = gridData[:,3:5].astype( 'float' )
    flagsObs = gridData[:,5].astype( 'bool' )
    passage = gridData[:,6].astype( 'int' )
    numObs = gridData[:,7].astype( 'int' )


    if args.plot or args.write : #Write parset files
        gridCoord0 = co.SkyCoord( radec , location=LOC_NANCAY , unit=( u.hourangle , u.deg ) , frame='icrs' )
        
        
    for date in dateData : #Each date of obs being planned
        print('Planning observations for', date, 'date[0]', date[0])
        tps , pts = FindPointing( flagsObs , date[0].split('-') , date[1:] )
        
        if len( tps ) == 0 :
            print("No pointings scheduled.")
            continue
        
        if args.write : #Write parset files
            st = StartTimeComputing( tps )
            actObsNum = GetNumObs( numObs )
            for t , p in zip( st , pts ) :
                print('gridCoord0[p]', gridCoord0[ p ])
                ParsetUserWriting( actObsNum  , t , p , gridCoord0[ p ] )
                numObs[ p ] = actObsNum
                actObsNum += 1

        if args.flag :
            TODAY = datetime.datetime.today().strftime('%Y%m%d')
            np.savetxt( 'Extended_grid_'+TODAY+'.csv', np.column_stack( (numpt,radec,sizes,flagsObs.astype('int'),passage,numObs) ) , fmt='%s' )
