#!/usr/bin/env python3
#
# Steven J. Gibbons
# NGI, Sognsveien 72, Oslo
# 2023/05/28
#
# VIScERaL_GD.py
#
# Vectorized Iterative SeismiC Event RelAtive Location - Gradient Descent.
#
# Solves for relative locations of seismic events and/or
# the slowness vectors which describe the speed and direction
# with which the seismic waves leave the source area.
# The user only ever needs relate to the geographical coordinates
# of the seismic events even though it is the Cartesian position
# vectors that are used in the solution.
# We need to specify three different input files.
#
# (1) --abslocflagsfile (file containing the absolute locations of events)
# (2) --slovecsflagsfile (file containing the slowness vectors)
# (3) --dtfile (file containing the differential times)
#
# We take a starting set of locations and slowness vectors (possibly with
# some random perturbations from the starting values) and we repeatedly
# evaluate a cost-function and its gradients moving the events and/or
# slowness vectors until a certain number of iterations is reached for
# which the cost-function does not decrease.
#
# If all the slowness vectors are fixed, the seismic events will tend to
# converge to the same locations on each run.
#
# If all the seismic events are fixed, the slowness vectors will tend to
# converge to the same states on each run.
#
# If both seismic events and slowness vectors are allowed to move
# then there is a fundamental non-uniqueness in the solutions.
# (If the velocities are faster then the events move further apart.
# If the velocities are slow then the events will move closer together.)
# However, repeat calculations with perturbed starting conditions will
# likely demonstrate significant relations between different events and
# different slowness vectors that become apparent in a statistic sense.
# If you are running a calculation with both variable events and slowness
# vectors it is recommended a three stage process:
#
#  (1) First solve event locations for the most realistic set of slowness
#       vectors at your disposal.
#        The program statphase2slowvec
#        https://github.com/stevenjgibbons/statphase2slowvec
#         will generate an appropriate file based upon the ak135 model.
#      Each line from this file needs to have a letter "F" (fixed) or
#       "S" (solve) appended.
#
#  (2) Fix the event locations according to the output obtained in
#       step 1 and then solve for the slowness vectors
#        which best fit these event locations.
#
#  (3) Solve for both event locations and slowness vectors using
#       starting conditions obtained from step (2).
#
# There are only two other mandatory parameters:
#
#   --reflat  (a reference latitude within the source area)
#   --reflon  (a reference longitude within the source area)
#
#  (Note that reflat and reflon are usually the coordinates with which the
#   statphase2slowvec program has calculated the slowness vectors)
#
# All other parameters are optional with default values set
#
#   --delkm   (step length = learning rate for distance, def 0.1 km)
#   --delslow (step length = learning rate for slowness, def 0.1 s/km)
#
#   --delkmmax (the largest distance in km an event can be moved
#                in a single iteration, def 0.005 km)
#   --maxdistkm (the largest distance in km an event can be moved
#                 from its original location in the whole run, def 999.9 km)
#   --randomizeloc (dist within which location can be moved randomly
#                    at start, def = 0.0 km)
#
#   --delslomax (distance in s/km a slowness can be moved in a single
#                 iteration, def 0.0005 s/km )
#   --maxsdist (the largest distance in s/km a slowness vector can be moved
#            from its original specification in the whole run, def 999.9 s/km)
#   --randomizeslo (dist within which slowness vector can be moved randomly
#                    at start, def = 0.0 s/km)
#
#   --writehistory (specify to write out the locations of every variable
#                    at every iteration. In general you will not want to
#                     do this as it generates huge ASCII files which you will
#                      not want unless you want to visualize the evolution
#                       of the solution. Use sparingly! def = False )
#
#   --allowmissing (Skip any lines in dtfile containing station/phases
#                    or events not found. def = False, but you will often
#                     want to set it to True.)
#
#   --FixEvent0 (even if the first line in the list of events ends in "S"
#                 we shift the whole set of events such that this event
#                   stays in the same location, def = False )
#                    
#   --locoutfile (locations output file, def = "locations.txt")
#   --slooutfile (slowness output file, def = "slowness_vectors.txt")
#   --nrmoutfile (norms output file, def = "normValues.txt" )
#
#   --blocksize (Number of iterations at which we solve either for
#                 events or slowness vectors, assuming we are solving for
#                  both. def = 0 - i.e. we solve for both at each iteration.
#                   I am not convinced there is any point in setting > 0.)
#
#   --numrandom (At each iteration, for each unknown, we pick a random number
#                 of measurements from which we estimate the gradient.
#                  def = 20. If there are fewer observations than this for
#                   an unknown event or slowness vector it will just keep
#                    picking the observations randomly this number of times.)
#
#   --maxiter (The maximum number of iterations, def = 10000 )
#   --maxpositiveslopes (The maximum number of times that the median slope
#                         of the median residual norm can be >= zero
#                           before we end the procedure. def = 100 but you
#                             will likely want to increase this to a few
#                               hundred. Do not make this value too large
#                                 or the solution may just evolve to something
#                                   far from the starting value.)
#
# Formats of the input files:
#
# (1) --abslocflagsfile (file containing the absolute locations of events)
# ------------------------------------------------------------------------
# Time_in_UTC_format      latitude   longitude   event_code  Fixed/Solve/Ignore
# ------------------------------------------------------------------------
# e.g.
# 2007-08-15T07:59:59.936 67.93590352 25.83491289 H01  F
# 2007-08-17T11:00:00.380  67.93590352  25.83491289  H06   S
# 2007-12-01T00:00:00.000  67.94        25.835       HXX   I
# ------------------------------------------------------------------------
# (2) --slovecsflagsfile (file containing the slowness vectors)
# ------------------------------------------------------------------------
# station phase   statlat    statlon   reflat     reflon      Sx            Sy       flag
# ------------------------------------------------------------------------
#  e.g.
# ARE0  P1        69.53490   25.50580  67.93590   25.83491   -0.00888570    0.12337091 F
# ARE0  S1        69.53490   25.50580  67.93590   25.83491   -0.01594664    0.22140645 F
# KEV   P1        69.75530   27.00670  67.93590   25.83491    0.02687697    0.12073203 F
# KEV   S1        69.75530   27.00670  67.93590   25.83491    0.04823417    0.21666916 F
# SGF   P1        67.44211   26.52611  67.93590   25.83491    0.08181667   -0.15176232 F
# SGF   S1        67.44211   26.52611  67.93590   25.83491    0.13714934   -0.25439928 F
# LP34  P1        67.26574   28.12528  67.93590   25.83491    0.13871543   -0.10238026 F
# LP34  S1        67.26574   28.12528  67.93590   25.83491    0.23252879   -0.17162010 F
# LP53  P1        68.08434   27.18877  67.93590   25.83491    0.16493699    0.05021592 F
# LP53  S1        68.08434   27.18877  67.93590   25.83491    0.27648398    0.08417697 F
# LP61  P1        67.91408   23.93216  67.93590   25.83491   -0.17239066   -0.00260066 F
# LP61  S1        67.91408   23.93216  67.93590   25.83491   -0.28897858   -0.00435949 F
# ------------------------------------------------------------------------
# (3) --dtfile (file containing the differential times)
# ------------------------------------------------------------------------
# ev1  ev2   starting_time_ev1_template   CC_time_for_ev2    station phase  ignored
# ------------------------------------------------------------------------
#  e.g.
# H01  H02  2007-08-15T08:00:19.109  2007-08-15T12:00:19.299  LP34   P1  0.850
# H01  H02  2007-08-15T08:00:33.100  2007-08-15T12:00:33.286  LP34   S1  0.798
# H01  H02  2007-08-15T08:00:08.467  2007-08-15T12:00:08.709  LP53   P1  0.926
# H01  H02  2007-08-15T08:00:15.228  2007-08-15T12:00:15.476  LP53   S1  0.813
# ------------------------------------------------------------------------
#
try:
    import os
    import sys
    import argparse
    import statistics
    from scipy import stats
    import random
    import numpy as np
    import geographiclib
    from geographiclib.geodesic import Geodesic
    import math
    import obspy
    from obspy import UTCDateTime
except ImportError as ie:
    miss_mod = ie.args[0].split()[3]
    print("\nThe Python module '" + miss_mod + "' is required.")
    print("Please install it and run again.\n")
    exit(1)

#==========================================================================
def vector_dxdy( dxarr, dyarr, cfarr, wgtarr ):
    xderivs = np.asarray( dxarr )
    yderivs = np.asarray( dyarr )
    costfns = np.asarray( cfarr )
    weights = np.asarray( wgtarr )
    ntotal  = len( xderivs )
    xymags  = np.zeros_like( xderivs )
    for i in range( ntotal ):
        x = xderivs[i]
        y = yderivs[i]
        z = np.sqrt( x*x + y*y )
        xymags[i] = z

    magmedian = np.median( xymags )
    maglimit  = 2.0 * magmedian
    magoutvec = 0.5 * magmedian

    indices = np.argwhere( xymags < maglimit )

    xdersel = xderivs[indices]
    ydersel = yderivs[indices]
    costsel = costfns[indices]
    wgtssel = weights[indices]
    nval    = len( wgtssel )

    costfv  = np.mean( costsel )
    xnew    = np.dot( xdersel.reshape(1,nval), wgtssel.reshape(nval,1) ).item()
    ynew    = np.dot( ydersel.reshape(1,nval), wgtssel.reshape(nval,1) ).item()
    vmagnew = np.sqrt( xnew*xnew + ynew*ynew )
    ratio   = magoutvec/vmagnew
    #
    # Note that we only want to scale by this ratio if
    # the required size of the output vector is smaller than this!
    # (i.e. we do not want to increase the size of the output vector!)
    #
    if ( ratio < 1.0 ):
        dbydx   = xnew*ratio
        dbydy   = ynew*ratio
    else:
        dbydx   = xnew
        dbydy   = ynew

    return dbydx, dbydy, costfv

#==========================================================================
#
# We intend to increase sx by dsx and sy by dsy
# However, we do not want magnitude of (dsx,dsy) to be greater than
# the specified percentage of the magnitude of (sx,sy) so we return
# dsxscale, dsyscale which have the same ratio as dsx and dsy
#
#==========================================================================
def slowinclimit( sx, sy, dsx, dsy, percent):
    dsxscale = dsx
    dsyscale = dsy
    sxsymag  = math.sqrt( sx*sx + sy*sy )
    dsxsymag = math.sqrt( dsx*dsx + dsy*dsy )
    maglimit = 0.01 * percent * sxsymag
    if ( dsxsymag > maglimit ):
        ratio    = maglimit / dsxsymag
        dsxscale = dsx * ratio
        dsyscale = dsy * ratio
    return dsxscale, dsyscale

#==========================================================================
def sxsyazivel( sx, sy):
    s   = np.sqrt( sx*sx + sy*sy )
    slo = 0.000001
    if ( np.abs( s ) < slo ):
        s = slo
    vel = 1.0/s
    azrad = math.atan2( sx, sy )
    azdeg = math.degrees( azrad )
    if ( azdeg < 0.0 ):
        azdeg += 360.0
    return azdeg, vel

#==========================================================================
def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

#==========================================================================
def file_exist(file):
    """Check if a file (with full path) exist"""
    if not os.path.isfile(file):
        print("File: ",file," does not exist. Bailing out ...")
        exit()

#==========================================================================
#
# calcTotalNormSetOfSloVecs returns the sum of the norms of all cost
# functions pertaining to this set of slowness vectors
# 
#==========================================================================
def calcTotalNormSetOfSloVecs( SetOfSlowVecs, XYOrigs, SloVecFlags ):
    Total    = 0.0
    numterms = 0
    for sv in range(0,len(SetOfSlowVecs) ):
        navail = len( SetOfSlowVecs[sv].obsPairs )
        for i in range(0,navail):
            Theta, Psi = Theta_ijab( SetOfSlowVecs[sv].obsPairs[i],
                                      XYOrigs, SloVecFlags )
            Total     += Psi
            numterms  += 1

    return Total, numterms

def calcMedianNormSetOfSloVecs( SetOfSlowVecs, XYOrigs, SloVecFlags ):
    psivals  = []
    numterms = 0
    for sv in range(0,len(SetOfSlowVecs) ):
        navail = len( SetOfSlowVecs[sv].obsPairs )
        for i in range(0,navail):
            Theta, Psi = Theta_ijab( SetOfSlowVecs[sv].obsPairs[i],
                                      XYOrigs, SloVecFlags )
            psivals.append( Psi )
            numterms  += 1

    medval = statistics.median( psivals )
    return medval, numterms

#==========================================================================
#
# calcTotalNormSetOfLocs returns the sum of the norms of all cost
# functions pertaining to this set of locations
# 
#==========================================================================
def calcTotalNormSetOfLocs( SetOfLocs, XYOrigs, SloVecFlags ):
    Total    = 0.0
    numterms = 0
    for loc in range(0,len(SetOfLocs) ):
        navail = len( SetOfLocs[loc].obsPairs )
        for i in range(0,navail):
            Theta, Psi = Theta_ijab( SetOfLocs[loc].obsPairs[i],
                                      XYOrigs, SloVecFlags )
            Total     += Psi
            numterms  += 1

    return Total, numterms

def calcMedianNormSetOfLocs( SetOfLocs, XYOrigs, SloVecFlags ):
    psivals  = []
    numterms = 0
    for loc in range(0,len(SetOfLocs) ):
        navail = len( SetOfLocs[loc].obsPairs )
        for i in range(0,navail):
            Theta, Psi = Theta_ijab( SetOfLocs[loc].obsPairs[i],
                                      XYOrigs, SloVecFlags )
            psivals.append( Psi )
            numterms  += 1

    medval = statistics.median( psivals )
    return medval, numterms

#==========================================================================
class llLocation:
    def __init__( self, lat, lon ):
        self.lat = lat
        self.lon = lon

#==========================================================================
class xyOrigin:
    def __init__( self, absloc, refloc, origintime, evid, ndt, flag ):
        self.refloc     = refloc
        self.absloc     = absloc
        self.abslocorig = absloc
        self.evid       = evid
        self.origintime = origintime
        self.ndt        = ndt
        self.flag       = flag
        distkm          = dist_between_locs_km( refloc, absloc )
        azifromref      = receiver_to_source_backazimuth( refloc, absloc )
        azirad          = math.radians( azifromref )
        self.x          = distkm * math.sin( azirad )
        self.y          = distkm * math.cos( azirad )
        self.history    = []

    def incxyonly( self, incxkm, incykm ):
        self.x = self.x + incxkm
        self.y = self.y + incykm
        # Note that incxyonly updates x and y without
        # updating absloc which means that absloc will be
        # inconsistent with x and y after this operation.
        # To rectify this requires updateabsloc.
        # incxy performs both incxyonly and updateabsloc

    # To be used if incxyonly has changed the x and y values
    def updateabsloc( self ):
        azim_rad = math.atan2( self.x, self.y )
        azim_deg = math.degrees( azim_rad )
        geod = Geodesic.WGS84
        reflat = self.refloc.lat
        reflon = self.refloc.lon
        distkm = math.sqrt( self.x * self.x + self.y * self.y )
        g = geod.Direct( reflat, reflon, azim_deg, 1000.0 * distkm )
        self.absloc.lat = g['lat2']
        self.absloc.lon = g['lon2']

    def incxy( self, incxkm, incykm ):
        self.x = self.x + incxkm
        self.y = self.y + incykm
        azim_rad = math.atan2( self.x, self.y )
        azim_deg = math.degrees( azim_rad )
        geod = Geodesic.WGS84
        reflat = self.refloc.lat
        reflon = self.refloc.lon
        distkm = math.sqrt( self.x * self.x + self.y * self.y )
        g = geod.Direct( reflat, reflon, azim_deg, 1000.0 * distkm )
        self.absloc.lat = g['lat2']
        self.absloc.lon = g['lon2']

    def incxyLimited( self, incxkm, incykm, maxldistkm ):
        provx  = self.x + incxkm
        provy  = self.y + incykm
        azim_rad = math.atan2( provx, provy )
        azim_deg = math.degrees( azim_rad )
        geod = Geodesic.WGS84
        reflat = self.refloc.lat
        reflon = self.refloc.lon
        distkm = math.sqrt( provx * provx + provy * provy )
        g = geod.Direct( reflat, reflon, azim_deg, 1000.0 * distkm )
        provloc = llLocation( g['lat2'], g['lon2'] )
        distfromorigloc = dist_between_locs_km( self.abslocorig, provloc )
        if ( distfromorigloc < maxldistkm ):
            self.x = self.x + incxkm
            self.y = self.y + incykm
            self.absloc.lat = g['lat2']
            self.absloc.lon = g['lon2']

    def xyLoc2llLoc( self ):
        azim_rad = math.atan2( self.x, self.y )
        azim_deg = math.degrees( azim_rad )
        geod = Geodesic.WGS84
        reflat = self.refloc.lat
        reflon = self.refloc.lon
        distkm = math.sqrt( self.x * self.x + self.y * self.y )
        g = geod.Direct( reflat, reflon, azim_deg, 1000.0 * distkm )
        return llLocation( g['lat2'], g['lon2'] )

    def writehistory( self ):
        lenhist = len( self.history )
        if ( lenhist > 0 ):
            increment = int( lenhist / 1000 )
            if ( increment < 1 ):
                increment = 1
            filename = "orig_" + self.evid + "_history.txt" 
            exists = os.path.isfile( filename )
            if exists:
                os.remove( filename )
            f = open( filename, "w" )
            for i in range( 0, lenhist, increment ):
                historyvec = self.history[i]
                lat = historyvec[0]
                lon = historyvec[1]
                x   = historyvec[2]
                y   = historyvec[3]
                dis = historyvec[4]
                costfun  = historyvec[5]
                line  = "{:05d}".format(       i ).rjust(6)   + " "
                line += "{:.6f}".format(     lat ).rjust(10)  + " "
                line += "{:.6f}".format(     lon ).rjust(11)  + " "
                line += "{:.6f}".format(       x ).rjust(12)  + " "
                line += "{:.6f}".format(       y ).rjust(12)  + " "
                line += "{:.6f}".format(     dis ).rjust(12)  + " "
                line += "{:.6f}".format( costfun ).rjust(12)  + "\n"
                print ( line )
                f.write( line )

            f.close()

        print ( self.evid, lenhist )

# Here we want to read in an xyOrigin,
# a distance range in km (randomizeloc) and
# a maximum distance (maxldistkm)
# and we pick a random dx between -randomizeloc and +randomizeloc and
# a random dy between -randomizeloc and +randomizeloc
# and we relocate the vector by (dx, dy) provided that the new location
# is not greater than maxldistkm from the original location.
#==========================================================================
def addRandomness2xyOrigin( xyoriginstance, randomizeloc, maxldistkm ):
    nval = 1000
    ix  = random.randint( -nval, nval )
    rx  = float( ix )/float( nval )
    dx  = randomizeloc * rx
    iy  = random.randint( -nval, nval )
    ry  = float( iy )/float( nval )
    dy  = randomizeloc * ry
    xyoriginstance.incxyLimited( dx, dy, maxldistkm )
    return xyoriginstance

# And here exactly the same for a SloVecFlag object.
# a distance range randomizeslo and max. displacement maxsdist
# and we pick a random dsx between -randomizeslo and +randomizeslo and
# a random dsy between -randomizeslo and +randomizeslo
# and we move the slowness vector by (dsx,dsy) provided that we do not 
# shift it more than maxsdist from it's original vector specification.
#==========================================================================
def addRandomness2SloVec( slovecinstance, randomizeslo, maxsdist ):
    nval = 1000
    isx  = random.randint( -nval, nval )
    rsx  = float( isx )/float( nval )
    dsx  = randomizeslo * rsx
    isy  = random.randint( -nval, nval )
    rsy  = float( isy )/float( nval )
    dsy  = randomizeslo * rsy
    slovecinstance.incSxSyLimited( dsx, dsy, maxsdist )
    return slovecinstance

#==========================================================================
def dist_between_locs_km( loc1, loc2 ):
    geod = Geodesic.WGS84
    g = geod.Inverse( loc1.lat, loc1.lon, loc2.lat, loc2.lon )
    return 0.001 * g['s12']

#==========================================================================
def dist_between_locs_deg( loc1, loc2 ):
    geod = Geodesic.WGS84
    g = geod.Inverse( loc1.lat, loc1.lon, loc2.lat, loc2.lon )
    return g['a12']

#==========================================================================
def source_to_receiver_azimuth( rloc, sloc ):
    geod = Geodesic.WGS84
    g = geod.Inverse( sloc.lat, sloc.lon, rloc.lat, rloc.lon )
    azimuth = g['azi1']
    if ( azimuth < 0.0 ):
        azimuth = azimuth + 360.0
    return azimuth

#==========================================================================
def receiver_to_source_backazimuth( rloc, sloc ):
    geod = Geodesic.WGS84
    g = geod.Inverse( rloc.lat, rloc.lon, sloc.lat, sloc.lon )
    backazimuth = g['azi1']
    if ( backazimuth < 0.0 ):
        backazimuth = backazimuth + 360.0
    return backazimuth

#==========================================================================
def new_location_azi_distkm( loc1, azi, distkm ):
    geod = Geodesic.WGS84
    g = geod.Direct( loc1.lat, loc1.lon, azi, 1000.0 * distkm )
    return llLocation( g['lat2'], g['lon2'] )

#==========================================================================
#
# If we have slowness vector islo and slowness vector jslo
# and evinta and evintb Then Theta_ij_ab returns
# valid = 0 if we do not have corresponding observations
# and valid = 1 if we do and it returns Theta(i,j,a,b) 
#
def Theta_ij_ab( islo, jslo, evinta, evintb, display,
                 XYOrigs, SloVecFlags, RelTimeObs ):
    valid = 0
    Theta = 0.0
    if ( islo == jslo ):
        return valid, Theta
    if ( evinta == evintb ):
        return valid, Theta
    obsi = -1
    for i in range(0,len(RelTimeObs)):
        if ( RelTimeObs[i].evinda == evinta and
             RelTimeObs[i].evindb == evintb and
             RelTimeObs[i].svind  == islo ):
            obsi = i
    if ( obsi == -1 ):
        return valid, Theta
    obsj = -1
    for j in range(0,len(RelTimeObs)):
        if ( RelTimeObs[j].evinda == evinta and
             RelTimeObs[j].evindb == evintb and
             RelTimeObs[j].svind  == jslo ):
            obsj = j
    if ( obsj == -1 ):
        return valid, Theta
    tib = RelTimeObs[obsi].abstb
    tjb = RelTimeObs[obsj].abstb
    tia = RelTimeObs[obsi].absta
    tja = RelTimeObs[obsj].absta
    Tijab  = (tib-tjb) - (tia-tja)
    svecx  = SloVecFlags[islo].sx - SloVecFlags[jslo].sx
    svecy  = SloVecFlags[islo].sy - SloVecFlags[jslo].sy
    rvecx  = XYOrigs[evintb].x    - XYOrigs[evinta].x
    rvecy  = XYOrigs[evintb].y    - XYOrigs[evinta].y
    valid  = 1
    Theta  = Tijab + svecx*rvecx + svecy*rvecy
    if ( display ):
        vecterm = svecx*rvecx + svecy*rvecy
        string1 = "{:.8f}".format(  Tijab  ).rjust(13)
        string2 = "{:.8f}".format(  vecterm).rjust(13)
        string3 = "{:.8f}".format(  Theta  ).rjust(13)
        print (islo,jslo,evinta,evintb,string1,string2,string3)
    return valid, Theta

#==========================================================================
def read_abslocflagsfile( filename, refloc, randomizeloc, maxldistkm ):
    file_exist( filename )
    infile = open( filename, 'r' )
    XYOrigs      = []
    numLoc2Solve = 0
    for line in infile:
        words = line.split()
        ISOtimestring = words[0]
        t             = UTCDateTime( ISOtimestring )
        lat           = float( words[1] )
        lon           = float( words[2] )
        evloc         = llLocation( lat, lon )
        evid          = words[3]
        flag          = words[4]
        if ( flag != "F" and flag != "S" and flag != "I" ):
           print ("Need to specify F, S, or I as fifth column")
           exit()
        # Check that we do not already have event evid
        for i in range(0,len(XYOrigs)):
            if ( XYOrigs[i].evid == evid ):
                print ("Already read in event with id= ", evid )
                exit()
        ndt       = 0
        newOrig = xyOrigin( evloc, refloc, t, evid, ndt, flag )
        # We only add the event if we have not set the flag to "Ignore"
        if ( flag != "I" ):
            if ( flag == "S" ):
                numLoc2Solve += 1
                print ('before ', newOrig.x,' , ',newOrig.y )
                newOrig = addRandomness2xyOrigin( newOrig, randomizeloc, maxldistkm )
                print ('after  ', newOrig.x,' , ',newOrig.y )
            XYOrigs.append( newOrig )
    infile.close()
    return XYOrigs, numLoc2Solve

#==========================================================================
# Note that the SloVecFlag class is different to SloVec in that
# it is not constructed from geographical coordinates and
# velocity models. It only needs statname phase sx sy and S/F flag
# ndt is the number of observations associated with that slowness vector
# However, we take with us statlat, statlon, reflat, reflon
# for future book-keeping.
#
class SloVecFlag:
    def __init__( self, statname, phasename, sx, sy, flag, ndt,
                        statlat, statlon, reflat, reflon ):
        self.statname  = statname
        self.phasename = phasename
        self.statlat   = statlat
        self.statlon   = statlon
        self.reflat    = reflat
        self.reflon    = reflon
        self.origsx    = sx
        self.origsy    = sy
        self.sx        = sx
        self.sy        = sy
        self.mag       = math.sqrt( sx*sx + sy*sy )
        self.flag      = flag
        self.ndt       = ndt
        self.history   = []
 
    def incSxSy( self, incsxsperkm, incsysperkm ):
        self.sx  = self.sx + incsxsperkm
        self.sy  = self.sy + incsysperkm
        self.mag = math.sqrt( self.sx*self.sx + self.sy*self.sy )

    def incSxSyLimited( self, incsxsperkm, incsysperkm, maxsdist ):
        provsx   = self.sx + incsxsperkm
        provsy   = self.sy + incsysperkm
        distsx   = provsx  - self.origsx
        distsy   = provsy  - self.origsy
        distance = math.sqrt( distsx*distsx + distsy*distsy )
        if ( distance < maxsdist ):
            self.sx  = provsx
            self.sy  = provsy
            self.mag = math.sqrt( self.sx*self.sx + self.sy*self.sy )

    def scale( self, scalingfac ):
        self.sx  = self.sx * scalingfac
        self.sy  = self.sy * scalingfac
        self.mag = math.sqrt( self.sx*self.sx + self.sy*self.sy )

    def writehistory( self ):
        lenhist = len( self.history )
        if ( lenhist > 0 ):
            increment = int( lenhist / 1000 )
            if ( increment < 1 ):
                increment = 1
            filename = "sv_" + self.statname + "_" + self.phasename + "_history.txt"
            exists = os.path.isfile( filename )
            if exists:
                os.remove( filename )
            f = open( filename, "w" )
            for i in range( 0, lenhist, increment ):
                historyvec = self.history[i]
                sx  = historyvec[0]
                sy  = historyvec[1]
                dis = historyvec[2]
                costfun  = historyvec[3]
                line  = "{:05d}".format(       i ).rjust(6)   + " "
                line += "{:.7f}".format(      sx ).rjust(13)  + " "
                line += "{:.7f}".format(      sy ).rjust(13)  + " "
                line += "{:.6f}".format(     dis ).rjust(12)  + " "
                line += "{:.6f}".format( costfun ).rjust(12)  + "\n"
                print ( line )
                f.write( line )

            f.close()
        print ( self.statname, self.phasename, lenhist )
        


#==========================================================================
def read_slovecsflagsfile( filename, randomizeslo, maxsdist ):
    file_exist( filename )
    infile = open( filename, 'r' )
    SloVecFlags  = []
    numSlo2Solve = 0
    for line in infile:
        words = line.split()
        statname  = words[0]
        phasename = words[1]
        statlat   = float( words[2] )
        statlon   = float( words[3] )
        reflat    = float( words[4] )
        reflon    = float( words[5] )
        # Word 2 is the station latitude
        # Word 3 is the station longitude
        # Word 4 is the reference latitude
        # Word 5 is the reference latitude
        sx        = float( words[6] )
        sy        = float( words[7] )
        flag      = words[8]
        if ( flag != "F" and flag != "S" and flag != "I" ):
            print ("Need to specify F, S, or I as eighth column")
            exit()
        # Check that we do not already have event evid
        for i in range(0,len(SloVecFlags)):
            if ( SloVecFlags[i].statname  == statname and
                 SloVecFlags[i].phasename == phasename  ):
                print ("Already read in slowness vector with " )
                print ("Station name = ", stationname )
                print ("Phase name   = ", phasename )
                exit()
        # ndt is zero to start with when we read in the file
        ndt       = 0
        newSloVecFlag = SloVecFlag( statname, phasename, sx, sy, flag, ndt,
                                    statlat, statlon, reflat, reflon )
        # We only add the event if we have not set the flag to "Ignore"
        if ( flag != "I" ):
            if ( flag == "S" ):
                numSlo2Solve += 1
                print ('Before ', newSloVecFlag.sx,' , ',newSloVecFlag.sy )
                newSloVecFlag = addRandomness2SloVec( newSloVecFlag, randomizeslo, maxsdist )
                print ('After  ', newSloVecFlag.sx,' , ',newSloVecFlag.sy )
            SloVecFlags.append( newSloVecFlag )
    infile.close()
    return SloVecFlags, numSlo2Solve

#==========================================================================
def setUpUnknownRelLocCalcs( numLoc2Solve, XYOrigs, SloVecFlags, RelTimeObs ):
    print ("In setUpUnknownRelLocCalcs")
    print ("numLoc2Solve = ", numLoc2Solve )
    UnknownLocs = []
    numLocs     = len(XYOrigs)
    for iLoc in range(0,numLocs):
        flag = XYOrigs[iLoc].flag
        if ( flag == "S" ):
            newUnknownLoc = unknownRelLoc( iLoc, XYOrigs, 
                                           SloVecFlags, RelTimeObs )
            UnknownLocs.append( newUnknownLoc )

    return UnknownLocs
        
#==========================================================================
class unknownRelLoc:
    def __init__( self, evind, XYOrigs, SloVecFlags, RelTimeObs ):
        self.evind    = evind
        self.lochist  = XYOrigs[evind]
        self.obsPairs = []
        print ("New unknownRelLoc for event ", evind )
        #
        # We now need to loop around all of the actual observations
        # and find all pairs which have evind
        #
        for iobs in range(0,len(RelTimeObs)-1):
            # RelTimeObs[iobs].display()
            if ( RelTimeObs[iobs].evindb == self.evind ):
                RelTimeObs[iobs].swap()
            if ( RelTimeObs[iobs].evinda == self.evind ):
                evinda = evind
                evindb = RelTimeObs[iobs].evindb
                for jobs in range(iobs+1,len(RelTimeObs)):
                    if ( RelTimeObs[jobs].evindb == evinda ):
                        RelTimeObs[jobs].swap()
                    if ( RelTimeObs[jobs].evinda == evinda     and
                         RelTimeObs[jobs].evindb == evindb     and
                         RelTimeObs[jobs].svind != RelTimeObs[iobs].svind ):
                        newObsPair = obsPair( RelTimeObs[iobs],
                                              RelTimeObs[jobs]  )
                        self.obsPairs.append( newObsPair )
        
#==========================================================================
def setUpUnknownSloVecCalcs( numSlo2Solve, 
                             XYOrigs, SloVecFlags, RelTimeObs ):
    UnknownSloVecs = []
    numSloVecs     = len(SloVecFlags)
    for iSloVec in range(0,numSloVecs):
        flag = SloVecFlags[iSloVec].flag
        if ( flag == "S" ):
            newUnknownSloVec = unknownSloVec( iSloVec, XYOrigs,
                                           SloVecFlags, RelTimeObs )
            numobspairs = len(newUnknownSloVec.obsPairs)
            if ( numobspairs == 0 ):
                print ("Slowness vector = ",iSloVec)
                print ("No observation pairs found: ",numobspairs)
                exit()
        
            UnknownSloVecs.append( newUnknownSloVec )

    return UnknownSloVecs
        
#==========================================================================
class unknownSloVec:
    def __init__( self, svind, XYOrigs, SloVecFlags, RelTimeObs ):
        self.svind    = svind
        self.svhist   = SloVecFlags[svind]
        self.obsPairs = []
        print ("New unknownSloVec for slowness vector ", svind )
        #
        # We now need to loop around all of the actual observations
        # and find all pairs which have svind
        #
        # print ("len ", len(RelTimeObs) )
        # print ("Starting here")
        # print (SloVecFlags[svind].sx)
        # print (SloVecFlags[svind].sy)
        # for iobs in range(0,len(RelTimeObs)-1):
        for iobs in range(0,len(RelTimeObs)-1):
            # RelTimeObs[iobs].display()
            if ( RelTimeObs[iobs].svind == self.svind ):
                # print ("iobs", iobs)
                # RelTimeObs[iobs].display()
                evinda = RelTimeObs[iobs].evinda
                evindb = RelTimeObs[iobs].evindb
                # for jobs in range(iobs+1,len(RelTimeObs)):
                for jobs in range(0,len(RelTimeObs)):
                    # RelTimeObs[jobs].display()
                    # print ("ab ",evinda,evindb)
                    # RelTimeObs[jobs].display()
                    # RelTimeObs[jobs].display()
                    if ( RelTimeObs[jobs].evindb == evinda ):
                        RelTimeObs[jobs].swap()
                    if ( RelTimeObs[jobs].evinda == evinda     and
                         RelTimeObs[jobs].evindb == evindb     and
                         RelTimeObs[jobs].svind != RelTimeObs[iobs].svind ):
                        newObsPair = obsPair( RelTimeObs[iobs],
                                              RelTimeObs[jobs]  )
                        self.obsPairs.append( newObsPair )
                        # print ( len( self.obsPairs ) )

#==========================================================================
class obsPair:
    def __init__( self, RelObsi, RelObsj ):
        self.evinda = RelObsi.evinda
        self.evindb = RelObsi.evindb
        if ( self.evinda != RelObsj.evinda ):
            print ( 'obsPair problem with event indices' )
            print ( 'You should not have got here' )
            exit()
        if ( self.evindb != RelObsj.evindb ):
            print ( 'obsPair problem with event indices' )
            print ( 'You should not have got here' )
            exit()
        self.svindi = RelObsi.svind
        self.svindj = RelObsj.svind
        tib         = RelObsi.abstb
        tjb         = RelObsj.abstb
        tia         = RelObsi.absta
        tja         = RelObsj.absta
        # print("generating ijab ", self.svindi, self.svindj, self.evinda, self.evindb )
        self.Tijab  = (tib-tjb) - (tia-tja)

#==========================================================================
def Theta_ijab( oneObsPair, XYOrigs, SloVecFlags ):
    islo   = oneObsPair.svindi
    jslo   = oneObsPair.svindj
    evinta = oneObsPair.evinda
    evintb = oneObsPair.evindb
    svecx  = SloVecFlags[islo].sx - SloVecFlags[jslo].sx
    svecy  = SloVecFlags[islo].sy - SloVecFlags[jslo].sy
    rvecx  = XYOrigs[evintb].x    - XYOrigs[evinta].x
    rvecy  = XYOrigs[evintb].y    - XYOrigs[evinta].y
    Theta  = oneObsPair.Tijab + svecx*rvecx + svecy*rvecy
    Psi    = 0.5*Theta*Theta
    return Theta, Psi

#==========================================================================
def dPsi_ijab_dsi( oneObsPair, XYOrigs, SloVecFlags ):
    islo   = oneObsPair.svindi
    jslo   = oneObsPair.svindj
    evinta = oneObsPair.evinda
    evintb = oneObsPair.evindb
    svecx  = SloVecFlags[islo].sx - SloVecFlags[jslo].sx
    svecy  = SloVecFlags[islo].sy - SloVecFlags[jslo].sy
    rvecx  = XYOrigs[evintb].x    - XYOrigs[evinta].x
    rvecy  = XYOrigs[evintb].y    - XYOrigs[evinta].y
    Theta  = oneObsPair.Tijab + svecx*rvecx + svecy*rvecy
    dmagd  = math.sqrt( svecx*svecx + svecy*svecy )
    dmag1  = SloVecFlags[islo].mag
    dmag2  = SloVecFlags[jslo].mag
    weight = dmagd/(dmag1+dmag2)
    Psi    = 0.5*Theta*Theta
    dbysx  = Theta * rvecx
    dbysy  = Theta * rvecy
    return dbysx, dbysy, Theta, Psi, weight

#==========================================================================
def dPsi_ijab_dsj( oneObsPair, XYOrigs, SloVecFlags ):
    islo   = oneObsPair.svindi
    jslo   = oneObsPair.svindj
    evinta = oneObsPair.evinda
    evintb = oneObsPair.evindb
    svecx  = SloVecFlags[islo].sx - SloVecFlags[jslo].sx
    svecy  = SloVecFlags[islo].sy - SloVecFlags[jslo].sy
    rvecx  = XYOrigs[evintb].x    - XYOrigs[evinta].x
    rvecy  = XYOrigs[evintb].y    - XYOrigs[evinta].y
    Theta  = oneObsPair.Tijab + svecx*rvecx + svecy*rvecy
    dmagd  = math.sqrt( svecx*svecx + svecy*svecy )
    dmag1  = SloVecFlags[islo].mag
    dmag2  = SloVecFlags[jslo].mag
    weight = dmagd/(dmag1+dmag2)
    Psi    = 0.5*Theta*Theta
    dbysx  = Theta * rvecx
    dbysy  = Theta * rvecy
    return -dbysx, -dbysy, Theta, Psi, weight

#==========================================================================
def dPsi_ijab_dra( oneObsPair, XYOrigs, SloVecFlags ):
    islo   = oneObsPair.svindi
    jslo   = oneObsPair.svindj
    evinta = oneObsPair.evinda
    evintb = oneObsPair.evindb
    svecx  = SloVecFlags[islo].sx - SloVecFlags[jslo].sx
    svecy  = SloVecFlags[islo].sy - SloVecFlags[jslo].sy
    rvecx  = XYOrigs[evintb].x    - XYOrigs[evinta].x
    rvecy  = XYOrigs[evintb].y    - XYOrigs[evinta].y
    Theta  = oneObsPair.Tijab + svecx*rvecx + svecy*rvecy
    dmagd  = math.sqrt( svecx*svecx + svecy*svecy )
    dmag1  = SloVecFlags[islo].mag
    dmag2  = SloVecFlags[jslo].mag
    weight = dmagd/(dmag1+dmag2)
    Psi    = 0.5*Theta*Theta
    dbyrx  = Theta * svecx
    dbyry  = Theta * svecy
    return -dbyrx, -dbyry, Theta, Psi, weight

#==========================================================================
def dPsi_ijab_drb( oneObsPair, XYOrigs, SloVecFlags ):
    islo   = oneObsPair.svindi
    jslo   = oneObsPair.svindj
    evinta = oneObsPair.evinda
    evintb = oneObsPair.evindb
    svecx  = SloVecFlags[islo].sx - SloVecFlags[jslo].sx
    svecy  = SloVecFlags[islo].sy - SloVecFlags[jslo].sy
    rvecx  = XYOrigs[evintb].x    - XYOrigs[evinta].x
    rvecy  = XYOrigs[evintb].y    - XYOrigs[evinta].y
    Theta  = oneObsPair.Tijab + svecx*rvecx + svecy*rvecy
    dmagd  = math.sqrt( svecx*svecx + svecy*svecy )
    dmag1  = SloVecFlags[islo].mag
    dmag2  = SloVecFlags[jslo].mag
    weight = dmagd/(dmag1+dmag2)
    Psi    = 0.5*Theta*Theta
    dbyrx  = Theta * svecx
    dbyry  = Theta * svecy
    return dbyrx, dbyry, Theta, Psi, weight

#==========================================================================
class RelTObservation:
    def __init__( self, evinda, evindb, svind, absta, abstb, ccval ):
        self.evinda = evinda
        self.evindb = evindb
        self.svind  = svind
        self.absta  = absta
        self.abstb  = abstb
        self.ccval  = ccval

    def swap( self ):
        newevinda   = self.evindb
        newevindb   = self.evinda
        newabsta    = self.abstb
        newabstb    = self.absta
        self.evinda = newevinda
        self.evindb = newevindb
        self.absta  = newabsta
        self.abstb  = newabstb

    def display( self ):
        a = self.evinda
        b = self.evindb
        s = self.svind
        print ("Events ", a, b, " slovec ", s )

#==========================================================================
def newObs( singleObs, RelTimeObs ):
    isNewObs = True
    numCurrObs = len( RelTimeObs )
    for i in range( 0, numCurrObs ):
        if ( singleObs.evinda == RelTimeObs[i].evinda     and
             singleObs.evindb == RelTimeObs[i].evindb     and
             singleObs.svind  == RelTimeObs[i].svind  ):
            isNewObs = False
        if ( singleObs.evinda == RelTimeObs[i].evindb     and
             singleObs.evindb == RelTimeObs[i].evinda     and
             singleObs.svind  == RelTimeObs[i].svind  ):
            isNewObs = False
    # print ('hello ',singleObs.evinda, singleObs.evindb,isNewObs )
    return isNewObs

#==========================================================================
def modifySloVecsOnly( numrand, numSlo2Solve, XYOrigs, SloVecFlags, RelTimeObs,
                    maxiter, delslow, 
                    delslomax, maxsdist, maxpositiveslopes ):
    status = 0
    if ( numSlo2Solve < 1 ):
        print ("Called modifySloVecsOnly but numSlo2Solve = ", numSlo2Solve)
        return status, SloVecFlags, normValues

    print ("in modifySloVecsOnly")
    print ("len(RelTimeObs) = ", len(RelTimeObs) )
    unknownSloVecs = setUpUnknownSloVecCalcs( numSlo2Solve,
                             XYOrigs, SloVecFlags, RelTimeObs )
    print ("Number of unknown slowness vectors = ", len(unknownSloVecs) )
    print ("obspairs = ", len(unknownSloVecs[0].obsPairs) )

    SDistMoved = 99999.9
    Iterations = 0

#    scaleslow  = 1.0/float( len(SloVecFlags) )
#    slowmag0   = 0.0
#    avgsx0     = 0.0
#    avgsy0     = 0.0
#    for i in range( 0, len(SloVecFlags) ):
#        slowmag0 += scaleslow*SloVecFlags[i].mag
#        avgsx0   += scaleslow*SloVecFlags[i].sx
#        avgsy0   += scaleslow*SloVecFlags[i].sy

    normValues    = []
    percentnorm   = []
    runningmedian = []
    numpositiveslopes = 0
    prevslope     = 10.0
    while ( Iterations < maxiter ):

        Iterations += 1
        slovecvecs = []
        for unknownsv in range(0,len(unknownSloVecs) ):

            dsxarr = []
            dsyarr = []
            spiarr = []
            wgtarr = []
            #
            # Choose numrand random numbers between 0 and 
            # len(unknownSloVecs[unknownsv].obsPairs)
            #
            navail = len(unknownSloVecs[unknownsv].obsPairs) - 1
            for j in range(0,numrand):
                i      = random.randint(0,navail)
                dbysx,dbysy,Theta,Psi,weight = dPsi_ijab_dsi(
                                   unknownSloVecs[unknownsv].obsPairs[i],
                                   XYOrigs, SloVecFlags )
                dsxarr.append( dbysx )
                dsyarr.append( dbysy )
                spiarr.append( Psi )
                wgtarr.append( weight )

            dbydx, dbydy, costfv = vector_dxdy( dsxarr, dsyarr, spiarr, wgtarr )
            delslx = dbydx * delslow
            delsly = dbydy * delslow
            if ( delslx < -delslomax ):
                delslx = -delslomax
            if ( delslx >  delslomax ):
                delslx =  delslomax
            if ( delsly < -delslomax ):
                delsly = -delslomax
            if ( delsly >  delslomax ):
                delsly =  delslomax
            medspi = costfv
            newsshiftvec = [ delslx, delsly, medspi ]
            slovecvecs.append( newsshiftvec )

        #
        # Now we have a "shift" vector for every unknown slovec
        #
        percent    = 1.0
        SDistPrev  = SDistMoved
        SDistMoved = 0.0
        TotalSpi   = 0.0
        for unknownsv in range(0,len(unknownSloVecs) ):
            svind  = unknownSloVecs[unknownsv].svind
            delslx = -slovecvecs[unknownsv][0]
            delsly = -slovecvecs[unknownsv][1]
            medspi =  slovecvecs[unknownsv][2]
            sdist  = np.sqrt( delslx*delslx + delsly*delsly )
            SDistMoved += sdist
            TotalSpi   += medspi
            slowvec     = SloVecFlags[ svind ]
            sx          = slowvec.sx
            sy          = slowvec.sy
            print (sx, sy, delslx, delsly )
            delslx, delsly = slowinclimit( sx, sy, delslx, delsly, percent)
            SloVecFlags[ svind ].incSxSyLimited( delslx, delsly, maxsdist )

            slowvec     = SloVecFlags[ svind ]
            historyvec  = [ slowvec.sx, slowvec.sy, sdist, medspi ]
            SloVecFlags[ svind ].history.append( historyvec )
            azi, vel    = sxsyazivel( slowvec.sx, slowvec.sy)
            sxstring    = "{:.8f}".format(  slowvec.sx ).rjust(12)
            systring    = "{:.8f}".format(  slowvec.sy ).rjust(12)
            azistring   = "{:.3f}".format(  azi        ).rjust(7)
            velstring   = "{:.3f}".format(  vel        ).rjust(9)
            print (Iterations,svind,sxstring,systring,azistring,velstring)
   
        total2, numterms = calcMedianNormSetOfSloVecs( unknownSloVecs, XYOrigs, SloVecFlags )
        print ("SDistmoved ",SDistMoved," s/km SPI ", TotalSpi,
               " Ratio ", SDistMoved/SDistPrev )

        # TOTAL     = total2/float( numterms )
        # Total is a median value so do not divide by numterms
        TOTAL     = total2
        if ( Iterations == 1 ):
            initnorm   = TOTAL
            percentval = 100.0
        else:
            percentval  = 100.0*TOTAL/initnorm

        percentnorm.append( percentval )
        medval    = statistics.median( percentnorm[-25:] )
        runningmedian.append( medval )
        y         = np.array( runningmedian[-25:] )
        x         = np.arange(0,len(y))
        slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
        # if ( slope > 0.0  and  prevslope < 0.0 ):
        if ( Iterations == 1 ):
            slope = 0.0
        if ( slope > 0.0  ):
            numpositiveslopes += 1
        prevslope = slope

        if ( numpositiveslopes > maxpositiveslopes ):
            return status, SloVecFlags, normValues

        line = "ITERATION "
        line += "{:05d}".format( Iterations ).rjust(8)
        line += " TOTAL " + "{:.9f}".format( TOTAL ).rjust(15)
        line += " percent " + "{:.9f}".format( percentval ).rjust(15)
        line += " medval " + "{:.9f}".format( medval ).rjust(15)
        line += " slope " + "{:.9f}".format( slope ).rjust(15)
        line += "{:05d}".format( numpositiveslopes ).rjust(8)
        line += "\n"
        print( line )

        normvec = [ TOTAL, percentval, medval, slope, numpositiveslopes ]
        normValues.append( normvec )


    return status, SloVecFlags, normValues

#==========================================================================
def modifyLocsOnly( numrand, numLoc2Solve, XYOrigs, SloVecFlags, RelTimeObs,
                    maxiter, delkm, 
                    delkmmax, maxldistkm, FixEvent0, maxpositiveslopes ):
    status = 0
    if ( numLoc2Solve < 1 ):
        print ("Called modifyLocsOnly but numLoc2Solve = ", numLoc2Solve )
        return status, XYOrigs, normValues
    print ("in modifyLocsOnly")
    print ("len(RelTimeObs) = ", len(RelTimeObs) )
    unknownLocs = setUpUnknownRelLocCalcs( numLoc2Solve,
                     XYOrigs, SloVecFlags, RelTimeObs )
    print ("Number of unknown locations = ", len(unknownLocs) )
  
    DistMoved  = 99999.9
    Iterations = 0
    rx0        = XYOrigs[0].x
    ry0        = XYOrigs[0].y

    normValues    = []
    percentnorm   = []
    runningmedian = []
    numpositiveslopes = 0
    prevslope     = 10.0
    while ( Iterations < maxiter ):

        Iterations += 1
        rellocvecs = []
        for unknownev in range(0,len(unknownLocs) ):

            ddxarr = []
            ddyarr = []
            psiarr = []
            wgtarr = []
            #
            # Choose numrand random numbers between 0 and 
            # len(unknownLocs[unknownev].obsPairs)
            #
            navail = len(unknownLocs[unknownev].obsPairs) - 1
            for j in range(0,numrand):
                i      = random.randint(0,navail)
                ddx,ddy,Theta,Psi,weight = dPsi_ijab_dra( 
                                       unknownLocs[unknownev].obsPairs[i],
                                       XYOrigs, SloVecFlags )
                ddxarr.append( ddx )
                ddyarr.append( ddy )
                psiarr.append( Psi )
                wgtarr.append( weight )

            dbydx, dbydy, costfv = vector_dxdy( ddxarr, ddyarr, psiarr, wgtarr )
            delkmx = dbydx * delkm
            delkmy = dbydy * delkm
            if ( delkmx < -delkmmax ):
                delkmx = -delkmmax
            if ( delkmx >  delkmmax ):
                delkmx =  delkmmax
            if ( delkmy < -delkmmax ):
                delkmy = -delkmmax
            if ( delkmy >  delkmmax ):
                delkmy =  delkmmax
            medpsi = costfv
            newrellocvec = [ delkmx, delkmy, medpsi ]
            rellocvecs.append( newrellocvec )
            # print ( "Here ", unknownev, delkmx, delkmy )

        #
        # Now we have a relocation vector for every unknown location
        #
        DistPrev  = DistMoved
        DistMoved = 0.0
        TotalPsi  = 0.0
        for unknownev in range(0,len(unknownLocs) ):
            evind = unknownLocs[unknownev].evind
            delkmx = -rellocvecs[unknownev][0]
            delkmy = -rellocvecs[unknownev][1]
            medpsi =  rellocvecs[unknownev][2]
            dist   = np.sqrt( delkmx*delkmx + delkmy*delkmy )
            DistMoved += dist
            TotalPsi  += medpsi
            XYOrigs[ evind ].incxyLimited( delkmx, delkmy, maxldistkm )
            
            # print ( "There", unknownev, delkmx, delkmy )
            location = XYOrigs[ evind ].absloc
            historyvec  = [ location.lat, location.lon,
                            XYOrigs[ evind ].x, XYOrigs[ evind ].y,
                            dist, medpsi ]
            XYOrigs[ evind ].history.append( historyvec )
            latstring = "{:.5f}".format(  location.lat ).rjust(9) 
            lonstring = "{:.5f}".format(  location.lon ).rjust(9) 
            print (Iterations,evind,latstring,lonstring)
        
        #
        # Now shift all events laterally such that event0 is
        # back where it started - if we so chose.
        #
        if ( FixEvent0 == True ):
            diffx0 = rx0 - XYOrigs[ 0 ].x
            diffy0 = ry0 - XYOrigs[ 0 ].y
            for i in range(0,len(XYOrigs) ):
                XYOrigs[ i ].incxyLimited( diffx0, diffy0, maxldistkm )

        print ("Distmoved ",DistMoved," km PSI ", TotalPsi,
               " Ratio ", DistMoved/DistPrev )

        totalloc, numtermloc = calcMedianNormSetOfLocs( unknownLocs, XYOrigs, SloVecFlags )

        # Total is a median value so do not divide by numterms
        TOTAL     = totalloc
        if ( Iterations == 1 ):
            initnorm   = TOTAL
            percentval = 100.0
        else:
            percentval  = 100.0*TOTAL/initnorm
 
        percentnorm.append( percentval )
        medval    = statistics.median( percentnorm[-25:] )
        runningmedian.append( medval )
        y         = np.array( runningmedian[-25:] )
        x         = np.arange(0,len(y))
        slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
        # if ( slope > 0.0  and  prevslope < 0.0 ):
        if ( Iterations == 1 ):
            slope = 0.0
        if ( slope > 0.0 ):
            numpositiveslopes += 1
        prevslope = slope

        if ( numpositiveslopes > maxpositiveslopes ):
            return status, XYOrigs, normValues

        line = "ITERATION "
        line += "{:05d}".format( Iterations ).rjust(8)
        line += " TOTAL " + "{:.9f}".format( TOTAL ).rjust(15)
        line += " percent " + "{:.9f}".format( percentval ).rjust(15)
        line += " medval " + "{:.9f}".format( medval ).rjust(15)
        line += " slope " + "{:.9f}".format( slope ).rjust(15)
        line += "{:05d}".format( numpositiveslopes ).rjust(8)
        line += "\n"
        print( line )

        normvec = [ TOTAL, percentval, medval, slope, numpositiveslopes ]
        normValues.append( normvec )

        # print ( "rellocvecs" )
        # print (  rellocvecs )
 
    return status, XYOrigs, normValues

#==========================================================================
def modifyLocsAndSloVecs( numrand, numLoc2Solve, numSlo2Solve,
                          XYOrigs, SloVecFlags, RelTimeObs,
                          maxiter, delkm, delslow, delkmmax,
                          delslomax, maxsdist, maxldistkm, FixEvent0,
                          maxpositiveslopes, blocksize ):
    status = 0
    if ( numSlo2Solve < 1 ):
        print ("Called modifyLocsAndSloVecs but numSlo2Solve = ", numSlo2Solve)
        return status, XYOrigs, SloVecFlags, normValues
    if ( numLoc2Solve < 1 ):
        print ("Called modifyLocsAndSloVecs but numLoc2Solve = ", numLoc2Solve)
        return status, XYOrigs, SloVecFlags, normValues

    unknownSloVecs = setUpUnknownSloVecCalcs( numSlo2Solve,
                             XYOrigs, SloVecFlags, RelTimeObs )
    unknownLocs = setUpUnknownRelLocCalcs( numLoc2Solve,
                     XYOrigs, SloVecFlags, RelTimeObs )
    SDistMoved = 99999.9
    DistMoved  = 99999.9
    Iterations = 0
    rx0        = XYOrigs[0].x
    ry0        = XYOrigs[0].y
    if ( blocksize < 0 ):
        blocksize = 0
    whichblock = 0

#    scaleslow  = 1.0/float( len(SloVecFlags) )
#    slowmag0   = 0.0
#    avgsx0     = 0.0
#    avgsy0     = 0.0
#    for i in range( 0, len(SloVecFlags) ):
#        slowmag0 += scaleslow*SloVecFlags[i].mag
#        avgsx0   += scaleslow*SloVecFlags[i].sx
#        avgsy0   += scaleslow*SloVecFlags[i].sy

    normValues    = []
    percentnorm   = []
    runningmedian = []
    numpositiveslopes = 0
    prevslope     = 10.0
    while ( Iterations < maxiter ):

        Iterations += 1
        MoveLocs    = False
        MoveSloVecs = False
        if ( blocksize == 0 ):
            MoveLocs    = True
            MoveSloVecs = True
        else:
            remainder = Iterations % blocksize
            if ( remainder == 0 ):
                whichblock = (whichblock+1) % 2
            if ( whichblock == 0 ):
                MoveLocs    = True
            if ( whichblock == 1 ):
                MoveSloVecs = True

        rellocvecs = []
        slovecvecs = []
        if ( MoveLocs ):
            for unknownev in range(0,len(unknownLocs) ):
    
                ddxarr = []
                ddyarr = []
                psiarr = []
                wgtarr = []
                #
                # Choose numrand random numbers between 0 and 
                # len(unknownLocs[unknownev].obsPairs)
                #
                navail = len(unknownLocs[unknownev].obsPairs) - 1
                for j in range(0,numrand):
                    i      = random.randint(0,navail)
                    ddx,ddy,Theta,Psi,weight = dPsi_ijab_dra( 
                                           unknownLocs[unknownev].obsPairs[i],
                                           XYOrigs, SloVecFlags )
                    ddxarr.append( ddx )
                    ddyarr.append( ddy )
                    psiarr.append( Psi )
                    wgtarr.append( weight )
    
                dbydx, dbydy, costfv = vector_dxdy( ddxarr, ddyarr, psiarr, wgtarr )
                delkmx = dbydx * delkm
                delkmy = dbydy * delkm
                if ( delkmx < -delkmmax ):
                    delkmx = -delkmmax
                if ( delkmx >  delkmmax ):
                    delkmx =  delkmmax
                if ( delkmy < -delkmmax ):
                    delkmy = -delkmmax
                if ( delkmy >  delkmmax ):
                    delkmy =  delkmmax
                medpsi = costfv
                newrellocvec = [ delkmx, delkmy, medpsi ]
                rellocvecs.append( newrellocvec )
                # print ( "Here ", unknownev, delkmx, delkmy )

        if ( MoveSloVecs ):
            for unknownsv in range(0,len(unknownSloVecs) ):
    
                dsxarr = []
                dsyarr = []
                spiarr = []
                wgtarr = []
                #
                # Choose numrand random numbers between 0 and 
                # len(unknownSloVecs[unknownsv].obsPairs)
                #
                navail = len(unknownSloVecs[unknownsv].obsPairs) - 1
                for j in range(0,numrand):
                    i      = random.randint(0,navail)
                    dbysx,dbysy,Theta,Psi,weight = dPsi_ijab_dsi(
                                       unknownSloVecs[unknownsv].obsPairs[i],
                                       XYOrigs, SloVecFlags )
                    dsxarr.append( dbysx )
                    dsyarr.append( dbysy )
                    spiarr.append( Psi )
                    wgtarr.append( weight )
    
                dbydx, dbydy, costfv = vector_dxdy( dsxarr, dsyarr, spiarr, wgtarr )
                delslx = dbydx * delslow
                delsly = dbydy * delslow
                if ( delslx < -delslomax ):
                    delslx = -delslomax
                if ( delslx >  delslomax ):
                    delslx =  delslomax
                if ( delsly < -delslomax ):
                    delsly = -delslomax
                if ( delsly >  delslomax ):
                    delsly =  delslomax
                medspi = costfv
                newsshiftvec = [ delslx, delsly, medspi ]
                slovecvecs.append( newsshiftvec )

        #
        # Now we have a relocation vector for every unknown location
        #
        TotalPsi  = 0.0
        Ratio     = 999.9
        if ( MoveLocs ):
            DistPrev  = DistMoved
            DistMoved = 0.0
            for unknownev in range(0,len(unknownLocs) ):
                evind = unknownLocs[unknownev].evind
                delkmx = -rellocvecs[unknownev][0]
                delkmy = -rellocvecs[unknownev][1]
                medpsi =  rellocvecs[unknownev][2]
                dist   = np.sqrt( delkmx*delkmx + delkmy*delkmy )
                DistMoved += dist
                TotalPsi  += medpsi
                XYOrigs[ evind ].incxyLimited( delkmx, delkmy, maxldistkm )
    
                # print ( "There", unknownev, delkmx, delkmy )
                location = XYOrigs[ evind ].absloc
                historyvec  = [ location.lat, location.lon,
                                XYOrigs[ evind ].x, XYOrigs[ evind ].y,
                                dist, medpsi ]
                XYOrigs[ evind ].history.append( historyvec )
                latstring = "{:.5f}".format(  location.lat ).rjust(9)
                lonstring = "{:.5f}".format(  location.lon ).rjust(9)
                print (Iterations,evind,latstring,lonstring)
                Ratio = DistMoved/DistPrev

        print ("Distmoved ",DistMoved," km PSI ", TotalPsi,
               " Ratio ", Ratio )
        #
        # Now shift all events laterally such that event0 is
        # back where it started - if we so chose.
        #
        if ( MoveLocs ):
            if ( FixEvent0 == True ):
                diffx0 = rx0 - XYOrigs[ 0 ].x
                diffy0 = ry0 - XYOrigs[ 0 ].y
                for i in range(0,len(XYOrigs) ):
                    XYOrigs[ i ].incxyLimited( diffx0, diffy0, maxldistkm )
        #
        # Now we have a "shift" vector for every unknown slovec
        #
        percent    = 1.0
        TotalSpi   = 0.0
        Ratio      = 999.9
        if ( MoveSloVecs ):
            SDistPrev  = SDistMoved
            SDistMoved = 0.0
            for unknownsv in range(0,len(unknownSloVecs) ):
                svind  = unknownSloVecs[unknownsv].svind
                delslx = -slovecvecs[unknownsv][0]
                delsly = -slovecvecs[unknownsv][1]
                medspi =  slovecvecs[unknownsv][2]
                sdist  =  np.sqrt( delslx*delslx + delsly*delsly )
                SDistMoved += np.sqrt( delslx*delslx + delsly*delsly )
                TotalSpi   += medspi
                # print ("SG HERE ",delslx, delsly)
                slowvec     = SloVecFlags[ svind ]
                sx          = slowvec.sx
                sy          = slowvec.sy
                print (sx, sy, delslx, delsly )
                delslx, delsly = slowinclimit( sx, sy, delslx, delsly, percent)
                SloVecFlags[ svind ].incSxSyLimited( delslx, delsly, maxsdist )
    
                slowvec     = SloVecFlags[ svind ]
                historyvec  = [ slowvec.sx, slowvec.sy, sdist, medspi ]
                SloVecFlags[ svind ].history.append( historyvec )
                azi, vel    = sxsyazivel( slowvec.sx, slowvec.sy)
                sxstring    = "{:.8f}".format(  slowvec.sx ).rjust(12)
                systring    = "{:.8f}".format(  slowvec.sy ).rjust(12)
                azistring   = "{:.3f}".format(  azi        ).rjust(7)
                velstring   = "{:.3f}".format(  vel        ).rjust(9)
                print (Iterations,svind,sxstring,systring,azistring,velstring)
                Ratio       = SDistMoved/SDistPrev

        # total2, numterms = calcMedianNormSetOfSloVecs( unknownSloVecs, XYOrigs, SloVecFlags )
        totalloc, numtermloc = calcMedianNormSetOfLocs( unknownLocs, XYOrigs, SloVecFlags )

        print ("SDistmoved ",SDistMoved," s/km SPI ", TotalSpi,
               " Ratio ", Ratio, " T ",
               totalloc, numtermloc )

        # Total is a median value so do not divide by numterms
        # TOTAL     = totalloc/float( numtermloc )
        TOTAL     = totalloc
        if ( Iterations == 1 ):
            initnorm   = TOTAL
            percentval = 100.0
        else:
            percentval  = 100.0*TOTAL/initnorm

        percentnorm.append( percentval )
        medval    = statistics.median( percentnorm[-25:] )
        runningmedian.append( medval )
        y         = np.array( runningmedian[-25:] )
        x         = np.arange(0,len(y))
        slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
        # if ( slope > 0.0  and  prevslope < 0.0 ):
        if ( Iterations == 1 ):
            slope = 0.0
        if ( slope > 0.0 ):
            numpositiveslopes += 1
        prevslope = slope

        if ( numpositiveslopes > maxpositiveslopes ):
            return status, XYOrigs, SloVecFlags, normValues

        if ( MoveLocs ):
            if ( MoveSloVecs ):
                MoveString = "LS"
            else:
                MoveString = "L_"
        else:
            if ( MoveSloVecs ):
                MoveString = "_S"
            else:
                MoveString = "__"

        line = "ITERATION "
        line += "{:05d}".format( Iterations ).rjust(8)
        line += " TOTAL " + "{:.9f}".format( TOTAL ).rjust(15)
        line += " percent " + "{:.9f}".format( percentval ).rjust(15)
        line += " medval " + "{:.9f}".format( medval ).rjust(15)
        line += " slope " + "{:.9f}".format( slope ).rjust(15)
        line += "{:05d}".format( numpositiveslopes ).rjust(8)
        line += " " + MoveString + "\n"
        print( line )

        normvec = [ TOTAL, percentval, medval, slope, numpositiveslopes ]
        normValues.append( normvec )

    return status, XYOrigs, SloVecFlags, normValues

#==========================================================================
def read_dtfile( filename, allowmissing, XYOrigs, SloVecFlags ):
    file_exist( filename )
    infile = open( filename, 'r' )
    RelTimeObs = []
    numevents  = len( XYOrigs )
    numphases  = len( SloVecFlags )
    print ("In read_dtfile numevents = ", numevents )
    print ("In read_dtfile numphases = ", numphases )
    for line in infile:
        useme  = True
        words  = line.split()
        evidi  = words[0]
        evinda =  -1
        for i in range (0,len(XYOrigs)):
            if ( XYOrigs[i].evid == evidi ):
                evinda = i
        if ( evinda == -1 ):
            useme = False
        evidj  = words[1]
        evindb =  -1
        for j in range (0,len(XYOrigs)):
            if ( XYOrigs[j].evid == evidj ):
                evindb = j
        if ( evindb == -1 ):
            useme = False
        if ( evinda == evindb ):
            useme = False
        station = words[4]
        phase   = words[5]
        svind   =  -1
        for k in range (0,len(SloVecFlags)):
            if ( SloVecFlags[k].statname == station   and
                 SloVecFlags[k].phasename ==  phase ):
                svind = k
        if ( svind == -1 ):
           useme = False
        #
        # It looks like we now have a line we want to read.
        #
        if ( useme ):
            ISOtimestring = words[2]
            absta         = UTCDateTime( ISOtimestring )
            ISOtimestring = words[3]
            abstb         = UTCDateTime( ISOtimestring )
            ccval         = float( words[6] )
            newRelTObs = RelTObservation( evinda, evindb, svind,
                                          absta, abstb, ccval)
            isNewObs = newObs( newRelTObs, RelTimeObs )
            if ( isNewObs == True ):
                XYOrigs[evinda].ndt    += 1
                XYOrigs[evindb].ndt    += 1
                SloVecFlags[svind].ndt += 1
                RelTimeObs.append( newRelTObs )
        else:
            # Need to check that if allowmissing is false
            # that we stop if we have an unfamiliar
            # event or slowness vector
            if ( allowmissing == False ):
                if ( evinda == -1 or evindb == -1 or svind == -1 ):
                    print ("The following line contains an unknown")
                    print (line)
                    exit()

    infile.close()
    print ("leaving ", len(RelTimeObs) )
    # exit()
    return RelTimeObs

#==========================================================================
def write_rellocsfile( filename, XYOrigs ):
    exists = os.path.isfile( filename )
    if exists:
        os.remove( filename )
    f = open( filename, "w" )
    norigs = len( XYOrigs )
    for i in range( 0, norigs ):
        orig1   = XYOrigs[i]
        line    = str(orig1.origintime) + " "
        evid    = orig1.evid
        flag    = orig1.flag
        ndt     = orig1.ndt
        refloc  = orig1.refloc
        absloc  = orig1.absloc
        abslat  = absloc.lat
        abslon  = absloc.lon
        reflat  = refloc.lat
        reflon  = refloc.lon
        xval    = orig1.x
        yval    = orig1.y
        line   += "{:.5f}".format(  abslat ).rjust(9)  + " "
        line   += "{:.5f}".format(  abslon ).rjust(10) + " "
        line   += evid                                 + " "
        line   += flag                                 + " "
        line   += "{:.5f}".format(  reflat ).rjust(9)  + " "
        line   += "{:.5f}".format(  reflon ).rjust(10) + " "
        line   += "{:.8f}".format(    xval ).rjust(13) + " "
        line   += "{:.8f}".format(    yval ).rjust(13) + " "
        line   += "{:05d}".format(     ndt ).rjust(6)  + "\n"
        print ( line )
        f.write( line )

    f.close()

#==========================================================================
def write_nrmvaluesfile( filename, normValues ):
    exists = os.path.isfile( filename )
    if exists:
        os.remove( filename )
    f = open( filename, "w" )
    niter = len( normValues )
    for i in range( 0, niter ):
        nrmvec     = normValues[i]
        TOTAL        = nrmvec[0]
        percentval   = nrmvec[1]
        medval       = nrmvec[2]
        slope        = nrmvec[3]
        numposslopes = nrmvec[4]
        line = "IT "
        line += "{:06d}".format( i ).rjust(8)
        line += " T " + "{:.6e}".format( TOTAL ).rjust(15)
        line += " pc " + "{:.9f}".format( percentval ).rjust(15)
        line += " med " + "{:.9f}".format( medval ).rjust(15)
        line += " slope " + "{:.9f}".format( slope ).rjust(15)
        line += "{:05d}".format( numposslopes ).rjust(8)
        line += "\n"
        print ( line )
        f.write( line )

    f.close()

#==========================================================================
def write_slovecflagsfile( filename, SloVecFlags ):
    exists = os.path.isfile( filename )
    if exists:
        os.remove( filename )
    f = open( filename, "w" )
    nvecs = len( SloVecFlags )
    for i in range( 0, nvecs ):
        vec1 = SloVecFlags[i]
        Station = vec1.statname
        line    = Station.ljust(5) + " "
        Phase   = vec1.phasename
        line   += Phase.ljust(8) + " "
        statlat = vec1.statlat
        statlon = vec1.statlon
        reflat  = vec1.reflat
        reflon  = vec1.reflon
        sx      = vec1.sx
        sy      = vec1.sy
        origsx  = vec1.origsx
        origsy  = vec1.origsy
        flag    = vec1.flag
        ndt     = vec1.ndt
        line   += "{:.5f}".format( statlat ).rjust(9)  + " "
        line   += "{:.5f}".format( statlon ).rjust(10) + " "
        line   += "{:.5f}".format(  reflat ).rjust(9)  + " "
        line   += "{:.5f}".format(  reflon ).rjust(10) + " "
        line   += "{:.8f}".format(      sx ).rjust(13) + " "
        line   += "{:.8f}".format(      sy ).rjust(13) + " "
        line   += flag                                 + " "
        line   += "{:.8f}".format(  origsx ).rjust(13) + " "
        line   += "{:.8f}".format(  origsy ).rjust(13) + " "
        line   += "{:05d}".format(     ndt ).rjust(6)  + "\n"
        print ( line )
        f.write( line )

    f.close()

#
scriptname = sys.argv[0]
numarg     = len(sys.argv) - 1
text       = 'Specify '
text      += '--reflat [reflat] '
text      += '--reflon [reflon] '
text      += '--abslocfile [absolutelocationsfile] '
text      += '--locoutfile       [location_outfile] '
text      += '--slooutfile       [slowness_outfile] '
parser     = argparse.ArgumentParser( description = text )
parser.add_argument("--reflat", help="Reference latitude", default=None, required=True )
parser.add_argument("--reflon", help="Reference longitude", default=None, required=True )
parser.add_argument("--delkm", help="step length (km)", default=0.1000, required=False )
parser.add_argument("--delkmmax", help="max dist move (km)", default=0.0050, required=False )
parser.add_argument("--delslomax", help="max slowness move (s/km)", default=0.0005, required=False )
parser.add_argument("--delslow", help="step length (slowness)", default=0.10, required=False )
parser.add_argument("--maxsdist", help="max. dist allowed in slowness from original value", default=999.9, required=False )
parser.add_argument("--maxldistkm", help="max. dist in km allowed in location from original value", default=999.9, required=False )
parser.add_argument("--randomizeslo", help="dist within which to randomly shift slovecs", default=0.0, required=False )
parser.add_argument("--randomizeloc", help="dist within which to randomly shift locations", default=0.0, required=False )
parser.add_argument("--maxpositiveslopes", help="max. number of iterations norm can grow", default=100, required=False )
parser.add_argument("--blocksize", help="number of iterations to spend on either locations or slowness vectors", default=0, required=False )
parser.add_argument("--maxiter", help="max. iterations", default=10000, required=False )
parser.add_argument("--numrandom", help="num random picks", default=20, required=False )
parser.add_argument("--abslocflagsfile", help="absolute locations and flags file", default=None, required=True )
parser.add_argument("--slovecsflagsfile", help="slowness vectors and flags file", default=None, required=True )
parser.add_argument("--dtfile", help="differential times file", default=None, required=True )
parser.add_argument("--FixEvent0", type=str2bool, nargs='?',
                        const=True, default=False,
                        help="shift all events at each iteration to ensure event0 stays fixed",
                        required=False )
parser.add_argument("--writehistory", type=str2bool, nargs='?',
                        const=True, default=False,
                        help="write history file for each unknown",
                        required=False )
parser.add_argument("--allowmissing", type=str2bool, nargs='?',
                        const=True, default=False,
                        help="skip lines in dtfile if stati/phase not found",
                        required=False )
parser.add_argument("--locoutfile", help="location output file", default="locations.txt", required=False )
parser.add_argument("--slooutfile", help="slowness output file", default="slowness_vectors.txt", required=False )
parser.add_argument("--nrmoutfile", help="norms output file", default="normValues.txt", required=False )

args = parser.parse_args()

reflat           = float( args.reflat  )
reflon           = float( args.reflon  )
delkm            = float( args.delkm   )
delkmmax         = float( args.delkmmax   )
delslomax        = float( args.delslomax   )
delslow          = float( args.delslow )
maxsdist         = float( args.maxsdist )
maxldistkm       = float( args.maxldistkm )
randomizeslo     = float( args.randomizeslo )
randomizeloc     = float( args.randomizeloc )
maxpositiveslopes    = int( args.maxpositiveslopes )
maxiter          = int( args.maxiter  )
numrand          = int( args.numrandom  )
blocksize        = int( args.blocksize  )
abslocflagsfile  = args.abslocflagsfile
slovecsflagsfile = args.slovecsflagsfile
dtfile           = args.dtfile          
locoutfile       = args.locoutfile
slooutfile       = args.slooutfile
nrmoutfile       = args.nrmoutfile
allowmissing     = args.allowmissing
writehistory     = args.writehistory
FixEvent0        = args.FixEvent0

print ("allowmissing ", allowmissing )

refloc = llLocation( reflat, reflon )

XYOrigs, numLoc2Solve     = read_abslocflagsfile( abslocflagsfile, refloc, randomizeloc, maxldistkm )
SloVecFlags, numSlo2Solve = read_slovecsflagsfile( slovecsflagsfile, randomizeslo, maxsdist )
RelTimeObs      = read_dtfile( dtfile, allowmissing, XYOrigs, SloVecFlags )
numevents   = len( XYOrigs )
numslovecs  = len( SloVecFlags )
numdtobs    = len( RelTimeObs )
print ("Number of events = ", numevents )
print ("Number of slowness vectors = ", numslovecs )
print ("Number of dt observations = ", numdtobs )
print ("Number of events to solve = ", numLoc2Solve )
print ("Number of slowness vectors to solve = ", numSlo2Solve )
for i in range(0,numdtobs):
    print( i, RelTimeObs[i].evinda, str( RelTimeObs[i].absta ),
              RelTimeObs[i].evindb, str( RelTimeObs[i].abstb ),
              RelTimeObs[i].svind )
for j in range(0,numevents):
    print( XYOrigs[j].evid, XYOrigs[j].ndt )
for j in range(0,numslovecs):
    print( SloVecFlags[j].statname, SloVecFlags[j].phasename, 
           SloVecFlags[j].ndt )

if ( numLoc2Solve == 0 and numSlo2Solve == 0 ):
    print ("You have no unknowns. Nothing to solve for!")
    exit()

if ( numLoc2Solve > 0 and numSlo2Solve == 0 ):
    status, XYOrigs, normValues = modifyLocsOnly( numrand, numLoc2Solve,
                                      XYOrigs, SloVecFlags, RelTimeObs,
                                      maxiter, delkm, 
                                      delkmmax, maxldistkm, FixEvent0,
                                      maxpositiveslopes )
    print ("status = ", status )
    write_rellocsfile( locoutfile , XYOrigs )
    write_slovecflagsfile( slooutfile , SloVecFlags )
    write_nrmvaluesfile( nrmoutfile , normValues )
    if ( writehistory == True ):
        for i in range(0,numevents):
            if ( len( XYOrigs[i].history ) > 1 ):
                XYOrigs[i].writehistory()
    exit()

if ( numLoc2Solve == 0 and numSlo2Solve > 0 ):
    status, SloVecFlags, normValues = modifySloVecsOnly( numrand, numSlo2Solve,
                                      XYOrigs, SloVecFlags, RelTimeObs,
                                      maxiter, delslow, 
                                      delslomax, maxsdist,
                                      maxpositiveslopes )
    print ("status = ", status )
    write_rellocsfile( locoutfile , XYOrigs )
    write_slovecflagsfile( slooutfile , SloVecFlags )
    write_nrmvaluesfile( nrmoutfile , normValues )
    if ( writehistory == True ):
        for i in range(0,numslovecs):
            if ( len( SloVecFlags[i].history ) > 1 ):
                SloVecFlags[i].writehistory()
    exit()

#
# We can now assume that both numLoc2Solve and numSlo2Solve
# are greater than zero
#
status, XYOrigs, SloVecFlags, normValues = modifyLocsAndSloVecs( numrand,
                      numLoc2Solve, numSlo2Solve,
                      XYOrigs, SloVecFlags, RelTimeObs,
                      maxiter, delkm, delslow, 
                      delkmmax, delslomax, maxsdist, maxldistkm, FixEvent0,
                      maxpositiveslopes, blocksize )
print ("status = ", status )
write_rellocsfile( locoutfile , XYOrigs )
write_slovecflagsfile( slooutfile , SloVecFlags )
write_nrmvaluesfile( nrmoutfile , normValues )
if ( writehistory == True ):
    print ("hello steve")
    for i in range(0,numslovecs):
        print (len( SloVecFlags[i].history ) )
        if ( len( SloVecFlags[i].history ) > 0 ):
            SloVecFlags[i].writehistory()
    for i in range(0,numevents):
        if ( len( XYOrigs[i].history ) > 0 ):
            XYOrigs[i].writehistory()

exit()
