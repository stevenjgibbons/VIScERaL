#!/usr/bin/env python3
#
# Steven J. Gibbons
# NGI, Sognsveien 72, Oslo
# 2021/06/14
# Find the theoretical slowness vectors from a reference
# location to specified stations for specified phases.
#
# Obtained from https://github.com/stevenjgibbons/statphase2slowvec 2023-06-17
# We just include this file here for convenience as it is used to generate
# the initial slowness vector files for the program VIScERaL.
#
try:
    import os
    import sys
    import argparse
    import numpy as np
    import geographiclib
    from geographiclib.geodesic import Geodesic
    import math
    import obspy
    from obspy.taup import TauPyModel
except ImportError as ie:
    miss_mod = ie.args[0].split()[3]
    print("\nThe Python module '" + miss_mod + "' is required.")
    print("Please install it and run again.\n")
    exit(1)

#==========================================================================
def file_exist(file):
    """Check if a file (with full path) exist"""
    if not os.path.isfile(file):
        print("File: ",file," does not exist. Bailing out ...")
        exit()

#==========================================================================
class llLocation:
    def __init__( self, lat, lon ):
        self.lat = lat
        self.lon = lon

#==========================================================================
class hSlownessVector:
    def __init__( self, sx, sy ):
        self.sx = sx
        self.sy = sy

#==========================================================================
class SloVec:
    def __init__( self, statname, phasename, refloc, statloc ):
        self.statname  = statname
        self.phasename = phasename
        self.refloc    = refloc
        self.statloc   = statloc
        svec           = src_rec_hSlownessVector( statloc, refloc, phasename )
        self.sx        = svec.sx
        self.sy        = svec.sy

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
def src_rec_hSlownessVector( rloc, sloc, phase ):
    ddeg = dist_between_locs_deg( rloc, sloc )
    model = TauPyModel( model = "ak135" )
    if ( phase == "P1" ):
        arrivals = model.get_travel_times( source_depth_in_km = 0.0,
                                           distance_in_degree = ddeg )
        ind = -1
        numarrivals = len( arrivals )
        for iarr in range( 0, numarrivals ):
            if ( ind == -1 and arrivals[iarr].name[0] == 'P' ):
                ind = iarr

    elif ( phase == "S1" ):
        arrivals = model.get_travel_times( source_depth_in_km = 0.0,
                                           distance_in_degree = ddeg )
        ind = -1
        numarrivals = len( arrivals )
        for iarr in range( 0, numarrivals ):
            if ( ind == -1 and arrivals[iarr].name[0] == 'S' ):
                ind = iarr
    else:
        arrivals = model.get_travel_times( source_depth_in_km = 0.0,
                                           distance_in_degree = ddeg,
                                           phase_list = [ phase ] )
        numarrivals = len( arrivals )
        if ( numarrivals == 0 ):
          print ('No arrivals for phase = ', phase)
          exit()
        ind = 0

    incidence = arrivals[ind].incident_angle
    inc_rad   = math.radians( incidence )
    if ( phase == "P" ):
      locvel    = 5.80
    elif ( phase == "P1" ):
      locvel    = 5.80
    elif ( phase == "Pn" ):
      locvel    = 5.80
    elif ( phase == "Pg" ):
      locvel    = 5.80
    elif ( phase == "S" ):
      locvel    = 3.46
    elif ( phase == "S1" ):
      locvel    = 3.46
    elif ( phase == "Sn" ):
      locvel    = 3.46
    elif ( phase == "Sg" ):
      locvel    = 3.46
    else:
      print ('Phase ',phase,' not familiar')
      exit()

    appvel    = locvel/math.sin( inc_rad )
    hSlow     = 1.0/appvel

    #
    azimuth_deg = source_to_receiver_azimuth( rloc, sloc )
    # print ( azimuth_deg )
    azimuth_rad = math.radians( azimuth_deg )
    sx = hSlow * math.sin( azimuth_rad )
    sy = hSlow * math.cos( azimuth_rad )
    # print (incidence, appvel)
    return hSlownessVector( sx, sy )

#==========================================================================
def write_slowvecs( filename, SlowVecs ):
    exists = os.path.isfile( filename )
    if exists:
        os.remove( filename )
    f = open( filename, "w" )
    nvecs = len( SlowVecs )
    for i in range( 0, nvecs ):
        vec1 = SlowVecs[i]
        Station = vec1.statname
        line    = Station.ljust(5) + " "
        Phase   = vec1.phasename
        line   += Phase.ljust(8) + " "
        refloc  = vec1.refloc
        statloc = vec1.statloc
        statlat = statloc.lat
        statlon = statloc.lon
        reflat  = refloc.lat
        reflon  = refloc.lon
        sx      = vec1.sx
        sy      = vec1.sy
        line   += "{:.5f}".format( statlat ).rjust(9)  + " "
        line   += "{:.5f}".format( statlon ).rjust(10) + " "
        line   += "{:.5f}".format(  reflat ).rjust(9)  + " "
        line   += "{:.5f}".format(  reflon ).rjust(10) + " "
        line   += "{:.8f}".format(      sx ).rjust(13) + " "
        line   += "{:.8f}".format(      sy ).rjust(13) + "\n"
        print ( line )
        f.write( line )

    f.close()

#==========================================================================
def read_phaselistfile( filename, refloc ):
    file_exist( filename )
    infile = open( filename, 'r' )
    SlowVecs = []
    for line in infile:
        words = line.split()
        statname  = words[0]
        phasename = words[1]
        # Check that we have not already read in this phase and station
        # combination
        for i in range(0,len(SlowVecs)):
            if ( SlowVecs[i].statname == statname and
                  SlowVecs[i].phasename == phasename ):
                print ("Already read in vectors with " )
                print (" station name = ", statname    )
                print (" phase name   = ", phasename    )
                exit()
        lat = float( words[2] )
        lon = float( words[3] )
        statloc = llLocation( lat, lon )
        newSloVec = SloVec( statname, phasename, refloc, statloc )
        SlowVecs.append( newSloVec )
    infile.close()
    return SlowVecs
#

#
scriptname = sys.argv[0]
numarg     = len(sys.argv) - 1
text       = 'Specify '
text      += '--reflat [reflat] '
text      += '--reflon [reflon] '
text      += '--phaselistfile [phaselistfile] '
text      += '--outfile       [outfile] '
parser     = argparse.ArgumentParser( description = text )
parser.add_argument("--reflat", help="Reference latitude", default=None, required=True )
parser.add_argument("--reflon", help="Reference longitude", default=None, required=True )
parser.add_argument("--phaselistfile", help="Phase list file", default=None, required=True )
parser.add_argument("--outfile", help="output file", default=None, required=True )

args = parser.parse_args()

reflat        = float( args.reflat )
reflon        = float( args.reflon )
phaselistfile = args.phaselistfile
outfile       = args.outfile

refloc = llLocation( reflat, reflon )
SlowVecs = read_phaselistfile( phaselistfile, refloc )
write_slowvecs( outfile, SlowVecs )
exit ()
nvecs = len( SlowVecs )
for i in range( 0, nvecs ):
    vec1 = SlowVecs[i]
    print ( vec1.sx, vec1.sy, vec1.statname, vec1.phasename )
