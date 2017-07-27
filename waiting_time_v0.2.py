
#!/usr/bin/env python

"""
    Purpose: compute time delay before observability
    Authors: A. Rossi, G. Stratta :)
    version 0.2
"""

print '------ Compute waiting time'
print '------ input is ra(hours) deg(deg) day  month  year  ut(hours)  '

#Import packages
import sys
import os
from astropy.coordinates import ICRS, Galactic
from astropy import units as u
#import scipy as sc
##from scipy.interpolate import interp1d
import numpy as np
from numpy import pi, arccos, arcsin, tan, sin , cos

# Default values , i.e Paranal
min_a=40. # degrees. Is the minimum observale elevation for big a big telescope

#TZp=-3.
#LTZp=-75.
#LAp=-24.6253   ##24:37:38
#LOp=289.5972   ##-70.4041666666667  ##70deg24'15''

#LActaN=28.7619444
#LOctaN=342.11

#LActaS=-24.6272222
#LOctaS=289.5958

print ''
print 'Lat and Long of famous sites: '
print 'Paranal: -24.6253 289.5972'
print 'CTA North: 28.7619 342.1100'
print 'CTA South: -24.6272 289.5958'


#LActaN=28.7619444
#LOctaN=342.11

#LActaS=-24.6272222
#LOctaS=289.5958'

tel=raw_input('Lat and Long of Observatory [e.g. Paranal -24.6253 289.5972]:') or '-24.6253 289.5972'
LAp=float(tel[0:7])
LOp=float(tel[8:16])

# ask for alpha,dec in hours,deg
#ra = float(sys.argv[1])
#dec = float(sys.argv[2])

ra_deg=raw_input('RA in deg [e.g. 49.81]:') or '49.81'
dec_deg=raw_input('Dec in deg [e.g. -24.75]:') or '-24.75'
ra=float(ra_deg)*24.0/360.0
dec=float(dec_deg)

#ask for UT of trigger in hours as ut(hours),d,m,y
##jd = float(sys.argv[3])
#dd=int(sys.argv[3])
#mm=int(sys.argv[4])
#yy=int(sys.argv[5])
#ut=float(sys.argv[6])
date=raw_input('date of observtion as YYYY-MM-DD [e.g. 2010-10-01]]:') or '2010-10-01'
yy=int(date[0:4])
mm=int(date[5:7])
dd=int(date[-2])

utall=raw_input('time of observation in hh:mm:ss UT [e.g. 11:22:33]:') or '11:22:33'
uts=float(utall[-2:])
utm=float(utall[3:5])
uth=float(utall[0:2])
ut=(((uts/60.)+utm)/60.)+uth

# compute julian date
#ddx=float(float(ut)/24.0)
#dd=dd+ddx
if mm<3 :
    yy=yy-1
    mm=mm+12

aa=int(yy/100.)
bb=2-aa+(aa/4.)
jd0=int(365.25*(yy+4716))+int(30.6001*(mm+1))+dd+bb-1524.5

#find sideral time at 0UT, GMST in hours
# http://aa.usno.navy.mil/faq/docs/GAST.php
h0=ut-0.5
jd=jd0+h0/24.
djd0=jd0-2451545.0
djd=jd - 2451545.0
T=(djd)/36525

stUT0=6.697374558 + 0.06570982441908*djd0 + 1.00273790935*h0 + 0.000026*T**2
if stUT0<0: stUT0=stUT0-int(stUT0/24.-1)*24.
else: stUT0=stUT0-int(stUT0/24.)*24.

# use default
#LA=LAp
#LO=LOp

LA=LAp
LO=LOp

# S1 is sideral time at UT trigger time
S1=stUT0+LOp/15.
if S1>24: S1=S1-int(S1/24)*24.
else: S1=S1

# Compute S2, i.e. sideral time when object comes visible with elevation min_a
# take the rising time, so minus, and convert to hours
h=(-arccos( -tan(dec*pi/180.)*tan(LA*pi/180.) + sin(min_a*pi/180.)/(cos(dec*pi/180.)*cos(LA*pi/180.)) ) ) *180./pi /15.
print h
S2=ra+h

# Compute S3, i.e. sideral time when object goes below elevation min_a
# take the setting time, so minus, and convert to hours
h2=(arccos( -tan(dec*pi/180.)*tan(LA*pi/180.) + sin(min_a*pi/180.)/(cos(dec*pi/180.)*cos(LA*pi/180.)) ) ) *180./pi /15.
print h2
S3=ra+h2

if S3-S1<0: S2=S2+24; S3=S3+24


#Waiting time for raising is
WTr=S2-S1  # in hours
#Waiting time for raising is
WTs=S3-S1  # in hours

print " "
print "JD at trigger is", jd,"days"
print "GMST at trigger is", stUT0,"hours"
print "hour angle at trigger is", S1-ra,"hours"
print "hour angle when gets visible is", h,"hours"
print "Sideral time at trigger is", S1,"hours"
print "Sideral time when gets visible is", S2,"hours"
print "Sideral time when is setting below ",min_a," degrees", S3,"hours"
print " "

print "Waiting time for raising is ",WTr,"hours"
print "Waiting time for setting is ",WTs,"hours"

print " "



#
## find number of days since vernal equinox (do better for bisestile)
#ml=(1,2,3,4,5,6,7,8,9,10,11,12)
#dml=(0,31,28,31,30,31,30,31,31,30,31,30,31)
#i=1
#ndm=0
#while (i<=m):
#	#print i
#	ndm=dml[i-1]+ndm
#	i=i+1
#vdm=31+28+21 # vernal equinox
#n=d+ndm-vdm  # days since vernal equinox
#
#print n,ndm
