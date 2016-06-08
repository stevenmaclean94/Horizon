function [lat, lon, alt] = eci2lla2(JD, eci)

% ; NAME:
% ;	ECI_TO_LLA
% ;
% ; PURPOSE:
% ;	Converts position vectors in ECI coordinates into Lat/Lon/Alt.
% ;
% ; CATEGORY:
% ;	Utility
% ;
% ; CALLING SEQUENCE:
% ;	ECI_TO_LLA, year, doy, utc, eci, lla
% ; 
% ; INPUTS:
% ;	year	Year, yyyy, longword integer
% ;	doy	Day of Year, ddd, longword integer
% ;	utc	Coordinated Universal Time of day in seconds, floating point
% ;       eci	ECI position vector, x, y, x, in km.
% ;
% ; OUTPUTS:
% ;	lla	latitude, longitude, altitude (degrees, degrees, km)
% ;
% ; KEYWORDS:
% ;	None
% ;
% ; COMMON BLOCKS:
% ;	None.
% ;
% ; PROCEDURE:
% ;	Transform Earth-Centered-Inertial position vector into
% ;	Geodetic latitude, longitude, and altitude above the surface.
% ;	Uses SUNCOR to find Greenwich sidereal time (GST), the angle between
% ;	the Greenwich meridian and the vernal equinox.  
% ;	Uses oblate spheroid approximation to shape of the Earth for altitude
% ;	and geodetic latitude calculation (ref.: W.J. Larson & J.R. Wertz,
% ;	Space Mission Analysis and Design, p. 809)
% ;	Arrays of vectors are OK!
% ;	
% ; ROUTINES USED:
% ;	SUNCOR - calculates coordinates of sun and Greenwich sidereal time
% ;
% ; MODIFICATION HISTORY:
% ;       Stan Solomon, 3/00
% ;
% ;-

% pro eci_to_lla, year, doy, utc, eci, lla

% ; f = Earth oblateness flattening factor, re = equatorial radius:
% f = 1./298.257D
% Using spherical Earth
f = 0;
re = 6378.137;

% ; Get Greenwich sidereal time:
J0 = floor(JD + 0.5) - 0.5;
UT = (JD - J0)*24.0;
T0 = (J0 - 2451545.0)/36525.0;
c = [100.4606184, 36000.77004, 0.000387933, -2.583e-8, 360.98564724];
g0 = c(1) + c(2)*T0 + c(3)*T0^2 + c(4)*T0^3;

g0 = g0 - 360.0*floor(g0/360.0);

gst = g0 + c(5)*UT/24.0;
gst = gst - 360.0*floor(gst/360.0);
gst = gst*pi/180;
% ; Calculate length of position vector:
rs = norm(eci);
%rs = sqrt(eci[0,*]^2+eci[1,*]^2+eci[2,*]^2)

% ; Calculate normalized position vector:
rn = eci/rs;
%rnx=eci[0,*]/rs
%rny=eci[1,*]/rs
%rnz=eci[2,*]/rs

% ; Calculate declination, geodetic latitude and altitude above oblate spheroid:
dec = asin(rn(3));
% dec = asin(rnz)
lat = atan( tan(dec)/(1-f)^2 );
alt = re * (rs/re-(1-f)/(sqrt(1-f*(2-f)*(cos(dec))^2)));

% ; Calculate  right ascension and geocentric longitude of satellite:
ra = atan2(rn(2),rn(1));
% ra = atan(rny,rnx)
lon = atan2( sin(ra-gst), cos(ra-gst) );

% ; Convert radians into degrees:
lat = lat * 180/pi;
lon = lon * 180/pi;
