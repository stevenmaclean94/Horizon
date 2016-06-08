function [sdec, srasn, gst] = suncor(idate, ut);
% ;+
% ; NAME:
% ;	SUNCOR
% ;
% ; PURPOSE:
% ;	Generates solar coordinates
% ;
% ; CAtEGORY:
% ;	utility
% ;
% ; CALLING SEQUENCE:  
% ;	SUNCOR, idate, ut, SDEC, SRASN, gst
% ;
% ; INPutS:
% ;	idate	Date in yyyyddd, longword integer
% ;	ut	time of day in seconds, utC, floating point
% ;
% ; OutPutS:  
% ;	sdec	solar declination, radians, floating point
% ;	srasn	solar right ascension, radians, floating point
% ;	gst	Greenwich sidereal time, radians, floating point
% ;
% ; COMMON BLOCKS:
% ;	None.
% ;
% ; PROCEDURE:
% ;	For a date and time or array of dates and times, the spherical
% ;	coordinates of the sun in the Earth-Centered Inertial (ECI) system
% ;	are returned as right ascension and declination.  Greenwich sidereal
% ;	time (the angle between the Greenwich meridian and the vernal equinox)
% ;	is also returned.
% ;	Dates prior to 2000 in yyddd are also accepted for backward
% ;	compatibility but all dates in an array must be in the same format.
% ;	Will not work properly after year 2100 due to lack of leap year.
% ;
% ; REFERENCE:
% ;	C.t. Russell, Geophysical Coordinate transforms.
% ;
% ; MODIFICAtION HIStORY:
% ;	~1983	Version F.1	Stan Solomon	Coded in Fortran
% ;	11/97	Version F.2	John Fulmer	Accept dates in yyyyddd format	
% ;	3/98	Version 1.0	Stan Solomon	Made into an IDL procedure
% ;	2/00	Version 1.1	Stan Solomon	Made double precision
% ;
% ;+

% PRO SUNCOR, idate, ut, SDEC, SRASN, gst

fday = ut/86400;
iyr = idate/1000;
iday = idate - iyr*1000;

if (iyr(1) >= 100)
    iyr = iyr - 1900;
end

dj = 365*iyr + (iyr-1)/4 + iday + fday - 0.5;
t = dj/36525;
vl = mod((279.696678+.9856473354*dj), 360);
gst = mod(279.696678+.9856473354*dj+360.*fday+180, 360*pi/180 );% MOD 360. * !PI/180.
g = mod(358.475845+.985600267*dj, 360*pi/180);% MOD 360. * !PI/180.
slong = vl+(1.91946-.004789*t)*sin(g)+.020094*sin(2.*g);
obliq = (23.45229-0.0130125*t)* pi/180; % *!PI/180.
slp = (slong-.005686)*pi/180; % * !PI/180.
sind = sin(obliq)*sin(slp);
cosd = sqrt(1.-sind^2);
sdec =atan(sind/cosd);
srasn = pi-atan2(1./tan(obliq)*sind/cosd,-cos(slp)/cosd);
