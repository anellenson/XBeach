function [x_frf, y_frf] = latLon2FRF(lat,lon);

% convert lat, lon to FRF the way kent's program does it.  Not quite right, but at least its consistent!
% Stolen from Kent's program "frf_coords" aka "Zeiss2LatLon"
%
%
% [x_frf, y_frf] = latLon2FRF(lat,lon);
%
%
% Input
%    lat, lon are lattitude and longitude in dd.dd form 
%		(decimal degrees!)
%    NOTE - longitudes should be negative, of order -75.75 degrees
%		versus +285 degrees
% 
% Output
%    x_frf, y_frf are frf plane coords
%
% check gps base station coordinates
% gps (lat,lon) =(36.105573810,-75.450431529) % ddmmss

% use Kent's origin and rotation
frf_lat = dms2deg(36, 10, 39.368875); 
frf_lon = dms2deg(-75, 44, 58.882497); 
frf_rot = -18.14182902;   % degrees
frf_rot = frf_rot*pi/180;   % radians



% conversion for degrees lat,lon to meters AT the frf lat, lon
m_lat = 364039.0*0.3048;
m_lon = 295100.0*0.3048;

% get vector from pier origin (meters)
y_frf = (lat - frf_lat)*m_lat;
x_frf = (lon - frf_lon)*m_lon;

% get azimuth of vector, zero is east (pos. x')
az = atan2(y_frf,x_frf); % radians

% get radius (meters)
r = sqrt((x_frf.^2) + (y_frf.^2));

% rotate to frf
x_frf = r.*cos(az+frf_rot);
y_frf = r.*sin(az+frf_rot);

function dms = dms2deg( deg, min, sec )

dms = abs(deg) + min/60 + sec/3600;
dms = dms * sign(deg);



%key Duck coordinate 
%comment  Convert lat/long to FRF coordinates 
