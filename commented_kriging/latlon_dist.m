function dist = latlon_dist_rs(lon1,lat1,lon2,lat2)
% function dist = latlon_dist_rs(lon1,lat1,lon2,lat2)
%
% Calculates the distance between two points defined by lat/lon
% coordinates. A spherical Earth geoid is used.

% Earth radius
a = 6370; % km

% degrees -> radians
lat1 = lat1/180*pi;
lon1 = lon1/180*pi;
lat2 = lat2/180*pi;
lon2 = lon2/180*pi;

if(lat1 == lat2 & lon1 == lon2)
    dist = 0;
else
    % Positional vectors on a unit sphere
    r1 = [cos(lon1).*cos(lat1); sin(lon1).*cos(lat1); sin(lat1)];
    r2 = [cos(lon2).*cos(lat2), sin(lon2).*cos(lat2), sin(lat2)];

    dot_product = r2*r1;
    alpha = acos(dot_product);
    
    % Angle between these vectors
    dist = a * alpha; % in km
end;
