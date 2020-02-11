% krigpinta.m
% Determination of snow reference value by kriging interpolation
% JP 20.6.2006

function [pinta,errvar] = krigpinta(Z,lon,lat,Dist,...
					c0,c1,c2,n,m)

% Input:
% Z: snow depth or swe, [cm] or [mm], respectively
% lon: row vector on longitudes of reference stations
% lat: row vector on latitudes of reference stations
% Dist: distances between reference stations [km]
% c0,c1,c2: correlation function parameters
% n: longitude under investigation
% m: latitude under investigation

warning('off', 'MATLAB:rankDeficientMatrix');
warning('off', 'MATLAB:nearlySingularMatrix');
warning('off', 'MATLAB:singularMatrix');

warning off all;    % No more warnings...

dim=length(lat);
C=c1*exp(Dist.*c2)+c0;
C((C<0))=0;
C(dim+1,:)=ones(1,dim);
C(:,dim+1)=ones(dim+1,1);
C(dim+1,dim+1)=0;

% Distances from the current point to all reference points.
d=latlon_dist(lon,lat,n,m); % Spherical Earth
D=c1*exp(d*c2)+c0;

D((D < 0)) = 0;

D(length(D)+1)=1;

w = C\D';

pinta=transpose(Z)*w(1:dim);
        
errvar=C(1,1)-dot(w,D);
       
return       
