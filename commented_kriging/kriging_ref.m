% kriging_ref.m
% 
% Determination of spatial correlation parameters for kriging interpolation
% from observations representing a single day
%
% JP, 20.6.2006
%

function [c0,c1,c2,Dist] = kriging_ref(stations,lon,lat)

% Input:
% stations: row vector of snow depth [cm] or snow water equivalent
%           [mm] observed at weather stations or snow courses
% lon: row vector on longitudes of observations 
% lat: row vector on latitudes of observations 

X=lon';
Y=lat';
Z=stations';

classes=20; % distance classes

% Geoid:
%geoid = almanac('earth','geoid','kilometers');

dim=length(X);
Dist=zeros(dim,dim); % setting of upper triangle matrix of distances between stations

if(dim == 0)
    c0 = 0;
    c1 = 0;
    c2 = 0;
    Dist = NaN;
    return;
end;

%kkk = 1;
for iii=1:dim
    for jjj=iii:dim
	Dist(iii,jjj)=latlon_dist(X(iii),Y(iii),X(jjj),Y(jjj));
    end
end

if(dim>=classes)
    
[cov_eta, korr_eta, kaeta, frekv] = covnkorr(Dist,Z,classes);

Dist=Dist+transpose(triu(Dist,1));
[c0,c1,c2]=covarcoef(cov_eta,kaeta);

else
    Dist=Dist+transpose(triu(Dist,1));
    c0=1.04;
    c1=305;
    c2=-0.0008;
    
end

%figure; plot(kaeta, c1*exp(kaeta*c2)+c0)

return
