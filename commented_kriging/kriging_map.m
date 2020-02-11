% kriging_map.m
% 
% Main function to calculate kriging map
%
% MT, 10.10.2011
%
%[ivalue,ivar] = kriging_map(station_data,latitude_matrix,longitude_matrix,interact);
%
%station_data is a n x 3 matrix where first column is latitude, second
%column longitude and third column data value 
%
%latitude_matrix is the matrix containing latitudes of pixel centers
%
%longitude_matrix is the matrix containing longitudes of pixel centers
%
%Latitudes and longitudes are in degrees. If value is NaN no interpolation
%estimation is calculated for that pixel.
%
%If Interact is set to true waitbar is prompted, if false then no waitbar
%is available


function [ivalue,ivar] = kriging_map(station_data,latitude_matrix,longitude_matrix,interact);

%initialize result variables
ivalue=0*latitude_matrix;
ivar=ivalue;

%filter NaN-observations from station data

inds=find(~isnan(station_data(:,3)));

filtered_data=station_data(inds,:);

%Calculate parameters for ordinary kriging
if(interact==true)
h=waitbar(0.5,'Calculating kriging parameters');
end
[c0,c1,c2,Dist] = kriging_ref(filtered_data(:,3)', filtered_data(:,2)', filtered_data(:,1)');

if(interact==true)
   close(h); 
end
%find indices that are not NaN
pixels=find(~isnan(latitude_matrix));

if(interact==true)
h=waitbar(0,'Calculating kriging interpolation');
end

%Do the kriging for all pixels
for px=1:length(pixels),
 
cur=pixels(px);

[ivalue(cur), ivar(cur)] = krigpinta(filtered_data(:,3),filtered_data(:,2)', filtered_data(:,1)',Dist,c0,c1,c2,longitude_matrix(cur),latitude_matrix(cur));
   
if (interact==true)
waitbar(px/length(pixels));
end

end

if(interact==true)
   close(h); 
end

end
