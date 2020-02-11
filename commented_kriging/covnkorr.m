% covnkorr.m
% function returns covariance correlation fitted to data

function [cov_eta, korr_eta, kaeta, frekv] = covnkorr(D,Z,luokkamaara);

if nargin<3
    luokkamaara = 20;
end

maxdist=max(max(D));

sdnumber=length(D);

%the distances between stations are divided into n number of categories
%distjako = [maxdist/luokkamaara: maxdist/luokkamaara : maxdist];
%figure; plot(distjako,1:length(distjako),'x')
%hold on

%--------------
dists = D(find(D > 0));
[N,X] = hist(dists,length(dists));

scaled = 20*N/sum(N);
a = 0;
b = 0;
for i = 1:length(N)
    if(a >= 1)
	a = 0;
	b = b + 1;
	x(b) = X(i);
    end
    a = a + scaled(i);
end
x(length(x)+1) = max(max((D)));
distjako = x;

%--------------------------------

for k=1:luokkamaara
    if k==1
        [indx,indy] = find(0<D & D<=distjako(k)); %x row ja y column
  
        indeksi=find(0<D & D<=distjako(k));                        
        kaeta(k) = mean(D(indeksi));
    else
        [indx,indy] = find(D<=distjako(k) & D>distjako(k-1));
     
        indeksi=find(D<=distjako(k) & D>distjako(k-1));
        kaeta(k) = mean(D(indeksi));
    end    
    maara(k) = length(indx); %number of data points in each category
    if maara(k)>1
	Smiinush(k) = std(Z(indx));
	Splush(k) = std(Z(indy));
	ero_nelio(k) = sum((Z(indx) - Z(indy)).^2); %squared difference for pairs of data
	
	semivar_eta(1,k) = (ero_nelio(k))/(2*maara(k));
	tulot(k) = sum(Z(indx).*Z(indy));
	tulot2(k) = sum(Z(indx)).* sum(Z(indy));
	Mh(k) = mean(Z(indx));
	
	cov_eta(k)=(tulot(k)/maara(k) - tulot2(k)/(maara(k)^2));
	if(Smiinush(k) == 0 | Splush(k) == 0)
	    korr_eta(1,k) = 1;
	else
	    korr_eta(1,k) = cov_eta(k) / (Smiinush(k)*Splush(k));
	end
    else
	cov_eta(k) = -2;
	korr_eta(k)= -2;
	kaeta(k) = -2;
	semivar_eta(k) = -2;
	Sh1(1,k)=-2;
	Sh2(1,k)=-2;
    end
end

frekv=maara;


