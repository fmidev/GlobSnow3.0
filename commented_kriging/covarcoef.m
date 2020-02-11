%function returns the coefficients of the covariance curve
%y(h)=c0+c1*exp(-h*c2)
%kertoimet
%
% function [coef(1),coef(2),coef2(2)]=covarcoef(c,d)
function [c0,c1,c2]=covarcoef(c,d)

% Change 22.10.2004/JP; select only positive values
i=1;
ind_apu=[];

if c(length(c))>=0 %to ensure that no negative values are included
    c(length(c)+1)=-0.1;
end
    
while c(i)>=0
    ind_apu = [ind_apu i];
    i=i+1;
end

c_apu=c(ind_apu);
c=[c_apu]; % first negative value to zero
d=d(1:length(c)); % including same number of distances!

t=transpose(d);
y=transpose(c);

%coefficient c2 is obtained by fitting covariance to the data
%
%exponent function y(h)=a*exp(h*c2)
r=length(y);

if(r == 0)
    c0 = 0;
    c1 = 0;
    c2 = 0;
else

    a=1;
    y2=log(y(1:r/a)/max(c));
    E2=[zeros(size(t(1:r/a))) t(1:r/a)];

    coef2=E2\y2;

    E=[ones(size(t)) exp(-t*abs(coef2(2)))];
    coef=E\y;

    c0=coef(1);
    c1=coef(2);
    c2=real(coef2(2));
end