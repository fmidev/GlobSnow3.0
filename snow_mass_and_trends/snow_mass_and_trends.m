% snow_mass_and_trends.m
% JP, 5.1.2010
% includes Matlab routines to calculate snow mass and its trend 
% 
% 
clear all

load bias_GSv3_kriging_March %bias correction field and its variance
Vq2 = bias; %kriging-tulos (masked off areas as NaN)
Vq2_var = bias_var; % vastaava varianssi

SWE=zeros(39,721,721);
vuodet=[1980:1:2018]; %years of investigation 

% routine to load SWE data set (given in EaseGrid, 721 x 721)
for z=1:39 
    vuosi = 1979 + z; % 1980 - 2018
    v=num2str(vuosi);
    eval(['load ' v '03_northern_hemisphere_swe_0.25grid_monthly_param.mat'])
    SWE(z,:,:)=swe_ave_month; % value -2 in GlobSnow product is mountain
    SWE_corrected(z,:,:)=swe_ave_month-Vq2; % Bias correction!
end

load easelat % EaseGrid latitudes (721 x 721 array)
load easelon % muuttuja longitudes ((721 x 721 array))
lat=easelat;
lon=easelon;

[rows,cols] = find(lat>40); %hemispheric region of interest

for z=1:39
    SWE_vuosi(:,:) = SWE(z,:,:);
    SWE_vuosi_corrected(:,:) = SWE_corrected(z,:,:);
    ind_anal = find(lat>40 & SWE_vuosi>0 & isfinite(SWE_vuosi_corrected)==1); %masking of coasts and mountains
    NH_ala(z) = 625*length(ind_anal)/1e6;
    SWE_ave_March(z) = mean(SWE_vuosi(ind_anal));
    Snow_Mass_March(z) = 625e3 * length(ind_anal) * SWE_ave_March(z) /1e9;
    SWE_ave_March_corrected(z) = mean(SWE_vuosi_corrected(ind_anal));
    Snow_Mass_March_corrected(z) = 625e3 * length(ind_anal) * SWE_ave_March_corrected(z)/1e9;
    %variances
    VAR_SWE_ave_March_corrected(z) = (1/(length(ind_anal))^2) * sum(Vq2_var(ind_anal));
    VAR_Snow_Mass_March_corrected(z) = (625e3*length(ind_anal)/1e9)^2 * VAR_SWE_ave_March_corrected(z);
    
    % Continents
    indEUR = find(lat>40 & (lon>-10 | lon<-170) & SWE_vuosi>0 & isfinite(SWE_vuosi_corrected)==1);
    EUR_ala(z) = 625*length(indEUR)/1e6; % Snow area in Eurasia
    indNA = find(lat>40 & (lon<-10 & lon>-170) & SWE_vuosi>0 & isfinite(SWE_vuosi_corrected)==1);
    NA_ala(z) = 625*length(indNA)/1e6; % Snow area in North America
    SWE_ave_EUR_corrected(z) = mean(SWE_vuosi_corrected(indEUR));
    SWE_ave_EUR(z) = mean(SWE_vuosi(indEUR));
    Snow_Mass_EUR_corrected(z) = 625e3 * length(indEUR) * SWE_ave_EUR_corrected(z)/1e9;
    Snow_Mass_EUR(z) = 625e3 * length(indEUR) * SWE_ave_EUR(z)/1e9;
    SWE_ave_NA_corrected(z) = mean(SWE_vuosi_corrected(indNA));
    SWE_ave_NA(z) = mean(SWE_vuosi(indNA));
    Snow_Mass_NA_corrected(z) = 625e3 * length(indNA) * SWE_ave_NA_corrected(z)/1e9;
    Snow_Mass_NA(z) = 625e3 * length(indNA) * SWE_ave_NA(z)/1e9;
    %variances
    VAR_SWE_ave_EUR_corrected(z) = (1/(length(indEUR))^2) * sum(Vq2_var(indEUR));
    VAR_Snow_Mass_EUR_corrected(z) = (625e3*length(indEUR)/1e9)^2 * VAR_SWE_ave_EUR_corrected(z);
    VAR_SWE_ave_NA_corrected(z) = (1/(length(indNA))^2) * sum(Vq2_var(indNA));
    VAR_Snow_Mass_NA_corrected(z) = (625e3*length(indNA)/1e9)^2 * VAR_SWE_ave_NA_corrected(z);

end

% trends
XX = [[1:39]' ones(39,1)];
[Btrend,BINTtrend,Rtrend,RINTtrend,STATStrend] = regress(Snow_Mass_March_corrected',XX);
[Btrend_EUR,BINTtrend_EUR,Rtrend_EUR,RINTtrend_EUR,STATStrend_EUR] = regress(Snow_Mass_EUR_corrected',XX);
[Btrend_NA,BINTtrend_NA,Rtrend_NA,RINTtrend_NA,STATStrend_NA] = regress(Snow_Mass_NA_corrected',XX);

[sweBtrend,sweBINTtrend,sweRtrend,sweRINTtrend,sweSTATStrend] = regress(SWE_ave_March_corrected',XX);
[sweBtrend_EUR,sweBINTtrend_EUR,sweRtrend_EUR,sweRINTtrend_EUR,sweSTATStrend_EUR] = regress(SWE_ave_EUR_corrected',XX);
[sweBtrend_NA,sweBINTtrend_NA,sweRtrend_NA,sweRINTtrend_NA,sweSTATStrend_NA] = regress(SWE_ave_NA_corrected',XX);

Y = XX*Btrend; % trend!
Y_EUR = XX*Btrend_EUR;
Y_NA = XX*Btrend_NA;

sweY = XX*sweBtrend; % trend!
sweY_EUR = XX*sweBtrend_EUR;
sweY_NA = XX*sweBtrend_NA;


% decadal trends
decadal_trend = Btrend(1)*10;
decadal_trend_EUR = Btrend_EUR(1)*10;
decadal_trend_NA = Btrend_NA(1)*10;

swe_decadal_trend = sweBtrend(1)*10;
swedecadal_trend_EUR = sweBtrend_EUR(1)*10;
swe_decadal_trend_NA = sweBtrend_NA(1)*10;

limits = 10*[Btrend(1)-BINTtrend(1,1) Btrend(1) - BINTtrend(1,2)]; 
limits_EUR = 10*[Btrend_EUR(1)-BINTtrend_EUR(1,1) Btrend_EUR(1) - BINTtrend_EUR(1,2)];
limits_NA = 10*[Btrend_NA(1)-BINTtrend_NA(1,1) Btrend_NA(1) - BINTtrend_NA(1,2)];

swe_limits = 10*[sweBtrend(1)-sweBINTtrend(1,1) sweBtrend(1) - sweBINTtrend(1,2)]; 
swe_limits_EUR = 10*[sweBtrend_EUR(1)-sweBINTtrend_EUR(1,1) sweBtrend_EUR(1) - sweBINTtrend_EUR(1,2)]; 
swe_limits_NA = 10*[sweBtrend_NA(1)-sweBINTtrend_NA(1,1) sweBtrend_NA(1) - sweBINTtrend_NA(1,2)];


%STATStrend(3)
%STATStrend_EUR(3)
%STATStrend_NA(3)

%sweSTATStrend(3)
%sweSTATStrend_EUR(3)
%sweSTATStrend_NA(3)

%min, ave, max snow areas (Northern hemisphere, Eurasia, North America)
NH_alat = [min(NH_ala) mean(NH_ala) max(NH_ala)]; 
EUR_alat = [min(EUR_ala) mean(EUR_ala) max(EUR_ala)];
NA_alat = [min(NA_ala) mean(NA_ala) max(NA_ala)];


%trend maps:
ind_trendi = find(SWE_kartta>5); % only for areas with SWE>5mm
ind_ei_trendi = find(SWE_kartta<=5); % excluded

Trendi_kartta = NaN * ones(721,721); % initializatio of map
Trendi_kartta_corr = NaN * ones(721,721);

for k=1:length(ind_trendi)
        [B,BINT,R,RINT,STATS] = regress(SWE(:,ind_trendi(k)),XX);
        [Bcorr,BINTcorr,Rcorr,RINTcorr,STATScorr] = regress(SWE_corrected(:,ind_trendi(k)),XX);
        Trendi_kartta(ind_trendi(k)) = B(1)*10; %decadal trend!!!
end

% Routine to calculate snow mass in mountains (only applicable to MERRA2, Crocus and Brown):
% The GlobSnow-product includes the mask for mountains (grid cell value= -2) 

% mountains:
ind_vuoret=find(lat>40 & SWE_vuosi_GS==-2); % index to mountains
vuoret_ala(z) = 625*length(ind_vuoret)/1e6; % total area of mountains (hemispheric north of LAT 40)
SWE_ave_vuoret(z) = mean(SWE_vuosi(ind_vuoret)); % mean SWE on mountains, when SWE_vuosi includes
% average SWE over 1980-2018 by MERRA2, Crocus or Brown
Snow_Mass_vuoret(z) = 625e3 * length(ind_vuoret) * SWE_ave_vuoret(z) /1e9; % mean snow mass on mountains
    
indEUR = find(lat>40 & (lon>-10 | lon<-170) & SWE_vuosi_GS>0 & isfinite(SWE_vuosi_GS_corrected)==1);
EUR_ala(z) = 625*length(indEUR)/1e6; % area of Eurasian mountains

indNA = find(lat>40 & (lon<-10 & lon>-170) & SWE_vuosi_GS>0 & isfinite(SWE_vuosi_GS_corrected)==1);
NA_ala(z) = 625*length(indNA)/1e6; % area of North American mountains

SWE_ave_EUR(z) = mean(SWE_vuosi(indEUR));
Snow_Mass_EUR(z) = 625e3 * length(indEUR) * SWE_ave_EUR(z)/1e9;
SWE_ave_NA(z) = mean(SWE_vuosi(indNA));
Snow_Mass_NA(z) = 625e3 * length(indNA) * SWE_ave_NA(z)/1e9;


