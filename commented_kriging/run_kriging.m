%select dataset, load it
dname='GS_CP1_helmikuu';

fname=[dname '_lumilinja_biakset.mat'];

data=load(fname);

%select North America, Eurasia 1 and 2 pieces from data by longitude

na_index=find(data.longi<0);
eur1_index=find((data.longi>=0)&(data.longi<=65));
eur2_index=find(data.longi>55);

na_data=[data.lati(na_index);data.longi(na_index);data.BIAS_lumilinja(na_index)]';
eur1_data=[data.lati(eur1_index);data.longi(eur1_index);data.BIAS_lumilinja(eur1_index)]';
eur2_data=[data.lati(eur2_index);data.longi(eur2_index);data.BIAS_lumilinja(eur2_index)]';

%calculate kriging interpolation for North America and save results

na_k_string=['[' dname '_na_bias,' dname '_na_bias_var]=kriging_map(na_data,data.YI,data.XI,true);'];
eval(na_k_string);

save ([dname '_na_bias_map.mat'],[dname '_na_bias'],[dname '_na_bias_var']);

%calculate kriging interpolation for Eurasia and save results

eur1_k_string=['[' dname '_eur1_bias,' dname '_eur1_bias_var]=kriging_map(eur1_data,data.YI,data.XI,true);'];
eur2_k_string=['[' dname '_eur2_bias,' dname '_eur2_bias_var]=kriging_map(eur2_data,data.YI,data.XI,true);'];

eval(eur1_k_string);
eval(eur2_k_string);

save ([dname '_eur_bias_map.mat'],[dname '_eur1_bias'],[dname '_eur1_bias_var'],[dname '_eur2_bias'],[dname '_eur2_bias_var']);
