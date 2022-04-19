function [ pf, pe, pn, EY, NY, FY, GY, xgY, xnY, xfY]=aux_calib_Prod(MOM, pg, paramss,list, poll)

read_in_pars_calib;


pf=(1/MOM.FG)^(1/eppse)*pg;
pe=(pf^(1-eppse)+pg^(1-eppse))^(1/(1-eppse)); 
pn=((1-deltay^eppsy.*pe.^(1-eppsy))./(1-deltay)^(eppsy)).^(1/(1-eppsy)); % definition prices and numeraire

EY = 1/(deltay+(1-deltay)*(pe/pn*(1-deltay)/deltay)^(eppsy-1));
NY = (pe/pn*(1-deltay)/deltay)^eppsy*EY;
FY = EY/(1+(pf/pg)^(eppse-1))^(eppse/(eppse-1));
GY = FY/MOM.FG;

%G2=(pf/pg)^eppse*F;
% from production function intermediate goods and machine demand
% follos machines => together with Y determines C

%- production function
An=1;
LfAfY = FY*(alphaf*pf*(1-tauf))^(-alphaf/(1-alphaf));
LgAgY = GY*(alphag*pg)^(-alphag/(1-alphag));
LnAnY = NY*(alphan*pn)^(-alphan/(1-alphan));
LnY   = LnAnY/An; 
xnY      = (alphan*pn).^(1/(1-alphan)).*LnAnY;
xfY      = (alphaf*pf.*(1-tauf)).^(1/(1-alphaf)).*LfAfY;
xgY      = (alphag*pg).^(1/(1-alphag)).*LgAgY;




end