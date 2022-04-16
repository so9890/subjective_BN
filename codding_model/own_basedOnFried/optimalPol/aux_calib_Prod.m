function [Y, pf, pe, pn, E, N, F, G, xg, xn, xf, C]=aux_calib_Prod(MOM, pg, Y, C, paramss,list, poll)

read_in_pars_calib;


% Y=MOM.Y;
An=MOM.An;
pf=(1/MOM.FG)^(1/eppse)*pg;
pe=(pf^(1-eppse)+pg^(1-eppse))^(1/(1-eppse)); 
pn=((1-deltay^eppsy.*pe.^(1-eppsy))./(1-deltay)^(eppsy)).^(1/(1-eppsy)); % definition prices and numeraire

E = Y/(deltay+(1-deltay)*(pe/pn*(1-deltay)/deltay)^(eppsy-1));
N = (pe/pn*(1-deltay)/deltay)^eppsy*E;
F = E/(1+(pf/pg)^(eppse-1))^(eppse/(eppse-1));
G = F/MOM.FG;
%G2=(pf/pg)^eppse*F;
% from production function intermediate goods and machine demand
% follos machines => together with Y determines C

%- production function
LfAf = F*(alphaf*pf*(1-tauf))^(-alphaf/(1-alphaf));
LgAg = G*(alphag*pg)^(-alphag/(1-alphag));
LnAn = N*(alphan*pn)^(-alphan/(1-alphan));
Ln   = LnAn/An; 
xn      = (alphan*pn).^(1/(1-alphan)).*LnAn;
xf      = (alphaf*pf.*(1-tauf)).^(1/(1-alphaf)).*LfAf;
xg      = (alphag*pg).^(1/(1-alphag)).*LgAg;

% price labour input neutral sector
hhnhln  = MOM.hhnhln;
thetan  = % follows from optimality neutral sector
wln     = pn.^(1/(1-alphan)).*(1-alphan).*alphan.^(alphan/(1-alphan).*An); 
wh      = thetan*(hlnhhn).^(1-thetan).wln; % from optimality labour input producers fossil, and demand labour fossil
wl      = (1-thetan)*(hhnhln).^(thetan).*wln;
% households
muu=C^(-thetaa);


end