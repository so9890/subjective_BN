function f = calibProd(x, MOM, list, paramss, polhelp)
% function determines pn, pg, omega and delta in equilibrium! 

% parameters
% read_in_pars_calib;
eppse=paramss(list.paramsdir=='eppse');
eppsy=paramss(list.paramsdir=='eppsy');
tauf =0; %polhelp(list.poldir=='tauf');
% vars
pn = exp(x(list.prod=='pn'));
pg = exp(x(list.prod=='pg'));
% omegaa = exp(x(list.prod=='omegaa'));
deltay = 1/(1+exp(x(list.prod=='deltay')));

%1) perceived as if Af, Ag, An in 2015-2019 are parameters but from here back
% out Af0, Ag0, An0
%1
% final output and demand final sector
    Y   = MOM.Y; 

pf  = MOM.FG^(-1/eppse)*pg-tauf; % optimality energy
pe  = ((pf+tauf).^(1-eppse)+pg.^(1-eppse)).^(1/(1-eppse)); % definition prices
E   = MOM.EpeY*Y/pe; 
N   = (pe/pn)^eppsy*(1-deltay)/deltay*E; % optimality final good
F   = E* (1+(1/MOM.FG)^((eppse-1)/eppse))^(-eppse/(eppse-1)); % uses optimality energy producers
    G   = F/MOM.FG; 
Yout = (deltay^(1/eppsy)*E^((eppsy-1)/eppsy)+(1-deltay)^(1/eppsy)*N^((eppsy-1)/eppsy))^((eppsy)/(eppsy-1)); 


% omissions
% f(1) = omegaa - MOM.grosemissionsUS2019/F;
f(1) = Yout-Y; 
f(2) = (pf+tauf)*F+pn*N+pg*G-Y; %market clearing demand
f(3)=  1-(deltay*pe^(1-eppsy)+(1-deltay)*pn^(1-eppsy))^(1/(1-eppsy));
end