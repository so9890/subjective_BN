function f = calibProd_GE(x, MOM, list, paramss)
% function determines pn, pg, omega and delta in equilibrium! 

% parameters
% read_in_pars_calib;
eppse=paramss(list.paramsdir=='eppse');
eppsy=paramss(list.paramsdir=='eppsy');

% vars
pn = exp(x(list.prod=='pn'));
pg = exp(x(list.prod=='pg'));
omegaa = exp(x(list.prod=='omegaa'));
deltay = 1/(1+exp(x(list.prod=='deltay')));

%1) perceived as if Af, Ag, An in 2015-2019 are parameters but from here back
% out Af0, Ag0, An0
%1
% final output and demand final sector
    Y   = MOM.Y; 

    FE  =(1-(MOM.GE)^((eppse-1)/eppse)).^(eppse/(eppse-1)); 
    FG  = FE/MOM.GE;
pf  =FG^(-1/eppse)*pg; % optimality energy
pe  = (pf.^(1-eppse)+pg.^(1-eppse)).^(1/(1-eppse)); % definition prices
E   = MOM.EpeY*Y/pe; 
N   = (pe/pn)^eppsy*((1-deltay)/deltay)*E; % optimality final good
F   = E*FE; % uses optimality energy producers
G   = E*MOM.GE; 
%
Yout = (deltay^(1/eppsy)*E^((eppsy-1)/eppsy)+(1-deltay)^(1/eppsy)*N^((eppsy-1)/eppsy))^((eppsy)/(eppsy-1)); 
%Yout = (deltay*E^((eppsy-1)/eppsy)+(1-deltay)*N^((eppsy-1)/eppsy))^((eppsy)/(eppsy-1)); 


% omissions
f(1) = omegaa - MOM.emissionsUS2019/F;
f(2) = Yout-Y; 
f(3) = pf*F+pn*N+pg*G-Y; %market clearing demand
f(4)=  1-(deltay*pe^(1-eppsy)+(1-deltay)*pn^(1-eppsy))^(1/(1-eppsy));
end