function f = calibProd(x, MOM, list, paramss)
% function which determines prices, omega and delta in equilibrium! 
% from here follow:

% labour good producers
% AfLf    = F/(alphaf*pf*(1-tauf))^(alphaf/(1-alphaf)); % production 
% AgLg    = G/(alphag*pg)^(alphag/(1-alphag)); % production 
% AnLn    = N/(alphan*pn)^(alphan/(1-alphan)); % production 
% xn      = (alphan*pn).^(1/(1-alphan)).*AnLn; % machine demand
% xf      = (alphaf*pf.*(1-tauf)).^(1/(1-alphaf)).*AfLf;
% xg      = (alphag*pg).^(1/(1-alphag)).*AgLg;
% 
% % labour input good demand: so that producers demand any amount of Ln
% LnwlnD     = pn.^(1/(1-alphan)).*(1-alphan).*alphan.^(alphan/(1-alphan)).*AnLn; % price labour input neutral sector
% LgwlgD     = pg.^(1/(1-alphag)).*(1-alphag).*alphag.^(alphag/(1-alphag)).*AgLg;
% LfwlfD     = (1-alphaf)*alphaf^(alphaf/(1-alphaf)).*((1-tauf).*pf).^(1/(1-alphaf)).*AfLf; 

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

pf  = MOM.FG^(-1/eppse)*pg; % optimality energy
pe  = (pf.^(1-eppse)+pg.^(1-eppse)).^(1/(1-eppse)); % definition prices
E   = MOM.EpeY*Y/pe; 
N   = (pe/pn)^eppsy*(1-deltay)/deltay*E; % optimality final good
F   = E* (1+(1/MOM.FG)^((eppse-1)/eppse))^(-eppse/(eppse-1)); % uses optimality energy producers
    G   = F/MOM.FG; 
Yout = (deltay^(1/eppsy)*E^((eppsy-1)/eppsy)+(1-deltay)^(1/eppsy)*N^((eppsy-1)/eppsy))^((eppsy)/(eppsy-1)); 


% omissions
f(1) = omegaa - MOM.emissionsUS2019/F;
f(2) = Yout-Y; 
f(3) = pf*F+pn*N+pg*G-Y; %market clearing demand
f(4)=  1-(deltay*pe^(1-eppsy)+(1-deltay)*pn^(1-eppsy))^(1/(1-eppsy));
end