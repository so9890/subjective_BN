function [ C, Lnwln, Lgwlg, Lfwlf, pf, F, pn, pg, pee, E, Y, N, G, xn, xg, xf, ...
            AfLf, AgLg, AnLn, omegaa, deltay]=resProd(list, trProd, MOM, paramss, poll, quoi)

% parameters
if quoi== 'calib'
    read_in_pars_calib;
else
    params=paramss;
    pol=poll;
    read_in_params;
    read_in_pol;
end

pn=trProd(list.prod=='pn');
pg=trProd(list.prod=='pg');
% omegaa = trProd(list.prod=='omegaa');
deltay = trProd(list.prod=='deltay');

Y   = MOM.Y;
pf  = MOM.FG^(-1/eppse)*pg- tauf; % optimality energy
pee  = ((pf+tauf).^(1-eppse)+pg.^(1-eppse)).^(1/(1-eppse)); 
E   = MOM.EpeY*Y/pee; 
N   = (pee/pn)^eppsy*(1-deltay)/deltay*E; % optimality final good
F   = E* (1+(1/MOM.FG)^((eppse-1)/eppse))^(-eppse/(eppse-1)); % uses optimality energy producers
G   = F/MOM.FG; 
xn=pn*alphan*N;
xg=pg*(1+taus)*alphag*G;
xf=pf*alphaf*F;

C= Y-xn-xf-xg-MOM.Debt; 

AfLf    = F/(alphaf*pf)^(alphaf/(1-alphaf)); % production 
AgLg    = G/(alphag*pg)^(alphag/(1-alphag)); % production 
AnLn    = N/(alphan*pn)^(alphan/(1-alphan)); % production 

Lnwln   = pn.^(1/(1-alphan)).*(1-alphan).*alphan.^(alphan/(1-alphan)).*AnLn; % price labour input neutral sector
Lgwlg   = (pg.*(1+taus)).^(1/(1-alphag)).*(1-alphag).*alphag.^(alphag/(1-alphag)).*AgLg;
Lfwlf   = (1-alphaf)*alphaf^(alphaf/(1-alphaf)).*(pf).^(1/(1-alphaf)).*AfLf; 

% omega
omegaa = MOM.grosemissionsUS2019/F;
end