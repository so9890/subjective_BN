function [ C, Lnwln, Lgwlg, Lfwlf, pf, F, pn, pg, pee, E, Y, N, G, xn, xg, xf, ...
            AfLf, AgLg, AnLn, omegaa, deltay]=resProd_GE(list, trProd, MOM, paramss, poll, quoi)

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
omegaa = trProd(list.prod=='omegaa');
deltay = trProd(list.prod=='deltay');

Y = MOM.Y;
    FE  =(1-(MOM.GE)^((eppse-1)/eppse)).^(eppse/(eppse-1)); 
    FG  = FE/MOM.GE;
pf  =FG^(-1/eppse)*pg; % optimality energy
pee  = (pf.^(1-eppse)+pg.^(1-eppse)).^(1/(1-eppse)); % definition prices
E   = MOM.EpeY*Y/pee; 
N   = (pee/pn)^eppsy*((1-deltay)/deltay)*E; % optimality final good
F   = E*FE; % uses optimality energy producers
G   = E*MOM.GE; 
xn=pn*alphan*N;
xg=pg*alphag*G;
xf=pf*(1-tauf)*alphaf*F;

C= Y-xn-xf-xg; 

AfLf    = F/(alphaf*pf*(1-tauf))^(alphaf/(1-alphaf)); % production 
AgLg    = G/(alphag*pg)^(alphag/(1-alphag)); % production 
AnLn    = N/(alphan*pn)^(alphan/(1-alphan)); % production 

Lnwln   = pn.^(1/(1-alphan)).*(1-alphan).*alphan.^(alphan/(1-alphan)).*AnLn; % price labour input neutral sector
Lgwlg   = pg.^(1/(1-alphag)).*(1-alphag).*alphag.^(alphag/(1-alphag)).*AgLg;
Lfwlf   = (1-alphaf)*alphaf^(alphaf/(1-alphaf)).*((1-tauf).*pf).^(1/(1-alphaf)).*AfLf; 
end