function [ Lnwln, Lgwlg, pf, F, pee, E, Y, N, G, xn, xg,  ...
            AgLg, AnLn]=resProdTarget(list, pn, pf, pg, params, Ems)
% params
read_in_params;
%read_in_pol;

% vars: VECTORS
F       = (Ems+deltaa)/omegaa;
G       = F.*(pf./pg).^(eppse);
E       = F.*(1+(pf./pg).^(eppse-1)).^(eppse/(eppse-1));
pee     = (pf.^(1-eppse)+pg.^(1-eppse)).^(1/(1-eppse)); 
N       = (pee./pn).^eppsy*(1-deltay)./deltay.*E; % optimality final good
Y       = (deltay.^(1/eppsy).*E.^((eppsy-1)/eppsy)+(1-deltay).^(1/eppsy)*N.^((eppsy-1)/eppsy)).^(eppsy/(eppsy-1)); % final output production


xn=pn.*alphan.*N; % from optimality and production function
xg=pg.*alphag.*G;
% production 
AgLg    = G./(alphag*pg).^(alphag/(1-alphag)); % production 
AnLn    = N./(alphan*pn).^(alphan/(1-alphan)); % production 

Lnwln   = pn.^(1/(1-alphan)).*(1-alphan).*alphan.^(alphan/(1-alphan)).*AnLn; % price labour input neutral sector
Lgwlg   = pg.^(1/(1-alphag)).*(1-alphag).*alphag.^(alphag/(1-alphag)).*AgLg;



end