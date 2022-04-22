
function [A_lag, Lg, Ln, Lf, muu, E, SGov, N, Y,wln, wlg, wlf, xn, xg, xf ] ...
= auxiliary_stuff(params, list, pol, targets,  C, hhg, hhf, hhn, hlg, hln, hlf, F, G, wh, hh, hl, wl, An, Ag, Af, pn, pe, pf, pg , Ag_lag, An_lag, Af_lag, indic, thetag)

if indic.calib==1 
    paramss=params;
    poll=pol;
    read_in_pars_calib;
else
    read_in_params;
    read_in_pol;
end 

%- definitions
%A       = max([Af, Ag, An]')'; 
% Af_lag  = laggs(list.laggs=='Af_lag'); % shift Af backwards
% Ag_lag  = laggs(list.laggs=='Ag_lag');
% An_lag  = laggs(list.laggs=='An_lag');
A_lag   = (rhof*Af_lag+rhon*An_lag+rhog*Ag_lag)/(rhof+rhon+rhog);

Lg      = hhg.^thetag.*hlg.^(1-thetag);
Ln      = hhn.^thetan.*hln.^(1-thetan);
Lf      = hhf.^thetaf.*hlf.^(1-thetaf);

muu      = C.^(-thetaa); % same equation in case thetaa == 1
E       = (F.^((eppse-1)/eppse)+G.^((eppse-1)/eppse)).^(eppse/(eppse-1));

SGov    = zh*(wh.*hh*eh-lambdaa.*(wh.*hh*eh).^(1-taul))...
            +zl*(wl.*hl*el-lambdaa.*(wl.*hl*el).^(1-taul))...
            +tauf.*omegaa*pf.*F;
            % subsidies and profits and wages scientists cancel
N       =  (1-deltay)/deltay.*(pe./pn)^(eppsy).*E; % demand N final good producers 
Y       =  (deltay^(1/eppsy).*E.^((eppsy-1)/eppsy)+(1-deltay)^(1/eppsy).*N.^((eppsy-1)/eppsy)).^(eppsy/(eppsy-1)); % production function Y 

wln     = pn.^(1/(1-alphan)).*(1-alphan).*alphan.^(alphan/(1-alphan)).*An; % price labour input neutral sector
wlg     = pg.^(1/(1-alphag)).*(1-alphag).*alphag.^(alphag/(1-alphag)).*Ag;
wlf     = (1-alphaf)*alphaf^(alphaf/(1-alphaf)).*((1-tauf).*pf).^(1/(1-alphaf)).*Af; 

xn      = (alphan*pn).^(1/(1-alphan)).*Ln.*An;
xf      = (alphaf*pf.*(1-tauf)).^(1/(1-alphaf)).*Lf.*Af;
xg      = (alphag*pg).^(1/(1-alphag)).*Lg.*Ag;
end

