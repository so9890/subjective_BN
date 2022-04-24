function [muu, pf, Yout, pe, E, N, F, G, omegaa, Af, An, Ag, Lg, Ln, Lf, xf, xg, xn, ...
   SGov, hhD, hlD, hln, hlg, hlf, wlg, wln, wlf,...
  wh, wl] = aux_calibFinal(pn, hh, hl, deltay, eh, el, lambdaa, thetaf, thetag, thetan, C, pg, hhn, hhf, hhg, MOM, list, paramss, poll)
% parameters
read_in_pars_calib;


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
omegaa = MOM.emissionsUS2019/F;

% labour good producers
AfLf    = F/(alphaf*pf*(1-tauf))^(alphaf/(1-alphaf)); % production 
AgLg    = G/(alphag*pg)^(alphag/(1-alphag)); % production 
AnLn    = N/(alphan*pn)^(alphan/(1-alphan)); % production 
xn      = (alphan*pn).^(1/(1-alphan)).*AnLn; % machine demand
xf      = (alphaf*pf.*(1-tauf)).^(1/(1-alphaf)).*AfLf;
xg      = (alphag*pg).^(1/(1-alphag)).*AgLg;

% labour input good demand
LnwlnD     = pn.^(1/(1-alphan)).*(1-alphan).*alphan.^(alphan/(1-alphan)).*AnLn; % price labour input neutral sector
LgwlgD     = pg.^(1/(1-alphag)).*(1-alphag).*alphag.^(alphag/(1-alphag)).*AgLg;
LfwlfD     = (1-alphaf)*alphaf^(alphaf/(1-alphaf)).*((1-tauf).*pf).^(1/(1-alphaf)).*AfLf; 

% labour input
% optimality divided
hln = hhn*(1-thetan)/(thetan)*MOM.whwl; % hln
hlf = hhf*(1-thetaf)/(thetaf)*MOM.whwl; % hlf
hlg = hhg*(1-thetag)/(thetag)*MOM.whwl; % hlg 

% production: Contains info on whwl
Lg      = hhg.^thetag.*hlg.^(1-thetag); % supply
Ln      = hhn.^thetan.*hln.^(1-thetan);
Lf      = hhf.^thetaf.*hlf.^(1-thetaf);

Af = AfLf/Lf;
Ag = AgLg/Lg; 
An = AnLn/Ln; 

wln = LnwlnD/Ln; % demand= supply
wlg = LgwlgD/Lg;
wlf = LfwlfD/Lf;

wl = wlg/(MOM.whwl*thetag+(1-thetag));  %(1-thetag)*Lg*wlg/hlg; superfluous! instead use definition wlf 
wh = MOM.whwl*wl;

% low skill demand => determined by optimal divide and high skill demand

% total skill demand
hhD  = (hhn+hhf+hhg)/(zh*eh); % high skill market clearing
hlD  = (hln+hlf+hlg)/(zl*el); % low skill market clearing


%- Household  optimality

% budget (assume skill market clearing)
%CDemand = zh*wh*eh*hhD+ zl*wl*el*hlD+tauf*pf*omegaa*F; 

muu = C^(-thetaa);

%- research follows residually: 
%  Alag, Anlag, etc and scientistsm wages scientists
% A0   = (rhof*Af0+rhon*An0+rhog*Ag0)/(rhof+rhon+rhog);
% 
% wsf = (gammaa*etaa*(A0./Af0).^phii.*sf.^(etaa-1).*pf.*(1-tauf).*F*(1-alphaf))./(rhof^etaa.*Af./Af0); 
% wsg = (gammaa*etaa*(A0./Ag0).^phii.*sg.^(etaa-1).*pg.*G*(1-alphag))./(rhog^etaa.*(1-taus)*Ag./Ag0);
% wsn = (gammaa*etaa*(A0./An0).^phii.*sn.^(etaa-1).*pn.*N*(1-alphan))./(rhon^etaa.*An./An0);

SGov    = zh*(wh.*hh*eh-lambdaa.*(wh.*hh*eh).^(1-taul))...
            +zl*(wl.*hl*el-lambdaa.*(wl.*hl*el).^(1-taul))...
            +tauf.*omegaa*pf.*F;
end 