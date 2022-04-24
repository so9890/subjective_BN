function [hln, hlf, hlg,  hh, hl, Lg, Ln, Lf, pf, pe, pn, G, E, FF, N, Y, xn, xf, xg, wh, wl, whg, whn]...
=aux_calib2(MOM, deltay, hhn, hhg, hhf,zh, zl, el, eh, alphag, alphaf, alphan,  thetag, thetan, thetaf, eppsy, eppse, Ag, An, Af, pg, tauf)

hln = hhn*(1-thetan)/(thetan)*MOM.whwl; % hln
hlf = hhf*(1-thetaf)/(thetaf)*MOM.whwl; % hlf
hlg = hhg*(1-thetag)/(thetag)*MOM.whwl; % hlg 
%eh  = el/eleh;
hh  = (hhn+hhf+hhg)/(zh*eh); % high skill market clearing
hl  = (hln+hlf+hlg)/(zl*el); % low skill market clearing

Lg = hhg.^thetag.*hlg.^(1-thetag);
Ln = hhn.^thetan.*hln.^(1-thetan);
Lf = hhf.^thetaf.*hlf.^(1-thetaf);

% prices
pf=(1/MOM.FG)^(1/eppse)*pg;
pe=(pf^(1-eppse)+pg^(1-eppse))^(1/(1-eppse)); 
pn=((1-deltay.*pe.^(1-eppsy))./(1-deltay)).^(1/(1-eppsy)); % definition prices and numeraire


% output
G = (Ag.*Lg).*(pg.*alphag).^(alphag./(1-alphag));
FF = ((1-tauf)*alphaf*pf).^(alphaf/(1-alphaf))*(Af.*Lf); 
E = (FF.^((eppse-1)/eppse)+G.^((eppse-1)/eppse)).^(eppse/(eppse-1));
N = (An.*Ln).*(pn.*alphan).^(alphan./(1-alphan));
Y = (deltay^(1/eppsy).*E.^((eppsy-1)/eppsy)+(1-deltay)^(1/eppsy).*N.^((eppsy-1)/eppsy)).^(eppsy/(eppsy-1)); 

% machines
xn      = (alphan*pn)^(1/(1-alphan))*Ln*An;
xf      = (alphaf*pf.*(1-tauf))^(1/(1-alphaf))*Lf*Af;
xg      = (alphag*pg)^(1/(1-alphag))*Lg*Ag;

% wages 
wln     = pn*(1-alphan)*N/Ln; %pn.^(1/(1-alphan)).*(1-alphan).*alphan.^(alphan/(1-alphan)).*An; % price labour input neutral sector
wlg     = pg*(1-alphag)*G/Lg;%pg.^(1/(1-alphag)).*(1-alphag).*alphag.^(alphag/(1-alphag)).*Ag;
wlf     = pf*(1-tauf)*(1-alphaf)*FF/Lf; %(1-alphaf)*alphaf^(alphaf/(1-alphaf)).*((1-tauf).*pf).^(1/(1-alphaf)).*Af; 

wh      = thetaf*(hlf./hhf).^(1-thetaf).*wlf; % from optimality labour input producers fossil, and demand labour fossil
wl      = (1-thetaf)*(hhf./hlf).^(thetaf).*wlf;
whg     = thetag*(hlg./hhg).^(1-thetag).*wlg;
whn     = thetan*(hln./hhn).^(1-thetan).*wln;
end