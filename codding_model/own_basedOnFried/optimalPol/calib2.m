function f= calib2(x, MOM, list, paramss, poll, thetag, eleh)

% in this function el eh are not correctly determined! 

read_in_pars_calib;

% calibration to 2019: Af, Ag, An => 2019 values
pg = exp(x(list.calib2=='pg'));
hhn = exp(x(list.calib2=='hhn'));
hhg = exp(x(list.calib2=='hhg'));
hhf = exp(x(list.calib2=='hhf'));
thetan = 1/(1+exp(x(list.calib2=='thetan')));
thetaf = 1/(1+exp(x(list.calib2=='thetaf')));
Af = exp(x(list.calib2=='Af'));
Ag = exp(x(list.calib2=='Ag'));
An = exp(x(list.calib2=='An'));
el = exp(x(list.calib2=='el'));
eh = exp(x(list.calib2=='eh'));
C  = exp(x(list.calib2=='C'));
omegaa = exp(x(list.calib2=='omegaa'));
lambdaa = exp(x(list.calib2=='lambdaa'));
deltay = 1/(1+exp(x(list.calib2=='deltay')));

%-remaining variables (auxiliary)
[hln, hlf, hlg, hh, hl, Lg, Ln, Lf, pf, pe, pn, G, E, FF, N, Y, xn, xf, xg, wh, wl, whg, whn]...
    =aux_calib2(MOM, deltay, hhn, hhg, hhf,zh, zl, el, eh, alphag, alphaf, alphan, thetag, thetan, thetaf, eppsy, eppse, Ag, An, Af, pg, tauf);

wln     = pn.^(1/(1-alphan)).*(1-alphan).*alphan.^(alphan/(1-alphan)).*An; % price labour input neutral sector
wlg     = pg.^(1/(1-alphag)).*(1-alphag).*alphag.^(alphag/(1-alphag)).*Ag;

% target equations
f(1) = MOM.FG-FF/G; %Af/Ag
f(2) = E*pe/Y - MOM.EpeY; % market share Epe = determines deltay
f(3) = Y - MOM.Y; % scales model!=> Af
% need two more equations to determine thetaf and thetan: assume equal 
f(4) = thetaf-thetan;
f(5) = (hhn+hhf)/(hhn+hln+hhf+hlf)-MOM.sharehighnongreen;
f(6) = Ag/An-MOM.AgAn; %(xn+xg+xf)/Y-(1-MOM.labourshare);% income share labour versus machines => 0.66; el
%f(6) = wlg/wln-MOM.wlgwln; %  Income share low skill already determined; scale!
f(7) = Y-xn-xf-xg-C;  %=> pg
f(8) = (1-alphaf)*(1-tauf)*pf*FF-(hhf)*wh/thetaf; % labour demand => determines hhf
f(9) = (pn*N*(1-alphan))-hhn*wh/thetan; % labour demand: hhn
f(10) = (pg*G*(1-alphag))-hhg*wh/thetag; % labour demand: hhg
% consumption: to ensure positive! (includes Gov spending)
f(11) = -C+ zh*wh*eh*hh+zl*wl*el*hl+tauf*pf*omegaa*FF; % => C
f(12) = omegaa - MOM.emissionsUS2019/FF; %=> omegaa
% GOV budget
f(13) = - MOM.Debt + zh*(wh.*eh*hh-lambdaa.*(wh.*eh*hh).^(1-taul))...
             +zl*(wl.*el*hl-lambdaa.*(wl.*el*hl).^(1-taul))+tauf.*pf.*omegaa*FF; % lambdaa
% wage premium
f(14) = whg/wl-MOM.whwl; %=> determines An/Af
f(15) = hh/hl - (MOM.whwl/el*eh)^((1-taul)/(taul+sigmaa)); % skill supply as fcn of wage ratio

end