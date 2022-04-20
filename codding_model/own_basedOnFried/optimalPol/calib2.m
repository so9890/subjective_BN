function f= calib2(x, MOM, list, paramss, poll, thetag, eleh, hhgHH, hlgHL)


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
C  = exp(x(list.calib2=='C'));
omegaa = exp(x(list.calib2=='omegaa'));

%-remaining variables (auxiliary)
[hln, hlf, hlg, eh, hh, hl, Lg, Ln, Lf, pf, pe, pn, G, E, FF, N, Y, xn, xf, xg, wh, wl, whg, whn]...
    =aux_calib2(MOM, deltay, hhn, hhg, hhf,zh, zl, el, eleh, alphag, alphaf, alphan,  thetag, thetan, thetaf, eppsy, eppse, Ag, An, Af, pg, tauf);

% target equations
f(1) = MOM.FG-FF/G; %Af
f(2) = E*pe/Y - MOM.EpeY; % Ag
f(3) = An - 1; % An
% f(4) = hhgHH - hhg/(zh*eh*hh); % as targets for thetan and thetaf
% f(5) = hlgHL - hlg/(zl*el*hl);
%f(4)= hhg/(zh*eh*hh) - MOM.hg_total*(1+(1-zh)/zh*eleh*hl/hh)/(1+hlg/hhg); %=> determines share of high skill in green
%f(5)= hlg/(zl*el*hl) - MOM.hg_total*(1+zh/(1-zh)/eleh*hh/hl)/(1+hhg/hlg); %=> determines share of high skill in green
f(4) = wh-whg; 
f(5) = whg-whn; %MOM.hg_total-((hhg+hlg)/(eh*zh*hh+el*zl*hl)); => conflicts with output targets!
f(6) = hh/hl - (MOM.whwl*eh/el)^((1-taul)/(taul+sigmaa)); % el
f(7) = Y+xn+xf+xg-C;  %=> pg
f(8) = (1-alphaf)*(1-tauf)*pf*FF-(hhf)/thetaf; % labour demand
f(9) = (pn*N*(1-alphan))-hhn/thetan; % laboru demand
f(10) = (pg*G*(1-alphag))-hhg/(thetag); % labour demand
% consumption: to ensure positive!
f(11) = -C+ zh*wh*eh*hh+zl*wl*el*hl;
f(12) = omegaa - MOM.emissionsUS2019/FF; 
end