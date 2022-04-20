function f= calibjoint(x, MOM, list, paramss, poll, thetag, eleh)


read_in_pars_calib;

% calibration to 2019: Af, Ag, An => 2019 values
pg = exp(x(list.calibjo=='pg'));
hhn = exp(x(list.calibjo=='hhn'));
hhg = exp(x(list.calibjo=='hhg'));
hhf = exp(x(list.calibjo=='hhf'));
thetan = 1/(1+exp(x(list.calibjo=='thetan')));
thetaf = 1/(1+exp(x(list.calibjo=='thetaf')));
Af = exp(x(list.calibjo=='Af'));
Ag = exp(x(list.calibjo=='Ag'));
An = exp(x(list.calibjo=='An'));
el = exp(x(list.calibjo=='el'));
C  = exp(x(list.calibjo=='C'));
omegaa = exp(x(list.calibjo=='omegaa'));

% scientists
Af_lag  = exp(x(list.calibjo=='Af_lag'));
Ag_lag  = exp(x(list.calibjo=='Ag_lag'));
An_lag  = exp(x(list.calibjo=='An_lag'));
sf      = exp(x(list.calibjo=='sf'));
sg      = exp(x(list.calibjo=='sg'));
sn      = exp(x(list.calibjo=='sn'));

%-remaining variables (auxiliary)
[hln, hlf, hlg, eh, hh, hl, Lg, Ln, Lf, pf, pe, pn, G, E, FF, N, Y, xn, xf, xg, wh, wl, whg, whn]...
    =aux_calibjo(MOM, deltay, hhn, hhg, hhf,zh, zl, el, eleh, alphag, alphaf, alphan,  thetag, thetan, thetaf, eppsy, eppse, Ag, An, Af, pg, tauf);

A_lag  = (rhof*Af_lag+rhon*An_lag+rhog*Ag_lag)/(rhof+rhon+rhog);

wsf = (gammaa*etaa*(A_lag./Af_lag).^phii.*sf.^(etaa-1).*pf.*(1-tauf).*FF*(1-alphaf))./(rhof^etaa.*Af./Af_lag); 
wsg = (gammaa*etaa*(A_lag./Ag_lag).^phii.*sg.^(etaa-1).*pg.*G*(1-alphag))./(rhog^etaa.*(1-taus)*Ag./Ag_lag);
wsn = (gammaa*etaa*(A_lag./An_lag).^phii.*sn.^(etaa-1).*pn.*N*(1-alphan))./(rhon^etaa.*An./An_lag);
%SGov    = zh*(wh.*eh*hh-lambdaa.*(wh.*eh*hh).^(1-taul))...
%             +zl*(wl.*el*hl-lambdaa.*(wl.*el*hl).^(1-taul))...
%             +tauf.*pf.*omegaa*FF;

% target equations
f(1) = MOM.FG-FF/G; %Af
f(2) = E*pe/Y - MOM.EpeY; % Ag
f(3) = An - MOM.An; % An
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


% target equations
% f(1) = SGov; % lambdaa
% LOM => lagged technology
f(13) = Af- Af_lag*(1+gammaa*(sf/rhof)^etaa*(A_lag/Af_lag)^phii);
f(14) = Ag- Ag_lag*(1+gammaa*(sg/rhog)^etaa*(A_lag/Ag_lag)^phii);
f(15) = An- An_lag*(1+gammaa*(sn/rhon)^etaa*(A_lag/An_lag)^phii);
% wages scientists
f(16) = wsf-wsg;
f(17) = wsn-wsg;
% market clearing scientists
f(18) = -S+sf+sn+sg; 
end
