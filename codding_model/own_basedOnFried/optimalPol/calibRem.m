function f= calibRem(x, MOM, list, params, pol, targets, Af, An, Ag,  hhn, hhg, hhf, pg)
% this function backs out missing functions:
% 

read_in_params;
read_in_pol;

% calibration to 2019 (lag = 2010-2014)
%lambdaa = exp(x(list.calib3=='lambdaa'));
Af_lag  = exp(x(list.calib3=='Af_lag'));
Ag_lag  = exp(x(list.calib3=='Ag_lag'));
An_lag  = exp(x(list.calib3=='An_lag'));

sf      = exp(x(list.calib3=='sf'));
sg      = exp(x(list.calib3=='sg'));
sn      = exp(x(list.calib3=='sn'));
ws      =   exp(x(list.calib3=='ws'));
%-remaining variables (auxiliary)
[hln, hlf, hlg, eh, hh, hl, Lg, Ln, Lf, pf, pe, pn, G, E, FF, N, Y, xn, xf, xg, wh, wl, whg, whn]...
    =aux_calib2(MOM, deltay, hhn, hhg, hhf,zh, zl, el, el*eh, alphag, alphaf, alphan, thetag, thetan, thetaf, eppsy, eppse, Ag, An, Af, pg, tauf);
A_lag  = (rhof*Af_lag+rhon*An_lag+rhog*Ag_lag)/(rhof+rhon+rhog);

wsf = (gammaa*etaa*(A_lag./Af_lag).^phii.*sf.^(etaa-1).*pf.*(1-tauf).*FF*(1-alphaf))./(rhof^etaa.*Af./Af_lag); 
wsg = (gammaa*etaa*(A_lag./Ag_lag).^phii.*sg.^(etaa-1).*pg.*G*(1-alphag))./(rhog^etaa.*(1-taus)*Ag./Ag_lag);
wsn = (gammaa*etaa*(A_lag./An_lag).^phii.*sn.^(etaa-1).*pn.*N*(1-alphan))./(rhon^etaa.*An./An_lag);
%SGov    = zh*(wh.*eh*hh-lambdaa.*(wh.*eh*hh).^(1-taul))...
%             +zl*(wl.*el*hl-lambdaa.*(wl.*el*hl).^(1-taul))...
%             +tauf.*pf.*omegaa*FF;

% target equations
% f(1) = SGov; % lambdaa
% LOM => lagged technology
f(1) = Af- Af_lag*(1+gammaa*(sf/rhof)^etaa*(A_lag/Af_lag)^phii);
f(2) = Ag- Ag_lag*(1+gammaa*(sg/rhog)^etaa*(A_lag/Ag_lag)^phii);
f(3) = An- An_lag*(1+gammaa*(sn/rhon)^etaa*(A_lag/An_lag)^phii);
% wages scientists
f(4) = wsf-wsg;
f(5) = wsn-wsg;
% market clearing scientists
f(6) = -S+sf+sn+sg; 
f(7) = ws-wsg;
end
