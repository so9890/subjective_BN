function f= rd_target(x,  list,  params, laggs, Ems, T, indexx)
% function to find pf, tauf, taul and research allocation 
% under emission target
read_in_params;

%- variables
sff     = exp(x(list.targetALL=='sff'));
sg      = exp(x(list.targetALL=='sg'));
sn      = exp(x(list.targetALL=='sn'));
%ws     = exp(x(list.targetALL=='ws'));
taul    = (1-exp(x(list.targetALL=='taul')));
pf      = exp(x(list.targetALL=='pf'));
taus    = x(list.targetALL=='taus'); 

% read in variables determined in other loop
[pn, pf, pg, F, G, N, Af, Ag, An, tauf] = polExp(pf, params, list, taul, T, Ems, indexx); 

%- initial condition
Af_lag=laggs(list.init=='Af0');
An_lag=laggs(list.init=='An0');
Ag_lag=laggs(list.init=='Ag0');

A_lag  = (rhof*Af_lag+rhon*An_lag+rhog*Ag_lag)/(rhof+rhon+rhog);

wsf = (gammaa*etaa*(A_lag./Af_lag).^phii.*sff.^(etaa-1).*pf.*(1-tauf).*F*(1-alphaf)*Af_lag)./(rhof^etaa.*Af); 
wsg = (gammaa*etaa*(A_lag./Ag_lag).^phii.*sg.^(etaa-1).*pg.*G*(1-alphag)*Ag_lag)./(rhog^etaa.*(1-taus)*Ag);
wsn = (gammaa*etaa*(A_lag./An_lag).^phii.*sn.^(etaa-1).*pn.*N*(1-alphan)*An_lag)./(rhon^etaa.*An);

%% target equations
q=0;
% LOM => lagged technology
q=q+1;
f((q-1)*T+1:T*q) = Af- Af_lag*(1+gammaa*(sff/rhof)^etaa*(A_lag/Af_lag)^phii);
q=q+1;
f((q-1)*T+1:T*q) = Ag- Ag_lag*(1+gammaa*(sg/rhog)^etaa*(A_lag/Ag_lag)^phii);
q=q+1;
f((q-1)*T+1:T*q) = An- An_lag*(1+gammaa*(sn/rhon)^etaa*(A_lag/An_lag)^phii);
% wages scientists
q=q+1;
f((q-1)*T+1:T*q) = wsf-wsg;
q=q+1;
f((q-1)*T+1:T*q) = wsn-wsg;
% market clearing scientists
q=q+1;
f((q-1)*T+1:T*q) = -S+sff+sn+sg; 
% f(7) = ws-wsg;
end