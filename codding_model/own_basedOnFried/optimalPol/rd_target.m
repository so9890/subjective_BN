function f= rd_target(x,  list,  params, laggs, Ems, T, indexx)
% function to find pf, tauf, taul and research allocation 
% under emission target
read_in_params;
% eppsilonf=1e-3;
phif = 12;
phig = 12;
phin = 12;
Sf   = S;
Sn   = S;
Sg   = S;

%- variables
gammasf     = x(list.targetALL=='gammasf').^2;
gammasn     = x(list.targetALL=='gammasn').^2;
gammasg     = x(list.targetALL=='gammasg').^2;
sff     = Sf/(1+exp(x(list.targetALL=='sff')));
sn     = Sn/(1+exp(x(list.targetALL=='sn')));
sg     = Sg/(1+exp(x(list.targetALL=='sg')));
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

% wsf = (gammaa*etaa*(A_lag./Af_lag).^phii.*(sff+eppsilonf).^(etaa-1).*pf.*(1-tauf).*F*(1-alphaf)*Af_lag)./(rhof^etaa.*Af); 
% wsg = (gammaa*etaa*(A_lag./Ag_lag).^phii.*(sg+eppsilonf).^(etaa-1).*pg.*G*(1-alphag)*Ag_lag)./(rhog^etaa.*(1-taus)*Ag);
% wsn = (gammaa*etaa*(A_lag./An_lag).^phii.*(sn+eppsilonf).^(etaa-1).*pn.*N*(1-alphan)*An_lag)./(rhon^etaa.*An);


%% target equations
q=0;
% LOM => lagged technology
q=q+1;
f((q-1)*T+1:T*q) = Af- Af_lag.*(1+gammaa*(sff/rhof).^etaa*(A_lag./Af_lag).^phii);
q=q+1;
f((q-1)*T+1:T*q) = Ag- Ag_lag.*(1+gammaa*(sg/rhog).^etaa*(A_lag./Ag_lag).^phii);
q=q+1;
f((q-1)*T+1:T*q) = An- An_lag.*(1+gammaa*(sn/rhon).^etaa*(A_lag./An_lag).^phii);
% scientists

q=q+1;
f((q-1)*T+1:T*q) = sff - ((phif+gammasf).*rhof^etaa.*Af./Af_lag.*(Af_lag./A_lag).^phii./(gammaa*etaa*pf.*(1-tauf).*F*(1-alphaf))).^(1/(etaa-1));
q=q+1;
f((q-1)*T+1:T*q) = sg  - ((phig+gammasg).*rhog^etaa.*Ag./Ag_lag.*(Ag_lag./A_lag).^phii.*(1-taus)./(gammaa*etaa*pg.*G*(1-alphag))).^(1/(etaa-1));
q=q+1;
f((q-1)*T+1:T*q) = sn  - ((phin+gammasn).*rhon^etaa.*An./An_lag.*(An_lag./A_lag).^phii./(gammaa*etaa*pn.*N*(1-alphan))).^(1/(etaa-1));

% wages scientists
q=q+1;
f((q-1)*T+1:T*q) = gammasf.*(Sf-sff);
q=q+1;
f((q-1)*T+1:T*q) =  gammasn.*(Sn-sn);
% market clearing scientists
q=q+1;
f((q-1)*T+1:T*q) =  gammasg.*(Sg-sg);
% f(7) = ws-wsg;
end