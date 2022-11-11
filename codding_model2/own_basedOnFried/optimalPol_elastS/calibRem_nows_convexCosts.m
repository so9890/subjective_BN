function [c, ceq]= calibRem_nows_convexCosts(x, MOM, list, trProd, paramss, poll, Af, An, Ag)
% research with convex costs of research so that inner solution exists
% ohne market clearing 
read_in_pars_calib

% calibration to 2019 (lag = 2010-2014)
Af_lag  = exp(x(list.calib3=='Af_lag'));
Ag_lag  = exp(x(list.calib3=='Ag_lag'));
An_lag  = exp(x(list.calib3=='An_lag'));
% phis    = exp(x(list.calib3=='phis'));
sff     = exp(x(list.calib3=='sff'));
sg      = exp(x(list.calib3=='sg'));
sn      = exp(x(list.calib3=='sn'));
wsn      = exp(x(list.calib3=='wsn'));
wse      = exp(x(list.calib3=='wse'));

% sff      = S/(1+exp(x(list.calib3=='sff')));
% sg      = S/(1+exp(x(list.calib3=='sg')));
% sn      = S/(1+exp(x(list.calib3=='sn')));
% gammasg      = (x(list.calib3=='gammasg'))^2;
% gammasf      = (x(list.calib3=='gammasf'))^2;
% gammasn      = (x(list.calib3=='gammasn'))^2;


A_lag  = (rhof*Af_lag+rhon*An_lag+rhog*Ag_lag)/(rhof+rhon+rhog);

% required vars from previous calibration steps
[ ~, ~, ~, ~, pf, FF, pn, pg, ~, ~, ~, N, G,...
    ~, ~, ~, ~, ~, ~, ~, ~]=resProd(list, trProd, MOM , paramss, poll, 'calib'); 
   
% wsf = (gammaa*etaa*(A_lag./Af_lag).^phii.*sff.^(etaa-1).*pf.*(1-tauf).*FF*(1-alphaf)*Af_lag)./(rhof^etaa.*Af); 
% wsg = (gammaa*etaa*(A_lag./Ag_lag).^phii.*sg.^(etaa-1).*pg.*G*(1-alphag)*Ag_lag)./(rhog^etaa.*(1-taus)*Ag);
% wsn = (gammaa*etaa*(A_lag./An_lag).^phii.*sn.^(etaa-1).*pn.*N*(1-alphan)*An_lag)./(rhon^etaa.*An);

%% target equations
q=0;
% LOM => lagged technology
q=q+1;
ceq(q) = Af- Af_lag*(1+gammaa*(sff/rhof)^etaa*(A_lag/Af_lag)^phii);
q=q+1;
ceq(q) = Ag- Ag_lag*(1+gammaa*(sg/rhog)^etaa*(A_lag/Ag_lag)^phii);
q=q+1;
ceq(q) = An- An_lag*(1+gammaa*(sn/rhon)^etaa*(A_lag/An_lag)^phii);
% scientist demand
q=q+1;
ceq(q)= wse - (gammaa*etaa*(A_lag./Af_lag).^phii.*sff.^(etaa-1).*pf.*(1-tauf).*FF.*(1-alphaf).*Af_lag)./(rhof^etaa.*Af); 
%8
q=q+1;
ceq(q)= wse - (gammaa*etaa*(A_lag./Ag_lag).^phii.*sg.^(etaa-1).*pg.*G.*(1-alphag).*Ag_lag)./(rhog^etaa.*(1-taus)*Ag);
%9
q=q+1;
ceq(q)= wsn - (gammaa*etaa*(A_lag./An_lag).^phii.*sn.^(etaa-1).*pn.*N.*(1-alphan).*An_lag)./(rhon^etaa.*An);
% KT
% q=q+1;
% f(q)= gammasf.*(S-sff);
% q=q+1;
% f(q)= gammasn.*(S-sn);
% q=q+1;
% f(q)= gammasg.*(S-sg);
q=q+1;
ceq(q)= sn-S*2/3;
q=q+1;
ceq(q)= sff+sg-S*1/3;

c=[];
end