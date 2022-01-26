function levels= solution_SS(Ac, Ad, params, pols, list)

% function uses analytically derived BGP/initial SS equations to get model
% variables determined by initial values and parameters, and policy

% input:
% Ac, Ad:   initial technology
% params:   numeric vector of parameters
% pols:     numeric vector of policy variables

% output:
% levels: numeric vector of model variables on BGP, first period
%% preparations 

% read in required variables/parameters from input
thetac=params(list.params=='thetac');
thetad=params(list.params=='thetad');
sigmaa=params(list.params=='sigmaa');
zetaa=params(list.params=='zetaa');
eppsilon=params(list.params=='eppsilon');
alphaa=params(list.params=='alphaa');
psii=params(list.params=='psii');

tauul=pols(list.pol=='tauul');
vc=pols(list.pol=='vc');
vd=pols(list.pol=='vd');

% auxiliary variables/parameters
chic = (thetac/(zetaa*(1-thetac)))^(thetac)*((1-thetac)*(thetad-zetaa*(1-thetad)))/(thetad*(1-thetac)-thetac*(1-thetad));
chid = (thetad/(zetaa*(1-thetad)))^(thetad)*(1-((1-thetac)*(thetad-zetaa*(1-thetad)))/(thetad*(1-thetac)-thetac*(1-thetad)));

%% solutions

% sector good prices
pc = ((chic/chid*Ac/Ad)^((1-alphaa)*(1-eppsilon)/(alphaa+eppsilon*(1-alphaa)))+1)^(-1/(1-eppsilon));
pd = pc*(chic/chid*Ac/Ad)^((1-alphaa)/(alphaa+eppsilon*(1-alphaa)));

% labour input prices
pcL = (1-alphaa)*(alphaa/psii)^(alphaa/(1-alphaa))*pc^(1/(1-alphaa))*Ac;
pdL = (1-alphaa)*(alphaa/psii)^(alphaa/(1-alphaa))*pd^(1/(1-alphaa))*Ad;

% skill specific wage
wl = pdL*zetaa^(-thetad)*thetad^thetad*(1-thetad)^(1-thetad);
wh = zetaa*wl;

% household
H = (1-tauul)^(1/(1+sigmaa)); 
c = lambdaa*(H*wl)^(1-tauul);

% goods market clearing
Y = c;

% hl as defined in paper
hl= (alphaa/psii)^(-alphaa/(1-alphaa))*...
    ((pc^(alphaa/(1-alphaa))*chic*Ac)^((eppsilon-1)/eppsilon)+...
      (pd^(alphaa/(1-alphaa))*chid*Ad)^((eppsilon-1)/eppsilon)...
      )^(-eppsilon/(eppsilon-1))*lambdaa*(H*wl)^(1-tauul);

% high skill supply
hh = 1/zetaa*(H-hl);

% Labour input good production 
Lc=chic*hl;
Ld=chid*hl;

% skill inputs
llc = Lc*(thetac/(zetaa*(1-thetac)))^(-thetac);
lld = Ld*(thetad/(zetaa*(1-thetad)))^(-thetad);

lhc = thetac/(zetaa*(1-thetac))*llc;
lhd = thetad/(zetaa*(1-thetad))*lld;

% sector output
yc = 
% technology in next period
Acp = (1+vc)*Ac;
Adp = (1+vd)*Ad; 

end
