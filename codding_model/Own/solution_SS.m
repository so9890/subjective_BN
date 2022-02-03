function levels= solution_SS(x, symsparams, pols, list, vars_tosolve, init, indic)

% function uses analytically derived BGP/initial SS equations to get model
% variables determined by initial values and parameters, and policy

% input:
% Ac, Ad:       initial technology (can be numeric or symbolic)
% symsparams:   numeric/symbolic vector of parameters
% pols:         numeric/symbolic vector of policy variables
% vars_tosolve: symbolic vector of controls and predetermined variables
%               determines order of solutions in output (levels).    
% init:         initial values for skill supply

% output:
% levels:  vector of model variables in first period given initial
%          conditions. Ordered as order in vars_tosolve
%% preparations 

Ac=x(list.x=='Ac');
Ad=x(list.x=='Ad');

% read in required symbolic variables/parameters from input
% if instead of symsparams a numeric vector is provided 
% these are numbers

thetac=symsparams(list.params=='thetac');
thetad=symsparams(list.params=='thetad');
sigmaa=symsparams(list.params=='sigmaa');
zetaa=symsparams(list.params=='zetaa');
eppsilon=symsparams(list.params=='eppsilon');
alphaa=symsparams(list.params=='alphaa');
psii=symsparams(list.params=='psii');

tauul=pols(list.pol=='tauul');
lambdaa=pols(list.pol=='lambdaa');
vc=pols(list.pol=='vc');
vd=pols(list.pol=='vd');

% read in labour variables (guesses)
H  = init(list.hours=='H');
hl = init(list.hours=='hl');
hh = init(list.hours=='hh');

% auxiliary variables/parameters
% chic = (thetac/(zetaa*(1-thetac)))^(thetac)*...
%         ((1-thetac)*(thetad-zetaa*(1-thetad)))...
%         /(thetad*(1-thetac)-thetac*(1-thetad));
% chid = (thetad/(zetaa*(1-thetad)))^(thetad)*(1-((1-thetac)*(thetad-zetaa*(1-thetad)))/(thetad*(1-thetac)-thetac*(1-thetad)));

gammad = (thetad/(zetaa*(1-thetad)))^thetad;
gammac = (thetac/(zetaa*(1-thetac)))^thetac;
chii   = (1-thetad)*(1-thetac)/(thetac*(1-thetad)-thetad*(1-thetac));

labc   = H+thetad/(1-thetad)*hl;
labd   = 1/(1-thetac)*hl-H;

%% solutions

% labour input good
Lc= gammac*chii*labc;
Ld= gammad*chii*labd;

% sector good prices
pd = ((gammad/gammac*Ad/Ac*labd/labc)^((1-alphaa)*(1-eppsilon)/(alphaa+eppsilon*(1-alphaa)))+1)^(-1/(1-eppsilon));
pc = pd*(gammad/gammac*Ad/Ac*labd/labc)^((1-alphaa)/(alphaa+eppsilon*(1-alphaa)));

% labour input prices
pcL = (1-alphaa)*(alphaa/psii)^(alphaa/(1-alphaa))*pc^(1/(1-alphaa))*Ac;
pdL = (1-alphaa)*(alphaa/psii)^(alphaa/(1-alphaa))*pd^(1/(1-alphaa))*Ad;

% skill specific wage
wl = pdL*zetaa^(-thetad)*thetad^thetad*(1-thetad)^(1-thetad);
wh = zetaa*wl;

% household
%H = (1-tauul)^(1/(1+sigmaa)); 
c = lambdaa*(H*wl)^(1-tauul);

% skill inputs
llc = Lc/gammac;
lld = Ld/gammad;

lhc = gammac^(1/thetac)*llc;
lhd = gammad^(1/thetad)*lld;

% sector output
yc = (alphaa/psii)^(alphaa/(1-alphaa))*Ac*Lc;

yd= (pc/pd)^eppsilon*yc;

% technology in next period
Acp = (1+vc)*Ac;
Adp = (1+vd)*Ad; 

% machines
xd = (alphaa/psii*pd)^(1/(1-alphaa))*Ad*Ld;
xc = (alphaa/psii*pc)^(1/(1-alphaa))*Ac*Lc;

% gov budget
G = (H*wl-lambdaa*(H*wl)^(1-tauul));

% goods market clearing
if indic.fullDisposal==1
    Y = c+psii*(xd+xc);
else
    Y = c+psii*(xd+xc)+G;
end

%% summarise
levels=eval(vars_tosolve);
end
