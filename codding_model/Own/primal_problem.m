%%CONTINTINUE .W
%% Primal approach
% goal: find optimal allocation (c, llc, lld, lhc, lhd, yc, yd )
% Ad, Ac as given 

% variables over which to optimise

% lagrange multiplier as syms
syms mu_target mu_budget mu_opt_final mu_rc mu_imp real % exogenous emission target

% symbolic vector of variables over which to optimise
symms.optimPA=[c, llc, lld, lhc, lhd, yc, yd, mu_target, mu_budget, mu_opt_final, mu_rc, mu_imp ];

% auxiliary variables
hl = llc+lld;       % market clearing low skill
hh = lhc+lhd;       % market clearing high skill
H  = hh*zetaa+hl;   % definition H

Ld = lld^thetad*lhd^(1-thetad); % labour production dirty
Lc = llc^thetac*lhc^(1-thetac); % labour production clean
xd = (alphaa/psii*pd)^(1/(1-alphaa))*Ad*Ld; % demand machines dirty
xc = (alphaa/psii*pc)^(1/(1-alphaa))*Ac*Lc; % demand machines clean

Y = (yc^((eppsilon-1)/eppsilon)+yd^((eppsilon-1)/eppsilon))^(eppsilon/(eppsilon-1)); % final good production

%- prices
pc = yc/(Ac*Lc)^((1-alphaa)/alphaa)*(psii/alphaa); % clean output
pd = yd/(Ad*Ld)^((1-alphaa)/alphaa)*(psii/alphaa); % dirty output

pcL = (1-alphaa)*(alphaa/psii)^(alphaa/(1-alphaa))*pc^(1/(1-alphaa))*Ac; % clean labour demand
pdL = (1-alphaa)*(alphaa/psii)^(alphaa/(1-alphaa))*pd^(1/(1-alphaa))*Ad; % dirty labour demand

wl = pdL*zetaa^(-thetad)*thetad^thetad*(1-thetad)^(1-thetad); % FOCs dirty labour production (2 eq.)
wh = pcL*zetaa^(1-thetac)*thetac^thetac*(1-thetac)^(1-thetac); % FOCs clean labour production (2 eq.)

%- policy
tauul = 1-H^(1+sigmaa);              % Foc H
lambdaa = (H*wl-G)/(H*wl)^(1-tauul); % gov budget

%- Gov problem
U = log(c)-(hl+zetaa*hh)^(1+sigmaa)/(1+sigmaa);
W = U;

% constraints
imp = c*Muc+Muhh/zetaa*H^(-sigmaa);
rc = Y-(c+psii*(xd+xc)+G);
opt_fin = yd-(pc/pd)^eppsilon*yc;
targets = yd-(deltaa+E)/kappaa; 
budget = 0-G;

% objective function 
Obj_ram = W -mu_target*targets...
            -mu_budget*budget...
            -mu_opt_final*opt_fin...
            -mu_rc*rc...
            -mu_imp*imp; 

% can now take derivatives if static, or do value function iteration!



