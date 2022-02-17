function [symms, Obj_ramPA ]=primal_problem(y, x, list, symms, E)
%% Primal approach
% goal: find optimal allocation (c, llc, lld, lhc, lhd, yc, yd )
% Ad, Ac as symbols

% read in variables  over which to optimise and paramters
c=y(list.y=='c');
llc=y(list.y=='llc');
lld=y(list.y=='lld');
lhc=y(list.y=='lhc');
lhd=y(list.y=='lhd');
yc=y(list.y=='yc');
yd=y(list.y=='yd');
H=y(list.y=='H');

Ac=x(list.x=='Ac');
Ad=x(list.x=='Ad');

% read in required symbolic variables/parameters from input
% if instead of symms.params a numeric vector is provided 
% these are numbers

thetac=symms.params(list.params=='thetac');
thetad=symms.params(list.params=='thetad');
sigmaa=symms.params(list.params=='sigmaa');
zetaa=symms.params(list.params=='zetaa');
eppsilon=symms.params(list.params=='eppsilon');
alphaa=symms.params(list.params=='alphaa');
psii=symms.params(list.params=='psii');
G=symms.params(list.params=='G');
betaa=symms.params(list.params=='betaa');

%gammaa=symms.params(list.params=='gammaa');
%etaa= symms.params(list.params=='etaa');

% targets
deltaa =symms.targets(list.targets=='deltaa');
kappaa =symms.targets(list.targets=='kappaa');

Muc =symms.marginals(list.marginals=='Muc');
Muhh =symms.marginals(list.marginals=='Muhh');
%Muhl =symms.marginals(list.marginals=='Muhl');

% lagrange multiplier as syms
syms mu_target mu_budget mu_opt_final mu_rc mu_imp mu_defH real % exogenous emission target

% symbolic vector of variables over which to optimise
symms.optimPA=[c, llc, lld, lhc, lhd, yc, yd, H, mu_target, mu_budget, mu_opt_final, mu_rc, mu_imp, mu_defH ];
symms.ramsey_mu=[mu_target, mu_budget, mu_opt_final, mu_rc, mu_imp, mu_defH];
list.ramsey_mu=string(symms.ramsey_mu);

% auxiliary variables
hl = llc+lld;       % market clearing low skill
hh = lhc+lhd;       % market clearing high skill
%;   % definition H

Ld = lld^thetad*lhd^(1-thetad); % labour production dirty
Lc = llc^thetac*lhc^(1-thetac); % labour production clean

Y = (yc^((eppsilon-1)/eppsilon)+yd^((eppsilon-1)/eppsilon))^(eppsilon/(eppsilon-1)); % final good production

%- prices
pc = yc/(Ac*Lc)^((1-alphaa)/alphaa)*(psii/alphaa); % clean output
pd = yd/(Ad*Ld)^((1-alphaa)/alphaa)*(psii/alphaa); % dirty output

pcL = (1-alphaa)*(alphaa/psii)^(alphaa/(1-alphaa))*pc^(1/(1-alphaa))*Ac; % clean labour demand
pdL = (1-alphaa)*(alphaa/psii)^(alphaa/(1-alphaa))*pd^(1/(1-alphaa))*Ad; % dirty labour demand

wl = pdL*zetaa^(-thetad)*thetad^thetad*(1-thetad)^(1-thetad); % FOCs dirty labour production (2 eq.)
wh = pcL*zetaa^(1-thetac)*thetac^thetac*(1-thetac)^(1-thetac); % FOCs clean labour production (2 eq.)

xd = (alphaa/psii*pd)^(1/(1-alphaa))*Ad*Ld; % demand machines dirty
xc = (alphaa/psii*pc)^(1/(1-alphaa))*Ac*Lc; % demand machines clean

%- policy

tauul = 1-H^(1+sigmaa);              % Foc H
lambdaa = (H*wl-G)/(H*wl)^(1-tauul); % gov budget => already determined!

%- Gov problem
U = log(c)-(hl+zetaa*hh)^(1+sigmaa)/(1+sigmaa);
W = U;

% constraints
imp = c*Muc+Muhh/zetaa*H^(-sigmaa);
rc = Y-(c+psii*(xd+xc)+G);
opt_fin = yd-(pc/pd)^eppsilon*yc;
targets = yd-(deltaa+E)/kappaa; 
defH = H  -( hh*zetaa+hl);

% objective function 
Obj_ramPA = W -mu_target*targets...
            -mu_opt_final*opt_fin...
            -mu_rc*rc...
            -mu_imp*imp...
            -mu_defH*defH; 
           % -mu_budget*budget...


% dynamic: generate vectors of variables for 30 periods


  

