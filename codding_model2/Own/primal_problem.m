function [symms, list, Obj_ramPA ]=primal_problem(y, x, list, symms, E)
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
Hbar=symms.params(list.params=='Hbar');

%gammaa=symms.params(list.params=='gammaa');
%etaa= symms.params(list.params=='etaa');

% targets
deltaa =symms.targets(list.targets=='deltaa');
kappaa =symms.targets(list.targets=='kappaa');

Muc =symms.marginals(list.marginals=='Muc');
Muhh =symms.marginals(list.marginals=='Muhh');
%Muhl =symms.marginals(list.marginals=='Muhl');

% lagrange multiplier as syms
syms mu_target mu_opt_final mu_rc mu_defH kt_lab real % exogenous emission target

% symbolic vector of variables over which to optimise
symms.optimPA=[c, llc, lld, lhc, lhd, yc, yd, H, mu_target, mu_opt_final, mu_rc, mu_defH, kt_lab ];
symms.ramsey_mu=[mu_target, mu_opt_final, mu_rc, mu_defH, kt_lab];
list.ramsey_mu=string(symms.ramsey_mu);

% replace list optim for primal approach
list.optim = string(symms.optimPA);

%% change the following to be more flexible in terms of utility
% auxiliary variables
hl = llc+lld;       % market clearing low skill
hh = lhc+lhd;       % market clearing high skill
%;   % definition H

Ld = lld^thetad*lhd^(1-thetad); % labour production dirty
Lc = llc^thetac*lhc^(1-thetac); % labour production clean

%- prices
pc = yc/(Ac*Lc)^((1-alphaa)/alphaa)*(psii/alphaa); % clean output
pd = yd/(Ad*Ld)^((1-alphaa)/alphaa)*(psii/alphaa); % dirty output

pcL = (1-alphaa)*(alphaa/psii)^(alphaa/(1-alphaa))*pc^(1/(1-alphaa))*Ac; % clean labour demand
pdL = (1-alphaa)*(alphaa/psii)^(alphaa/(1-alphaa))*pd^(1/(1-alphaa))*Ad; % dirty labour demand

wl = (llc/lhc)^(1-thetac)*thetac*pcL; % FOCs dirty labour production (1 eq.)
wh = (lhc/llc)^(thetac)*(1-thetac)*pcL; % FOCs clean labour production (2 eq.)

xd = (alphaa/psii*pd)^(1/(1-alphaa))*Ad*Ld; % demand machines dirty
xc = (alphaa/psii*pc)^(1/(1-alphaa))*Ac*Lc; % demand machines clean

Y = (yc^((eppsilon-1)/eppsilon)+yd^((eppsilon-1)/eppsilon))^(eppsilon/(eppsilon-1)); % final good production
%- policy

tauul = 1-H^(1+sigmaa);     % Foc H => follows from budget and foc consumption, foc hl 
                            % => equivalent to IMP
                            % only remains to ensure wages satisfy foc
                            % labour supply
lambdaa = (H*wl-G)/((H*wl)^(1-tauul)); % gov budget => already determined!

%- Gov problem
U = log(c)-(hl+zetaa*hh)^(1+sigmaa)/(1+sigmaa);
W = U;

% constraints
imp         = c*Muc+Muhh/zetaa*H^(-sigmaa);
rc          = Y-(c+psii*(xd+xc)+G); % superfluous as gov budget is determined and hh budget ensured to hold
opt_fin     = yd-(pc/pd)^eppsilon*yc;
targets     = yd-(deltaa+E)/kappaa; 
defH        = H  -( hh*zetaa+hl);
free_lab    = wh/wl-zetaa;
% Kuhn tucker on labour supply
lab_con     = Hbar -(hl+hh);

% objective function 
Obj_ramPA = W -mu_target*targets...
            -mu_opt_final*opt_fin...  
            -mu_free_lab*free_lab... 
            -mu_defH*defH...
            -mu_rc*rc...
            -kt_lab*lab_con; 
           % -mu_budget*budget...% holds by walras law
           % -mu_imp*imp...    % holds by definition of taul     


% dynamic: generate vectors of variables for 30 periods


  

