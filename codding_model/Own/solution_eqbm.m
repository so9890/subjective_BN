function levels= solution_eqbm(x, symsparams, pols, list, vars_tosolve)

% function uses analytically derived equilibrium to use in initial period
% as guess
% THIS VERSION: with no disposal of government revenues

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
gammaa=symsparams(list.params=='gammaa');
etaa= symsparams(list.params=='etaa');
G= symsparams(list.params=='G');


tauul=pols(list.pol=='tauul');
%lambdaa=pols(list.pol=='lambdaa');
vc=pols(list.pol=='vc');
vd=pols(list.pol=='vd');

%% read in competitive eqbm solution
solution_eqbm_syms;

%% summarise
levels=eval(vars_tosolve);
end
