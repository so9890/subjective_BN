function [params, targets,  pol, list, symms, Ems ]=get_params(F0, T, indic, lengthh)

% function to read in parameter values and to numerically calculate initial 
% period values of endogenous variables.

% input
% symmparams:   model parameters
% f:            model in symbolic variables
% pol:          symbolic vector of policy variables
% x:            symbolic vector of states


% output
% params:       numeric vector of calibrated parameters
% pols_num:     numeric vector of policy
% x_init:       initial conditions
% vars_tosolve: ordered list of variables as they enter in model function 

%% symbolic vector and list
syms sigmaa...      % 1/sigmaa = Frisch elasticity of labour
     thetaa...      % courvature consumption utility
     betaa...       % discount factor 5 years
     zh ...         % share high skilled
     zl ...         % share low skilled
     barHl...       % time endowment low 
     barHh...       % time endowment high  
     alphaf ...     % machine share fossil
     alphan ...     % machine share neutral
     alphag ...     % machine share green
     thetaf ...     % high skill weight fossil sector
     thetan ...     % high skill weight neutral sector
     thetag ...     % high skill weight green sector
     eppsy ...      % elasticity final good production
     eppse ...      % elasticity energy production
     deltay ...     % weight on energy in final production
     gammaa ...     % quality research
     etaa ...       % returns to scale scientists
     rhof ...       % number tasks scientists fossil
     rhon ...       % number tasks scientists neutral
     rhog ...       % number tasks scientists green
     phii ...       % spill over curvature
     Ems ...        % vector of net emission targets
     deltaa ...     % regeneration rate nature
     omegaa ...     % emission share of dirty output
     S ...          % number scientists
     Af0 ...        % initial technology level fossil
     Ag0 ...        % initial technology level green
     An0 ...        % initial technology level neutral     
    real 
 
syms taul ...       % income tax progressivity
     taus ...       % subsidy on green research
     tauf ...       % sales tax fossil
     lambdaa ...    % scale tax scheme
     real
 
symms.params = [sigmaa, thetaa, betaa, zh, zl , barHl, barHh, alphaf, alphan, alphag,...
                thetaf, thetan, thetag, eppsy, eppse, deltay, ...
                gammaa, etaa, rhof, rhon, rhog, phii, S, Af0, An0, Ag0];   
list.params  = string(symms.params);

symms.targets = [deltaa, omegaa];
list.targets  = string(symms.targets);

symms.pol     = [taul, taus, tauf, lambdaa];
list.pol      = string(symms.pol);

%% Calibration 

sigmaa   = 1/0.75;      % from Chetty et al 

if indic.util== 0
    thetaa   = 1;
else
    thetaa   = 2; % look up in Boppart
end

betaa    = (.985)^5;  % Barrage, but here for 5 years
barHl    = 1;
barHh    = barHl; 
zl       = 0.4;
zh       = 1-zl; 

eppse    = 1.5;            % Fried
eppsy    = 0.05;           % Fried
deltay   = 1.44e-38;       % Fried
alphaf   = 1-0.28;         % Fried: fossil has a higher labour share!
alphag   = 1-0.09;         % Fried
alphan   = 1-0.64;         % Fried
 
thetan   = 0.5;
thetag   = 0.6;
thetaf   = thetag*0.5;

gammaa   = 3.96;           % Fried
etaa     = 0.79; 
rhof     = 0.01;
rhon     = 1; 
rhog     = 0.01;
phii     = 0.5;            % Fried
S        = 0.01; 

Af0     = 1; % match output levels 
Ag0     = 0.5;
An0     = 1; 

%- emissions
[deltaa, omegaa, Ems]= calibration_emissions(F0,T, lengthh); 

%- policies
taul    = 0.181;
taus    = 0; 
tauf    = 0; 
lambdaa = 1;

% saving to vector
params  = eval(symms.params);
targets = eval(symms.targets);
pol     = eval(symms.pol);
end