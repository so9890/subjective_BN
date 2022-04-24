function [params, targets,  pol, list, symms, Ems ]=get_params( T, indic, lengthh)

% function to read in parameter values and to calibrate direct parameters
% calls on calibration_matching and calibration_emissions 

% input
% symmparams:   model parameters
% f:            model in symbolic variables
% pol:          symbolic vector of policy variables
% x:            symbolic vector of states


% output
% params:       numeric vector of calibrated parameters and initial conditions
% pols_num:     numeric vector of policy
% vars_tosolve: ordered list of variables as they enter in model function 

%% symbolic vector and list
syms sigmaa...      % 1/sigmaa = Frisch elasticity of labour
     thetaa...      % courvature consumption utility
     betaa...       % discount factor 5 years
     zh ...         % share high skilled
     zl ...         % share low skilled
     eh ...         % effective labour productivity high skill
     el ...         % effective productivity low skill
     chii ...       % disutility labour
     upbarH...        % time endowment 
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
 
symms.params = [sigmaa, thetaa, betaa, zh, zl, el, eh, chii, upbarH, alphaf, alphan, alphag,...
                thetaf, thetan, thetag, eppsy, eppse, deltay, ...
                gammaa, etaa, rhof, rhon, rhog, phii, S, Af0, An0, Ag0];   
list.params  = string(symms.params);

symms.targets = [deltaa, omegaa];
list.targets  = string(symms.targets);

symms.pol     = [taul, taus, tauf, lambdaa];
list.pol      = string(symms.pol);

% parameters directly calibrated
symms.paramsdir = [sigmaa, thetaa, betaa, zh, zl , upbarH, alphaf, alphan, alphag,...
                eppsy, eppse, S, ...
                gammaa, etaa, rhof, rhon, rhog, phii];   
list.paramsdir  = string(symms.paramsdir);

symms.poldir     = [taul, taus, tauf];
list.poldir      = string(symms.poldir);
%% Calibration 

sigmaa   = 1/0.75;      % from Chetty et al 

if indic.util== 0
    thetaa   = 1;
else
    thetaa   = 2; % look up in Boppart
end

betaa    = (.985)^5;  % Barrage, but here for 5 years
upbarH     = 1;
zh       = 0.3169;       % Slavik paper! to be updated
zl       = 1-zh; 

eppse    = 1.5;            % Fried
eppsy    = 0.05;           % Fried
%deltay   = 1.44e-38;       % Fried
alphaf   = 1-0.28;         % Fried: fossil has a higher labour share!
alphag   = 1-0.09;         % Fried
alphan   = 1-0.64;         % Fried
 
gammaa   = 3.96;           % Fried
etaa     = 0.79; 
rhof     = 0.01;
rhon     = 1; 
rhog     = 0.01;
phii     = 0.5;            % Fried
S        = 0.01; 

%- policies
taul    = 0.181;
taus    = 0; 
tauf    = 0; 

%% - indirect calibration 
%-- get moments
MOM = calibration_moments();
MOM.lowskill = 0.4;
MOM.AgAn =1.5; 
% thetan   = 0.5;
% thetag   = 0.6;
% thetaf   = thetag*0.5;
% Af0     = 1.877; % Fried 
% Ag0     = 0.9196;
% An0     = 1; 
% lambdaa

%% - emissions
[deltaa, Ems, MOM]= calibration_emissions(T, lengthh, MOM); 
% -omegaa follows in main calibration
%% -others
parsHelp = eval(symms.paramsdir);
polhelp= eval(symms.poldir);
targetsHelp = eval(symms.targets(list.targets~='omegaa'));
%%
[An0, Af0, Ag0, thetaf, thetan, thetag, ...
    lambdaa, omegaa]= calibration_matching(MOM, symms, list, parsHelp, polhelp, targetsHelp);


% save
params  = eval(symms.params);

end