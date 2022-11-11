function [params, Sparams,  pol, init201014, init201519, list, symms, Ems,  Sall, x0LF, MOM , indexx]=get_params( T, indic, lengthh)

% function to read in parameter values and to calibrate direct parameters
% calls on calibration_matching and calibration_emissions 

% input
% symmparams:   model parameters
% f:            model in symbolic variables
% pol:          symbolic vector of policy variables
% x:            symbolic vector of states


% output
% params:       numeric vector of calibrated parameters and initial conditions
% pol:          numeric vector of policy
% targets:      numeric vector of emission targets
% init:         numeric vector of initial conditions
% Sall:         structure of all model variables in baseyear
% Ems:           numeric vector of emission targets
% x0LF:         numeric vector of LF solution in baseyear ordered as in list.choice (used for LF solution) 
% 
%% symbolic vector and list
syms sigmaa...      % 1/sigmaa = Frisch elasticity of labour
     thetaa...      % courvature consumption utility
     betaa...       % discount factor 5 years
     zh ...         % share high skilled
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
     sigmaas ...    % Frisch elasticity labour supply scientists
     chiis ...      % disutility of labour scientists
     upbS ...       % upper bound scientists
     Af0 ...        % initial technology level fossil
     Ag0 ...        % initial technology level green
     An0 ...        % initial technology level neutral 
     real 
%      phis ...       % cost of scientists fixed
   
 
syms taul ...       % income tax progressivity
     taus ...       % subsidy on green research
     tauf ...       % sales tax fossil
     lambdaa ...    % scale tax scheme
     real
 
symms.params = [sigmaa,sigmaas, chiis, thetaa, betaa, zh, chii, upbarH, alphaf, alphan, alphag,...
                thetaf, thetan, thetag, eppsy, eppse, deltay, upbS, ...
                gammaa, etaa, rhof, rhon, rhog, phii, deltaa, omegaa];   
list.params  = string(symms.params);

symms.init   = [Af0, An0, Ag0];
list.init    = string(symms.init);

symms.pol     = [taul, taus, tauf, lambdaa];
list.pol      = string(symms.pol);

% parameters directly calibrated
symms.paramsdir = [sigmaa, thetaa, betaa, upbarH, alphaf, alphan, alphag,...
                eppsy, eppse, sigmaas, upbS, ...
                etaa, phii,  rhof, rhon, rhog, deltaa];   
list.paramsdir  = string(symms.paramsdir);

symms.poldir     = [taul, taus, tauf];
list.poldir      = string(symms.poldir);
%% Calibration 

sigmaa   = 1/0.75;      % from Chetty et al 
sigmaas  = sigmaa; 
if indic.util== 0
    thetaa   = 1;
else
    thetaa   = 2; % look up in Boppart
end

betaa    = (.985)^5;  % Barrage, but here for 5 years
upbarH     = 1;

eppse    = 1.5;            % Fried
eppsy    = 0.05;           % Fried
alphaf   = 1-0.28;         % Fried: fossil has a higher labour share!
alphag   = 1-0.09;         % Fried
alphan   = 1-0.64;         % Fried
 
% gammaa   = 3.96;           % Fried
if indic.spillovers==0
    etaa     = 0.79; 
else
    etaa     = 1.2; % positive spillovers=> to accomodate zero scientists in competitive eqbm 
end
rhof     = 0.01;
rhon     = 1; 
rhog     = 0.01;
% rhof     = 0.01;
% rhon     = 1; 
% rhog     = 0.01;
phii     = 0.5;            % Fried
upbS     = 0.01;

% phis = 1; % scaled
%- policies
taul    = 0.181;
taus    = 0; 
tauf    = 0; 

%% - indirect calibration 
%-- get moments
MOM = calibration_moments();
% MOM.S = 0.01; % from fried: Supply scientists in base year
MOM.growth = (1.02)^5 -1; %5 year grwoth rate

%% - emissions
[deltaa, Ems, MOM]= calibration_emissions(T, lengthh, MOM); 
% -omegaa follows in main calibration
%% save directly calibrated variables
parsHelp = eval(symms.paramsdir);
polhelp= eval(symms.poldir);
%%
[x0LF, ~, ~, ~, Sall, ~,  init201014 , ~, init201519, Sparams, ~, params, pol, symms, MOM,indexx, list]...
    = calibration_matching(MOM, symms, list, parsHelp, polhelp);


% %% - parameters for no-skill version
% Sparams_noskill=Sparams;
% Sparams_noskill.thetaf=0.5;
% Sparams_noskill.thetag=0.5;
% Sparams_noskill.thetan=0.5;
% Sparams_noskill.zh=0.5;
% 
% params_noskill = params;
% params_noskill(list.params=='thetan')=0.5;
% params_noskill(list.params=='thetaf')=0.5;
% params_noskill(list.params=='thetag')=0.5;
% params_noskill(list.params=='zh')=0.5;


end