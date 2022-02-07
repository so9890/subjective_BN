function [params, pols_num, model_pars, E]=params_bgp_rep_agent(symsparams, f, pol, indic, T, tauul_opt)

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

%% symbolic solution: choice variables as a function of Ad, Ac, policy and parameters
% uses analytically derived equations

% some mistake; check later
%analy_solution=solution_SS(x, symsparams, pol, list, [y, xp]);

%% Calibration 
sigmaa   = 1/0.75;      % from Chetty et al 

if indic.util== 0
    gammaa   = 1;
    etaa     = 0;       % not relevant
else
    gammaa   = 2;
    etaa     = 4;       % exponent on leisure
end

if isfield(indic, 'zetaa')
    zetaa= indic.zetaa;
else
    zetaa    = 1.4;         % matches skill premium; with zeta==1 there is no 
                            % difference in skills from a hh perspective
end
                            
eppsilon = 0.4;         % elasticity of substitution clean and dirty production
alphaa   = 1/3;         % income share capital
psii     = alphaa^2;    % cost of machine production following AA12
thetac   = 0.7;         % high skill labour share clean sector
thetad   = thetac*0.8;  % high skill labour share dirty sector
Uppsilon = 0.08;        % sum of growth rates; used as an upper bound
betaa    = 0.999;       % matches time preference in AA12 (rho= 0.001; betaa=exp(-rho))
deltaa   = 20;           % regeneration rate of nature
kappaa   = 0.04;        % share of dirty output which translates into emissions 

E        = [repmat(30,1, 50-20-1), zeros(1,T-30+1)]; % from t=1 to t=29 E=30; from t=30 to t=T E=0

% policy variables

if isnan(tauul_opt) 
if ~isfield(indic,'tauul_ex')
    tauul=0.181;                % progressivity; taken from HSV
else
    tauul   = indic.tauul_ex;        
end

elseif ~isnan(tauul_opt)
    tauul = tauul_opt;
end

lambdaa = 1;             % as if not there

if indic.het_growth==1       % heterogeneous growth, there should be no structural change!
    vd      = 0.07;          % growth dirty sector
    vc      = Uppsilon-vd;   % growth clean sector
    
elseif indic.het_growth==0  % sectors grow at a equal rate
    
    vd      = 0.07;         % growth clean sector
    vc      = vd; 
    
end
% numeric vector of parameter values
params=eval(symsparams);

pols_num=eval(pol);

% substitute in model
model_pars=subs(f, [symsparams, pol], [params, pols_num]); % as a function of initial values Ad, Ac

end