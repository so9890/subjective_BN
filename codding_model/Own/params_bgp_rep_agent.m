function [params, pols_num, model_pars]=params_bgp_rep_agent(symsparams, f, pol, indic)

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
zetaa    = 1.4;         % matches skill premium; with zeta==1 there is no 
                        % difference in skills from a hh perspective
eppsilon = 0.4;           % elasticity of substitution clean and dirty production
alphaa   = 1/3;         % income share capital
psii     = alphaa^2;    % cost of machine production following AA12
thetac   = 0.7;         % high skill labour share clean sector
thetad   = thetac*0.8;  % high skill labour share dirty sector
Uppsilon = 0.08;         % sum of growth rates; used as an upper bound
betaa    = 0.999;       % matches time preference in AA12 (rho= 0.001; betaa=exp(-rho))

% numeric vector of parameter values
params=eval(symsparams);

% policy variables
tauul   = indic.tauul_ex;        % progressivity; taken from HSV
% have to ensure that tauul>(1-alphaa)*(1-eppsilon) in the case of
% substitutes
if eppsilon>1
   if tauul<(1-alphaa)*(1-eppsilon)
       display ('*********the program is not stable for that range, with goods being substitutes******************')
   end
elseif eppsilon<1
    if tauul>(1-alphaa)*(1-eppsilon)
       display ('*********the program is not stable for that range, with goods being complements******************')
    end
end

lambdaa = 1;             % as if not there

if indic.het_growth==1       % heterogeneous growth, there should be no structural change!
    vd      = 0.07;          % growth dirty sector
    vc      = Uppsilon-vd;   % growth clean sector
    
elseif indic.het_growth==0  % sectors grow at a equal rate
    
    vd      = 0.07;         % growth clean sector
    vc      = vd; 
    
end
pols_num=eval(pol);

% substitute in model
model_pars=subs(f, [symsparams, pol], [params, pols_num]); % as a function of initial values Ad, Ac

end