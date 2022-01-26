function [params, analy_solution]=params_bgp_rep_agent(symsparams, f, pol, xp, y, x )

% function to read in parameter values and balanced growth path
% retrieves solution to analytically solvable model as symbolic vector
% which are then passed as initial values for numeric solution to first
% period solution. 

% input
% symmparams: model parameters
% f: model in absolute variables
% pol: symbolic vector of policy variables
% y, xp, x: controls, predetermined variables, and states

% output
% params:         numeric vector of calibrated parameters
% analy_solution: solution of choice and predetermined variables as a
%                 function of parameters, initial conditions, and policy

%% symbolic solution: choice variables as a function of Ad, Ac, policy and parameters
% uses analytically derived equations

analy_solution=solution_SS(x, symsparams, pol, list, [y, xp]);

%% Calibration 
sigmaa   = 1/0.75;      % from Chetty et al 
zetaa    = 1.4;         % matches skill premium; with zeta==1 there is no 
                        % difference in skills from a hh perspective
eppsilon = 2;           % elasticity of substitution clean and dirty production
alphaa   = 1/3;         % income share capital
psii     = alphaa^2;    % cost of machine production following AA12
thetac   = 0.6;         % high skill labour share clean sector
thetad   = thetac*0.7;  % high skill labour share dirty sector
Uppsilon = 0.1;         % sum of growth rates; used as an upper bound
betaa    = 0.999;       % matches time preference in AA12 (rho= 0.001; betaa=exp(-rho))

% numeric vector of parameter values
params=eval(symsparams);

% policy variables
tauul   = 0.181;        % progressivity; taken from HSV
lambdaa = 1;            % as if not there
vc      = 0.01;         % growth clean sector
vd      = Uppsilon-vc;  % growth dirty sector

pols_num=eval(pol);

%% numeric solution 
% uses analytical results as initial values despite modelling differences

% initial conditions
Ad=8;
Ac=4;
x_init=eval(x);

% Test if analytic solution to SS sets model equations numerically to zero
model_num=subs(f, [symsparams, pol, x], [params, pols_num, x_init]); % as a function of initial values Ad, Ac
num_solution=eval(subs(analy_solution, [symsparams, pol, x], [params, pols_num, x_init]));

modelSS= eval(subs(model_num, [y, xp], num_solution)); % important to have y and x in order as passes to solution_SS!

if max(modelSS)>1.e-14
    error('SS not found')
end

%%
% write model as function 
matlabFunction(model_num, 'vars', {[y, xp]}, 'File', 'model_fun' ); % vars: variables for which to solve model
                                                                    % File: filename of function 
% get variables for which to solve model                                                                    
vars_tosolve=symvar(model_num);


% initial values: use analytically derived solution as starting values
Acp=(1+vc)*Ac_init;
Adp=(1+vd)*Ad_init; 
Lc=
x0=[]
modFF = @(x)model_num(x);
options = optimoptions('fsolve', 'MaxFunEvals',8e5, 'MaxIter', 3e5, 'Algorithm', 'levenberg-marquardt');%, 'TolFun', 10e-8);%, );%, 'Display', 'Iter', );

[y, fval, exitflag]=fsolve(modFF, guess_trans, options);
end