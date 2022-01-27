function [params, pols, analy_solution]=params_bgp_rep_agent(symsparams, f, pol, xp, y, x )

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

% some mistake; check later
%analy_solution=solution_SS(x, symsparams, pol, list, [y, xp]);

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
% uses analytical functions as initial values (despite modelling
% differences)

%-- evaluate model at initial conditions, parameters, and policy
%   instruments. Then, extract variables for which to solve model and
%   write model as a function. 

% initial conditions
Ad=8;
Ac=4;
x_init=eval(x);

% substitute numeric values
model_num=subs(f, [symsparams, pol, x], [params, pols_num, x_init]); % as a function of initial values Ad, Ac

% variables for which to solve model numerically
vars_tosolve=symvar(model_num);

% write model as function 
matlabFunction(model_num, 'vars', {vars_tosolve}, 'File', 'model_tosolve' );

%-- initial values: use analytically derived solution as starting values
%   and guesses for hours supplied

H=(1-tauul)^(1/(1+sigmaa));
% guess distribution of H on hl and hh
hh=H/4; % one fourth in high skill hours
hl=H-zetaa*hh;

% vector of guesses and list
init=[H, hh, hl];
list.init=["H", "hh", "hl"];

%-- solve function using fsolve
%   pass numeric versions of initial values, parameters, policy to
%   analytic (half) solution
%   results ordered as in vars_tosolve (= order of variables in model_tosolve)

guess=solution_SS(x_init, params, pols_num, list, vars_tosolve, init);

modFF = @(x)model_tosolve(x);
options = optimoptions('fsolve', 'MaxFunEvals',8e5, 'MaxIter', 3e5, 'TolFun', 10e-15);%, 'Algorithm', 'levenberg-marquardt');%, );%, );%, 'Display', 'Iter', );

[vars_solution, fval, exitflag]=fsolve(modFF, guess, options);

%% evaluate vectors of y and xp
endo_vars_num=eval(subs([y, xp], vars_tosolve, vars_solution));
y=endo_vars_num(1:21);
x=endo_vars_num(22:23);
end