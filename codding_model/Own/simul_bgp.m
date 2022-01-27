function [ybgp, xpbgp]= simul_bgp(hours, list, x, x_init, params, pols_num, model_pars, symsparams, pol)

% function to calculate 
% endogenous variables given calibration, policies, and initial conditions.

% input
% list:         structure of lists of model parameters and variables
% hours:        starting value for H, hh, hl to solve for initial guess

% output
% ybgp:          numeric controls in current iteration
% xpbgp:         numeric predetermined variables in current iteration

%% write new function for each new initial condition
% substitute numeric values
model_num=subs(model_pars, x, x_init); % as a function of initial values Ad, Ac

% variables for which to solve model numerically
vars_tosolve=symvar(model_num);

% write model as function 
matlabFunction(model_num, 'vars', {vars_tosolve}, 'File', 'model_tosolve' );


%% solve function using fsolve
%   pass numeric versions of initial values at beginning of period,
%   parameters, policy to
%   solution of model as function of hours supplied
%   results ordered as in vars_tosolve (= order of variables in model_tosolve)

guess=solution_SS(x_init, params, pols_num, list, vars_tosolve, hours);

modFF = @(x)model_tosolve(x);
options = optimoptions('fsolve', 'MaxFunEvals',8e5, 'MaxIter', 3e5, 'TolFun', 10e-14);%, 'Algorithm', 'levenberg-marquardt');%, );%, );%, 'Display', 'Iter', );

[vars_solution, fval, exitflag]=fsolve(modFF, guess, options);

%% save numeric vectors of y and xp

ybgp=zeros(size(list.y));
xpbgp=zeros(size(list.xp));

for i =list.y
   ybgp(list.y==i) =vars_solution(vars_tosolve==i);
end

for i = list.xp
    xpbgp(list.xp==i)= vars_solution(vars_tosolve==i);
end
