function [params, Xss, Yss]=params_bgp_rep_agent(symsparams, f )

%function to read in parameter values and balanced growth path

% input
% symmparams: model parameters
% f: model in absolute variables

%% parameters 
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

% numerical vector of parameter values
params=eval(symsparams);

%% BGP
% policy variables
tauul   = 0.181;        % progressivity; taken from HSV
lambdaa = 1;            % as if not there
vc      = 0.01;         % growth clean sector
vd      = Uppsilon-vc;  % growth dirty sector

pols=eval(pol); 
% analytically derived variable values on BGP
% write in terms of growth rates, then simulate by determining initial
% variables

% use ana and solve model as function of Ad and Ac (initial conditions)

model=subs(f, [symsparams, pol], [params, pols]); % can be solved numerically,
                                                  % as a function of initial values Ad, Ac
%symvar(model)
%S = solve(model, [y xp]); 

% numerically solve for SS
Ad_init=8;
Ac_init=4;
model_num= subs(model, [Ad Ac] , [Ad_init, Ac_init]);

% write model as function 
matlabFunction(model_num, 'vars', {[y, xp]}, 'File', 'model_fun' ); % vars: variables for which to solve model
                                                                    % File: filename of function 
% get variables for which to solve model                                                                    
vars_tosolve=symvar(model_num);

% initial values: use analytically derived
Acp=(1+vc)*Ac_init;
Adp=(1+vd)*Ad_init; 
Lc=
x0=[]
modFF = @(x)model_num(x);
options = optimoptions('fsolve', 'MaxFunEvals',8e5, 'MaxIter', 3e5, 'Algorithm', 'levenberg-marquardt');%, 'TolFun', 10e-8);%, );%, 'Display', 'Iter', );

[y, fval, exitflag]=fsolve(modFF, guess_trans, options);
end