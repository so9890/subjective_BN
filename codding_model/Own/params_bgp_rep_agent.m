function [params, Xss, Yss]=params_bgp_rep_agent(symsparams, f )

%function to read in parameter values and balanced growth path

% input
% symmparams: model parameters
% f: model in absolute variables

%% parameters 
sigmaa   = 1/0.75;      % from Chetty et al 
zetaa    = 1.4;         % matches skill premium    
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
vd     = Uppsilon-vc;  % growth dirty sector

pols=eval(pol); 
% analytically derived variable values on BGP
% write in terms of growth rates, then simulate by determining initial
% variables

%- use symbolic toolbox to solve
%S = solve(f, [y xp]); % f= model, params, Ad, Ac taken as given 

% use ana and solve model as function of Ad and Ac (initial conditions)

cell.y=arrayfun(@char, list.y, 'uniform', 0);
cell.x=arrayfun(@char, list.x, 'uniform', 0);
cell.yp=arrayfun(@char, list.yp, 'uniform', 0);
cell.xp=arrayfun(@char, list.xp, 'uniform', 0);

%modelSS=subs(model.f, [cell.x, cell.xp, cell.yp, cell.y], [Xss', Xss', Yss', Yss']);

model=subs(f, [symsparams, pol], [params, pols]);
symvar(model)

S = solve(model, [y xp]); 
end