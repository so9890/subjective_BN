function [ model, model_param, varsModel, paramsModel]=model_eq_Ramsey(Obj, symmsoptim, symmsparams, params, symmssolve, file_name, x_init, x, non_pol_num, non_pol_sym)

% function to differentiate the Objective functions 
% and to return a model file of a numeric function that can be solved using fsolve. 
% Note: that this specific model cannot be solved analytically. 

% input
% Obj:          objective function to be differentiated
% symmsoptim:   symbolic vector of variables wrt which Obj is to be
%               differentiated: lagrange on target and tauul
% symmsparams:  symbolic vector of parameters in model
% params:       numerical vector of parameter values, same order as
%               symmsparam
% listparams:   list of parameters as in symmsparams
% listoptim:    list of choice variables wrt which Obj is differentiated
% file_name:    name of file to which model should be written
% symmssolve:   symbolic vector for which to solve equilibrium

% output
% 

%% get model and write model file to be passed to fsolve

% differentiate wrt each variables separately to abvoid replacement and
% unreadable output

% initialise model by vector of eqbm variables
model=symmssolve'; 

for i=1:length(symmsoptim)
model(i)=jacobian(Obj, symmsoptim(i)); % should give derivative=0 if not present
end

% retrieve variables in model for substitution
model_vars=symvar(model)';
paramsModel=symmsparams(ismember(symmsparams,model_vars)); 

%-- replace parameters by values

model_param=subs(model',[symmsparams(ismember(symmsparams,model_vars)), x, non_pol_sym], [params(ismember(symmsparams,model_vars)), x_init, non_pol_num]);


%-- write parameterised model to file
 % second argument ensures same ordering of input args as in symmsoptim
 % the second output gives the vector or input arguments!
 
% the inputs are given by symmsoptim: the non-transformed variables!
matlabFunction(model_param, 'vars', {symmssolve(ismember(symmssolve,model_vars))}, 'File', file_name );

% get vector of variables in model to check if correct
varsModel=symvar(model_param);
end
