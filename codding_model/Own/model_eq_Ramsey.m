function [indexx, model, model_param_trans, varsModel, paramsModel]=model_eq_Ramsey(Obj, symmsoptim, symmsparams, params, symmssolve, file_name, x_init, x, non_pol_num, non_pol_sym, list)

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
% symmssolve:   variables for which to solve model
% file_name:    name of file to which model should be written
% x_init:       initial values of state variables
% x:            symbolic vector of state variables
% non_pol_num:  numeric vector of potential policy variables which are treated as
%               parameters in this model version
% non_pol_sym:  as above but symbolic
% list:         structure of variable lists, here: needed to have a list
%               for which to solve the system for transformation 

% output
% model:        symbolic model: jacobian of ramsey objective function 
% model_param_trans:  model with substituted numeric parameter values; and
%                     variables transformed to unbounded variables
% varsModel:    variables for which to solve Ramsey problem
% paramsModel:  parameters in model

%% get model and write model file to be passed to fsolve

% differentiate wrt each variables separately to abvoid replacement and
% unreadable output

% initialise model by vector of eqbm variables
model=transpose(symmssolve); 

for i=1:length(symmsoptim)
model(i)=jacobian(Obj, symmsoptim(i)); % should give derivative=0 if not present
end


% retrieve variables in model for substitution
model_vars=transpose(symvar(model));
paramsModel=symmsparams(ismember(symmsparams,model_vars)); 


%-- replace parameters by values

model_param=subs(transpose(model),[symmsparams(ismember(symmsparams,model_vars)), x, non_pol_sym], [params(ismember(symmsparams,model_vars)), x_init, non_pol_num]);

%-- transform remaining variables to undbounded ones

% get symbolic vector of transformed variables, as they should occur in
% model to be passed to numeric
% indexx: structure of indices for how to transform variables in ramsey model
% lagrange multis in ramsey are not transformed!

indexx.oneabove=(list.optim=='tauul');

%transformation: how variables should occur in model

trans_optim=symmssolve;
trans_optim(indexx.oneabove)=(1./(1+exp(symmssolve(indexx.oneabove))));

% replace untransformed with transformed variables

model_param_trans=subs(model_param, symmssolve, trans_optim);

%-- write parameterised model to file
 % second argument ensures same ordering of input args as in symmsoptim
 % the second output gives the vector or input arguments!
 
% the inputs are given by symmsoptim: the non-transformed variables!
matlabFunction(model_param_trans, 'vars', {symmssolve(ismember(symmssolve,model_vars))}, 'File', file_name );

% get vector of variables in model to check if correct
varsModel=sort(symvar(model_param_trans));
end
