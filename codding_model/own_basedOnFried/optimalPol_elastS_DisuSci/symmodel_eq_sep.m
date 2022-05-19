% function so derive problem symbolically by taking derivatives
% parameters already numeric
function [indexx, model]=symmodel_eq_sep(OB_RAM,  x, params,  Ftarget, file_name, list, indic, indexx)

% function to differentiate the Objective functions 
% and to return a model file of a numeric function that can be solved using fsolve. 
% Note: that this specific model cannot be solved analytically. 

% input
% OB_RAM:       objective function to be differentiated
% params:       numeric vector of parameter values, same order as
%               symmsparam
% file_name:    name of file to which model should be written
% x:            symbolic vector of model variables
% list:         structure of variable lists, here: needed to have a list
%               for which to solve the system for transformation 


% output
% model:        symbolic model: jacobian of ramsey objective function 
% indexx:       structure of indices now containing index of symbolic model
% symms:        structure of symbolic vectors now containing variables
%               ordered as occuring in model

%% get model and write model file to be passed to fsolve

% differentiate wrt each variables separately to abvoid replacement and
% unreadable output

% initialise model by vector of eqbm variables
model=x; 

for i=1:length(x)
model(i)=jacobian(OB_RAM, x(i)); % should give derivative=0 if not present
end

% replace derivative wrt kt multiplier with complementary slackness
% condition (economic time endowment)
% 
%  model(contains(list.optALL,'KT_hh'))=x(contains(list.optALL,'KT_hh')).*(params(list.params=='upbarH')-x(contains(list.optALL,'HH')));
% model(contains(list.optALL,'KT_hl'))=x(contains(list.optALL,'KT_hl')).*(params(list.params=='upbarH')-x(contains(list.optALL,'HL')));
%  model(contains(list.optALL,'KT_S'))=x(contains(list.optALL,'KT_S')).*(params(list.params=='upbarH')-x(startsWith(list.optALL,'S')));
% % 

if indic.target==1
    model(contains(list.optALL,'mu_target'))=x(contains(list.optALL,'mu_target')).*(Ftarget-x(contains(list.optALL,'F')));
end

% retrieve variables in model for substitution
model_vars=transpose(symvar(model));
if sum(~ismember(x,model_vars))>0
    error('there are variables in the choice set (i.e. in x) which are not in the model ')
end

%-- transform variables to undbounded ones

% get symbolic vector of transformed variables, 
% AS THEY SHOULD OCCUR IN THE MODEL FOR SUBSTITUTION
% lagrange multis in ramsey are not transformed!
upbarH=params(list.params=='upbarH');
out_trans=exp(x);
if indic.noskill==0
    out_trans(contains(list.optALL,'HL'))=upbarH./(1+exp(x(contains(list.optALL,'HL'))));
 out_trans(contains(list.optALL,'HH'))=upbarH./(1+exp(x(contains(list.optALL,'HH'))));
else
 out_trans(contains(list.optALL,'H'))=upbarH./(1+exp(x(contains(list.optALL,'H'))));
end

if params(list.params=='etaa')>=1
    if indic.target == 0
     out_trans(startsWith(list.optALL,'sff'))=upbarH./(1+exp(x(startsWith(list.optALL,'sff'))));
      out_trans(startsWith(list.optALL,'sn'))=upbarH./(1+exp(x(startsWith(list.optALL,'sn'))));
       out_trans(startsWith(list.optALL,'sg'))=upbarH./(1+exp(x(startsWith(list.optALL,'sg'))));
    else
        out_trans(startsWith(list.optALL,'sg')) = (x(startsWith(list.optALL,'sg'))).^2;
          out_trans(startsWith(list.optALL,'sn')) = (x(startsWith(list.optALL,'sn'))).^2;
            out_trans(startsWith(list.optALL,'sff')) = (x(startsWith(list.optALL,'sff'))).^2;
    end
else
   out_trans(startsWith(list.optALL,'sff'))=upbarH./(1+exp(x(startsWith(list.optALL,'sff'))));
   out_trans(startsWith(list.optALL,'sn'))=upbarH./(1+exp(x(startsWith(list.optALL,'sn'))));
   out_trans(startsWith(list.optALL,'sg'))=upbarH./(1+exp(x(startsWith(list.optALL,'sg'))));
end

out_trans(contains(list.optALL, 'mu')|contains(list.optALL, 'KT'))=x(contains(list.optALL, 'mu')|contains(list.optALL, 'KT'));
% if indic.taus==1
%     out_trans((find(list.opt=='sg')-1)*T+1:find(list.opt=='sg')*T) = (x((find(list.opt=='sg')-1)*T+1:find(list.opt=='sg')*T)).^2;
% end
if indic.target==1
    out_trans(startsWith(list.optALL,'F'))=Ftarget./(1+exp(x(startsWith(list.optALL,'F'))));
end
    % replace untransformed with transformed variables

model_trans=subs(model, x, out_trans);

%-- write model to file
 % the second output gives the vector or input arguments!
 
% the inputs are given by x: the non-transformed variables!
% if ~isfile(sprintf('%s.m', file_name))
    matlabFunction(model_trans, 'vars', {x(ismember(x,model_vars))}, 'File', file_name );
% end
end
