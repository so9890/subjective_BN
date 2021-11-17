%% Model file of the simple, static model, 
% eppsilons =1 => no low skill in sustainable sector
% no innovation 
% no monopolistic competition
% no labour sector choice

% provide all variables as a function of ps, pn
clear, clc

cd('/home/sonja/Documents/projects/Overconsumption/codding_model')
% variables
syms ps pn real % all other variables are functions of ps, pn

% parameters
syms eppsilons eppsilonn omegas alphaa psii B lambdaa As An real

% vector of variables for which to solve model
symm.variables=[ps]; %(only as function of pn)
symm.params= [eppsilons, eppsilonn, omegas, alphaa, psii, B, lambdaa, As, An];
symm.params_noncomp=[eppsilons, eppsilonn, alphaa, psii, lambdaa, As, An];

%% 
% -- auxiliary equations

% write down already solved equations, that is variables as function of
% pn and parameters
pnL=(1-alphaa)*(alphaa/psii)^(alphaa/(1-alphaa))*An*pn^(1/(1-alphaa)); 
psL=(1-alphaa)*(alphaa/psii)^(alphaa/(1-alphaa))*As*ps^(1/(1-alphaa)); 

% wages from labour firm optimality conditions and free movement
wh=psL; 
wl=pnL^(1/(1-eppsilonn))*(eppsilonn/wh)^(eppsilonn/(1-eppsilonn))*(1-eppsilonn); 

% knowing these prices, sustainable demand and labour supply follow from 
% hh side

csh=wh^2*(omegas/ps)^(1+omegas)*((1-omegas)/pn)^(1-omegas)/(1+wh^2*(omegas/ps)^(2*omegas)*((1-omegas)/pn)^(2*(1-omegas)))*B;
%cnh=wh^2*(omegas/ps)^(omegas)*((1-omegas)/pn)^(2-omegas)/(1+wh^2*(omegas/ps)^(2*omegas)*((1-omegas)/pn)^((1-2*omegas)))*B;
hh= wh*(omegas/ps)^(omegas)*((1-omegas)/pn)^(1-omegas)/(1+wh^2*(omegas/ps)^(2*omegas)*((1-omegas)/pn)^(2*(1-omegas)))*B;
% c=c+1;
% f(c)=Ch-csh^omegas*cnh^(1-omegas);

% low skill
csl=wl^2*(omegas/ps)^(1+omegas)*((1-omegas)/pn)^(1-omegas)/(1+wl^2*(omegas/ps)^(2*omegas)*((1-omegas)/pn)^(2*(1-omegas)))*B;
%cnl=wl^2*(omegas/ps)^(omegas)*((1-omegas)/pn)^(2-omegas)/(1+wl^2*(omegas/ps)^(2*omegas)*((1-omegas)/pn)^((1-2*omegas)))*B;
hl=wl*(omegas/ps)^(omegas)*((1-omegas)/pn)^(1-omegas)/(1+wl^2*(omegas/ps)^(2*omegas)*((1-omegas)/pn)^(2*(1-omegas)))*B;

% using labour market clearing and optimal labour input in labour firm:
% lls, lhs, lln, lhn follow

lhn= (pnL/wh*eppsilonn)^(1/(1-eppsilonn))*(hl);

lln=(pnL*(1-eppsilonn)/wl)^(1/(eppsilonn))*lhn;
lhs=hh-lhn;
lls=0;
%lls=hl-lln;

% production labour input good
Ls=lhs; %^eppsilons*lls^(1-eppsilons);
Ln=lhn^eppsilonn*lln^(1-eppsilonn);

% production final good
ys=(alphaa*ps/psii)^(alphaa/(1-alphaa))*As*Ls;


%-- model equations
% sustainable market clearing to pin down ps

c=0;
c=c+1;
f(c)=ys-(lambdaa*csh+(1-lambdaa)*csl);

%--end model block
fprintf('model equations %d', length(f))
fprintf('variables %d', length(symm.variables))

%-- Solve analytically
% S=solve(f==0, symm.variables);

%% - set parameter values
eppsilonn   = 0.5;    % high-skill labour share in unsustainable sector
eppsilons   = 1;    % high-skill labour share in sustainable labour input
omegas      = 0.5;    % preference for sustainable good
alphaa      = 0.5;    % labour share final good production
psii        = 1;      % cost machine production 
B           = 20;     % satiation point
lambdaa     = 0.5;    % share high skill workers
As          = 1;      % initial average productivity of machines in sector s
An          = 2;      % initial average productivity of machines in sector n

%% evaluate model
params      = eval(symm.params);
params_noncomp=eval(symm.params_noncomp);

model_param = subs(f, [symm.params, pn], [params, 1]); % model only dependent on ps
model_param_noncomp= subs(f, [symm.params_noncomp, pn, omegas], [params_noncomp, 1, 0.5]); % model only dependent on ps, omega and B

% solving model including omega and B as fixed
SNum=vpasolve(model_param, [symm.variables]); % numerical solution
S_numana= solve(model_param, symm.variables, 'ReturnConditions', true, 'Real', true)

% solution as function of omega and B
S_numana= solve(model_param_noncomp, symm.variables, 'ReturnConditions', true, 'Real', true)
