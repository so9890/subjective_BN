function [LF_SIM]=solve_LF_VECT(T, list, pol, params,symms, init201519, helper)
%test OPT policy result without target in competitive equilibrium
read_in_params;

%helper=load(sprintf('FB_LF_SIM_NOTARGET_spillover%d.mat', indic.spillovers));
varrs=helper.LF_SIM;
y=log(varrs);

% create new list
syms HL HH lambdaa real
symms.test= [symms.choice(list.choice~='hl'&list.choice~='hh'), HL, HH, lambdaa];
list.test=string(symms.test);

% read in results from social planner
HL= log((params(list.params=='upbarH')-varrs(list.allvars=='hl', :))./(varrs(list.allvars=='hl', :)))';
HH= log((params(list.params=='upbarH')-varrs(list.allvars=='hh', :))./(varrs(list.allvars=='hh', :)))';
hhf=y(list.allvars=='hhf', :)';
hhn=y(list.allvars=='hhn', :)'; 
hhg=y(list.allvars=='hhg', :)';
hln =y(list.allvars=='hln', :)';
hlg=y(list.allvars=='hlg', :)';
hlf=y(list.allvars=='hlf', :)';
C=y(list.allvars=='C', :)';
F=y(list.allvars=='F', :)';
G=y(list.allvars=='G', :)';
Af=y(list.allvars=='Af', :)';
Ag =y(list.allvars=='Ag', :)';
An =y(list.allvars=='An', :)';
sff =y(list.allvars=='sff', :)';
sg =y(list.allvars=='sg', :)';
sn =y(list.allvars=='sn', :)';
S =y(list.allvars=='S', :)';
gammalh =sqrt(varrs(list.allvars=='gammalh', :))';
gammall =sqrt(varrs(list.allvars=='gammall', :))';
wh =y(list.allvars=='wh', :)';
wl=y(list.allvars=='wl', :)';
ws=y(list.allvars=='ws', :)';
pg=y(list.allvars=='pg', :)';
pn=y(list.allvars=='pn', :)';
pee=y(list.allvars=='pee', :)';
pf=y(list.allvars=='pf', :)';
lambdaa=y(list.allvars=='lambdaa', :)';


x0=eval(symms.test);
x0=x0(:);
% f=laissez_faireVECT(sol, params, list, varrs, init201519,T);
% 
% if max(abs(f))>1e-9
%    % equations where results are off
%     ind=1:length(f);
%     pos=ind(abs(f)>1e-9);
%     eqset= floor(pos./T);
%     error('optimal policy does not solve laissez faire. For the relevant equations see eqset');
%     % to evaluate stuff
% end

lb=[];
ub=[];

objf=@(x)objectiveCALIBSCI(x);
constLF=@(x)laissez_faireVECT_fmincon(x, params, list, varrs, init201519,T);
options = optimset('algorithm','active-set','TolCon', 1e-11,'Tolfun',1e-26,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
[x,fval,exitflag,output,lambda] = fmincon(objf,x0,[],[],[],[],lb,ub,constLF,options);

% test solution to 
f=laissez_faireVECT(x, params, list, varrs, init201519,T);
if max(abs(f))>1e-9
    error('LF function does not solve')
end

% save results
LF_SIM=aux_solutionLF_VECT(x, pol, list, symms, varrs, params, T);
end
