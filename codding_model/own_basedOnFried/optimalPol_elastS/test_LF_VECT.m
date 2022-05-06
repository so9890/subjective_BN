function []=test_LF_VECT(T, list,  params,symms, init201519, helper)
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

% test solution to 
f=laissez_faireVECT(x0, params, list, varrs, init201519,T);

if max(abs(f))>1e-9
    error('LF function does not solve at 1e-9')
else
    fprintf('Solution solves LF problem')
end

end
