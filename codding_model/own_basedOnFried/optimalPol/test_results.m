%test SP result without target in competitive equilibrium

helper=load('OPT_notarget.mat');
varrs=helper.opt_all';
y=log(varrs);

% create new list
syms HL HH real
symms.test= [symms.choice(list.choice~='hl'&list.choice~='hh'), HL, HH];
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
gammalh =sqrt(varrs(list.allvars=='gammalh', :))';
gammall =sqrt(varrs(list.allvars=='gammall', :))';
wh =y(list.allvars=='wh', :)';
wl=y(list.allvars=='wl', :)';
ws=y(list.allvars=='ws', :)';
pg=y(list.allvars=='pg', :)';
pn=y(list.allvars=='pn', :)';
pee=y(list.allvars=='pee', :)';
pf=y(list.allvars=='pf', :)';

x=eval(symms.test);
x=x(:);
f=laissez_faireVECT(x, params, list, varrs, init201519,T);


% to evaluate stuff
list.testVEC=repmat(list.test, 12,1);
list.testVEC=list.testVEC(:);

% position where results are off
ind=1:length(f);
pos=ind(abs(f)>1e-7);