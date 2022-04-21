%function [An0, Af0, Ag0, thetaf, thetan, thetag, lambdaa]= calibration_matching(MOM, symms, list, parsHelp, polhelp, targets)
% function to match moments to model equations
% i.e. solving model plus additional equations for paramters

% To be chosen: 
% An0, Af0, Ag0 (2019)
% thetaf, thetan, thetag : share of high skill in labour input 
% lambdaa to have gov budget = 0 in baseyear

% data
% Energy consumption share, GDP, Fossil to Green output share

% numeraire: Consumption Good

% balanced budget, 
% skills: match Consoli; or wage premia (then includes)

%% initaliase stuff needed to check LF solves
syms Ag_lag Af_lag An_lag real
symms.laggs = [Ag_lag, Af_lag, An_lag];
laggs=[Ag0, Af0, An0];
list.laggs= string(symms.laggs);
syms hhf hhg hhn hln hlf hlg C F G Af Ag An hl hh sf sg sn wh wl ws pg pn pe pf gammalh gammall real

symms.choice = [hhf, hhg, hhn, hln, hlf, hlg, C, F, G, Af, Ag, An, hl, hh,  sf, sg, sn, wh, wl, ws, pg, pn, pe, pf, gammalh, gammall];
list.choice  = string(symms.choice);

%% solve for skill 

eleh0   = 0.43;
thetag0 = 0.5; % bounded above and below
hhgHH0  = 0.1;
hlgHL0  = 0.1;
x0=[log((1-thetag0)/thetag0),log(eleh0), log((1-hhgHH0)/hhgHH0), log((1-hlgHL0)/hlgHL0)]; %log((1-zh0)/zh0)];
skillf = @(x)aux_calib_skill(x, MOM, parsHelp,list, polhelp);
options = optimoptions('fsolve', 'TolFun', 10e-16, 'MaxFunEvals',8e3, 'MaxIter', 3e5);%, 'Algorithm', 'levenberg-marquardt');%, );%, );%, 'Display', 'Iter', );
[solNOM] = fsolve(skillf, x0, options);

thetag =1/(1+exp(solNOM(1))); 
eleh   = exp(solNOM(2));
hhgHH  = 1/(1+exp(solNOM(3))); % as targets for thetan and thetaf
hlgHL  = 1/(1+exp(solNOM(4)));

%% First calibration reduced model
%- variables
syms hhf hhg hhn C Af Ag An  ...
     pg thetan thetaf el omegaa lambdaa deltay real

symms.calib2 = [hhf hhg hhn C Af Ag An ...
     pg thetan thetaf el omegaa lambdaa deltay];

list.calib2  = string(symms.calib2);

%- initial values
Af     = 1.877; % Fried 
Ag     = 0.9196;
An     = 1; 
thetan = 0.5;
thetaf = 0.5;
lambdaa= 1;
omegaa = 90;
deltay = 0.3;
el     = 1; 
hhf    =.02; % hhf
hhg    =.04; % hhg
hhn    =.04; % hhn
pg     = 1; 
C      = 1; 
MOM.An = 1; 
MOM.SGov=0; 

x0= eval(symms.calib2);

% - transforming variables to unbounded variables
%-- index for transformation 
indexx.lab = boolean(zeros(size(list.calib2)));
indexx.exp = boolean(zeros(size(list.calib2)));
indexx.sqr = boolean(zeros(size(list.calib2)));
indexx.oneab = boolean(zeros(size(list.calib2)));

indexx.lab(list.calib2=='hl'| list.calib2=='hh')=1;
indexx.exp(list.calib2~='deltay' &list.calib2~='thetan' & list.calib2~='thetaf' )=1;
indexx.sqr(list.calib2=='gammall'| list.calib2=='gammalh')=1;
indexx.oneab(list.calib2=='deltay'| list.calib2=='thetan'| list.calib2=='thetaf') = 1;

%-- transform
guess_trans=trans_guess(indexx, x0, parsHelp, list.paramsdir);

%- test
tet=calib2(guess_trans, MOM, list, parsHelp, polhelp, thetag, eleh);

%% - solving model
% saved old solution to be passed to model, not yet transformed back
%soll=load('solutionCalib2.mat');
modFF = @(x)calib2(x, MOM, list, parsHelp, polhelp, thetag, eleh);

options = optimoptions('fsolve', 'TolFun', 10e-12, 'MaxFunEvals',8e5, 'MaxIter', 3e5, 'Algorithm', 'levenberg-marquardt');%'trust-region-dogleg');%'levenberg-marquardt');%, );%, );%, );%, 'Display', 'Iter', )
[sol, fval, exitf] = fsolve(modFF, guess_trans, options);
%[solNOM, fval, exitf] = fsolve(modFF, soll.sol(1:end-1), options);
%save('solutionCalib2', 'sol');
options = optimoptions('fsolve', 'TolFun', 10e-12, 'MaxFunEvals',8e5, 'MaxIter', 3e5);%, 'Algorithm', 'levenberg-marquardt');%'trust-region-dogleg');%'levenberg-marquardt');%, );%, );%, );%, 'Display', 'Iter', )
% [sol, fval, exitf] = fsolve(modFF, guess_trans, options);
[solNOM, fval, exitf] = fsolve(modFF, sol, options);

allo_trans=trans_allo_out(indexx, solNOM, parsHelp, list.paramsdir);

%% Second Calibration Research
syms Af_lag Ag_lag An_lag ...
     sf sg sn ws real
symms.calib3 = [Af_lag Ag_lag An_lag ...
     sf sg sn ws ];
list.calib3 =string(symms.calib3);

%-variables: in baseperiod 2015-2019 relevant to solve 3rd calibration step
Af     = allo_trans(list.calib2=='Af');
Ag     = allo_trans(list.calib2=='Ag');
An     = allo_trans(list.calib2=='An'); 
hhf    = allo_trans(list.calib2=='hhf');
hhg    = allo_trans(list.calib2=='hhg');
hhn    = allo_trans(list.calib2=='hhn');
pg     = allo_trans(list.calib2=='pg'); 
deltay    = allo_trans(list.calib2=='deltay');
thetaf     = allo_trans(list.calib2=='thetaf');
thetan     = allo_trans(list.calib2=='thetan'); 
omegaa     = allo_trans(list.calib2=='omegaa'); 
lambdaa     = allo_trans(list.calib2=='lambdaa'); 
C     = allo_trans(list.calib2=='C'); 
el     = allo_trans(list.calib2=='el'); 


[hln, hlf, hlg, eh, hh, hl, Lg, Ln, Lf, pf, pe, pn, G, E, F, N, Y, xn, xf, xg, wh, wl, whg, whn]...
    =aux_calib2(MOM,deltay, hhn, hhg, hhf,zh, zl, el, eleh, alphag, alphaf, alphan,  thetag, thetan, thetaf, eppsy, eppse, Ag, An, Af, pg, tauf);
%% remaining variables
sg  = 0.2;
sf  = 0.4;
sn  = 0.4;
ws = 1;
Af_lag = An;
An_lag = An;
Ag_lag = Ag;
x0= eval(symms.calib3);

Af0=1;Ag0=1;An0=1; % for evaluation only;
params = eval(symms.params);
pol    = eval(symms.pol);
targets = eval(symms.targets);

f= calibRem(log(x0), MOM, list, params, pol,targets, Af, An, Ag,  hhn, hhg, hhf, pg);

%% - solving model

modF3 = @(x)calibRem(x, MOM, list, params, pol,targets, Af, An, Ag,  hhn, hhg, hhf, pg);
options = optimoptions('fsolve', 'TolFun', 10e-12, 'MaxFunEvals',8e3, 'MaxIter', 3e5);%, 'Algorithm', 'levenberg-marquardt');%, );%, );%, 'Display', 'Iter', );

[sol3, fval, exitf] = fsolve(modF3, log(x0), options);
trans_sol3=exp(sol3);

Af0 = trans_sol3(list.calib3=='Af_lag');
Ag0 = trans_sol3(list.calib3=='Ag_lag');
An0 = trans_sol3(list.calib3=='An_lag'); 
%% save results
params = eval(symms.params);
% pol    = eval(symms.pol);
% targets = eval(symms.targets);
baseLF = eval(sort([symms.calib2, symms.calib3]));
list.baseLF =  string(sort([symms.calib2, symms.calib3]));

%% - test laissez faire
%- preparation for laissez faire function
gammalh = 0;
gammall = 0;
x0= eval(symms.choice);

indexx.lab = boolean(zeros(size(list.choice)));
indexx.exp = boolean(zeros(size(list.choice)));
indexx.sqr = boolean(zeros(size(list.choice)));
indexx.oneab = boolean(zeros(size(list.choice)));

indexx.lab(list.choice=='hl'| list.choice=='hh')=1;
indexx.exp(list.choice~='hl'& list.choice~='hh' & list.choice~='gammall'& list.choice~='gammalh' )=1;
indexx.sqr(list.choice=='gammall'| list.choice=='gammalh')=1;

guess_trans=trans_guess(indexx, x0, params, list.params);

f=laissez_faire(guess_trans, params, list, pol, laggs, targets);

modLF = @(x)laissez_faire(x, params, list, pol, laggs, targets);
options = optimoptions('fsolve', 'TolFun', 10e-12, 'MaxFunEvals',8e3, 'MaxIter', 3e5,'Algorithm', 'levenberg-marquardt');%, 'Algorithm', 'levenberg-marquardt');%, );%, );%, 'Display', 'Iter', );
[sol, fval, exitf] = fsolve(modLF, solsav, options);
%%
solsav=sol;
%% final calibration full model

%- choiceCALIB variabes
syms hhf hhg hhn C Af Ag An ...
     sf sg sn ws pg gammalh gammall deltay ...
     Af_lag Ag_lag An_lag thetan thetaf lambdaa el omegaa real

symms.choiceCALIB = [ hhf hhg hhn C Af Ag An ...
     sf sg sn ws pg gammalh gammall deltay ...
     Af_lag Ag_lag An_lag thetan thetaf lambdaa el omegaa];

list.choiceCALIB  = string(symms.choiceCALIB);

% initial guess => results from separate problems
Af     = allo_trans(list.calib2=='Af');
Ag     = allo_trans(list.calib2=='Ag');
An     = allo_trans(list.calib2=='An'); 
hhf    = allo_trans(list.calib2=='hhf');
hhg    = allo_trans(list.calib2=='hhg');
hhn    = allo_trans(list.calib2=='hhn');
pg     = allo_trans(list.calib2=='pg'); 
deltay = allo_trans(list.calib2=='deltay');
thetaf = allo_trans(list.calib2=='thetaf');
thetan = allo_trans(list.calib2=='thetan'); 
% thetag=params(list.params=='thetag');
omegaa = allo_trans(list.calib2=='omegaa'); 
lambdaa = allo_trans(list.calib2=='lambdaa'); 
C      = allo_trans(list.calib2=='C'); 
el     = allo_trans(list.calib2=='el'); 
Af_lag=trans_sol3(list.calib3=='Af_lag');
Ag_lag=trans_sol3(list.calib3=='Ag_lag');
An_lag=trans_sol3(list.calib3=='An_lag');
sf=trans_sol3(list.calib3=='sf');
sg=trans_sol3(list.calib3=='sg');
sn=trans_sol3(list.calib3=='sn');
ws = trans_sol3(list.calib3=='ws');
gammalh = 0;
gammall = 0;

x0=eval(symms.choiceCALIB);
% - transforming variables to unbounded variables
%-- index for transformation 
indexx.lab = boolean(zeros(size(list.choiceCALIB)));
indexx.exp = boolean(zeros(size(list.choiceCALIB)));
indexx.sqr = boolean(zeros(size(list.choiceCALIB)));
indexx.oneab = boolean(zeros(size(list.choiceCALIB)));

indexx.exp( list.choiceCALIB~='deltay'& list.choiceCALIB~='gammall'&...
    list.choiceCALIB~='gammalh'& list.choiceCALIB~='thetan' &  ...
    list.choiceCALIB~='thetaf' )=1;
indexx.sqr(list.choiceCALIB=='gammall'| list.choiceCALIB=='gammalh')=1;
indexx.oneab(list.choiceCALIB=='deltay'|list.choiceCALIB=='thetan'| list.choiceCALIB=='thetaf') = 1;

%-- transform
guess_trans=trans_guess(indexx, x0,parsHelp, list.paramsdir);

%- test
f=target_equ(guess_trans, MOM, parsHelp, list, polhelp, eleh, thetag);

%% -solve
modF4 = @(x)target_equ(x, MOM, parsHelp, list, polhelp, eleh, thetag);
options = optimoptions('fsolve', 'TolFun', 10e-6, 'MaxFunEvals',8e3, 'MaxIter', 3e5, 'Algorithm', 'levenberg-marquardt');%);%, );%, );%, 'Display', 'Iter', );
[sol4, fval, exitf] = fsolve(modF4, guess_trans, options);

%%
exitfl=0;
count=0;
countmax=100;
save=zeros(2,countmax);

while exitfl<=0 && count<countmax
    count=count+1;
    guess_trans=solNOM;
    save(1,count)=max(fval);
    save(2,count)=exitf;
end