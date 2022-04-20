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

%% solve for intermediate and final good producers, prices

% x0=log(50); % guess pg
% prod = @(x)init_calib(x, MOM, parsHelp,list, polhelp);
% options = optimoptions('fsolve', 'TolFun', 10e-12, 'MaxFunEvals',8e3, 'MaxIter', 3e5);%, 'Algorithm', 'levenberg-marquardt');%, );%, );%, 'Display', 'Iter', );
% [sol, fval, exitf] = fsolve(prod, x0, options);
% pg=exp(sol);
% [pf, pe, pn, EY, NY, FY, GY, xgY, xnY, xfY]=aux_calib_Prod(MOM, pg, parsHelp,list, polhelp);

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
thetaf = thetag*0.5;
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
%[solNOM, fval, exitf] = fsolve(modFF, sol, options);
%save('solutionCalib2', 'sol');
options = optimoptions('fsolve', 'TolFun', 10e-12, 'MaxFunEvals',8e5, 'MaxIter', 3e5);%, 'Algorithm', 'levenberg-marquardt');%'trust-region-dogleg');%'levenberg-marquardt');%, );%, );%, );%, 'Display', 'Iter', )
% [sol, fval, exitf] = fsolve(modFF, guess_trans, options);
[solNOM, fval, exitf] = fsolve(modFF, sol, options);

%%
% sol2STR= load('sol2');
% [sol2, fval, exitf] = fsolve(modFF, sol2STR.sol2, options);
% save('sol204', 'sol')
% save_sol=sol.*(1+rand(size(sol))*1e-3);
allo_trans=trans_allo_out(indexx, solNOM, parsHelp, list.paramsdir);

%% Second Calibration full model
syms Af_lag Ag_lag An_lag ...
     sf sg sn real
symms.calib3 = [Af_lag Ag_lag An_lag ...
     sf sg sn ];
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
% f= calibjoint(x, MOM, list, paramss, poll, thetag, eleh)
%% remaining variables
sg  = 0.2;
sf  = 0.4;
sn  = 0.4;
Af_lag = An;
An_lag = An;
Ag_lag = Ag;
x0= eval(symms.calib3);

% - transforming variables to unbounded variables
%-- index for transformation 
%indexx3.exp = boolean(zeros(size(list.calib3)));
%indexx3.exp(list.calib3~='thetan' & list.calib3~='thetaf' )=1;

%-- transform
%guess_trans=trans_guess(indexx3, x0, parsHelp, list.paramsdir);
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


%% final calibration full model

%- choiceCALIB variabes
syms hhf hhg hhn hln hlf hlg C F G Af Ag An ...
     hl hh sf sg sn wh wl ws pg pn pe pf gammalh gammall Y...
     Af0 Ag0 An0 thetan thetaf lambdaa eh el omegaa real

symms.choiceCALIB = [hhf, hhg, hhn, hln, hlf, hlg, C, F, G, ...
    Af, Ag, An, hl, hh,  sf, sg, sn, wh, wl, ws, pg, pn,...
    pe, pf, gammalh, gammall, Y, Af0, Ag0, An0, thetan, thetaf,...
    lambdaa, eh, el, omegaa];

list.choiceCALIB  = string(symms.choiceCALIB);

% - transforming variables to unbounded variables
%-- index for transformation 
indexx.lab = boolean(zeros(size(list.choiceCALIB)));
indexx.exp = boolean(zeros(size(list.choiceCALIB)));
indexx.sqr = boolean(zeros(size(list.choiceCALIB)));
indexx.oneab = boolean(zeros(size(list.choiceCALIB)));

indexx.lab(list.choiceCALIB=='hl'| list.choiceCALIB=='hh')=1;
indexx.exp(list.choiceCALIB~='hl'& list.choiceCALIB~='hh' & list.choiceCALIB~='gammall'&...
    list.choiceCALIB~='gammalh'& list.choiceCALIB~='thetan' & list.choiceCALIB~='thetag'& ...
    list.choiceCALIB~='thetaf' )=1;
indexx.sqr(list.choiceCALIB=='gammall'| list.choiceCALIB=='gammalh')=1;
indexx.oneab(list.choiceCALIB=='thetan'| list.choiceCALIB=='thetag'| list.choiceCALIB=='thetaf') = 1;

%-- transform
guess_trans=trans_guess(indexx, x0, paramss, list);

%- test
f=target_equ(guess_trans, MOM, paramss, list, poll);

exitfl=0;
count=0;
countmax=100;
save=zeros(2,countmax);

while exitfl<=0 && count<countmax
    count=count+1;
    [solNOM, fval, exitf] = fsolve(modFF, guess_trans, options);
    guess_trans=solNOM;
    save(1,count)=max(fval);
    save(2,count)=exitf;
end