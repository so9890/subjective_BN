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

%% initaliase stuff needed for diverse solution 

%- lags
syms Ag_lag Af_lag An_lag real
symms.laggs = [Ag_lag, Af_lag, An_lag];
list.laggs= string(symms.laggs);

%- all variables: to save base year variables!
syms chii hhf hhg hhn hln hlf hlg C F G N Y E Af Ag An hl hh sf sg sn wh wl ws pg pn pe pf gammalh gammall wlg wln wlf xn xg xf real
symms.allvars= [hhf, hhg, hhn, hln, hlf, hlg, C, F, G, N, Y, E, Af, Ag, An, hl, hh,  sf, sg, sn, wh, wl, ws, pg, pn, pe, pf,  wlg, wln, wlf, xn, xg, xf, gammalh, gammall];
list.allvars  = string(symms.allvars);

%- addiitonal symbolic variables for first calibration routine
syms thetan thetaf el eh omegaa lambdaa deltay real
symms.calib2 = [hhf hhg hhn C Af Ag An ...
     pg thetan thetaf el eh omegaa lambdaa deltay];
list.calib2  = string(symms.calib2);

%- choiceCALIB variabes: Full calibration
syms  thetag real
symms.choiceCALIB = [  hhf hhg hhn hlf hlg hln F G C Af Ag An hl hh ...
     sf sg sn ws pg wh wl pn pe pf gammalh gammall ...
     Af_lag Ag_lag An_lag  deltay thetan thetaf thetag...
     lambdaa el eh omegaa chii ];
list.choiceCALIB  = string(symms.choiceCALIB);

%- variables for laissez faire
symms.choice = [hhf, hhg, hhn, hln, hlf, hlg, C, F, G, Af, Ag, An, hl, hh,  sf, sg, sn, wh, wl, ws, pg, pn, pe, pf, gammalh, gammall];
list.choice  = string(symms.choice);

%- alternative calibration 
symms.calib = [pg pn hhn hhg hhf hh hl gammalh gammall C deltay thetan thetaf thetag ...
     lambdaa el eh chii];
 list.calib = string(symms.calib);
%% solve for skill : thetag fix and others as initial guess

eleh0   = 0.43;
thetag0 = 0.5; % bounded above and below
thetan0  = 0.2;
thetaf0  = 0.2;
x0=[log((1-thetag0)/thetag0),log(eleh0),log((1-thetaf0)/thetaf0),log((1-thetan0)/thetan0)]; %log((1-zh0)/zh0)];
skillf = @(x)aux_calib_skill(x, MOM, parsHelp,list, polhelp);
options = optimoptions('fsolve', 'TolFun', 10e-16, 'MaxFunEvals',8e3, 'MaxIter', 3e5);%, 'Algorithm', 'levenberg-marquardt');%, );%, );%, 'Display', 'Iter', );
[solNOM] = fsolve(skillf, x0, options);

% in these results thetag is excat, thetaf and thetan are only starting
% values
thetag =1/(1+exp(solNOM(1))); 
eleh   = exp(solNOM(2));
thetaf  = 1/(1+exp(solNOM(3))); % as targets for thetan and thetaf
thetan  = 1/(1+exp(solNOM(4)));

%% First calibration reduced model

%- initial values
Af     = 87; % Fried 
Ag     = 0.9196;
An     = 60; 
% thetan = 0.5; FOLLOW FROM RESULTS ABOVES 
% thetaf = 0.5;
lambdaa= 1;
omegaa = 90;
deltay = 0.15;
el     = 0.5; 
eh     = 0.5; 
hhf    =.02; % hhf
hhg    =.04; % hhg
hhn    =.04; % hhn
pg     = 1; 
C      = 1; 

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

solnew=guess_trans;
options = optimoptions('fsolve', 'TolFun', 10e-12, 'MaxFunEvals',8e5, 'MaxIter', 3e5, 'Algorithm', 'levenberg-marquardt');%'trust-region-dogleg');%'levenberg-marquardt');%, );%, );%, );%, 'Display', 'Iter', )
[sol, fval, exitf] = fsolve(modFF,solnew , options);
%%
solnew=sol;
options = optimoptions('fsolve', 'TolFun', 10e-12, 'MaxFunEvals',8e5, 'MaxIter', 3e5);%, 'Algorithm', 'levenberg-marquardt');%'trust-region-dogleg');%'levenberg-marquardt');%, );%, );%, );%, 'Display', 'Iter', )
[sol, fval, exitf] = fsolve(modFF,solnew , options);
%[solNOM, fval, exitf] = fsolve(modFF, soll.sol(1:end-1), options);
%save('solutionCalib2', 'sol');
%options = optimoptions('fsolve', 'TolFun', 10e-7, 'MaxFunEvals',8e5, 'MaxIter', 3e5);%, 'Algorithm', 'levenberg-marquardt');%'trust-region-dogleg');%'levenberg-marquardt');%, );%, );%, );%, 'Display', 'Iter', )
% [sol, fval, exitf] = fsolve(modFF, guess_trans, options);
%[solNOM, fval, exitf] = fsolve(modFF, sol, options);

allo_trans=trans_allo_out(indexx, sol, parsHelp, list.paramsdir);

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
eh     = allo_trans(list.calib2=='eh');

[hln, hlf, hlg,  hh, hl, Lg, Ln, Lf, pf, pe, pn, G, E, F, N, Y, xn, xf, xg, wh, wl, whg, whn]...
    =aux_calib2(MOM,deltay, hhn, hhg, hhf,zh, zl, el, eh, alphag, alphaf, alphan,  thetag, thetan, thetaf, eppsy, eppse, Ag, An, Af, pg, tauf);
%% remaining variables
sg  = 0.2;
sf  = 0.4;
sn  = 0.4;
ws = 1;
Af_lag = Af;
An_lag = An;
Ag_lag = Ag;
x0= eval(symms.calib3);

Af0=1;Ag0=1;An0=1; % for evaluation only;
chii=10;
params = eval(symms.params);
pol    = eval(symms.pol);
targets = eval(symms.targets);

f= calibRem(log(x0), MOM, list, params, pol,targets, Af, An, Ag,  hhn, hhg, hhf, pg);

% solving model

modF3 = @(x)calibRem(x, MOM, list, params, pol,targets, Af, An, Ag,  hhn, hhg, hhf, pg);
options = optimoptions('fsolve', 'TolFun', 10e-12, 'MaxFunEvals',8e3, 'MaxIter', 3e5);%, 'Algorithm', 'levenberg-marquardt');%, );%, );%, 'Display', 'Iter', );

[sol3, fval, exitf] = fsolve(modF3, log(x0), options);
trans_sol3=exp(sol3);

Af0 = trans_sol3(list.calib3=='Af_lag');
Ag0 = trans_sol3(list.calib3=='Ag_lag');
An0 = trans_sol3(list.calib3=='An_lag'); 

%% final calibration full model

% initial guess => results from separate problems

deltay = allo_trans(list.calib2=='deltay');
thetaf = allo_trans(list.calib2=='thetaf');
thetan = allo_trans(list.calib2=='thetan'); 
thetag = params(list.params=='thetag');
omegaa = allo_trans(list.calib2=='omegaa'); 
lambdaa = allo_trans(list.calib2=='lambdaa'); 
el     = allo_trans(list.calib2=='el'); 
eh     = allo_trans(list.calib2=='eh'); 
Af_lag=trans_sol3(list.calib3=='Af_lag');
Ag_lag=trans_sol3(list.calib3=='Ag_lag');
An_lag=trans_sol3(list.calib3=='An_lag');
 
% variables
Af     = allo_trans(list.calib2=='Af');
Ag     = allo_trans(list.calib2=='Ag');
An     = allo_trans(list.calib2=='An'); 
C      = allo_trans(list.calib2=='C'); 
hhf    = allo_trans(list.calib2=='hhf');
hhg    = allo_trans(list.calib2=='hhg');
hhn    = allo_trans(list.calib2=='hhn');
sf=trans_sol3(list.calib3=='sf');
sg=trans_sol3(list.calib3=='sg');
sn=trans_sol3(list.calib3=='sn');
ws = trans_sol3(list.calib3=='ws');
pg     = allo_trans(list.calib2=='pg'); 
gammalh = 0;
gammall = 0;
chii=10;

[hln, hlf, hlg, hh, hl, Lg, Ln, Lf, pf, pe, pn, G, E, F, N, Y, xn, xf, xg, wh, wl, whg, whn]...
    =aux_calib2(MOM,deltay, hhn, hhg, hhf,zh, zl, el, eh, alphag, alphaf, alphan, ...
                thetag, thetan, thetaf, eppsy, eppse, Ag, An, Af, pg, tauf);
while hh>=upbarH || hl>=upbarH
    hh=hh/10;
    hl=hl/10;
end

x0=eval(symms.choiceCALIB);

% - transforming variables to unbounded variables
%-- index for transformation 
indexx.lab = boolean(zeros(size(list.choiceCALIB)));
indexx.exp = boolean(zeros(size(list.choiceCALIB)));
indexx.sqr = boolean(zeros(size(list.choiceCALIB)));
indexx.oneab = boolean(zeros(size(list.choiceCALIB)));

indexx.exp( list.choiceCALIB~='deltay'& list.choiceCALIB~='gammall'&...
    list.choiceCALIB~='gammalh'& list.choiceCALIB~='thetan' &  ...
    list.choiceCALIB~='thetaf'& list.choiceCALIB~='thetag'...
    & list.choiceCALIB~='hl'& list.choiceCALIB~='hh')=1;
indexx.sqr(list.choiceCALIB=='gammall'| list.choiceCALIB=='gammalh')=1;
indexx.oneab(list.choiceCALIB=='deltay'|list.choiceCALIB=='thetan'|...
    list.choiceCALIB=='thetaf'| list.choiceCALIB=='thetag') = 1;
indexx.lab(list.choiceCALIB=='hl'| list.choiceCALIB=='hh')=1;

%-- transform
guess_trans=trans_guess(indexx, x0,parsHelp, list.paramsdir);

%- test
f=target_equFinal(guess_trans, MOM, parsHelp, list, polhelp, targets);
%% -solve
modF4 = @(x)target_equFinal(x, MOM, parsHelp, list, polhelp, targets);
% options = optimoptions('fsolve', 'TolFun', 10e-6, 'MaxFunEvals',8e3, 'MaxIter', 3e5, 'Algorithm', 'levenberg-marquardt');%);%, );%, );%, 'Display', 'Iter', );
% [sol4, fval, exitf] = fsolve(modF4, guess_trans, options);

soll=guess_trans;
count=0;
exitf=-1;
while exitf<=0 && count<=10
    options = optimoptions('fsolve', 'TolFun', 10e-12, 'MaxFunEvals',8e3, 'MaxIter', 3e5, 'Algorithm', 'levenberg-marquardt');
    [sol4, fval, exitf] = fsolve(modF4, soll, options);
    count=count+1;
    soll=sol4;
end

count=0;
tolfun= 10e-5; % function tolerance increase
while count<=10 
    if exitf>0 % that is: the code has solved
        tolfun= tolfun/10;
    else
        tolfun=tolfun*10;
    end
    soll=sol4;
    options = optimoptions('fsolve', 'TolFun', tolfun, 'MaxFunEvals',8e3, 'MaxIter', 3e5);%, 'Algorithm', 'levenberg-marquardt');
    [sol4, fval, exitf] = fsolve(modF4, soll, options);
    if exitf>=0 % bcs it solves at some point
    count=count+1;
    end
end
%save final round results
soll=sol4;
%tolfun=tolfun*10;
options = optimoptions('fsolve', 'TolFun', tolfun*10, 'MaxFunEvals',8e3, 'MaxIter', 3e5);%, 'Algorithm', 'levenberg-marquardt');
[sol4, fval, exitf] = fsolve(modF4, soll, options);
fprintf('full calibration solved at tolfun%d with exitflag %d', tolfun, exitf); 

save('soll', 'sol4');

%% save results
transsol=trans_allo_out(indexx, sol4, parsHelp, list.paramsdir);
[list.choiceCALIB; transsol]'

thetan=transsol(list.choiceCALIB=='thetan');
thetaf=transsol(list.choiceCALIB=='thetaf');
thetag=transsol(list.choiceCALIB=='thetag');

Ag0=transsol(list.choiceCALIB=='Ag_lag');
An0=transsol(list.choiceCALIB=='An_lag');
Af0=transsol(list.choiceCALIB=='Af_lag');
el=transsol(list.choiceCALIB=='el');
eh=transsol(list.choiceCALIB=='eh');
deltay=transsol(list.choiceCALIB=='deltay');
lambdaa=transsol(list.choiceCALIB=='lambdaa');
omegaa=transsol(list.choiceCALIB=='omegaa');


Ag=transsol(list.choiceCALIB=='Ag');
An=transsol(list.choiceCALIB=='An');
Af=transsol(list.choiceCALIB=='Af');
hhg=transsol(list.choiceCALIB=='hhg');
hhn=transsol(list.choiceCALIB=='hhn');
hhf=transsol(list.choiceCALIB=='hhf');
C=transsol(list.choiceCALIB=='C');

pg=transsol(list.choiceCALIB=='pg');

[hlnt, hlft, hlgt, hh, hl, Lg, Ln, Lf, pf, pe, pn, G, E, F, N, Y, xn, xf, xg, wh, wl, whg, whn]...
=aux_calib2(MOM, deltay, hhn, hhg, hhf,zh, zl, el, eh, alphag, alphaf,...
alphan,  thetag, thetan, thetaf, eppsy, eppse, Ag, An, Af, pg, tauf);

wln     = pn.^(1/(1-alphan)).*(1-alphan).*alphan.^(alphan/(1-alphan)).*An; % price labour input neutral sector
wlg     = pg.^(1/(1-alphag)).*(1-alphag).*alphag.^(alphag/(1-alphag)).*Ag;
wlf     = (1-alphaf)*alphaf^(alphaf/(1-alphaf)).*((1-tauf).*pf).^(1/(1-alphaf)).*Af; 

%% save stuff 
chii=10;
params = eval(symms.params);
pol    = eval(symms.pol);
targets = eval(symms.targets);

%- generate structures
cell_par=arrayfun(@char, symms.choiceCALIB, 'uniform', 0);
strLF=cell2struct(num2cell(transsol), cell_par, 2);

cell_par=arrayfun(@char, symms.params, 'uniform', 0);
strpar=cell2struct(num2cell(params), cell_par, 2);

%% Alternative
%- alternative calibration
while hh>=upbarH || hl>=upbarH
    hh=hh/10;
    hl=hl/10;
end

x0=eval(symms.calib);
%-- index for transformation 
indexx.lab = boolean(zeros(size(list.calib)));
indexx.exp = boolean(zeros(size(list.calib)));
indexx.sqr = boolean(zeros(size(list.calib)));
indexx.oneab = boolean(zeros(size(list.calib)));

indexx.exp( list.calib~='deltay'& list.calib~='gammall'&...
    list.calib~='gammalh'& list.calib~='thetan' &  ...
    list.calib~='thetaf'& list.calib~='thetag' & ...
    list.calib~='hl' & list.calib~='hh')=1;
indexx.sqr(list.calib=='gammall'| list.calib=='gammalh')=1;
indexx.oneab(list.calib=='deltay'|list.calib=='thetan'|...
    list.calib=='thetaf'| list.calib=='thetag') = 1;
indexx.lab(list.calib=='hl'| list.calib=='hh')=1;

%-- transform
guess_trans=trans_guess(indexx, x0,parsHelp, list.paramsdir);

% solve
f= calibration(guess_trans, MOM, list, parsHelp, polhelp);

%% -solve
MOM.lowskill=0.7;
modFfinal = @(x)calibration(x, MOM, list, parsHelp, polhelp);
options = optimoptions('fsolve', 'TolFun', 10e-6, 'MaxFunEvals',8e3, 'MaxIter', 3e5, 'Algorithm', 'levenberg-marquardt');%);%, );%, );%, 'Display', 'Iter', );
[solfin, fval, exitf] = fsolve(modFfinal, guess_trans, options);

soll=solfin;
count=0;
exitf=-1;
while exitf<=0 && count<=10
    [sol4, fval, exitf] = fsolve(modFfinal, soll, options);
    options = optimoptions('fsolve', 'TolFun', 10e-12, 'MaxFunEvals',8e3, 'MaxIter', 3e5, 'Algorithm', 'levenberg-marquardt');
    count=count+1;
    soll=sol4;
end

options = optimoptions('fsolve', 'TolFun', 10e-12, 'MaxFunEvals',8e3, 'MaxIter', 3e5); %, 'Algorithm', 'levenberg-marquardt');%);%, );%, );%, 'Display', 'Iter', );
[solfin, fval, exitf] = fsolve(modFfinal, soll, options);
%save('solutionfinal','solfin'); 
%solST=load('solutionfinal','solfin'); 
%% save results
transsol=trans_allo_out(indexx, solfin, parsHelp, list.paramsdir);
[list.calib; transsol]'
C=transsol(list.calib=='C');
pg=transsol(list.calib=='pg');
hhn=transsol(list.calib=='hhn');
hhg=transsol(list.calib=='hhg');
hhf=transsol(list.calib=='hhf');
gammalh=transsol(list.calib=='gammalh');
gammall=transsol(list.calib=='gammall');
thetan=transsol(list.calib=='thetan');
thetag=transsol(list.calib=='thetag');
thetaf=transsol(list.calib=='thetaf');
el=transsol(list.calib=='el');
eh=transsol(list.calib=='eh');
lambdaa=transsol(list.calib=='lambdaa');
deltay=transsol(list.calib=='deltay');
chii=transsol(list.calib=='chii');
hl=transsol(list.calib=='hl');
hh=transsol(list.calib=='hh');


[muu, pf, Yout, pe, E, N, F, G, omegaa, Af, An, Ag, Lg, Ln, Lf, xf, xg, xn, ...
   SGov, hhD, hlD, hln, hlg, hlf, wlg, wln, wlf,...
  wh, wl] = aux_calibFinal(pn, hh, hl, deltay, eh, el, lambdaa, thetaf, thetag, thetan, C, pg, hhn, hhf, hhg, MOM, list, parsHelp, polhelp);
%- generate structures
cell_par=arrayfun(@char, symms.allvars, 'uniform', 0);
strLF=cell2struct(num2cell(eval(symms.allvars)), cell_par, 2);

%% Test laissez Faire baseyear
%- initial condition as laggs
An_lag=params(list.params=='An0');
Ag_lag=params(list.params=='Ag0');
Af_lag=params(list.params=='Af0');
laggs=eval(symms.laggs);

%-- index for transformation 
indexx.lab = boolean(zeros(size(list.choice)));
indexx.exp = boolean(zeros(size(list.choice)));
indexx.sqr = boolean(zeros(size(list.choice)));
indexx.oneab = boolean(zeros(size(list.choice)));

indexx.lab(list.choice=='hl'| list.choice=='hh')=1;
indexx.exp(list.choice~='hl'& list.choice~='hh' & list.choice~='gammall'& list.choice~='gammalh' )=1;
indexx.sqr(list.choice=='gammall'| list.choice=='gammalh')=1;

x0= eval(symms.choice); 
guess_trans=trans_guess(indexx, x0, params, list.params);
f=laissez_faire(guess_trans, params, list, pol, laggs, targets);


