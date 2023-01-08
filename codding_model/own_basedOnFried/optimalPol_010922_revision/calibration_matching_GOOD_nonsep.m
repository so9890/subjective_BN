function [x0LF, SL, SP, SR, Sall, Sinit201014, init201014 , Sinit201519, init201519, Sparams, Spol, params, pol, symms, MOM, indexx, list]...
    = calibration_matching_GOOD_nonsep(MOM, symms, list, parsHelp, polhelp, indic)
% function to match moments to model equations
% i.e. solving model plus additional equations for paramters
% to this end, the model is solved 
% 1) from an aggregate pproduction side
% 2) the household probelm and labour producing firms
% 3) research side. 

% To be chosen: 
% An0, Af0, Ag0 (2019)
% thetaf, thetan, thetag : share of high skill in labour input 
% lambdaa to have gov budget = 0 in baseyear
% zh: rationalises wage premium together with thetag, thetaf, thetan
% deltay: weight on energy in production function 
% omegaa: emissions per fossil 

% output
% indexx: map of indices indicating variable transformation for code.
%         contains laissez-faire index and calibration index. 

if polhelp(list.poldir=='tauf')~=0 || polhelp(list.poldir=='taus')~=0 
error('tauf not 0 or taus not =0; have to change calibration and pay attention how env revs are recycled!')
end
%% initaliase stuff needed for diverse solution 

%- calibration research side
syms Af_lag Ag_lag An_lag ...
     sff sg sn  ws gammaa chiis real
symms.calib3 = [Af_lag Ag_lag An_lag ...
     sff sg sn ws chiis gammaa];
 
list.calib3 =string(symms.calib3);

%- all variables: to save base year variables!
syms muu chii hhf hhg hhn hln hlf hlg C F G N Y E Af Ag An hl hh sff sg sn ...
    wh wl pg pn pee pf gammalh gammall wlg wln wlf xn xg xf SGov Emnet A Ch Cl muuh muul...
     tauf taus taul lambdaa Ln Lg Lf SWF S gammas real
symms.allvars= [muu, hhf, hhg, hhn, hln, hlf, hlg, C, F, G, N, Y, E, Af, Ag, An, ...
    hl, hh,  sff, sg, sn, wh, wl, ws, pg, pn, pee, pf,  wlg, wln, wlf, xn, xg, xf, ...
    gammalh, gammall, SGov, Emnet, A, tauf, taus, taul, lambdaa, Ln, Lg, Lf, SWF, S ];
list.allvars  = string(symms.allvars);
symms.allvars_ineq=[muuh, muul, Ch, Cl,symms.allvars];
list.allvars_ineq=string(symms.allvars_ineq);

%- variables and index for laissez faire
symms.choice = [hhf, hhg, hhn, hln, hlf, hlg, C, F, G, Af, Ag, An, hl, hh,  sff, sg, sn,...
             wh, wl, pg, pn, pee, pf, gammalh, gammall, lambdaa, ws, S, gammas];
list.choice  = string(symms.choice);

indexxLF.lab = boolean(zeros(size(list.choice)));
indexxLF.exp = boolean(zeros(size(list.choice)));
indexxLF.sqr = boolean(zeros(size(list.choice)));
indexxLF.oneab = boolean(zeros(size(list.choice)));
indexxLF.Sc = boolean(zeros(size(list.choice)));


indexxLF.lab(list.choice=='hl'| list.choice=='hh')=1;
indexxLF.Sc( list.choice=='S')=1;
indexxLF.exp(  list.choice~='lambdaa'& list.choice~='hl'& list.choice~='hh' & list.choice~='gammas'&list.choice~='gammall'& list.choice~='gammalh'& list.choice~='S' )=1;
indexxLF.sqr(list.choice=='gammall'| list.choice=='gammalh'| list.choice=='gammas' )=1;

%- variables for version with inequality
symms.choice_ineq=[symms.choice(list.choice~='C'), Ch, Cl];
list.choice_ineq =string(symms.choice_ineq);

indexxLF_ineq.lab = boolean(zeros(size(list.choice_ineq)));
indexxLF_ineq.exp = boolean(zeros(size(list.choice_ineq)));
indexxLF_ineq.sqr = boolean(zeros(size(list.choice_ineq)));
indexxLF_ineq.oneab = boolean(zeros(size(list.choice_ineq)));

indexxLF_ineq.lab(list.choice_ineq=='hl'| list.choice_ineq=='hh'| list.choice_ineq=='S')=1;
indexxLF_ineq.exp(list.choice_ineq~='hl'& list.choice_ineq~='hh' & list.choice_ineq~='gammas'&list.choice_ineq~='gammall'& list.choice_ineq~='gammalh'& list.choice_ineq~='S' )=1;
indexxLF_ineq.sqr(list.choice_ineq=='gammall'| list.choice_ineq=='gammalh'| list.choice_ineq=='gammas' )=1;

%- calibration productivity
syms omegaa deltay real
symms.prod= [pn, pg, deltay, omegaa];
list.prod = string(symms.prod);


%- calibration labour firms and HH side
syms thetan thetaf thetag lambdaa chii zh real
symms.calib = [hhn hhg hhf hh hl gammalh gammall wh wl thetan thetaf thetag ...
     lambdaa zh chii];
 list.calib = string(symms.calib);
%-- indexx
indexxcalib.lab = boolean(zeros(size(list.calib)));
indexxcalib.exp = boolean(zeros(size(list.calib)));
indexxcalib.sqr = boolean(zeros(size(list.calib)));
indexxcalib.oneab = boolean(zeros(size(list.calib)));
indexxcalib.Sc = boolean(zeros(size(list.choice)));

indexxcalib.lab(list.calib=='hl'| list.calib=='hh')=1;
indexxcalib.exp(list.calib~='thetan' & list.calib~='thetaf' & list.calib~='thetag'& list.calib~='zh' ...
    & list.calib~='hl' & list.calib~='hh'& list.calib~='gammalh' & list.calib~='gammall')=1;
indexxcalib.sqr(list.calib=='gammall'| list.calib=='gammalh')=1;
indexxcalib.oneab(list.calib=='thetag'| list.calib=='thetan'| ...
    list.calib=='thetaf' | list.calib=='zh') = 1;

% save indices to map 
indexx = containers.Map({'LF', 'calib', 'LF_ineq'}, {indexxLF, indexxcalib, indexxLF_ineq});
 %% First calibration reduced model: ONLY producers' side

pn=log(1);
pg=log(1);
deltay =  log((1-0.4)/0.4);

x0=eval(symms.prod(symms.prod~='omegaa'));
%- solve
f = calibProd(x0, MOM, list, parsHelp);
prodf = @(x)calibProd(x,  MOM, list, parsHelp);
options = optimoptions('fsolve', 'TolFun', 10e-12, 'MaxFunEvals',8e3, 'MaxIter', 3e5, 'Algorithm', 'levenberg-marquardt');%, );%, );%, 'Display', 'Iter', );
[solProd, fval] = fsolve(prodf, x0, options);
trProd=exp(solProd);
trProd(list.prod=='deltay')=1/(1+trProd(list.prod=='deltay'));

% required for next functions
[C, Lnwln, Lgwlg, Lfwlf, pf, F, pn, pg, pee, E, Y, N, G, xn, xg, xf, ...
            AfLf, AgLg, AnLn, omegaa, deltay]=resProd(list, trProd, MOM, parsHelp, polhelp, 'calib');
trProd =eval(symms.prod);
%% Labour side

%- initial guess
hhn = 0.02;
hhg = 0.01;
hhf = 0.02;
zh=0.2;
wl=1;
wh=MOM.whwl*wl;
thetan=0.4;
thetaf=0.4;
thetag=0.5;
    hln = hhn*(1-thetan)/(thetan)*MOM.whwl; % hln
    hlf = hhf*(1-thetaf)/(thetaf)*MOM.whwl; % hlf
    hlg = hhg*(1-thetag)/(thetag)*MOM.whwl; % hlg 
hl = (hln+hlf+hlg)/((1-zh));
hh = (hhn+hhf+hhg)/zh;
gammall = 0;
gammalh = 0;
chii =10;
lambdaa=1;

x0=eval(symms.calib);

%- transform
guess_trans=trans_guess(indexxcalib, x0, parsHelp, list.paramsdir);

%test:
f =calibLabour(guess_trans,  MOM, C, Lnwln, Lgwlg, Lfwlf, pf, F, parsHelp, list, polhelp);

%- solve
Labf = @(x)calibLabour(x,  MOM, C, Lnwln, Lgwlg, Lfwlf, pf, F, parsHelp, list, polhelp);
options = optimoptions('fsolve', 'TolFun', 10e-12, 'MaxFunEvals',8e3, 'MaxIter', 3e5, 'Algorithm', 'levenberg-marquardt');%, );%, );%, 'Display', 'Iter', );
[solLab, fval] = fsolve(Labf, guess_trans, options);
trLab=trans_allo_out(indexxcalib, solLab, parsHelp, list.paramsdir, indic);


%% Research side
%- from previous
cell_par=arrayfun(@char, symms.calib, 'uniform', 0);
Sparams=cell2struct(num2cell(trLab), cell_par, 2);

hhn = trLab(list.calib=='hhn');
hhg = trLab(list.calib=='hhg');
hhf = trLab(list.calib=='hhf');
lambdaa = trLab(list.calib=='lambdaa');
hln = hhn*(1-Sparams.thetan)/(Sparams.thetan)*MOM.whwl; % hln
hlf = hhf*(1-Sparams.thetaf)/(Sparams.thetaf)*MOM.whwl; % hlf
hlg = hhg*(1-Sparams.thetag)/(Sparams.thetag)*MOM.whwl; % hlg 

Lg = hhg.^Sparams.thetag.*hlg.^(1-Sparams.thetag);
Ln = hhn.^Sparams.thetan.*hln.^(1-Sparams.thetan);
Lf = hhf.^Sparams.thetaf.*hlf.^(1-Sparams.thetaf);

% initial guess
% ss=load('calib_initRes');
%%
 sg  = log(0.1);
 ws  = log(3);
 gammaa =log(1);
 sff = log(0.1);
 
 sn  = log(MOM.targethour-exp(sg)-exp(sff)); 
    Af = AfLf/Lf;
    An = AnLn/Ln;
    Ag = AgLg/Lg;
 Af_lag = log(Af/(1.02));
 An_lag = log(An/(1.02));
 Ag_lag = log(Ag/(1.02));
 chiis =log(4);
% rhon=log(1);
% rhog=log(1);
% rhof=log(1);

x0 = eval(symms.calib3);

% x0=ss.x1;
%%

%-- find upper bound on growth rate
% xup=[log(Af), log(Ag), log(An), log(0.1)];
% 
% [ceq]= calibRem_upbar_gammaa(xup, MOM, list, parsHelp,polhelp,  Af, An, Ag);
% modg=@(x)calibRem_upbar_gammaa(x, MOM, list, parsHelp,polhelp,  Af, An, Ag);
% options = optimoptions('fsolve', 'TolFun', 10e-10, 'MaxFunEvals',8e4, 'MaxIter', 3e5, 'Display', 'Iter'); %,'Algorithm', 'levenberg-marquardt');%, );%, );%, 'Display', 'Iter', );
% [x, fval, exitf] = fsolve(modg, xup, options);
% gammaa=exp(x(4));
%%
%- test
%  f= calibRem_nows_fsolve(x, MOM, list, trProd, parsHelp, polhelp, Af, An, Ag);
% 
% % solving model
modF3 = @(x)calibRem_nows_fsolve_GOOD_nonsep(x, MOM, list, trProd, parsHelp, polhelp, Af, An, Ag, lambdaa); 
options = optimoptions('fsolve', 'TolFun', 10e-10, 'MaxFunEvals',8e4, 'MaxIter', 3e5, 'Display', 'Iter','Algorithm', 'levenberg-marquardt');%, );%, );%, 'Display', 'Iter', );
% rr=load(sprintf('calib_2209_sep'), 'x');
rr=load(sprintf('calib_0501_sep'), 'xSc');
[xSc, fval, exitf] = fsolve(modF3, [rr.xSc], options);
%[x, fval, exitf] = fsolve(modF3, x0, options);

save(sprintf('calib_0501_sep'), 'xSc');

%%

  x0=[solLab,xSc];
  list.calibBoth=[list.calib, list.calib3];
%% Jointly calibrate labor and scientists
f=calibLabourAnsSc(x0,  MOM, trProd,  parsHelp, list, polhelp);
modFJoint = @(x)calibLabourAnsSc(x, MOM, trProd,  parsHelp, list, polhelp);
options = optimoptions('fsolve', 'TolFun', 10e-10, 'MaxFunEvals',8e4, 'MaxIter', 3e5, 'Display', 'Iter','Algorithm', 'levenberg-marquardt');%, );%, );%, 'Display', 'Iter', );
[x, fval, exitf] = fsolve(modFJoint, x0, options);

% read in results
trLabJoint=trans_allo_out(indexxcalib, x(1:length(list.calib)), parsHelp, list.paramsdir, indic);
%   Af_lag  = exp(x(list.calib3=='Af_lag'));
%   Ag_lag  = exp(x(list.calib3=='Ag_lag'));
%   An_lag  = exp(x(list.calib3=='An_lag'));
%   chiis    = exp(x(list.calib3=='chiis'));
%   ws       = exp(x(list.calib3=='ws'));
%   gammaa   = exp(x(list.calib3=='gammaa'));
%   sff     = exp(x(list.calib3=='sff'));
%   sg      = exp(x(list.calib3=='sg'));
%   sn      = exp(x(list.calib3=='sn'));
resSciJoint=[exp(x(length(list.calib)+1:end))];
%%
 
% syms gammaa real 
% symms.calib3=[symms.calib3, gammaa];
% get calibrated parameters and policy
[Sparams, Spol, params, pol]=parsSol_GOOD(symms,trProd, trLabJoint, resSciJoint, parsHelp, list, polhelp);

%% save all results 
[x0LF, SL, SP, SR, Sall, Sinit201014, init201014 , Sinit201519, init201519]...
    =fsolution_GOOD(symms, trProd, trLabJoint, resSciJoint, parsHelp, list, polhelp, MOM, indic); 

%% Test if is calibration and baseline model solve LF in baseyear
guess_transLF=trans_guess(indexxLF, x0LF, params, list.params);
indic.noskill=0;
indic.ineq=0;
f=laissez_faire_nows_sep(guess_transLF, params, list, pol, init201014, indic);


if max(abs(f))>1e-9
    error('calibration is not a solution to LF')
else
    fprintf('Hurray!!! LF solves at baseline calibration!!!');
end
end